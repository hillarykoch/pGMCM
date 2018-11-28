module mcmc_functions

using Distributed
using Statistics
using Distributions
using LinearAlgebra
using StatsBase

export propose_prop
function propose_prop(prop, tune_var = 0.75)
    """
    * prop is an nw x nt array of mixing proportions, each with length-nm entries
    * tune_var is an adaptive tuning variance for the proposal
        (this should be updated so that each walker tracks its own tuning variance)
    * nm is the cluster number
    * nw is the number of walkers
    * nt is the number of temperatures
    * Possibly need to be storing the dimension of data as well
    """
    nw, nt = size(prop)
    nm = size(prop[1])[1]
    prop_out = copy(prop)

    # For each walker, at each temperature,
    #   make a new proposal for mixing proportions
    for i in 1:nw
        for j in 1:nt
            log_prop = log.(prop[i,j])
            log_proposal = rand(MvNormal(log_prop, Matrix{Float64}(I,nm,nm) * sqrt(tune_var)))
            proposal = exp.(log_proposal)
            prop_out[i,j] = proposal / sum(proposal)
        end
    end
    prop_out
end

export make_gibbs_update
function make_gibbs_update(dat, hyp, z, prop)
    """
    * dat is an n x dm data matrix (maybe need it to be data frame)
    * hyp is a tuple with
        hyp[1] = vector of kappa0, and sum(kappa0) = n
        hyp[2] = dm x nm array of mu0
        hyp[3] = dm x dm x nm array of Psi0
    * z is an int array of size nw x nt, where each entry is of length n
        classifying each observation for each walker and temperature
    * prop is an nw x nt array, where each entry is of length nm
    """
    nw, nt = size(prop)
    nm = size(prop[1])[1]
    kappa0, mu0, Psi0 = hyp
    n, dm = size(dat)

    for i in 1:nw
        for j in 1:nt
            """
            Count the number of observations in each class
                for the current walker, current temperature
            """
            rles = @> begin
                    z[i,j]
                    @>sort()
                    @>rle()
                end
            nz = zeros(nm)
            nz[rles[1]] = rles[2]

            """
            Compute the d-dimensional sample mean for each class
                for the current walker, current temperature
            """
            # xbar is an array of d dictionaries
            xbar = colwise(x -> tapply(x, z[i,j], mean), dat)

            """
            Draw NIW random variables
            If there are enough (dm) observations in the class, sample from the posterior
                Otherwise, just draw from the prior
            Store the NIW in an array of dictionaries
            """
            NIW = Array{Dict{String,Array{Float64,N} where N}}(undef,nw,nt,nm)
            for m in 1:nm
                if nz[m] >= dm
                    # Draw from the posterior (I don't have the additional ifelse that is in my R code here)
                    Sigma = rand(
                                InverseWishart(kappa0[m] + nz[m],
                                               Psi0[:,:,m] * kappa0[m] +
                                                   (Matrix(dat[z[i,j] .== m,:]) .- hcat(map(x -> get(x,m,0), xbar))')' *
                                                        (Matrix(dat[z[i,j] .== m,:]) .- hcat(map(x -> get(x,m,0), xbar))') +
                                                    (kappa0[m] * nz[m]) / (kappa0[m] + nz[m]) *
                                                    (map(x -> get(x,m,0), xbar) - mu0[:,m]) *  (map(x -> get(x,m,0), xbar) - mu0[:,m])'
                                )
                    mu = rand(
                            MvNormal(
                                (kappa0[m] * mu0[:,m] + nz[m] * map(x -> get(x,m,0), xbar)) / (kappa0[m] - nz[m]),
                                Sigma / (kappa0[m] + nz[m])
                            )
                    )
                    NIW[i,j,m] = Dict("mu" => mu, "Sigma" => Sigma)
                else
                    # Draw from the prior
                    Sigma = rand(
                                InverseWishart(kappa0[m],
                                               Psi0[:,:,m] * kappa0[m]
                                )
                    mu = rand(
                            MvNormal(
                                mu0[:,m],
                                Sigma / kappa0[m]
                            )
                    )
                    NIW[i,j,m] = Dict("mu" => mu, "Sigma" => Sigma)
                end
            end
        end
    end

    """
    Draw new cluster labels
    Store in an n x nw x nt array
    """
    zout = copy(z)
    for i in 1:nw
        for j in 1:nt
            distns = map(x -> MvNormal(x["mu"], x["Sigma"]), NIW[i,j,:])
            p = Array{Float64,2}(undef,n,nm)
            for m in 1:nm
                p[:,m] = pdf(distns[m], Matrix(dat)') * prop[i,j]
            end
            zout[i,j] = mapslices(x -> sample(1:nm, pweights(x)), p, dims = 2)[:,1]
        end
    end

    (zout, NIW)
end

export loglike
function loglike(param, dat)
    """
    Define a log-likelihood function to get passed around
    * param is a tuple containing mus, sigmas, prop, and z for a given walker, temperature
    * specifically, param = ( Dict(NIW), Array(prop), Array(z) )
    """
    NIW, prop, z = param
    nm = size(prop)
    n, dm = size(dat)

    # Compute the log-likelihood of the data given the NIW param and class labels
    ll_normal = zeros(nm)
    for m in 1:nm
        if sum(z .== m) >= dm
            ll_normal[m] =
                sum(
                    logpdf(
                            MvNormal(get(NIW[m], "mu", 0), get(NIW[m], "Sigma", 0)),
                            Matrix(dat[z .== m,:])'
                    )
                )
        else
            ll_normal[m] = 0
        end
    end

    # Compute the log-likelihood of the class labels given the NIW param and mixing proportions
    distns = map(x -> MvNormal(x["mu"], x["Sigma"]), NIW)
    p = Array{Float64,2}(undef,n,nm)
    for m in 1:nm
        p[:,m] = pdf(distns[m], Matrix(dat)') * prop[m]
    end

    # Convert density to probabilities
    prob = @> begin
            p'
            @>>mapslices(x -> all(x .== 0) ? rep(1, times = nm) : x, dims = 1)
            @>>mapslices(x -> x/sum(x), dims = 1)
        end

    mult_distns = mapslices(x -> Multinomial(1,x), prob, dims = 1)
    ll_mult = zeros(n)
    for i in 1:n
        occ = zeros(nm)
        occ[z[i]] = 1
        ll_mult[i] = logpdf(mult_distns[i], occ)
    end
    return sum(ll_mult) + sum(ll_normal)
end

export logprior
function logprior(prop, alpha)
    """
    Define a log-likelihood function to get passed around
    * alpha is an nm array of hyperparameters for the mixing proportions
    * prop is an nm array of a current sample of mixing proportions
    """
    return logpdf(Dirichlet(alpha), prop)
end

export lnlike_lnprior
function lnlike_lnprior(dat, param, alpha, ll, lp)
    """
    For the mixing proportions at a given walker, temperature
        evaluates log-likelihood and log prior at the proposal
    This is:
    1. log-likelihood of data given all param (except mixing prop)
    2. log-likelihood of labels given param and mixing prop
    3. Dirichlet prior on the mixing prop
    Returns Inf if log prior is Inf at the parameters
    """
    NIW, prop, z = param
    p = lp(prop, alpha)
    if p == -Inf
        (-Inf, -Inf)
    else
        (ll(param, dat), p)
    end
end

export lnlikes_lnpriors
function lnlikes_lnpriors(dat, param, alpha, ll, lp)
  """
  * dat is the n x dm data matrix (currently treating it like a DataFrame I think)
  * param is a multidimensional nt x nw array that contains the tuple ( Dict(NIW), Array(prop), Array(z) )
        for each walker at each temperature
  * alpha is an nm array of prior weights for the mixing proportions
  * ll = log-likelihood function
  * lp = log prior function
  """
    nw, nt = size(param)

    # If there is capacity to parallelize, then parallelize
    if length(workers()) > 1
        parr = Tuple{Array{Dict{String,Array{Float64,N} where N},3},Array{Float64,1},Array{Int64,1}}[]
        for j in 1:nw
            for i in 1:nt
                push!(parr, param[i,j])
            end
        end
        lls_lps = pmap(x -> lnlike_lnprior(dat, x, alpha, ll, lp), parr, batch_size=div(nw*nt,(4*length(workers()))))
        lls = reshape(Float64[l for (l,p) in lls_lps], (nw, nt))
        lps = reshape(Float64[p for (l,p) in lls_lps], (nw, nt))
    else
        lls = zeros(nw, nt)
        lps = zeros(nw, nt)
        for i in 1:nt
            for j in 1:nw
                lls, lps = lnlike_lnprior(dat, param[i,j], alpha, ll, lp)
            end
        end
    end
    (lls, lps)
end

export lhastings
function lhastings(prop, prop_star)
    """
    * prop, prop_star are each nw x nt arrays,
        with length-nm entries
    """
    for i in 1:nw
        for j in 1:nt

        end
    end
    prop_old_given_new = logpdf(Dirichlet(prop_star), prop)
    prop_new_given_old = logpdf(Dirichlet(prop), prop_star)

    return prop_old_given_new - prop_new_given_old
end

function make_update(dat, param, hyp, alpha, llps, lpps, qs, ll, lp, betas)
    nw, nt = size(param)
    NIW = map(x -> x[1], param)
    prop = map(x -> x[2], param)
    z = map(x -> x[3], param)

    """
    1. Make gibbs updates of mus, Sigmas, and cluster labels
        * z is an int nw x nt dimensional array, where each entry is of length n
            classifying each observation for each walker and temperature
        * NIW is an an nw x nt array of dictionaries, each with nm keys "mu" and "Sigma"
    2. Propose new mixing proportions
        * prop_star is an nw x nt array of proposed mixing proportions,
            where each entry is of length nm
        * tune_var needs to be updated eventually, right now it is static
    3. Compute the log-likelihood and log prior at the prop_star proposal
    4. Compute the hastings ratio for the proposed mixing proportions
    """
    z, NIW = make_gibbs_update(dat, hyp, z, prop)
    # Update param to reflect gibbs updates
    for i in 1:nt
        for j in 1:nw
            param[i,j][1] = NIW[i,j]
            param[i,j][3] = z
        end
    end

    prop_star = propose_prop(prop, tune_var = 0.75)
    # Propose a temporary version of param reflecting prop_star to pass to lnlikes_lnpriors
    param_temp = copy(param)
    for i in 1:nt
        for j in 1:nw
            param[i,j][2] = prop_star[i,j]
        end
    end
    prop_llps, prop_lpps = lnlikes_lnpriors(dat, param_star, alpha, ll, lp)


    """
    1. propose new ps, zs (the parameters)
    2. compute the log-likelihood and log prior at the proposal
    """
    prop_ps, zs = propose(ps, qs)
    prop_llps, prop_lpps = lnlikes_lnpriors(prop_ps, ll, lp)

    new_ps = zeros(size(ps)...)
    new_llps = zeros(size(llps)...)
    new_lpps = zeros(size(lpps)...)

    for i in 1:nw
        for j in 1:nt
          """
          For each temperature and walker, grab the corresponding parameter set
          Compute the difference between proposal and current likelihoods
            as well as the difference between priors
          And also a term that is probably the Hastings ratio part for asymmetric proposals?
          The target is augmented by the temperature
          Smaller beta -> Larger temperature -> flatter target surface
          """
            lpacc = betas[j]*(prop_llps[i,j] - llps[i,j]) + prop_lpps[i,j] - lpps[i,j] + nd*log(zs[i,j])

            # If log(uniform RV) < log(probability of acceptance)
            if log(rand()) < lpacc
                new_ps[:,i,j] = prop_ps[:,i,j]
                new_llps[i,j] = prop_llps[i,j]
                new_lpps[i,j] = prop_lpps[i,j]
            else
                new_ps[:,i,j] = ps[:,i,j]
                new_llps[i,j] = llps[i,j]
                new_lpps[i,j] = lpps[i,j]
            end
        end
    end

    (new_ps, new_llps, new_lpps)
end
