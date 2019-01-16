module mcmc_functions

using Distributed
using Statistics
using Distributions
using LinearAlgebra
using StatsBase
using Lazy
using RLEVectors
using ProgressMeter
using DataFrames

export tapply_mean
function tapply_mean(subs, val, sz=(maximum(subs),))
    A = zeros(eltype(val), sz...)
    counter = zeros(Int64, sz...)
    for i = 1:length(val)
        @inbounds A[subs[i]] += val[i]
        @inbounds counter[subs[i]] += 1
    end
    A./counter
end

export propose_prop
function propose_prop(prop; tune_var=0.05)
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

    NIW = Array{Dict{String,Array{Float64,N} where N}}(undef,nw,nt,nm)
    for i in 1:nw
        for j in 1:nt
            """
            Count the number of observations in each class
                for the current walker, current temperature
            """
            rles = @> begin
                    z[i,j]
                    @>sort()
                    @>StatsBase.rle()
                end
            nz = zeros(nm)
            nz[rles[1]] = rles[2]

            """
            Compute the d-dimensional sample mean for each class
                for the current walker, current temperature
            """
            # xbar is an array of d dictionaries
            xbar = DataFrames.colwise(x -> tapply_mean(z[i,j], x), dat)

            """
            Draw NIW random variables
            If there are enough (dm) observations in the class, sample from the posterior
                Otherwise, just draw from the prior
            Store the NIW in an array of dictionaries
            """
            for m in 1:nm
                if nz[m] >= dm
                    # Draw from the posterior (I don't have the additional ifelse that is in my R code here)
                    Sigma = rand(
                                InverseWishart(kappa0[m] + nz[m],
                                               round.(Psi0[:,:,m] * kappa0[m] +
                                                  (Matrix(dat[z[i,j] .== m,:]) .- map(x -> x[m], xbar)')' *
                                                       (Matrix(dat[z[i,j] .== m,:]) .- map(x -> x[m], xbar)') +
                                                       (kappa0[m] * nz[m]) / (kappa0[m] + nz[m]) *
                                                       (map(x -> x[m], xbar) - mu0[:,m]) *  (map(x -> x[m], xbar) - mu0[:,m])'; digits=6)
                            )
                    )
                    mu = rand(
                            MvNormal(
                                (kappa0[m] * mu0[:,m] + nz[m] * map(x -> x[m], xbar)) / (kappa0[m] - nz[m]),
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
                p[:,m] = pdf(distns[m], Matrix(dat)') * prop[i,j][m]
            end

            # WRITE SOMETHING THAT SAMPLES RANDOMLY IF WEIGHTS ARE EQUAL?
            # RIGHT NOW ITS BIASED TOWARDS SAMPLING JUST 1.
            # IDK IF IT REALLY MATTERS THOUGH
            zout[i,j] = mapslices(x -> sample(1:nm, pweights(x)), p, dims = 2)[:,1]
        end
    end

    (zout, NIW)
end

export loglike
function loglike(param::Tuple{Array{Dict{String,Array{Float64,N} where N},1},Array{Float64,1},Array{Int64,1}}, dat)
    """
    Define a log-likelihood function to get passed around
    * param is a tuple containing mus, sigmas, prop, and z for a given walker, temperature
    * specifically, param = ( Dict(NIW), Array(prop), Array(z) )
    """
    NIW, prop, z = param
    nm = size(prop, 1)
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
            ll_normal[m] = -Inf # used to be 0
        end
    end

    # Compute the log-likelihood of the class labels given the NIW param and mixing proportions
    distns = map(x -> MvNormal(get(x, "mu", 0), get(x, "Sigma", 0)), NIW)
    p = Array{Float64,2}(undef,n,nm)
    for m in 1:nm
        p[:,m] = pdf(distns[m], Matrix(dat)') * prop[m]
    end

    # Convert density to probabilities
    prob = @> begin
            p'
            @>exp.() # consider rounding "x" to avoid inexact error?
            #@>round.(digits=25)
            @>>mapslices(x -> all(x .== 0.0) ? rep(1, times = nm) : x, dims = 1)
            @>>mapslices(x -> x/sum(x), dims = 1)
    end

    z = Int64.(z) # is this my problem????
    mult_distns = mapslices(x -> Multinomial(1,x), prob, dims = 1)
    ll_mult = zeros(n)
    for i in 1:n
        occ = Int64.(zeros(nm))
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
        ## for debugging, delete later
        #NIW, prop, z = param
        #nm = size(prop, 1)
        #n, dm = size(dat)
#
        ## Compute the log-likelihood of the data given the NIW param and class labels
        #ll_normal = zeros(nm)
        #for m in 1:nm
            #if sum(z .== m) >= dm
                #ll_normal[m] =
                    #sum(
                        #logpdf(
                                #MvNormal(get(NIW[m], "mu", 0), get(NIW[m], "Sigma", 0)),
                                #Matrix(dat[z .== m,:])'
                        #)
                    #)
            #else
                #ll_normal[m] = -Inf # used to be 0
            #end
        #end
#
        ## Compute the log-likelihood of the class labels given the NIW param and mixing proportions
        #distns = map(x -> MvNormal(get(x, "mu", 0), get(x, "Sigma", 0)), NIW)
        #q = Array{Float64,2}(undef,n,nm)
        #for m in 1:nm
            #q[:,m] = pdf(distns[m], Matrix(dat)') * prop[m]
        #end
#
        ## Convert density to probabilities
        #prob = @> begin
                #q'
                #@>exp.() # consider rounding "x" to avoid inexact error?
                ##@>round.(digits=25)
                #@>>mapslices(x -> all(x .== 0.0) ? rep(1, times = nm) : x, dims = 1)
                #@>>mapslices(x -> x/sum(x), dims = 1)
        #end
#
        #z = Int64.(z) # is this my problem????
        #mult_distns = mapslices(x -> Multinomial(1,x), prob, dims = 1)
        #ll_mult = zeros(n)
        #for i in 1:n
            #occ = Int64.(zeros(nm))
            #occ[z[i]] = 1
            #ll_mult[i] = logpdf(mult_distns[i], occ)
        #end
        (ll(param, dat), p)
    end
end

export lnlikes_lnpriors
function lnlikes_lnpriors(dat, param, alpha, ll, lp)
  """
  * dat is the n x dm data matrix (currently treating it like a DataFrame I think)
  * param is a multidimensional nw x nt array that contains the tuple ( Dict(NIW), Array(prop), Array(z) )
        for each walker at each temperature
  * alpha is an nm array of prior weights for the mixing proportions
  * ll = log-likelihood function
  * lp = log prior function
  """
    nw, nt = size(param)

    # If there is capacity to parallelize, then parallelize
    if nworkers() > 1
        # unfolds param column-wise
        parr = reshape(param, nw*nt)
        lls_lps = pmap(x -> lnlike_lnprior(dat, x, alpha, ll, lp), parr, batch_size=max(div(nw*nt,4*nworkers()),1))
        lls = reshape(Float64[l for (l,p) in lls_lps], (nw, nt))
        lps = reshape(Float64[p for (l,p) in lls_lps], (nw, nt))
    else
        lls = zeros(nw, nt)
        lps = zeros(nw, nt)
        for i in 1:nw
            for j in 1:nt
                lls[i,j], lps[i,j] = lnlike_lnprior(dat, param[i,j], alpha, ll, lp)
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
    nw, nt = size(prop)
    prop_old_given_new = zeros(nw, nt)
    prop_new_given_old = zeros(nw, nt)
    for i in 1:nw
        for j in 1:nt
            prop_old_given_new[i,j] = logpdf(Dirichlet(prop_star[i,j]), prop[i,j])
            prop_new_given_old[i,j] = logpdf(Dirichlet(prop[i,j]), prop_star[i,j])
        end
    end

    return prop_old_given_new - prop_new_given_old
end

export make_mcmc_move
function make_mcmc_move(dat, param, hyp, alpha, ll, lp, betas)
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
    for i in 1:nw
        for j in 1:nt
            map!(x -> x, param[i,j][1], NIW[i,j,:])
            map!(x -> x, param[i,j][3], z[i,j])
        end
    end
    """
    * llps is an nw x nt array of log-likelihoods of given parameters
    * lpps is an nw x nt array of log priors of given parameters
    """
    llps, lpps = lnlikes_lnpriors(dat, param, alpha, ll, lp)

    # Propose a temporary version of param reflecting prop_star to pass to lnlikes_lnpriors
    prop_star = propose_prop(prop; tune_var=0.05)
    param_star = deepcopy(param)
    for i in 1:nw
        for j in 1:nt
            map!(x -> x, param_star[i,j][2], prop_star[i,j])
        end
    end
    prop_llps, prop_lpps = lnlikes_lnpriors(dat, param_star, alpha, ll, lp)

    # Account for asymmetry in mixing proportion proposal
    lhaste = lhastings(prop, prop_star)

    new_param = Array{Tuple{Array{Dict{String,Array{Float64,N} where N},1},
                            Array{Float64,1},
                            Array{Int64,1}},2}(undef, nw, nt)
    new_llps = zeros(size(llps)...)
    new_lpps = zeros(size(lpps)...)
    for i in 1:nw
        for j in 1:nt
            """
            For each temperature and walker, grab the corresponding parameter set
            Compute the difference between proposal and current likelihoods
              as well as the difference between priors,
              and also account for asymmatric proposals
            The target is augmented by the temperature (specifically, only the likelihood part)
            Smaller beta -> Larger temperature -> flatter target surface
            """
            lpacc = betas[j]*(prop_llps[i,j] - llps[i,j]) + prop_lpps[i,j] - lpps[i,j] + lhaste[i,j]

            if log(rand()) < lpacc
                # Make update
                new_param[i,j] = param_star[i,j]
                new_llps[i,j] = prop_llps[i,j]
                new_lpps[i,j] = prop_lpps[i,j]
            else
                # Only retain the Gibbs updates
                new_param[i,j] = param[i,j]
                new_llps[i,j] = llps[i,j]
                new_lpps[i,j] = lpps[i,j]
            end
        end
    end
    (new_param, new_llps, new_lpps)
end

export make_tswap
function make_tswap(param, lnlikes, lnpriors, betas)
    """
    * twaps counts the number of swaps that happen between various temps
    """
    nw, nt = size(param)
    new_param = copy(param)
    new_lnlikes = copy(lnlikes)
    new_lnpriors = copy(lnpriors)
    tswaps = zeros(Int, nt-1)

    for k in nt:-1:2 # From nt to 2, decreasing by 1
        for j in 1:nw
            # Randomly select another walker index
            jj = rand(1:nw)

            """
            log-likelihood for current walker, temperature
            and some other randomly selected walker, adjacent temperature
            """
            llhigh = new_lnlikes[j,k]
            lllow = new_lnlikes[jj, k-1]

            """
            Computing the difference in temperatures of two walkers selected
             (higher temp - lower temp)
            """
            diff_beta = betas[k] - betas[k-1]

            """
            probability of acceptance = [beta1-beta2] * (log-likelihood2 - log-likelihood1)
            This maintains detailed balance, and note that if diff_beta were too large,
              then we would almost never accept. That is why we propose swaps between consecutive temperatures
            """
            lpacc = diff_beta*(lllow - llhigh)

            # I think i am wasting some memory in this function with the update step
            if log(rand()) < lpacc
                """
                If accepted, swap the log-likelihood, log prior, and parameter estimates
                  between the chains at the two given temperatures
                """
                # Copy these so I don't write over them and lose them
                phigh = new_param[j,k]
                lphigh = new_lnpriors[j,k]

                # Make the swap
                new_param[j,k] = new_param[jj,k-1]
                new_lnlikes[j,k] = new_lnlikes[jj,k-1]
                new_lnpriors[j,k] = new_lnpriors[jj,k-1]

                new_param[jj,k-1] = phigh
                new_lnlikes[jj,k-1] = llhigh
                new_lnpriors[jj,k-1] = lphigh

                # Note that I have make a swap at this between temps k and k-1
                tswaps[k-1] += 1
            end # If not accepted, do nothing
        end
    end
    new_param, new_lnlikes, new_lnpriors, tswaps
end

"""
    tevolve(swapfraction, betas, tc)

Returns `new_betas`, which are evolved according to a rule to equalise
temperature swap rates over a timescale `tc`.

The rule is formulated in terms of a quantity `x[i]`, which is the logit
transform of `beta[i]` relative to its neighbours (recall that `beta` is the
*inverse* temperature, so decreases as `i` increases):

    x[i] = log(beta[i] - beta[i+1]) - log(beta[i-1] - beta[i])

`x` maps the range between `beta[i+1]` and `beta[i-1]` to `-Inf` to `Inf`.  We
evolve `x` via

    x_new[i] = x[i] + 2.0*(swapfraction[i] - swapfraction[i-1])/(tc*(swapfraction[i] + swapfraction[i-1]))

where `swapfraction[i]` measures the fraction of accepted swaps between
temperature `i` and temperature `i+1` (i.e. between chain `i` and the
next-highest temperature chain) and `tc` is a "time constant" that controls the
(exponential) rate of convergence of `x`.

This is similar to the evolution rule in [Vousden, Farr, & Mandel
(2016)](https://ui.adsabs.harvard.edu/#abs/2016MNRAS.455.1919V/abstract).

"""

export tevolve
function tevolve(swapfraction, betas, tc)
    """
    * tc is the "time constant" controlling the convergence rate
    * swapfraction is tswaps/nw (tswaps output from make_tswap)
    """
    new_betas = copy(betas)

    # This seems problematic if number of temperatures is < 3
    for i in 2:size(betas, 1)-1
        if swapfraction[i] == 0 && swapfraction[i-1] == 0
            # Do nothing since we don't know which way to move
        else
            x = log(betas[i] - betas[i+1]) - log(betas[i-1] - betas[i])
            xnew = x + 2.0*(swapfraction[i] - swapfraction[i-1])/(tc*(swapfraction[i] + swapfraction[i-1]))
            exnew = exp(xnew)
            new_betas[i] = (new_betas[i+1] + exnew*new_betas[i-1])/(1.0 + exnew)
        end
    end

    new_betas
end

export run_mcmc
function run_mcmc(dat, param, hyp, alpha, loglike, logprior, betas, nstep, burnin)
    """
      The function will return `(chain, chainloglike, chainlogprior, new_betas)` where
        each returned chain value has an extra dimension appended counting steps of the
        chain (so `chain` is of shape `(ndim, nwalkers, ntemp, nstep)`, for example).
      * dat is an n x nd array of observations
      * alpha is an nm array of hyperparameters for the mixing proportions
      * loglike = function giving the log-likelihood for parameters.
      * logprior = function giving the log-prior for parameters.
      * betas = an [ntemp,] array of inverse temperatures.
      * burnin = number of steps used to tune the temperature (to be thrown away)
      * nstep = the number of steps of the already tuned mcmc
      * nd = number of dimensions
      * nw = number of walkers
      * nt = number of temperatures
      * param contains current parameter estimates for each dimension
            across walkers and temperatures
      * the time constant is one quarter of burnin
    """

    # Instantiate each object to the correct size
    nw, nt = size(param)
    n, nd = size(dat)

    @assert size(betas,1)==nt "size(betas,1) must be equal to the number of temperatures"

    # compute the log likelihood and log prior at current parameters
    # lnlikes, lnpriors = lnlikes_lnpriors(dat, param, alpha, loglike, logprior)

    # The burn-in period
    @showprogress 1 "Computing for burn-in..." for i in 1:burnin
        # Propose new everything, accept or reject, store results
        param, lnlikes, lnpriors = make_mcmc_move(dat, param, hyp, alpha, loglike, logprior, betas)

        # Propose swaps between various chains at various temperatures
        # Also track (during this burnin period) how often and which walkers swap
        param, lnlikes, lnpriors, tswaps = make_tswap(param, lnlikes, lnpriors, betas)
        swapfraction = tswaps / nw
        tc = 10.0 + i/10.0 # tc > 10 always, but asymptotes to "current length / 10"

        # Adaptively tune the temperatures (looks a lot like formula for a spring)
        betas = tevolve(swapfraction, betas, tc)
    end

    """
    Pre-allocating memory, and setting up thinning (defaults to thin=1 though)
    Chain stores parameter estimates, chainlnlike the log-likelihoods, and
      chainlnprior the log priors for each dimension, walker, and temperature
    """
    chain = Array{
                    Tuple{Array{Dict{String,Array{Float64,N} where N},1},
                    Array{Float64,1},
                    Array{Int64,1}}
            }(undef, nw, nt, nstep)
    chainlnlike = zeros(nw, nt, nstep)
    chainlnprior = zeros(nw, nt, nstep)

    @showprogress 1 "Computing for main Markov chain..."  for i in 1:nstep
      """
      Alternate between usual MCMC stuff and swapping temperatures
      Temperatures no longer evolve (to satisfy detailed balance)
      May also want to adjust -- multiple MCMC updates in between each swap step
      """
        param, lnlikes, lnpriors = make_mcmc_move(dat, param, hyp, alpha, loglike, logprior, betas)
        param, lnlikes, lnpriors, _ = make_tswap(param, lnlikes, lnpriors, betas)

        chain[:,:,i] = param
        chainlnlike[:,:,i] = lnlikes
        chainlnprior[:,:,i] = lnpriors
    end

    chain, chainlnlike, chainlnprior, betas
end

end # End the module
