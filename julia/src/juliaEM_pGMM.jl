#using Distributions

function jfpGMM{T1<:Number, T2<:Number}(x::Array{Float64, 2},
  prop::Vector{Float64},
  mu::Array{Float64, 2},
  sigma::Array{Array, 1},
  k::Int64,
  df::T1,
  lambda::T2,
  itermax::Int64 = 300,
  tol::Float64 = 1e-06
  )
  # mu a k by d matrix
  # sigma a [d,d,k] array
  const n::Int64 = size(x, 1)
  const d::Int64 = size(x, 2)

  delta::Float64 = 1.0
  prop_old = copy(prop)
  mu_old = deepcopy(mu)
  sigma_old = deepcopy(sigma)

  stp::Int64 = 0
  pdf_est = Array{Float64,2}(n,k)
  prob0 = Array{Float64,2}(n,k)
  while stp < itermax
    # E step
    for i in 1:k
      tmp_mu = mu_old[i,:]
      tmp_sigma = sigma_old[i]
      pdf_est[:,i] = pdf(MvNormal(tmp_mu, PDMat(cholfact(Hermitian(tmp_sigma)))),x')
      prob0[:,i] = pdf_est[:,i].*prop_old[i]
    end
    h_est = prob0./sum(prob0,2)

    # M step
    mu_new = *(h_est',x)./sum(h_est,1)'
    sigma_new = deepcopy(sigma_old)
    for i in 1:k
      sigma_new[i] = *((x-reshape(kron(mu[i,:], rep(1,times = n)), (n,d)))',
                      diagm(h_est[:,i]),
                      (x-reshape(kron(mu[i,:], rep(1,times = n)), (n,d))))./sum(h_est[:,i])
    end
    prop_new = vec((sum(h_est,1).-lambda*df)./(n-k*lambda*df))

    if any(prop_new .< 1e-04)
      prop_new[prop_new .< 1e-04] = 0
    end
    prop_new = prop_new./sum(prop_new)

    delta = sum(abs(prop_new .- prop_old))
    if delta < tol
      break
    end

    if any(prop_new .== 0)
        keepidx = find(x -> x > 0, prop_new)
        k = length(keepidx)

        prop_old = prop_new[keepidx]
        mu_old = mu_new[keepidx,:]
        sigma_old = sigma_new[keepidx]
        pdf_est = pdf_est[:,keepidx]
        prob0 = h_est[:,keepidx]
        h_est = h_est[:,kepidx]
        delta = 1
      else
        prop_old = prop_new
        mu_old = mu_new
        sigma_old = sigma_new
    end

    if k <= 1
      break
    end
    stp += 1
  end
  tag = mapslices(indmax, h_est, 2)
  return Dict("k" => k,
              "prop" => prop_old,
              "mu" = mu_old,
              "sigma" => sigma_old,
              "pdf_est" => pdf_est,
              "ll" => sum(log(sum(prob0,1)))),
              "cluster" => tag)
end
