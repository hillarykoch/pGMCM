cd("/Users/hillarykoch/Box Sync/School/research - Qunhua/Project_2/pGMCM/julia/src")

using CSV
using Gadfly
using DataFrames
using CategoricalArrays
using Distributed
using LinearAlgebra
using Lazy
using Distributions
using RLEVectors
using Random
if Distributed.nworkers() > 1
        # everywhere macro loads onto all workers
        @everywhere include("mcmc_functions.jl")
        @everywhere using .mcmc_functions
else
        include("mcmc_functions.jl")
        using .mcmc_functions
end


# Truth:
# mu1 = (0,0), Sigma1 = ( (1, 0),
#                         (0, 1))
# mu2 = (4,4), Sigma2 = ( (1.9, 0.7),
#                         (0.7, 1.9))
# pi = (.4,.6)

nw = 2 # number of walkers per temp
nt = 3 # number of temperatures
dm = 2 # dimension of data
mu1 = [0,0]
mu2 = [4,4]
mu3 = [-2.5,6]
sig1 = Matrix{Float64}(I,dm,dm)*2
sig2 = Matrix{Float64}(I,dm,dm)*1.2 .+ .7
sig3 = reshape([.85, .6, .6,.95], (2,2))
prop = [.2,.45, .35]
nm = size(prop)[1] # Number of clusters
n = 1000

# set seed
Random.seed!(1234);

# Simulate the labels/data
zprop = @as zprop Multinomial(n, prop) begin
        @> rand(zprop)
        @> rep(collect(1:1:nm); each = zprop)
end;

dat = hcat(
        hcat(rand(MvNormal(mu1, sig1), sum(zprop .== 1)),
           rand(MvNormal(mu2, sig2), sum(zprop .== 2)),
           rand(MvNormal(mu3, sig3), sum(zprop .== 3)))',
        zprop
        );

df1 = @> begin
        DataFrame(dat, [:x, :y, :cluster])
        @> categorical!(:cluster)
end;
#Gadfly.plot(df1, x=:x, y=:y, color = :cluster, Geom.point)

# Choose initial values
mu0 = hcat(mu1, mu2, mu3)'
nu0 = kappa0 = round.(prop*n)
alpha = prop
Psi0 = reshape(hcat(sig1, sig2, sig3), (dm,dm,nm))
hyp = (kappa0, mu0, Psi0)

# Initial estimates at NIW, prop, and z for each walker and each temperature
param = Array{Tuple{Array{Dict{String,Array{Float64,N} where N},1},Array{Float64,1},Array{Int64,1}}}(undef, (nw, nt))
for i in 1:nw
        for j in 1:nt
                dictionary = Dict{String,Array{Float64,N} where N}[]
                for m in 1:nm
                        Sigma = rand(InverseWishart(500, 500 * 2 * Matrix{Float64}(I,dm,dm)))
                        mu = rand(MvNormal(mu0[m,:], Sigma / 500))
                        push!(dictionary, Dict("mu" => mu, "Sigma" => Sigma))
                end
                prop = rand(Dirichlet(alpha))
                z = rand(1:nm, n)
                param[i,j] = (dictionary, prop, z)
        end
end
z = map(x -> x[3], param)
prop = map(x -> x[2], param)

betas = collect(range(1, stop=0.001, length=nt))
nstep = 4000
burnin = 1000

ll = mcmc_functions.loglike
lp = mcmc_functions.logprior

chain, chainlnlike, chainlnprior, betas =
        mcmc_functions.run_mcmc(df1[[:x,:y]], param, hyp, alpha, ll, lp, betas, nstep, burnin);

# Get mu estimates
norm_chain = map(x -> x[1], chain[:,1,:])
mu_ests = [
        [
        @> begin
                norm_chain
                @>> map(x -> x[1])
                @>> map(x -> get(x, "mu", 0)[1])
                @> mean
        end
        ,
        @> begin
                norm_chain
                @>> map(x -> x[1])
                @>> map(x -> get(x, "mu", 0)[2])
                @> mean
        end
        ],
        [
        @> begin
                norm_chain
                @>> map(x -> x[2])
                @>> map(x -> get(x, "mu", 0)[1])
                @> mean
        end
        ,
        @> begin
                norm_chain
                @>> map(x -> x[2])
                @>> map(x -> get(x, "mu", 0)[2])
                @> mean
        end
        ],
        [
        @> begin
                norm_chain
                @>> map(x -> x[3])
                @>> map(x -> get(x, "mu", 0)[1])
                @> mean
        end
        ,
        @> begin
                norm_chain
                @>> map(x -> x[3])
                @>> map(x -> get(x, "mu", 0)[2])
                @> mean
        end
        ]]

# Plot mu chains
mudf = @> begin
        norm_chain
        @>> map(x -> x[1])
        @>> map(x -> get(x, "mu", 0)[1]) # Can change this one for other walkers
        @> adjoint
        @> DataFrame
end
Gadfly.plot(mudf, x=vcat(1:nrow(mudf)...), y=:x1, Geom.line)

# Get Sigma estimates
Sigma_ests = [
        [
        @> begin
                norm_chain
                @>> map(x -> x[1])
                @>> map(x -> get(x, "Sigma", 0)[1])
                @> mean
        end
        ,
        @> begin
                norm_chain
                @>> map(x -> x[1])
                @>> map(x -> get(x, "Sigma", 0)[3])
                @> mean
        end
        ,
        @> begin
                norm_chain
                @>> map(x -> x[1])
                @>> map(x -> get(x, "Sigma", 0)[4])
                @> mean
        end
        ],
        [
        @> begin
                norm_chain
                @>> map(x -> x[2])
                @>> map(x -> get(x, "Sigma", 0)[1])
                @> mean
        end
        ,
        @> begin
                norm_chain
                @>> map(x -> x[2])
                @>> map(x -> get(x, "Sigma", 0)[3])
                @> mean
        end
        ,
        @> begin
                norm_chain
                @>> map(x -> x[2])
                @>> map(x -> get(x, "Sigma", 0)[4])
                @> mean
        end
        ],
        [
        @> begin
                norm_chain
                @>> map(x -> x[3])
                @>> map(x -> get(x, "Sigma", 0)[1])
                @> mean
        end
        ,
        @> begin
                norm_chain
                @>> map(x -> x[3])
                @>> map(x -> get(x, "Sigma", 0)[3])
                @> mean
        end
        ,
        @> begin
                norm_chain
                @>> map(x -> x[3])
                @>> map(x -> get(x, "Sigma", 0)[4])
                @> mean
        end
        ]]

# Plot Sigma chains
Sigmadf = @> begin
        norm_chain
        @>> map(x -> x[1])
        @>> map(x -> get(x, "Sigma", 0)[1]) # Can change this one for other walkers
        @> adjoint
        @> DataFrame
end
Gadfly.plot(Sigmadf, x=vcat(1:nrow(Sigmadf)...), y=:x1, Geom.line)

# Get prop estimates
prop_chain = map(x -> x[2], chain[:,1,:])
prop_est = [
        @> begin
                prop_chain
                @>> map(x -> x[1])
                @> mean
        end
        ,
        @> begin
                prop_chain
                @>> map(x -> x[2])
                @> mean
        end
        ,
        @> begin
                prop_chain
                @>> map(x -> x[3])
                @> mean
        end]

# Plot prop chains
propdf = @> begin
        prop_chain
        @>> map(x -> x[1])
        @> adjoint
        @> DataFrame
end
Gadfly.plot(propdf, x=vcat(1:nrow(propdf)...), y=:x1, Geom.line)

# Get cluster labels
z_chain = map(x -> x[3], chain[:,1,:])
z_est = @> begin
                hcat(hcat(z_chain[1,:]...), hcat(z_chain[2,:]...))
                @>> mapslices(sort; dims = 2)
                @>> mapslices(StatsBase.rle; dims = 2)
                @>> map(x -> x[1][argmax(x[2])])
                @> reshape((n,))
end;

# Compute number of correct classifications
sum(z_est[collect(1:1:sum(zprop.==1))] .== 1) +
        sum(z_est[collect((sum(zprop.==1) + 1):1:(sum(zprop.==1) + sum(zprop .== 2)))] .== 2) +
        sum(z_est[collect((sum(zprop.==1) + sum(zprop.==2) + 1):1:n)] .== 3)

# Plot data with estimated classes
Gadfly.plot(df1, x=:x, y=:y, color = categorical(z_est), Geom.point)
