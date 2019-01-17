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
include("singlegibbs.jl")
using .singlegibbs



# Truth:
# mu1 = (0,0), Sigma1 = ( (1, 0),
#                         (0, 1))
# mu2 = (4,4), Sigma2 = ( (1.9, 0.7),
#                         (0.7, 1.9))
# pi = (.4,.6)

dm = 2 # dimension of data
mu1 = [0,0]
mu2 = [4,4]
sig1 = Matrix{Float64}(I,dm,dm)*2
sig2 = Matrix{Float64}(I,dm,dm)*1.2 .+ .7
prop = [.4,.6]
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
           rand(MvNormal(mu2, sig2), sum(zprop .== 2)))',
        zprop
        );

df1 = @> begin
        DataFrame(dat, [:x, :y, :cluster])
        @> categorical!(:cluster)
end;
#Gadfly.plot(df1, x=:x, y=:y, color = :cluster, Geom.point)

# Choose priors values
mu0 = hcat(mu1, mu2)'
nu0 = kappa0 = round.(prop*n)
alpha = prop
Psi0 = reshape(hcat(sig1, sig2), (dm,dm,nm))
hyp = (kappa0, mu0, Psi0)

# Starting values for NIW, prop, and z
dictionary = Dict{String,Array{Float64,N} where N}[]
for m in 1:nm
        Sigma = rand(InverseWishart(500, 500 * 2 * Matrix{Float64}(I,dm,dm)))
        mu = rand(MvNormal(mu0[m,:], Sigma / 500))
        push!(dictionary, Dict("mu" => mu, "Sigma" => Sigma))
end
prop = [.5, .5]
z = zprop
param = (dictionary, prop, z)

nstep = 100
burnin = 10
betas = [1, 2] # dont need this bc no temperatures being used

ll = loglike
lp = logprior

chain, chainlnlike, chainlnprior, betas =
        singlegibbs.run_mcmc(df1[[:x,:y]], param, hyp, alpha, ll, lp, betas, nstep, burnin)


norm_chain = map(x -> x[1], chain[:,1,:])

#mu1_chain = map(x -> get(x[1], "mu", 0), norm_chain)
#mu_ests = mapslices(mean, mu1_chain; dims = 2)
#mu1_est = mu_ests[:,1]
#mu2_est = mu_ests[:,2]

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
        ]]

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
        end]

propdf = @> begin
        prop_chain
        @>> map(x -> x[1])
        @> adjoint
        @> DataFrame
end
Gadfly.plot(propdf, x=vcat(1:nrow(propdf)...), y=:x1, Geom.line)

z_chain = map(x -> x[3], chain[:,1,:])
#z_est = @>begin
                #hcat(hcat(z_chain[1,:]...), hcat(z_chain[2,:]...))
                #@>> mapslices(sort; dims = 2)
                #@>> mapslices(StatsBase.rle; dims = 2)
                #@>> map(x -> argmax(x[2]))
#end
