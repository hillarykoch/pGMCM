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

# Choose initial values
mu0 = hcat(mu1, mu2)'
nu0 = kappa0 = round.(prop*n)
alpha = prop
Psi0 = reshape(hcat(sig1, sig2), (dm,dm,nm))
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
                prop = rand(Dirichlet([.5,.5]))
                z = rand(1:nm, n)
                param[i,j] = (dictionary, prop, z)
        end
end
z = map(x -> x[3], param)
prop = map(x -> x[2], param)

betas = collect(range(1, stop=0.001, length=nt))
nstep = 2000
burnin = 1000

ll = loglike
lp = logprior

chain, chainlnlike, chainlnprior, betas =
        mcmc_functions.run_mcmc(df1[[:x,:y]], param, hyp, alpha, ll, lp, betas, nstep, burnin)


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

z_chain = map(x -> x[3], chain[:,1,:])
z_est = @>begin
                hcat(hcat(z_chain[1,:]...),hcat(z_chain[2,:]...))
                @>> mapslices(sort; dims = 2)
                @>> mapslices(rle; dims = 2)
                @>> map(x -> argmax(x[2]))
end


# Iris data set analysis
"""
# Instantiate initial parameters, etc, so that I can test functions
#
#
## The data
#iris = CSV.read(joinpath(dirname(pathof(DataFrames)), "../test/data/iris.csv"));
#dat = iris[filter(x -> !(x in [:Species]), names(iris))]
#
## Indexing
#nw = 2 # walkers
#nt = 3 # temperatures
#nm = 3 # clusters
#n, dm = size(dat)
#
## Initial estimates at NIW, prop, and z for each walker and each temperature
#param = Array{Tuple{Array{Dict{String,Array{Float64,N} where N},1},Array{Float64,1},Array{Int64,1}}}(undef, (nw, nt))
#for i in 1:nw
        #for j in 1:nt
                #dictionary = Dict{String,Array{Float64,N} where N}[]
                #for m in 1:nm
                        #Sigma = rand(InverseWishart(100, Matrix{Float64}(I,dm,dm)))
                        #mu = rand(MvNormal(100*DataFrames.colwise(mean, dat), Sigma / 100))
                        #push!(dictionary, Dict("mu" => mu, "Sigma" => Sigma))
                #end
                #prop = rand(Dirichlet([.35,.3,.35]))
                #z = rand(1:nm, n)
                #param[i,j] = (dictionary, prop, z)
        #end
#end
#z = map(x -> x[3], param)
#prop = map(x -> x[2], param)
#
#
## hyperparameters
#kappa0 = [50, 50, 50]
#mu0 = rand(MvNormal(Matrix{Float64}(I,dm,dm)), nm)
#Psi0 = Array{Float64}(undef, (dm, dm, nm))
#for m in 1:nm
        #Psi0[:,:,m] = rand(InverseWishart(100, Matrix{Float64}(I,dm,dm)))
#end
#
#hyp = (kappa0, mu0, Psi0)
#alpha = [.4,.3,.3]
#betas = collect(range(1, stop=0, length=nt))
#nstep = 1000
#burnin = 100
#
#ll = loglike
#lp = logprior
#
#chain, chainlnlike, chainlnprior, betas =
        #mcmc_functions.run_mcmc(dat, param, hyp, alpha, loglike, logprior, betas, nstep, burnin)
"""
