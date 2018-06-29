using Distributions
using Plots
using Gadfly
using DataFrames
using DataFramesMeta
using RLEVectors
using Lazy

const n = 1500
const m = 3
const d = 2
prop = vcat(1/3,1/3,1/3)
sampsize = round.(Int, 1500.*prop)

mu = vcat(hcat(-1, 1), hcat(1, 1), hcat(0, -sqrt(2)))
sigma = Array[vcat(hcat(0.65, 0.7794), hcat(0.7794, 1.55)),
              vcat(hcat(0.65, -0.7794), hcat(-0.7794, 1.55)),
              diagm(vcat(2,0.2))]
x = Array{Float64,2}(sum(sampsize),2)

en = cumsum(sampsize)
strt = en .- sampsize .+ 1
for i in 1:size(mu,1)
  dist = MvNormal(mu[i,:], sigma[i])
  x[strt[i]:en[i],:] = rand(dist, sampsize[i])'
end

df = @> begin
  x
  @>>convert(DataFrame)
  @transform(group = rep(collect(1:3), each = sampsize))
end

Gadfly.plot(df, x="x1", y="x2", Geom.point, color = "group")
