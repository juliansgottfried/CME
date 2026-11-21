import DelimitedFiles as DF
import Plots

include("sde_functions_7.jl")

xbar = DF.readdlm("xbar.csv", ',')
S = size(xbar)[1]
κ = 117
π = 0.0455
R0s = 2 .^ (0:1:10)
η = 0

for i in eachindex(R0s)
    helper.runner(i, S, xbar, κ, π, R0s, η)
end
