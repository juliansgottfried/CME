import StochasticDiffEq as SDE
import DiffEqCallbacks as CB
import DelimitedFiles as DF
import Plots

include("sde_functions_8.jl")

xbar = DF.readdlm("xbar.csv", ',')
xbar = xbar[xbar .> 0.005]
Xbar = sum(xbar)
κ = 117
π = 0.0455
R0 = [0; 2 .^ (0:1:10)]
η = 10κ

i = 1
collector = helper.overbins(i, Xbar, xbar, κ, π, R0, η, 100)

Plots.histogram(collector[:, 3], color = :white, label = false)
Plots.histogram(collector[:, 4], color = :white, label = false)
Plots.histogram(collector[:, 4] ./ collector[:, 3], color = :white, label = false)

x = collector[:, 3:2:size(collector)[2]]
y = collector[:, 4:2:size(collector)[2]]
DF.writedlm("simulations/x_$i.csv", x, ',')
DF.writedlm("simulations/y_$i.csv", y, ',')
