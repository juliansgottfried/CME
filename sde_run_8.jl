import StochasticDiffEq as SDE
import DiffEqCallbacks as CB
import DelimitedFiles as DF
import Plots

include("sde_functions_8.jl")

xbar = DF.readdlm("xbar.csv", ',')
xbar = xbar[xbar .> 0.005]
Xbar = sum(xbar)
κ = 117

π = 0.045
R0 = 0.0001:((1000000 - 0.0001) / 29):1000000
η = 1000
# η = 103448.3
# η = 310344.8
# η = 827586.2

for i in 1:30
    collector = helper.overbins(i, Xbar, xbar, κ, π, R0, η, 42)

    #=
    Plots.plot(collector[:, 63], color = :black, label = false)
    Plots.plot(collector[:, 64], color = :black, label = false)
    Plots.plot(collector[:, 64] ./ collector[:, 63], color = :black, label = false, ylims = (0, 1))

    Plots.histogram(collector[:, 3], color = :white, label = false)
    Plots.histogram(collector[:, 4], color = :white, label = false)
    Plots.histogram(collector[:, 4] ./ collector[:, 3], color = :white, label = false)
    =#

    x = collector[:, 3:2:size(collector)[2]]
    y = collector[:, 4:2:size(collector)[2]]

    # DF.writedlm("simulations/x_$i.csv", x, ',')
    # DF.writedlm("simulations/y_$i.csv", y, ',')

    DF.writedlm("simulations_etaverysmall/z_$i.csv", y ./ x, ',')
end