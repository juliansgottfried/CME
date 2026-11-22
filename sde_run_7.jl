import StochasticDiffEq as SDE
import DiffEqCallbacks as CB
import DelimitedFiles as DF
import Plots

include("sde_functions_7.jl")

xbar = DF.readdlm("xbar.csv", ',')
S = size(xbar)[1]
κ = 117
π = 0.0455
R0s = 2 .^ (0:1:10)
η = κ

for i in eachindex(R0s) helper.runner(i, S, xbar, κ, π, R0s, η) end

p = (S, xbar, κ, π, R0s[1], η)

x0 = xbar
y0 = π * xbar
u0 = helper.makeu0(x0, y0)

prob = SDE.SDEProblem(helper.f!, helper.g!,
    u0, (0.0, 1.0), p, noise_rate_prototype = helper.setA(S))
sol = SDE.solve(prob, SDE.SOSRA(),
    callback = CB.PositiveDomain())

outx, outy = helper.getsol(sol, S)
timeidx = helper.gettime(sol, 42)

# Plots.plot(outx[timeidx, :], color = :black, label = false)
# Plots.plot(outy[timeidx, :], color = :black, label = false)
Plots.plot(outy[timeidx, :] ./ outx[timeidx, :], color = :black, label = false)
