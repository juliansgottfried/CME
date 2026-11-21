import StochasticDiffEq as SDE
import DiffEqCallbacks as CB
import DelimitedFiles as DF
import Plots

include("sde_functions_7.jl")

xbar = sort(DF.readdlm("xbar.csv", ','), dims = 1, rev = true)
xbar = xbar[1:20]
S = size(xbar)[1]
κ = 117
π = 0.05
R0 = κ
η = κ
p = (S, xbar, κ, π, R0, η)

x0 = xbar
y0 = zeros(Float64, S)
u0 = helper.makeu0(x0, y0)

prob = SDE.SDEProblem(helper.f!, helper.g!,
    u0, (0.0, 10.0), p, noise_rate_prototype = helper.setA(S))
sol = SDE.solve(prob, SDE.ImplicitEM(),
    callback = CB.PositiveDomain())
outx, outy = helper.getsol(sol, S)
timeidx = helper.gettime(sol, 1000)

Plots.plot(outx[timeidx, :], label = false)
Plots.plot(outy[timeidx, :], label = false)
