import StochasticDiffEq as SDE
import StatsBase
import Plots

include("sde_functions_2.jl")

S = 10

K = 10^6
m = 10K^-1
xbar = helper.fakexbar(S)

R0 = 0
η = 0.5
π = 0.01

p = (S, K, m, xbar, R0, η, π)

x0 = xbar
y0 = ones(Float64, S) * minimum(xbar) / 10
u0 = helper.makeu0(x0, y0)

prob = SDE.SDEProblem(helper.f!, helper.g!,
    u0, (0.0, 10^7), p, noise_rate_prototype = helper.setA(S))
sol = SDE.solve(prob, SDE.LambaEM(), maxiters = 10^7)

outx, outy = helper.getsol(sol, S)
timeidx = helper.gettime(sol)

Plots.plot(outx[timeidx, :])
Plots.plot(outy[timeidx, :])

Plots.plot(sum(eachcol(outx[timeidx, :])))

Plots.plot(xbar, StatsBase.mean(eachrow(outx)))
Plots.plot(xbar .* (1 .- xbar) ./ (K * m + 1), StatsBase.var(eachrow(outx)))
