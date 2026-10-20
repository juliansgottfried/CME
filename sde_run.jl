import StochasticDiffEq as SDE
import DiffEqCallbacks as CB
import StatsBase
import StatsPlots
import Plots

include("sde_functions.jl")

S = 100

K = 10^6
m = 10K^-1
xbar = helper.fakexbar(S)

R0 = 2
η = 0.5
π = 0.01

isy = false

p = (S, K, m, xbar, R0, η, π, isy)

x0 = xbar
y0 = ones(Float64, S) * minimum(xbar) / 10
u0 = helper.makeu0(x0, y0, isy)

prob = SDE.SDEProblem(helper.f!, helper.g!,
    u0, (0.0, 10^6), p, noise_rate_prototype = helper.setA(S))
sol = SDE.solve(prob, SDE.SOSRA(),
    isoutofdomain = (u, p, t) -> any(x -> x < 0, u))
# sol = SDE.solve(prob, SDE.SOSRA(),
#     callback = CB.PositiveDomain())

outx, outyorz = helper.getsol(sol, S)
timeidx = helper.gettime(sol)

StatsPlots.density(outyorz[timeidx, :], color = :black, label = false)

#=
Plots.plot(sol.t[timeidx], outx[timeidx, :], color = :black, label = false)
if isy
    Plots.plot(sol.t[timeidx], outyorz[timeidx, :], color = :black, label = false)
else
    Plots.plot(sol.t[timeidx], outx[timeidx, :] .* outyorz[timeidx, :], color = :black, label = false)
end
if isy
    Plots.plot(sol.t[timeidx], outyorz[timeidx, :] ./ outx[timeidx, :], color = :black, label = false)
else
    Plots.plot(sol.t[timeidx], outyorz[timeidx, :], color = :black, label = false)
end
=#
