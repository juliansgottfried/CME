include("sde_functions.jl")

import StochasticDiffEq as SDE
import Plots

S = 10

K = 10^6
m = 10K^-1
xbar = helper.fakexbar(S)

R0 = 10
η = 0.5
π = 0.01

p = (S, K, m, xbar, R0, η, π)

x0 = xbar
y0 = ones(Float64, S) * minimum(xbar) / 10
u0 = helper.makeu0(x0, y0)

prob = SDE.SplitSDEProblem(helper.f_stiff!, helper.f_non!, helper.g!,
    u0, (0.0, 100000.0), p, noise_rate_prototype = helper.setA(S))
sol = SDE.solve(prob, SDE.ImplicitEM(autodiff = SDE.AutoFiniteDiff()))
# sol = SDE.solve(prob, SDE.SKenCarp(autodiff = SDE.AutoFiniteDiff()))

outx, outz = helper.getsol(sol, S)
timeidx = helper.gettime(sol)

Plots.plot(sol.t[timeidx], outx[timeidx, S],
    color = :black, label = false)
Plots.plot(sol.t[timeidx], outz[timeidx, S],
   color = :black, label = false)
Plots.plot(sol.t[timeidx], (outx .* outz)[timeidx, S],
    color = :black, label = false)


xbar = 0.1
K = 10^6
m = 10K^-1
prob = SDE.SDEProblem(helper.testf, helper.testg,
    log(xbar / (1 - xbar)), 
    (0.0, 100000000000.0), 
    (K, m, xbar), 
    noise_rate_prototype = 1.0)
sol = SDE.solve(prob, SDE.SOSRA())
timeidx = helper.gettime(sol)
outx = 1 ./ (1 .+ exp.(-sol.u))
Plots.plot(sol.t[timeidx], outx[timeidx],
    color = :black, label = false)
StatsBase.mean(outx)
xbar * (1 - xbar) / (K * m + 1)
StatsBase.var(outx)
