include("sde_functions.jl")

import StochasticDiffEq as SDE
import Plots

S = 100

K = 10 ^ 5
m = 0.2
xbar = helper.fakexbar(S)

R0 = 1
η = 0.25
π = 0.01

p = (S, K, m, xbar, R0, η, π)

x0 = xbar
y0 = ones(Float64, S) / K
u0 = helper.makeu0(x0, y0)

prob = SDE.SDEProblem(helper.f!, helper.g!, 
    u0, (0.0, 1000.0), p, 
    noise_rate_prototype = helper.setA(S))

sol = SDE.solve(prob, SDE.SRA1())
# sol = SDE.solve(prob, SDE.SOSRA())

outx, outz = helper.getsol(sol, S)

Plots.plot(sol.t, outx,
    color = :black, label = false)
Plots.plot(sol.t, outz,
    color = :black, label = false)
Plots.plot(sol.t, outz .* outx,
    color = :black, label = false)
