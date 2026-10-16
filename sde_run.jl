include("sde_functions.jl")

import StochasticDiffEq as SDE
import Plots

S = 10

K = 10^6
m = 10K^-1
xbar = helper.fakexbar(S)

R0 = 10^5 * K^-1
η = 0.5
π = 0.00001

p = (S, K, m, xbar, R0, η, π)

x0 = xbar
y0 = ones(Float64, S) * minimum(xbar) / 10
u0 = helper.makeu0(x0, y0)

# prob = SDE.SDEProblem(helper.f!, helper.g!, 
#     u0, (0.0, 1000.0), p, noise_rate_prototype = helper.setA(S))
# sol = SDE.solve(prob, SDE.SOSRA())
# sol = SDE.solve(prob, SDE.SRA1())
# sol = SDE.solve(prob, SDE.LambaEM())
# sol = SDE.solve(prob, SDE.LambaEulerHeun())

prob = SDE.SplitSDEProblem(helper.f_stiff!, helper.f_non!, helper.g!,
    u0, (0.0, 1000.0), p, noise_rate_prototype = helper.setA(S))
# sol = SDE.solve(prob, SDE.SKenCarp(autodiff = SDE.AutoFiniteDiff()))
sol = SDE.solve(prob, SDE.ImplicitEM(autodiff = SDE.AutoFiniteDiff()))

outx, outz = helper.getsol(sol, S)
timeidx = helper.gettime(sol)

# Plots.plot(sol.t[timeidx], outx[timeidx, :],
#     color = :black, label = false)
Plots.plot(sol.t[timeidx], outz[timeidx, :],
   color = :black, label = false)
# Plots.plot(sol.t[timeidx], (outx .* outz)[timeidx, :],
#     color = :black, label = false)

