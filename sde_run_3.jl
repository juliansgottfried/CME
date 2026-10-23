import StochasticDiffEq as SDE
import DiffEqCallbacks as CB
import StatsBase
import Plots

include("sde_functions_3.jl")

xbar = 0.01
K = 10^7
κ = 134.8684
m = κ / K

R0 = 20
η = 0.5
π = 0.05

p = (K, m, xbar, R0, η, π)

x0 = xbar
y0 = π * xbar
u0 = [x0, y0]

prob = SDE.SDEProblem(helper.f!, helper.g!, u0, (0.0, 10^8), p, 
    noise_rate_prototype = ones(Float64, 2, 2), maxiters = 10^8)
sol = SDE.solve(prob, SDE.ImplicitEM(), callback = CB.PositiveDomain())
# sol = SDE.solve(prob, SDE.LambaEM(), callback = CB.PositiveDomain())

outx, outy = helper.getsol(sol)
timeidx = helper.gettime(sol, 1000)

# plot x and y
Plots.plot(outx[timeidx], label = false)
Plots.plot(outy[timeidx], label = false)

# relationship
# Plots.scatter(outx[timeidx], outy[timeidx], label = false, color = :black)

# check x distribution
xbar
StatsBase.mean(outx)
xbar * (1 - xbar) / (K * m + 1) 
StatsBase.var(outx)
