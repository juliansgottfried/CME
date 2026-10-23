import StochasticDiffEq as SDE
import DiffEqCallbacks as CB
import StatsBase
import Plots

include("sde_functions_4.jl")

xbar = sort(DF.readdlm("xbar.csv", ','), dims = 1, rev = true)
p1 = maximum(xbar)
p2 = sum(xbar) - p1

K = 10^7
κ = 134.8684
m = κ / K

R0 = 0
η = 0.5
π = 0.05

p = (m, p1, p2, η, π, R0, K)

u0 = [p1; p2; π * p1; π * p2]

prob = SDE.SDEProblem(helper.f!, helper.g!, 
    u0, (0.0, 10^6), p, 
    noise_rate_prototype = ones(Float64, 4, 4))
sol = SDE.solve(prob, SDE.ImplicitEM(), callback = CB.PositiveDomain())

outx1, outx2, outy1, outy2 = helper.getsol(sol)
timeidx = helper.gettime(sol, 1000)

# plot x
Plots.plot(outx1[timeidx], label = false, color = :black)
Plots.plot(outx2[timeidx], label = false, color = :black)

# plot y
Plots.plot(outy1[timeidx], label = false, color = :black)
Plots.plot(outy2[timeidx], label = false, color = :black)

# check x1 distribution
p1
StatsBase.mean(outx1)
p1 * (1 - p1) / (K * m + 1) 
StatsBase.var(outx1)

# check x2 distribution
p2
StatsBase.mean(outx2)
p2 * (1 - p2) / (K * m + 1) 
StatsBase.var(outx2)
