import StochasticDiffEq as SDE
import DiffEqCallbacks as CB
import StatsBase
import Plots

include("sde_functions_5.jl")

xbar = sort(DF.readdlm("xbar.csv", ','), dims = 1, rev = true)
pT = sum(xbar)
p1 = 0.01
f = p1 / pT

K = 10^7
κ = 134.8684
m = κ / K

R0 = 0.8
η = 0.5
π = 0.05

p = m, pT, f, η, π, R0, K

u0 = [pT; p1; pT * π; p1 * π]

prob = SDE.SDEProblem(helper.f!, helper.g!, 
    u0, (0.0, 10^8), p, 
    noise_rate_prototype = ones(Float64, 4, 4))
sol = SDE.solve(prob, SDE.ImplicitEM(), callback = CB.PositiveDomain(), maxiters = 10^8)

outX, outx, outY, outy = helper.getsol(sol)
timeidx = helper.gettime(sol, 1000)

# X and x
Plots.plot(outX[timeidx], label = false, color = :black)
Plots.plot(outx[timeidx], label = false, color = :black)

# Y and y
Plots.plot(outY[timeidx], label = false, color = :black)
Plots.plot(outy[timeidx], label = false, color = :black)

# X vs. y and x vs. y
Plots.scatter(outX[timeidx], outY[timeidx], label = false, color = :black)
Plots.scatter(outx[timeidx], outy[timeidx], label = false, color = :black)

# X vs. x and Y vs. y
Plots.scatter(outX[timeidx], outx[timeidx], label = false, color = :black)
Plots.scatter(outY[timeidx], outy[timeidx], label = false, color = :black)

#=
# check x1 distribution
pT
StatsBase.mean(outX)
pT * (1 - pT) / (K * m + 1) 
StatsBase.var(outX)

# check x2 distribution
p1
StatsBase.mean(outx)
p1 * (1 - p1) / (K * m + 1) 
StatsBase.var(outx)
=#