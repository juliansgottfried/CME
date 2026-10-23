import StochasticDiffEq as SDE
import DiffEqCallbacks as CB
import DelimitedFiles as DF
# import StatsBase
import Plots

include("sde_functions_6.jl")

xbar = sort(DF.readdlm("xbar.csv", ','), dims = 1, rev = true)
pT = sum(xbar)
p1 = 0.01
f = p1 / pT

K = 10^7
κ = 21.4
m = κ / K

η = 0.5
π = 0.05
R0 = 1

p = m, pT, f, η, π, R0, K

u0 = [pT; p1; pT * π; p1 * π]

prob = SDE.SDEProblem(helper.f!, helper.g!, 
    u0, (0.0, 10^8), p, 
    noise_rate_prototype = ones(Float64, 4, 4))
sol = SDE.solve(prob, SDE.ImplicitEM(), callback = CB.PositiveDomain(), maxiters = 10^8)

outX, outx, outY, outy = helper.getsol(sol)
timeidx = helper.gettime(sol, 1000)

# X and x series
Plots.plot(outX[timeidx], label = false, color = :black)
Plots.plot(outx[timeidx], label = false, color = :black)

# Y and y series
Plots.plot(outY[timeidx], label = false, color = :black)
Plots.plot(outy[timeidx], label = false, color = :black)

# X and x hists
Plots.histogram(outX[timeidx], label = false, linecolor = :white, color = :black)
Plots.histogram(outx[timeidx], label = false, linecolor = :white, color = :black)

# Y and y hists
Plots.histogram(outY[timeidx], label = false, linecolor = :white, color = :black)
Plots.histogram(outy[timeidx], label = false, linecolor = :white, color = :black)

# X vs. x and Y vs. y
Plots.scatter(outX[timeidx], outx[timeidx], label = false, color = :black)
Plots.scatter(outY[timeidx], outy[timeidx], label = false, color = :black)

# X vs. y and x vs. y
Plots.scatter(outX[timeidx], outY[timeidx], label = false, color = :black)
Plots.scatter(outx[timeidx], outy[timeidx], label = false, color = :black)
