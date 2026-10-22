import StochasticDiffEq as SDE
import StatsBase
import Plots
import DelimitedFiles as DF

include("sde_functions_2.jl")

xbar = sort(DF.readdlm("xbar.csv", ','), dims = 1, rev = true)

S = size(xbar)[1]

K = 10^7
κ = 134.8684
m = κ / K

R0 = 1
η = 0.5
π = 0.01

p = (S, K, m, xbar, R0, η, π)

x0 = xbar
y0 = zeros(Float64, S)
u0 = helper.makeu0(x0, y0)

prob = SDE.SDEProblem(helper.f!, helper.g!,
    u0, (0.0, 10^5), p, noise_rate_prototype = helper.setA(S))
# sol = SDE.solve(prob, SDE.LambaEM(), maxiters = 10^7)
sol = SDE.solve(prob, SDE.SRA1(), maxiters = 10^7)

outx, outy = helper.getsol(sol, S)
timeidx = helper.gettime(sol)

# plot x and y
Plots.plot(outx[timeidx, :], label = false)
Plots.plot(outy[timeidx, :], label = false)

# plot a particular x and y
idx = 50
Plots.plot(outx[timeidx, idx], label = false, color = :black)
Plots.plot(outy[timeidx, idx], label = false, color = :black)
Plots.scatter(outx[timeidx, idx], outy[timeidx, idx], 
    label = false, color = :black)

# check x distribution
Plots.scatter(xbar, 
    [StatsBase.mean(eachrow(outx[timeidx, :]))],
    color = :black, label = false)
Plots.plot!(xbar, xbar, color = :red, label = false)
Plots.scatter(xbar .* (1 .- xbar) ./ (K * m + 1), 
    [StatsBase.var(eachrow(outx[timeidx, :]))],
    color = :black, label = false)
Plots.plot!(xbar .* (1 .- xbar) ./ (K * m + 1), xbar .* (1 .- xbar) ./ (K * m + 1), color = :red, label = false)

DF.writedlm("outx.csv", outx, ',')
DF.writedlm("outy.csv", outy, ',')