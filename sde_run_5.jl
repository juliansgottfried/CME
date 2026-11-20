import StochasticDiffEq as SDE
import DiffEqCallbacks as CB
import DelimitedFiles as DF
import StatsBase
import Plots
import Distributions

include("sde_functions_5.jl")

xbar = sort(DF.readdlm("xbar.csv", ','), dims = 1, rev = true)
pT = sum(xbar)
p1 = StatsBase.maximum(xbar)
f = p1 / pT

K = 10^13
κ = 117.4659
m = κ / K

π = 0.0454943

η = m
R0 = 5m

p = m, pT, f, η, π, R0, K

u0 = [pT; p1; pT * π; p1 * π]

prob = SDE.SDEProblem(helper.f!, helper.g!, 
    u0, (0.0, 100), p,
    noise_rate_prototype = ones(Float64, 4, 4))
sol = SDE.solve(prob, SDE.ImplicitEM(), callback = CB.PositiveDomain(), maxiters = 10^8)

outX, outx, outY, outy = helper.getsol(sol)
timeidx = helper.gettime(sol, 1000)

a = R0
b = m + η - R0 * pT
c = -m * pT * π
EY = (-b+sqrt(b^2-4a*c))/2a

Ex = p1
Ex2 = p1 * (κ * p1 + 1) / (κ + 1)
Varx = Ex2 - Ex^2
Ey = p1 * (K * R0 * EY + κ * π) / (K * R0 * EY + κ + K * η)
Exy = (Ey * κ * p1 + p1 ^ 2 * κ * π + Ex2 * K * R0 * EY) / (2κ + K * η + K * R0 * EY)
Ey2 = (Ey * (κ * p1 * π + 1) + K * R0 * EY * Exy) / (κ + K * η + K * R0 * EY + 1)
Vary = Ey2 - Ey^2
Covxy = Exy - Ex*Ey

A = [-2κ; K * R0 * EY; K * R0 * EY; 0;; 0; -2κ - K * η; 0; K * R0 * EY;; 0; 0; -2κ - K * η; K * R0 * EY;; 0; 0; 0; -2κ - 2K * η]
b = -[p1 * (1 - p1); sqrt(p1 * (1 - p1) * Ey * (1 - Ey)); sqrt(p1 * (1 - p1) * Ey * (1 - Ey)); Ey * (1 - Ey)]
P = A \ b
P = [P[1]; P[3];; P[3]; P[4]]

approx = Distributions.MvNormal([Ex; Ey], [Varx; Covxy;; Covxy; Vary])
samples = Distributions.rand(approx, 5000)
xsample = samples[1, :]
ysample = samples[2, :]

Plots.scatter(outx[timeidx], outy[timeidx], alpha = 0.1, xlims = (0, 0.05), ylims = (0, 0.03), label = false, color = :black)
Plots.scatter(xsample, ysample,  alpha = 0.1, xlims = (0, 0.05), ylims = (0, 0.03), label = false, color = :black)

# X and Y series
Plots.plot(outX[timeidx], label = false, color = :black)
Plots.plot(outY[timeidx], label = false, color = :black)

# x and y series
Plots.plot(outx[timeidx], label = false, color = :black)
Plots.plot(outy[timeidx], label = false, color = :black)

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

# X vs. Y and x vs. y
Plots.scatter(outX[timeidx], outY[timeidx], label = false, color = :black)
Plots.scatter(outx[timeidx], outy[timeidx], label = false, color = :black)

# examine xbar
# Plots.histogram(log.(xbar), bins = 50, label = false, linecolor = :white, color = :black)

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



ret = helper.sim1xi(m, pT, 0.000463, η, π, R0, K)
ret[:, 1] = round.(ret[:, 1], digits = 3)
ret = ret[ret[:, 1] .∈ Ref([0; 0.019]), :]
# can do the above for each species with its list
# across a grid of R0 and eta