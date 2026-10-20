import StochasticDiffEq as SDE
# import DiffEqCallbacks as CB
import StatsPlots
import Plots
import ADTypes as AD

include("sde_functions.jl")

S = 10

K = 10^6
m = 10K^-1
xbar = helper.fakexbar(S)

R0 = 0
η = 0.5
π = 0.01

isy = false

p = (S, K, m, xbar, R0, η, π, isy)

x0 = xbar
y0 = ones(Float64, S) * minimum(xbar) / 10
u0 = helper.makeu0(x0, y0, isy)

prob = SDE.SDEProblem(helper.f!, helper.g!,
    u0, (0.0, 10^6), p, noise_rate_prototype = helper.setA(S))
# sol = SDE.solve(prob, SDE.SOSRA())
# isoutofdomain = (u, p, t) -> any(x -> x < 0, u))
# sol = SDE.solve(prob, SDE.SRA1())
sol = SDE.solve(prob, SDE.LambaEM())

# sol = SDE.solve(prob, SDE.SOSRA(),
#     callback = CB.PositiveDomain())


# prob = SDE.SplitSDEProblem(helper.f_stiff!, helper.f_non!, helper.g!,
#     u0, (0.0, 10^4), p, noise_rate_prototype = helper.setA(S))
# sol = SDE.solve(prob, SDE.ImplicitEM(autodiff = AD.AutoFiniteDiff()))

outx, outyorz = helper.getsol(sol, S)
timeidx = helper.gettime(sol)

Plots.plot(((outx[timeidx, :])))
Plots.plot(((outyorz[timeidx, :])))

Plots.plot(sum(eachcol(outx)))
Plots.plot(xbar, StatsBase.mean(eachrow(outx)))
Plots.plot(xbar .* (1 .- xbar) ./ (K * m + 1), StatsBase.var(eachrow(outx)))

Plots.plot([outx[timeidx, S], outyorz[timeidx, S]], label = false)

Plots.plot(sum(eachcol(outx)))

ally = sum(eachcol(outyorz))
scale = @. (R0 * ally + m * (1 - η) * π) / (R0 * ally + m * (1 - η) + η)
yhat = similar(outyorz)
for i in 1:S
    yhat[:, i] = outx[:, i] .* scale
end

Plots.plot([outx[timeidx, S], yhat[timeidx, S]], label = false)

1

# Plots.plot(outx[timeidx, :], outyorz[timeidx, :], color = :black, label = false)
# StatsPlots.density(outyorz[timeidx, :], color = :black, label = false)

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

DelimitedFiles.writedlm("outx_10.csv", outx, ',')
DelimitedFiles.writedlm("outy_10.csv", outyorz, ',')
