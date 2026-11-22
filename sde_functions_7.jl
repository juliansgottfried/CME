module helper

import StochasticDiffEq as SDE
import DelimitedFiles as DF

μx(x, xbar, κ) = @. κ * (xbar - x)
μy(x, y, Y, xbar, κ, π, R0, η) = @. κ * (xbar * π - y) - η * y + R0 * Y * (x - y)

σxx(x) = @. 2x * (1 - x)
σyy(y) = @. 2y * (1 - y)
σxy(x, y) = @. 2y * (1 - x)

function load(u, S)
    (u[1:2:2S], u[2:2:2S])
end

function f!(du, u, p, t)
    S, xbar, κ, π, R0, η = values(p)
    x, y = load(u, S)
    Y = sum(y)
    du[1:2:2S] .= μx(x, xbar, κ)
    du[2:2:2S] .= μy(x, y, Y, xbar, κ, π, R0, η)
end

function gloop!(du, S, s11, s22, s12)
    for i in 1:S
        idx = 2(i - 1)
        du[idx .+ (1:2), idx + 1] .= [s11[i]; s12[i]]
        du[idx .+ (1:2), idx + 2] .= [s12[i]; s22[i]]
        @views du[idx .+ (1:2), idx .+ (1:2)] .= real.(sqrt(du[idx .+ (1:2), idx .+ (1:2)]))
    end    
end

function g!(du, u, p, t)
    S, __ = values(p)
    x, y = load(u, S)
    gloop!(du, S, σxx(x), σyy(y), σxy(x, y))
end

function setA(S)
    A = zeros(Float64, 2S, 2S)
    for i in 1:S
        idx = 2(i - 1)
        A[idx .+ (1:2), idx .+ (1:2)] .= oneunit(Float64)
    end
    A
end

function makeu0(x0, y0)
    [x0 y0]'[:]
end

function getsol(sol, S)
    outu = (reduce(hcat, [sol.u[i] for i in eachindex(sol.t)])')
    outx = outu[:, 1:2:2S]
    outy = outu[:, 2:2:2S]

    # outx[outx .> 1] .= 1
    # outx[outx .< 0] .= 0
    # outy[outy .> outx] .= outx[outy .> outx]
    # outy[outy .< 0] .= 0

    (outx, outy)
end

function gettime(sol, J)
    T = size(sol.t)[1]
    stepper = div(T, J) == 0 ? 1 : div(T, J)
    min(div(T, 2), J):stepper:T
end

function runner(i, S, xbar, κ, π, R0s, η)
    p = (S, xbar, κ, π, R0s[i], η)

    x0 = xbar
    y0 = zeros(Float64, S)
    u0 = makeu0(x0, y0)

    prob = SDE.SDEProblem(helper.f!, helper.g!,
        u0, (0.0, 1.0), p, noise_rate_prototype = setA(S))
    sol = SDE.solve(prob, SDE.ImplicitEM())

    outx, outy = helper.getsol(sol, S)
    timeidx = helper.gettime(sol, 42)

    DF.writedlm("simulations/x_$i.csv", outx[timeidx, :], ',')
    DF.writedlm("simulations/y_$i.csv", outy[timeidx, :], ',')
end

end