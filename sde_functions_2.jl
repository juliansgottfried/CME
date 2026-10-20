module helper

μx(x, m, xbar) = @. m * (xbar - x)
μy(x, K, y, m, xbar, R0, η, π) = R0 .* sum(y) .* (x .- y) .+ m .* (1 - η) .* (xbar .* π .- y) .- η .* y

σxx(x, K) = @. 2x * (1 - x) / K
σyy(x, y, K, R0, η) = R0 .* sum(y) .* (x .- y) ./ K .+ (2 .* y .* (1 .- η ./ 2 .- y .* (1 .- η))) ./ K
σxy(x, y, K, η) = @. (2 - η) * y * (1 - x) / K

function load(u, S)
    @views (u[1:2:2S], u[2:2:2S])
end

function muop!(du, S, K, m, xbar, R0, η, π, x, y)
    du[1:2:2S] .= μx(x, m, xbar)
    du[2:2:2S] .= μy(x, K, y, m, xbar, R0, η, π)
end

function f!(du, u, p, t)
    S, K, m, xbar, R0, η, π = values(p)
    x, y = load(u, S)
    x .= max.(x, 0.0)
    y .= max.(y, 0.0)
    muop!(du, S, K, m, xbar, R0, η, π, x, y)
end

function sigmaloop!(du, S, s11, s22, s12)
    for i in 1:S
        idx = 2(i - 1)
        du[idx .+ (1:2), idx + 1] .= [s11[i]; s12[i]]
        du[idx .+ (1:2), idx + 2] .= [s12[i]; s22[i]]
    end    
end

function g!(du, u, p, t)
    S, K, R0, η, = values(p)
    x, y = load(u, S)
    x .= max.(x, 0.0)
    y .= max.(y, 0.0)
    sigmaloop!(du, S, sqrt.(max.(σxx(x, K), 0)), 
        sqrt.(max.(σyy(x, y, K, R0, η), 0)), 
        sqrt.(max.(σxy(x, y, K, η), 0)))
end

function setA(S)
    A = zeros(Float64, 2S, 2S)
    for i in 1:S
        idx = 2(i - 1)
        A[idx .+ (1:2), idx .+ (1:2)] .= 1
    end
    A
end

function makeu0(x0, y0)
    [x0 y0]'[:]
end

function getsol(sol, S)
    outu = max.(reduce(hcat, [sol.u[i] for i in eachindex(sol.t)])', 0)
    (outu[:, 1:2:2S], outu[:, 2:2:2S])
end

function gettime(sol)
    T = size(sol.t)[1]
    stepper = div(T, 1000) == 0 ? 1 : div(T, 1000)
    min(div(T, 2), 1000):stepper:T
end

function fakexbar(S)
    xbar = exp.((1:S) / S)
    xbar ./ sum(xbar)
end

function testf(u, p, t)
    K, m, xbar = values(p)
    m * (xbar - u)
end

function testg(u, p, t)
    K, m, xbar = values(p)
    sqrt(max(2u * (1 - u) / K, tol))
end

end