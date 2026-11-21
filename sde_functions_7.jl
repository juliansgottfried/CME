module helper

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
    end    
end

function g!(du, u, p, t)
    S, __ = values(p)
    x, y = load(u, S)
    gloop!(du, S, σxx(x), σyy(y), σxy(x, y))
    du .= real.(sqrt(du))
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
    (outu[:, 1:2:2S], outu[:, 2:2:2S])
end

function gettime(sol, J)
    T = size(sol.t)[1]
    stepper = div(T, J) == 0 ? 1 : div(T, J)
    min(div(T, 2), J):stepper:T
end

end