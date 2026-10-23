module helper

transf(x) = x
inv_transf(x) = x
prime1(x) = oneunit(Float64)
prime2(x) = zero(Float64)

#=
transf(x) = log(x)
inv_transf(x) = exp(x)
prime1(x) = 1 / x
prime2(x) = -1 / x^2
=#

μx(x, y, K, m, xbar, R0, η, π) = @. m * (xbar - x)
μy(x, y, K, m, xbar, R0, η, π) = R0 .* sum(y) .* (x .- y) .+ m .* (1 - η) .* (xbar .* π .- y) .- η .* y

σxx(x, y, K, m, xbar, R0, η, π) = @. 2x * (1 - x) / K
σyy(x, y, K, m, xbar, R0, η, π) = R0 .* sum(y) .* (x .- y) ./ K .+ (2 .* y .* (1 .- η ./ 2 .- y .* (1 .- η))) ./ K
σxy(x, y, K, m, xbar, R0, η, π) = @. (2 - η) * y * (1 - x) / K

#=
σxx(x, y, K, m, xbar, R0, η, π) = @. (2x * (1 - x) + m * (1 - 2x) * (xbar - x)) / K
σyy(x, y, K, m, xbar, R0, η, π) = R0 .* sum(y) .* (x .- y) ./ K .+ m .* xbar .* π .* (1 .- η) .* (1 .- 2y) ./ K .- y .* (m .* (1 .- η) .* (1 .- 2y) .- 2y .* η .+ η .+ 2y .- 2) ./ K
σxy(x, y, K, m, xbar, R0, η, π) = @. ((1 - η) * (1 - x) * (m * xbar * π - m * y + y) - y * (m * (xbar - x) + x - 1)) / K
=#

function load(u, S)
    (inv_transf.(u[1:2:2S]), inv_transf.(u[2:2:2S]))
end

function muop!(du, S, x, y, K, m, xbar, R0, η, π, sum1, sum2)
    du[1:2:2S] .= prime1.(x) .* μx(x, y, K, m, xbar, R0, η, π) .+ prime2.(x) .* sum1 ./ 2
    du[2:2:2S] .= prime1.(y) .* μy(x, y, K, m, xbar, R0, η, π) .+ prime2.(y) .* sum2 ./ 2
end

function f!(du, u, p, t)
    S, K, m, xbar, R0, η, π = values(p)
    x, y = load(u, S)
    
    #=
    x[x .< 0] .= 0
    x[x .> 1] .= 1

    y[y .< 0] .= 0
    y[y .> x] .= x[y .> x]
    =#

    # y[y .< exp(-200)] .= exp(-200)

    #=
    s11 = max.(σxx(x, y, K, m, xbar, R0, η, π), 0.0)
    s22 = max.(σyy(x, y, K, m, xbar, R0, η, π), 0.0)
    s12 = max.(σxy(x, y, K, m, xbar, R0, η, π), 0.0)
    =#

    s11 = σxx(x, y, K, m, xbar, R0, η, π)
    s22 = σyy(x, y, K, m, xbar, R0, η, π)
    s12 = σxy(x, y, K, m, xbar, R0, η, π)

    muop!(du, S, x, y, K, m, xbar, R0, η, π, s11 + s12, s12 + s22)
end

function sigmaloop!(du, S, x, y, s11, s22, s12)
    for i in 1:S
        idx = 2(i - 1)
        du[idx .+ (1:2), idx + 1] .= [prime1(x[i]) * s11[i]; prime1(y[i]) * s12[i]]
        du[idx .+ (1:2), idx + 2] .= [prime1(x[i]) * s12[i]; prime1(y[i]) * s22[i]]
    end    
end

function g!(du, u, p, t)
    S, K, m, xbar, R0, η, π = values(p)
    x, y = load(u, S)
    
    #=
    x[x .< 0] .= 0
    x[x .> 1] .= 1

    y[y .< 0] .= 0
    y[y .> x] .= x[y .> x]
    =#

    # y[y .< exp(-200)] .= exp(-200)

    #=
    sigmaloop!(du, S, x, y, sqrt.(σxx(x, y, K, m, xbar, R0, η, π)), 
        sqrt.(σyy(x, y, K, m, xbar, R0, η, π)), 
        sqrt.(σxy(x, y, K, m, xbar, R0, η, π)))
    =#
    
    sigmaloop!(du, S, x, y, 
        sqrt.(abs.(σxx(x, y, K, m, xbar, R0, η, π))), 
        sqrt.(abs.(σyy(x, y, K, m, xbar, R0, η, π))), 
        sqrt.(abs.(σxy(x, y, K, m, xbar, R0, η, π))))

    #=
    sigmaloop!(du, S, x, y, 
        sqrt.(max.(σxx(x, y, K, m, xbar, R0, η, π), 0.0)), 
        sqrt.(max.(σyy(x, y, K, m, xbar, R0, η, π), 0.0)), 
        sqrt.(max.(σxy(x, y, K, m, xbar, R0, η, π), 0.0)))
    =#
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
    transf.([x0 y0]'[:])
end

function getsol(sol, S)
    # outu = inv_transf.(reduce(hcat, [sol.u[i] for i in eachindex(sol.t)])')
    outu = (reduce(hcat, [sol.u[i] for i in eachindex(sol.t)])')
    (outu[:, 1:2:2S], outu[:, 2:2:2S])
end

function gettime(sol, J)
    T = size(sol.t)[1]
    stepper = div(T, J) == 0 ? 1 : div(T, J)
    min(div(T, 2), J):stepper:T
end

end