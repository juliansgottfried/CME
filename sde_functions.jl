module helper

transf(x) = x
inv_transf(x) = x
prime1(x) = 1
prime2(x) = 0

#=
transf(x) = log(1 - x)
inv_transf(x) = 1 - exp(x)
prime1(x) = -1 / (1 - x)
prime2(x) = -1 / (1 - x)^2
=#

#=
transf(x) = log(x / (1 - x))
inv_transf(x) = 1 / (1 + exp(-x))
prime1(x) = 1 / (x * (1 - x))
prime2(x) = (2x - 1) / (2x^2 * (1 - x)^2)
=#

p11(x, y, K, m, xbar, R0, η, π) = @. K * (1 - x) * (1 - η) * (m * xbar * π + (1 - m) * y)
p10(x, y, K, m, xbar, R0, η, π) = @. K * (1 - x) * (m * xbar * (1 - π * (1 - η)) + (1 - m) * (x - y * (1 - η)))
p01non(x, y, K, m, xbar, R0, η, π) = @. K * (x - y) * (1 - η) * (m * xbar * π + (1 - m) * y)
p01stiff(x, y, K, m, xbar, R0, η, π) = @. K * R0 * sum(y) * (x - y)
p01(x, y, K, m, xbar, R0, η, π) = @. p01non(x, y, K, m, xbar, R0, η, π) + p01stiff(x, y, K, m, xbar, R0, η, π)
p0_1(x, y, K, m, xbar, R0, η, π) = @. K * y * (m * xbar * (1 - π * (1 - η)) + (1 - m) * (x - y * (1 - η)))
p_10(x, y, K, m, xbar, R0, η, π) = @. K * (x - y) * (m * (1 - xbar) + (1 - m) * (1 - x))
p_1_1(x, y, K, m, xbar, R0, η, π) = @. K * y * (m * (1 - xbar) + (1 - m) * (1 - x))

Δx11(K) = 1 / K
Δx10(K) = 1 / K
Δx01(K) = 0
Δx0_1(K) = 0
Δx_10(K) = -1 / K
Δx_1_1(K) = -1 / K

Δy11(K) = 1 / K
Δy10(K) = 0
Δy01(K) = 1 / K
Δy0_1(K) = -1 / K
Δy_10(K) = 0
Δy_1_1(K) = -1 / K

Δz11(x, y, K) = @. (x - y) / (x * (K * x + 1))
Δz10(x, y, K) = @. - y / (x * (K * x + 1))
Δz01(x, y, K) = @. 1 / (K * x)
Δz0_1(x, y, K) = @. -1 / (K * x)
Δz_10(x, y, K) = @. y / (x * (K * x - 1))
Δz_1_1(x, y, K) = @. (y - x) / (x * (K * x - 1))

μx(x, y, K, m, xbar, R0, η, π) = @. p11(x, y, K, m, xbar, R0, η, π) * Δx11(K) +
    p10(x, y, K, m, xbar, R0, η, π) * Δx10(K) +
    p01(x, y, K, m, xbar, R0, η, π) * Δx01(K) +
    p0_1(x, y, K, m, xbar, R0, η, π) * Δx0_1(K) +
    p_10(x, y, K, m, xbar, R0, η, π) * Δx_10(K) +
    p_1_1(x, y, K, m, xbar, R0, η, π) * Δx_1_1(K)

μy_non(x, y, K, m, xbar, R0, η, π) = @. p11(x, y, K, m, xbar, R0, η, π) * Δy11(K) +
    p10(x, y, K, m, xbar, R0, η, π) * Δy10(K) +
    p01non(x, y, K, m, xbar, R0, η, π) * Δy01(K) +
    p0_1(x, y, K, m, xbar, R0, η, π) * Δy0_1(K) +
    p_10(x, y, K, m, xbar, R0, η, π) * Δy_10(K) +
    p_1_1(x, y, K, m, xbar, R0, η, π) * Δy_1_1(K)
μy_stiff(x, y, K, m, xbar, R0, η, π) = @. p01stiff(x, y, K, m, xbar, R0, η, π) * Δy01(K)

μz_non(x, y, K, m, xbar, R0, η, π) = @. p11(x, y, K, m, xbar, R0, η, π) * Δz11(x, y, K) +
    p10(x, y, K, m, xbar, R0, η, π) * Δz10(x, y, K) +
    p01non(x, y, K, m, xbar, R0, η, π) * Δz01(x, y, K) +
    p0_1(x, y, K, m, xbar, R0, η, π) * Δz0_1(x, y, K) +
    p_10(x, y, K, m, xbar, R0, η, π) * Δz_10(x, y, K) +
    p_1_1(x, y, K, m, xbar, R0, η, π) * Δz_1_1(x, y, K)
μz_stiff(x, y, K, m, xbar, R0, η, π) = @. p01stiff(x, y, K, m, xbar, R0, η, π) * Δz01(x, y, K)

σxx(x, y, K, m, xbar, R0, η, π) = @. sqrt(abs(p11(x, y, K, m, xbar, R0, η, π) * Δx11(K)^2 +
    p10(x, y, K, m, xbar, R0, η, π) * Δx10(K)^2 +
    p01(x, y, K, m, xbar, R0, η, π) * Δx01(K)^2 +
    p0_1(x, y, K, m, xbar, R0, η, π) * Δx0_1(K)^2 +
    p_10(x, y, K, m, xbar, R0, η, π) * Δx_10(K)^2 +
    p_1_1(x, y, K, m, xbar, R0, η, π) * Δx_1_1(K)^2))
σyy(x, y, K, m, xbar, R0, η, π) = @. sqrt(abs(p11(x, y, K, m, xbar, R0, η, π) * Δy11(K)^2 +
    p10(x, y, K, m, xbar, R0, η, π) * Δy10(K)^2 +
    p01(x, y, K, m, xbar, R0, η, π) * Δy01(K)^2 +
    p0_1(x, y, K, m, xbar, R0, η, π) * Δy0_1(K)^2 +
    p_10(x, y, K, m, xbar, R0, η, π) * Δy_10(K)^2 +
    p_1_1(x, y, K, m, xbar, R0, η, π) * Δy_1_1(K)^2))
σzz(x, y, K, m, xbar, R0, η, π) = @. sqrt(abs(p11(x, y, K, m, xbar, R0, η, π) * Δz11(x, y, K)^2 +
    p10(x, y, K, m, xbar, R0, η, π) * Δz10(x, y, K)^2 +
    p01(x, y, K, m, xbar, R0, η, π) * Δz01(x, y, K)^2 +
    p0_1(x, y, K, m, xbar, R0, η, π) * Δz0_1(x, y, K)^2 +
    p_10(x, y, K, m, xbar, R0, η, π) * Δz_10(x, y, K)^2 +
    p_1_1(x, y, K, m, xbar, R0, η, π) * Δz_1_1(x, y, K)^2))

σxy(x, y, K, m, xbar, R0, η, π) = @. sqrt(abs(p11(x, y, K, m, xbar, R0, η, π) * Δx11(K) * Δy11(K) +
    p10(x, y, K, m, xbar, R0, η, π) * Δx10(K) * Δy10(K) +
    p01(x, y, K, m, xbar, R0, η, π) * Δx01(K) * Δy01(K) +
    p0_1(x, y, K, m, xbar, R0, η, π) * Δx0_1(K) * Δy0_1(K) +
    p_10(x, y, K, m, xbar, R0, η, π) * Δx_10(K) * Δy_10(K) +
    p_1_1(x, y, K, m, xbar, R0, η, π) * Δx_1_1(K) * Δy_1_1(K)))
σxz(x, y, K, m, xbar, R0, η, π) = @. sqrt(abs(p11(x, y, K, m, xbar, R0, η, π) * Δx11(K) * Δz11(x, y, K) +
    p10(x, y, K, m, xbar, R0, η, π) * Δx10(K) * Δz10(x, y, K) +
    p01(x, y, K, m, xbar, R0, η, π) * Δx01(K) * Δz01(x, y, K) +
    p0_1(x, y, K, m, xbar, R0, η, π) * Δx0_1(K) * Δz0_1(x, y, K) +
    p_10(x, y, K, m, xbar, R0, η, π) * Δx_10(K) * Δz_10(x, y, K) +
    p_1_1(x, y, K, m, xbar, R0, η, π) * Δx_1_1(K) * Δz_1_1(x, y, K)))

function load(u, S)
    v = inv_transf.(u)
    (v[1:2:2S], v[2:2:2S])
end

function muop_non!(du, S, K, m, xbar, R0, η, π, isy, x, yorz, s11, s22, s12)
    diag_ssT_1 = s11 + s12
    diag_ssT_2 = s12 + s22

    if isy
        du[1:2:2S] .= @. prime1(x) * μx(x, yorz, K, m, xbar, R0, η, π) + prime2(x) * diag_ssT_1
        du[2:2:2S] .= @. prime1(yorz) * μy_non(x, yorz, K, m, xbar, R0, η, π) + prime2(yorz) * diag_ssT_2
    else
        du[1:2:2S] .= @. prime1(x) * μx(x, x * yorz, K, m, xbar, R0, η, π) + prime2(x) * diag_ssT_1
        du[2:2:2S] .= @. prime1(yorz) * μz_non(x, x * yorz, K, m, xbar, R0, η, π) + prime2(yorz) * diag_ssT_2
    end
end

function f_non!(du, u, p, t)
    S, K, m, xbar, R0, η, π, isy = values(p)

    x, yorz = load(u, S)

    if isy
        s11 = σxx(x, yorz, K, m, xbar, R0, η, π) .^2
        s22 = σyy(x, yorz, K, m, xbar, R0, η, π) .^2
        s12 = σxy(x, yorz, K, m, xbar, R0, η, π) .^2
    else
        s11 = σxx(x, x .* yorz, K, m, xbar, R0, η, π) .^2
        s22 = σzz(x, x .* yorz, K, m, xbar, R0, η, π) .^2
        s12 = σxz(x, x .* yorz, K, m, xbar, R0, η, π) .^2
    end

    muop_non!(du, S, K, m, xbar, R0, η, π, isy, x, yorz, s11, s22, s12)
end

function muop_stiff!(du, S, K, m, xbar, R0, η, π, isy, x, yorz, s22, s12)
    diag_ssT_2 = s12 + s22

    du[1:2:2S] .= 0

    if isy
        du[2:2:2S] .= @. prime1(yorz) * μy_stiff(x, yorz, K, m, xbar, R0, η, π) + prime2(yorz) * diag_ssT_2
    else
        du[2:2:2S] .= @. prime1(yorz) * μz_stiff(x, x * yorz, K, m, xbar, R0, η, π) + prime2(yorz) * diag_ssT_2
    end
end

function f_stiff!(du, u, p, t)
    S, K, m, xbar, R0, η, π, isy = values(p)

    x, yorz = load(u, S)

    if isy
        s22 = σyy(x, yorz, K, m, xbar, R0, η, π) .^2
        s12 = σxy(x, yorz, K, m, xbar, R0, η, π) .^2
    else
        s22 = σzz(x, x .* yorz, K, m, xbar, R0, η, π) .^2
        s12 = σxz(x, x .* yorz, K, m, xbar, R0, η, π) .^2
    end

    muop_stiff!(du, S, K, m, xbar, R0, η, π, isy, x, yorz, s22, s12)
end

function f!(du, u, p, t)
    f_non!(du, u, p, t) .+ f_stiff!(du, u, p, t)
end

function sigmaloop!(du, S, x, yorz, s11, s22, s12)
    for i in 1:S
        idx = 2(i - 1)
        du[idx .+ (1:2), idx + 1] .= @. [prime1(x[i]) * s11[i]; prime1(yorz[i]) * s12[i]]
        du[idx .+ (1:2), idx + 2] .= @. [prime1(x[i]) * s12[i]; prime1(yorz[i]) * s22[i]]
    end    
end

function g!(du, u, p, t)
    S, K, m, xbar, R0, η, π, isy = values(p)

    x, yorz = load(u, S)

    if isy
        s11 = σxx(x, yorz, K, m, xbar, R0, η, π)
        s22 = σyy(x, yorz, K, m, xbar, R0, η, π)
        s12 = σxy(x, yorz, K, m, xbar, R0, η, π)
    else
        s11 = σxx(x, x .* yorz, K, m, xbar, R0, η, π)
        s22 = σzz(x, x .* yorz, K, m, xbar, R0, η, π)
        s12 = σxz(x, x .* yorz, K, m, xbar, R0, η, π)
    end

    sigmaloop!(du, S, x, yorz, s11, s22, s12)
end

function setA(S)
    A = zeros(Float64, 2S, 2S)
    for i in 1:S
        idx = 2(i - 1)
        A[idx .+ (1:2), idx .+ (1:2)] .= 1
    end
    A
end

function makeu0(x0, y0, isy)
    if isy
        v0 = [x0 y0]'[:]
    else
        z0 = y0 ./ x0
        v0 = [x0 z0]'[:]
    end
    transf.(v0)
end

function getsol(sol, S)
    outu = reduce(hcat, [sol.u[i] for i in eachindex(sol.t)])'
    outv = inv_transf.(outu)
    (outv[:, 1:2:2S], outv[:, 2:2:2S])
end

function gettime(sol)
    T = size(sol.t)[1]
    stepper = div(T, 1000) == 0 ? 1 : div(T, 1000)
    div(T, 2):stepper:T
end

function fakexbar(S)
    xbar = exp.((1:S) / S)
    xbar ./ sum(xbar)
end

end