module helper

p11(x, z, K, m, xbar, R0, η, π) = @. K * (1 - x) * (1 - η) * (m * xbar * π + (1 - m) * x * z)
p10(x, z, K, m, xbar, R0, η, π) = @. K * (1 - x) * (m * xbar * (1 - π * (1 - η)) + (1 - m) * (x - x * z * (1 - η)))
p01(x, z, K, m, xbar, R0, η, π) = @. K * (x - x * z) * (1 - η) * (m * xbar * π + (1 - m) * x * z) + R0 * K * x * z * (x - x * z)
p0_1(x, z, K, m, xbar, R0, η, π) = @. K * x * z * (m * xbar * (1 - π * (1 - η)) + (1 - m) * (x - x * z * (1 - η)))
p_10(x, z, K, m, xbar, R0, η, π) = @. K * (x - x * z) * (m * (1 - xbar) + (1 - m) * (1 - x))
p_1_1(x, z, K, m, xbar, R0, η, π) = @. K * x * z * (m * (1 - xbar) + (1 - m) * (1 - x))

Δx11(K) = 1 / K
Δx10(K) = 1 / K
Δx01(K) = 0
Δx0_1(K) = 0
Δx_10(K) = -1 / K
Δx_1_1(K) = -1 / K

Δz11(x, z, K) = @. (x - x * z) / (x * (K * x + 1))
Δz10(x, z, K) = @. - (x * z) / (x * (K * x + 1))
Δz01(x, z, K) = @. 1 / (K * x)
Δz0_1(x, z, K) = @. -1 / (K * x)
Δz_10(x, z, K) = @. (x * z) / (x * (K * x - 1))
Δz_1_1(x, z, K) = @. (x * z - x) / (x * (K * x - 1))

μx(x, z, K, m, xbar, R0, η, π) = @. p11(x, z, K, m, xbar, R0, η, π) * Δx11(K) +
    p10(x, z, K, m, xbar, R0, η, π) * Δx10(K) +
    p01(x, z, K, m, xbar, R0, η, π) * Δx01(K) +
    p0_1(x, z, K, m, xbar, R0, η, π) * Δx0_1(K) +
    p_10(x, z, K, m, xbar, R0, η, π) * Δx_10(K) +
    p_1_1(x, z, K, m, xbar, R0, η, π) * Δx_1_1(K)
μz(x, z, K, m, xbar, R0, η, π) = @. p11(x, z, K, m, xbar, R0, η, π) * Δz11(x, z, K) +
    p10(x, z, K, m, xbar, R0, η, π) * Δz10(x, z, K) +
    p01(x, z, K, m, xbar, R0, η, π) * Δz01(x, z, K) +
    p0_1(x, z, K, m, xbar, R0, η, π) * Δz0_1(x, z, K) +
    p_10(x, z, K, m, xbar, R0, η, π) * Δz_10(x, z, K) +
    p_1_1(x, z, K, m, xbar, R0, η, π) * Δz_1_1(x, z, K)

σxx(x, z, K, m, xbar, R0, η, π) = @. p11(x, z, K, m, xbar, R0, η, π) * Δx11(K)^2 +
    p10(x, z, K, m, xbar, R0, η, π) * Δx10(K)^2 +
    p01(x, z, K, m, xbar, R0, η, π) * Δx01(K)^2 +
    p0_1(x, z, K, m, xbar, R0, η, π) * Δx0_1(K)^2 +
    p_10(x, z, K, m, xbar, R0, η, π) * Δx_10(K)^2 +
    p_1_1(x, z, K, m, xbar, R0, η, π) * Δx_1_1(K)^2
σzz(x, z, K, m, xbar, R0, η, π) = @. p11(x, z, K, m, xbar, R0, η, π) * Δz11(x, z, K)^2 +
    p10(x, z, K, m, xbar, R0, η, π) * Δz10(x, z, K)^2 +
    p01(x, z, K, m, xbar, R0, η, π) * Δz01(x, z, K)^2 +
    p0_1(x, z, K, m, xbar, R0, η, π) * Δz0_1(x, z, K)^2 +
    p_10(x, z, K, m, xbar, R0, η, π) * Δz_10(x, z, K)^2 +
    p_1_1(x, z, K, m, xbar, R0, η, π) * Δz_1_1(x, z, K)^2
σxz(x, z, K, m, xbar, R0, η, π) = @. p11(x, z, K, m, xbar, R0, η, π) * Δx11(K) * Δz11(x, z, K) +
    p10(x, z, K, m, xbar, R0, η, π) * Δx10(K) * Δz10(x, z, K) +
    p01(x, z, K, m, xbar, R0, η, π) * Δx01(K) * Δz01(x, z, K) +
    p0_1(x, z, K, m, xbar, R0, η, π) * Δx0_1(K) * Δz0_1(x, z, K) +
    p_10(x, z, K, m, xbar, R0, η, π) * Δx_10(K) * Δz_10(x, z, K) +
    p_1_1(x, z, K, m, xbar, R0, η, π) * Δx_1_1(K) * Δz_1_1(x, z, K)

function loadxz(u, S)
    v = 1 ./ (1 .+ exp.(-u))
    (v[1:2:2S], v[2:2:2S])
end

function muop!(du, S,K, m, xbar, R0, η, π, x, z, sxx, szz, sxz)
    diag_ssT_x = sxx + sxz
    diag_ssT_z = sxz + szz

    du[1:2:2S] .= @. (1 / (x * (1 - x))) * (μx(x, z, K, m, xbar, R0, η, π) + ((2x - 1) / (2x * (1 - x))) * diag_ssT_x)
    du[2:2:2S] .= @. (1 / (z * (1 - z))) * (μz(x, z, K, m, xbar, R0, η, π) + ((2z - 1) / (2z * (1 - z))) * diag_ssT_z)
end

function f!(du, u, p, t)
    S, K, m, xbar, R0, η, π = values(p)

    x, z = loadxz(u, S)

    sxx = σxx(x, z, K, m, xbar, R0, η, π) .^2
    szz = σzz(x, z, K, m, xbar, R0, η, π) .^2
    sxz = σxz(x, z, K, m, xbar, R0, η, π) .^2

    muop!(du, S,K, m, xbar, R0, η, π, x, z, sxx, szz, sxz)
end

function sigmaloop!(du, S, x, z, sxx, szz, sxz)
    for i in 1:S
        idx = 2(i - 1)
        du[idx .+ (1:2), idx + 1] .= @. [sxx[i] / (x[i] * (1 - x[i])); sxz[i] / (z[i] * (1 - z[i]))]
        du[idx .+ (1:2), idx + 2] .= @. [sxz[i] / (x[i] * (1 - x[i])); szz[i] / (z[i] * (1 - z[i]))]
    end
end

function g!(du, u, p, t)
    S, K, m, xbar, R0, η, π = values(p)

    x, z = loadxz(u, S)

    sxx = σxx(x, z, K, m, xbar, R0, η, π)
    szz = σzz(x, z, K, m, xbar, R0, η, π)
    sxz = σxz(x, z, K, m, xbar, R0, η, π)

    sigmaloop!(du, S, x, z, sxx, szz, sxz)
end

function setA(S)
    A = zeros(Float64, 2S, 2S)
    for i in 1:S
        idx = 2(i - 1)
        A[idx .+ (1:2), idx .+ (1:2)] .= 1
    end
    A
end

function fakexbar(S)
    xbar = exp.((1:S) ./ S)
    xbar ./ sum(xbar)
end

function makeu0(x0, y0)
    z0 = y0 ./ x0
    v0 = [x0 z0]'[:]
    log.(v0 ./ (1 .- v0))
end

function getsol(sol, S)
    outu = reduce(hcat, [sol.u[i] for i in eachindex(sol.t)])'
    outv = 1 ./ (1 .+ exp.(-outu))
    (outv[:, 1:2:2S], outv[:, 2:2:2S])
end

end