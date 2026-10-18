module helper

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
μy(x, y, K, m, xbar, R0, η, π) = @. p11(x, y, K, m, xbar, R0, η, π) * Δy11(K) +
    p10(x, y, K, m, xbar, R0, η, π) * Δy10(K) +
    p01(x, y, K, m, xbar, R0, η, π) * Δy01(K) +
    p0_1(x, y, K, m, xbar, R0, η, π) * Δy0_1(K) +
    p_10(x, y, K, m, xbar, R0, η, π) * Δy_10(K) +
    p_1_1(x, y, K, m, xbar, R0, η, π) * Δy_1_1(K)
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
    (u[1:2:2S], u[2:2:2S])
end

function load_transf(u, S)
    v = 1 ./ (1 .+ exp.(-u))
    (v[1:2:2S], v[2:2:2S])
end

function muop_non!(du, S, K, m, xbar, R0, η, π, x, z, sxx, szz, sxz)
    diag_ssT_x = sxx + sxz
    diag_ssT_z = sxz + szz

    du[1:2:2S] .= @. (1 / (x * (1 - x))) * (μx(x, x * z, K, m, xbar, R0, η, π) + ((2x - 1) / (2x * (1 - x))) * diag_ssT_x)
    du[2:2:2S] .= @. (1 / (z * (1 - z))) * (μz_non(x, x * z, K, m, xbar, R0, η, π) + ((2z - 1) / (2z * (1 - z))) * diag_ssT_z)
end

function f_non!(du, u, p, t)
    S, K, m, xbar, R0, η, π = values(p)

    x, z = load_transf(u, S)

    sxx = σxx(x, x .* z, K, m, xbar, R0, η, π) .^2
    szz = σzz(x, x .* z, K, m, xbar, R0, η, π) .^2
    sxz = σxz(x, x .* z, K, m, xbar, R0, η, π) .^2

    muop_non!(du, S, K, m, xbar, R0, η, π, x, z, sxx, szz, sxz)
end

function muop_stiff!(du, S, K, m, xbar, R0, η, π, x, z, szz, sxz)
    diag_ssT_z = sxz + szz

    du[1:2:2S] .= 0
    du[2:2:2S] .= @. (1 / (z * (1 - z))) * (μz_stiff(x, x * z, K, m, xbar, R0, η, π) + ((2z - 1) / (2z * (1 - z))) * diag_ssT_z)
end

function f_stiff!(du, u, p, t)
    S, K, m, xbar, R0, η, π = values(p)

    x, z = load_transf(u, S)

    szz = σzz(x, x .* z, K, m, xbar, R0, η, π) .^2
    sxz = σxz(x, x .* z, K, m, xbar, R0, η, π) .^2

    muop_stiff!(du, S, K, m, xbar, R0, η, π, x, z, szz, sxz)
end

function f!(du, u, p, t)
    f_non!(du, u, p, t) .+ f_stiff!(du, u, p, t)
end

function muop_xy!(du, S, K, m, xbar, R0, η, π, x, y)
    du[1:2:2S] .= @. μx(x, y, K, m, xbar, R0, η, π)
    du[2:2:2S] .= @. μy(x, y, K, m, xbar, R0, η, π)
end

function f_xy!(du, u, p, t)
    S, K, m, xbar, R0, η, π = values(p)

    x, y = load(u, S)

    muop_xy!(du, S, K, m, xbar, R0, η, π, x, y)
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

    x, z = load_transf(u, S)

    sxx = σxx(x, x .* z, K, m, xbar, R0, η, π)
    szz = σzz(x, x .* z, K, m, xbar, R0, η, π)
    sxz = σxz(x, x .* z, K, m, xbar, R0, η, π)

    sigmaloop!(du, S, x, z, sxx, szz, sxz)
end

function sigmaloop_xy!(du, S, sxx, syy, sxy)
    for i in 1:S
        idx = 2(i - 1)
        du[idx .+ (1:2), idx + 1] .= @. [sxx[i]; sxy[i]]
        du[idx .+ (1:2), idx + 2] .= @. [sxy[i]; syy[i]]
    end
end

function g_xy!(du, u, p, t)
    S, K, m, xbar, R0, η, π = values(p)

    x, y = load(u, S)

    sxx = σxx(x, y, K, m, xbar, R0, η, π)
    syy = σyy(x, y, K, m, xbar, R0, η, π)
    sxy = σxy(x, y, K, m, xbar, R0, η, π)

    sigmaloop_xy!(du, S, sxx, syy, sxy)
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
    xbar = exp.(10(1:S) / S)
    xbar ./ 10sum(xbar)
end

function makeu0_trans(x0, y0)
    z0 = y0 ./ x0
    v0 = [x0 z0]'[:]
    log.(v0 ./ (1 .- v0))
end

function makeu0(x0, y0)
    [x0 y0]'[:]
end

function getsol_trans(sol, S)
    outu = reduce(hcat, [sol.u[i] for i in eachindex(sol.t)])'
    outv = 1 ./ (1 .+ exp.(-outu))
    (outv[:, 1:2:2S], outv[:, 2:2:2S])
end

function getsol(sol, S)
    outu = reduce(hcat, [sol.u[i] for i in eachindex(sol.t)])'
    (outu[:, 1:2:2S], outu[:, 2:2:2S])
end

function gettime(sol)
    T = size(sol.t)[1]
    stepper = div(T, 1000) == 0 ? 1 : div(T, 1000)
    div(T, 2):stepper:T
end

function testf_untrans(u, p, t)
    K, m, xbar = values(p)
    m * (xbar - u)
end

function testf(u, p, t)
    K, m, xbar = values(p)
    v = 1 / (1 + exp(-u))
    (1 / (v * (1 - v))) * (m * (xbar - v) + ((2v - 1) / (2v * (1 - v))) * sqrt(2v * (1 - v) / K)^2)
end

function testg_untrans(u, p, t)
    K, m, xbar = values(p)
    sqrt(2u * (1 - u) / K)
end

function testg(u, p, t)
    K, m, xbar = values(p)
    v = 1 / (1 + exp(-u))
    (1 / (v * (1 - v))) * sqrt(2v * (1 - v) / K)
end

end