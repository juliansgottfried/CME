module helper

μx1(x1, x2, y1, y2, m, p1, p2, η, π, R0, K) = m * (p1 - x1)
μx2(x1, x2, y1, y2, m, p1, p2, η, π, R0, K) = m * (p2 - x2)

μy1(x1, x2, y1, y2, m, p1, p2, η, π, R0, K) = m * (1 - η) * (p1 * π - y1) - η * y1 + R0 * (y1 + y2) * (x1 - y1)
μy2(x1, x2, y1, y2, m, p1, p2, η, π, R0, K) = m * (1 - η) * (p2 * π - y2) - η * y2 + R0 * (y1 + y2) * (x2 - y2)

sigx1x1(x1, x2, y1, y2, m, p1, p2, η, π, R0, K) = (2 * (1 - x1) * x1 + m * (1 - 2x1) * (p1 - x1)) / K
sigx2x2(x1, x2, y1, y2, m, p1, p2, η, π, R0, K) = (2 * (1 - x2) * x2 + m * (1 - 2x2) * (p2 - x2)) / K

sigy1y1(x1, x2, y1, y2, m, p1, p2, η, π, R0, K) = (m * (1 - η) * (1 - 2y1) * (p1 * π - y1) + 2y1 * (1 - 0.5η - y1 * (1 - η)) + R0 * (y1 + y2) * (x1 - y1)) / K
sigy2y2(x1, x2, y1, y2, m, p1, p2, η, π, R0, K) = (m * (1 - η) * (1 - 2y2) * (p2 * π - y2) + 2y2 * (1 - 0.5η - y2 * (1 - η)) + R0 * (y1 + y2) * (x2 - y2)) / K

sigx1x2(x1, x2, y1, y2, m, p1, p2, η, π, R0, K) = -(m * (p1 * x2 + p2 * x1) + 2x1 * x2 * (1 - m)) / K
sigy1y2(x1, x2, y1, y2, m, p1, p2, η, π, R0, K) = -(1 - η) * (m * π * (p1 * y2 + p2 * y1) + 2y1 * y2 * (1 - m)) / K

sigx1y1(x1, x2, y1, y2, m, p1, p2, η, π, R0, K) = (m * p1 * π * (1 - η) * (1 - x1) - m * y1 * (1 - 2x1 + p1 - η * (1 - x1)) + (2 - η) * (1 - x1) * y1) / K
sigx2y2(x1, x2, y1, y2, m, p1, p2, η, π, R0, K) = (m * p2 * π * (1 - η) * (1 - x2) - m * y2 * (1 - 2x2 + p2 - η * (1 - x2)) + (2 - η) * (1 - x2) * y2) / K

sigx1y2(x1, x2, y1, y2, m, p1, p2, η, π, R0, K) = (m * y2 * ((2 - η) * x1 - p1) - m * p2 * π * (1 - η) * x1 - (2 - η) * x1 * y2) / K
sigx2y1(x1, x2, y1, y2, m, p1, p2, η, π, R0, K) = (m * y1 * ((2 - η) * x2 - p2) - m * p1 * π * (1 - η) * x2 - (2 - η) * x2 * y1) / K

function f!(du, u, p, t)
    m, p1, p2, η, π, R0, K = values(p)

    x1 = u[1]
    x2 = u[2]
    y1 = u[3]
    y2 = u[4]

    du[1] = μx1(x1, x2, y1, y2, m, p1, p2, η, π, R0, K)
    du[2] = μx2(x1, x2, y1, y2, m, p1, p2, η, π, R0, K)
    du[3] = μy1(x1, x2, y1, y2, m, p1, p2, η, π, R0, K)
    du[4] = μy2(x1, x2, y1, y2, m, p1, p2, η, π, R0, K)
end

function g!(du, u, p, t)
    m, p1, p2, η, π, R0, K = values(p)
    
    x1 = u[1]
    x2 = u[2]
    y1 = u[3]
    y2 = u[4]

    sx1x1 = sqrt(abs(sigx1x1(x1, x2, y1, y2, m, p1, p2, η, π, R0, K)))
    sy1y1 = sqrt(abs(sigy1y1(x1, x2, y1, y2, m, p1, p2, η, π, R0, K)))
    
    sx2x2 = sqrt(abs(sigx2x2(x1, x2, y1, y2, m, p1, p2, η, π, R0, K)))
    sy2y2 = sqrt(abs(sigy2y2(x1, x2, y1, y2, m, p1, p2, η, π, R0, K)))
   
    sx1x2 = sqrt(abs(sigx1x2(x1, x2, y1, y2, m, p1, p2, η, π, R0, K)))
    sy1y2 = sqrt(abs(sigy1y2(x1, x2, y1, y2, m, p1, p2, η, π, R0, K)))

    sx1y1 = sqrt(abs(sigx1y1(x1, x2, y1, y2, m, p1, p2, η, π, R0, K)))
    sx2y2 = sqrt(abs(sigx2y2(x1, x2, y1, y2, m, p1, p2, η, π, R0, K)))

    sx1y2 = sqrt(abs(sigx1y2(x1, x2, y1, y2, m, p1, p2, η, π, R0, K)))
    sx2y1 = sqrt(abs(sigx2y1(x1, x2, y1, y2, m, p1, p2, η, π, R0, K)))

    du[1, 1] = sx1x1
    du[1, 2] = sx1x2
    du[1, 3] = sx1y1
    du[1, 4] = sx1y2

    du[2, 1] = sx1x2
    du[2, 2] = sx2x2
    du[2, 3] = sx2y1
    du[2, 4] = sx2y2

    du[3, 1] = sx1y1
    du[3, 2] = sx2y1
    du[3, 3] = sy1y1
    du[3, 4] = sy1y2

    du[4, 1] = sx1y2
    du[4, 2] = sx2y2
    du[4, 3] = sy1y2
    du[4, 4] = sy2y2
end

function getsol(sol)
    outu = (reduce(hcat, [sol.u[i] for i in eachindex(sol.t)])')
    (outu[:, 1], outu[:, 2], outu[:, 3], outu[:, 4])
end

function gettime(sol, J)
    T = size(sol.t)[1]
    stepper = div(T, J) == 0 ? 1 : div(T, J)
    min(div(T, 2), J):stepper:T
    # 1:stepper:T
end

end