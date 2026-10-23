module helpTer

μX(X, x, Y, y, m, pT, f, η, π, R0, K) = m * (pT - X)
μx(X, x, Y, y, m, pT, f, η, π, R0, K) = m * (f * pT - x)

μY(X, x, Y, y, m, pT, f, η, π, R0, K) = m * (1 - η) * (pT * π - Y) - η * Y + R0 * Y * (X - Y)
μy(X, x, Y, y, m, pT, f, η, π, R0, K) = m * (1 - η) * (f * pT * π - y) - η * y + R0 * Y * (x - y)

sigXX(X, x, Y, y, m, pT, f, η, π, R0, K) = (2X * (1 - X) + m * (1 - 2X) * (pT - X)) / K
sigxx(X, x, Y, y, m, pT, f, η, π, R0, K) = (2x * (1 - x) + m * (1 - 2x) * (f * pT - x)) / K

sigYY(X, x, Y, y, m, pT, f, η, π, R0, K) = (2Y * (1 - η/2 - Y * (1 - η)) + m * (1 - η) * (1 - 2Y) * (pT * π - Y) + R0 * Y * (X - Y)) / K
sigyy(X, x, Y, y, m, pT, f, η, π, R0, K) = (2y * (1 - η/2 - y * (1 - η)) + m * (1 - η) * (1 - 2y) * (f * pT * π - y) + R0 * Y * (x - y)) / K

sigXx(X, x, Y, y, m, pT, f, η, π, R0, K) = -(m * pT * (x - f * (1 - X)) - x * ((1 - m) * (1 - 2X) + 1)) / K
sigYy(X, x, Y, y, m, pT, f, η, π, R0, K) = -(m * pT * π * (1 - η) * (y - f * (1 - Y)) - y * (1 - η) * ((1 - m) * (1 - 2Y) + 1) - R0 * Y * (x - y)) / K

sigXY(X, x, Y, y, m, pT, f, η, π, R0, K) = (m * ((1 - η) * (1 - X) * (pT * π - Y) - Y * (pT - X)) + (2 - η) * (1 - X) * Y) / K
sigxy(X, x, Y, y, m, pT, f, η, π, R0, K) = (m * ((1 - η) * (1 - x) * (f * pT * π - y) - y * (f * pT - x)) + (2 - η) * (1 - x) * y) / K

sigXy(X, x, Y, y, m, pT, f, η, π, R0, K) = (m * ((1 - η) * (1 - X) * (f * pT * π - y) - y * (pT - X)) + (2 - η) * (1 - X) * y) / K
sigxY(X, x, Y, y, m, pT, f, η, π, R0, K) = (m * (1 - η) * ((f * pT * π - y - x * (pT * π - Y) - Y * (f * pT - x)) + (2 - η) * (y - x * Y))) / K

function f!(du, u, p, t)
    m, pT, f, η, π, R0, K = values(p)

    X = u[1]
    x = u[2]
    Y = u[3]
    y = u[4]

    du[1] = μX(X, x, Y, y, m, pT, f, η, π, R0, K)
    du[2] = μx(X, x, Y, y, m, pT, f, η, π, R0, K)
    du[3] = μY(X, x, Y, y, m, pT, f, η, π, R0, K)
    du[4] = μy(X, x, Y, y, m, pT, f, η, π, R0, K)
end

function g!(du, u, p, t)
    m, pT, f, η, π, R0, K = values(p)
    
    X = u[1]
    x = u[2]
    Y = u[3]
    y = u[4]

    sXX = sqrt((sigXX(X, x, Y, y, m, pT, f, η, π, R0, K)))
    sYY = sqrt((sigYY(X, x, Y, y, m, pT, f, η, π, R0, K)))
    
    sxx = sqrt((sigxx(X, x, Y, y, m, pT, f, η, π, R0, K)))
    syy = sqrt((sigyy(X, x, Y, y, m, pT, f, η, π, R0, K)))
   
    sXx = sqrt((sigXx(X, x, Y, y, m, pT, f, η, π, R0, K)))
    sYy = sqrt((sigYy(X, x, Y, y, m, pT, f, η, π, R0, K)))

    sXY = sqrt((sigXY(X, x, Y, y, m, pT, f, η, π, R0, K)))
    sxy = sqrt((sigxy(X, x, Y, y, m, pT, f, η, π, R0, K)))

    sXy = sqrt((sigXy(X, x, Y, y, m, pT, f, η, π, R0, K)))
    sxY = sqrt((sigxY(X, x, Y, y, m, pT, f, η, π, R0, K)))

    du[1, 1] = sXX
    du[1, 2] = sXx
    du[1, 3] = sXY
    du[1, 4] = sXy

    du[2, 1] = sXx
    du[2, 2] = sxx
    du[2, 3] = sxY
    du[2, 4] = sxy

    du[3, 1] = sXY
    du[3, 2] = sxY
    du[3, 3] = sYY
    du[3, 4] = sYy

    du[4, 1] = sXy
    du[4, 2] = sxy
    du[4, 3] = sYy
    du[4, 4] = syy
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