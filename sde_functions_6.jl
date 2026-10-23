module helper

μX(X, x, Y, y, m, pT, f, η, π, R0, K) = 0
μx(X, x, Y, y, m, pT, f, η, π, R0, K) = 0

μY(X, x, Y, y, m, pT, f, η, π, R0, K) = - η * Y + R0 * Y * (X - Y)
μy(X, x, Y, y, m, pT, f, η, π, R0, K) = - η * y + R0 * Y * (x - y)

sigXX(X, x, Y, y, m, pT, f, η, π, R0, K) = 2X * (1 - X) / K
sigxx(X, x, Y, y, m, pT, f, η, π, R0, K) = 2x * (1 - x) / K

sigYY(X, x, Y, y, m, pT, f, η, π, R0, K) = (2Y * (1 - η/2 - Y * (1 - η)) + R0 * Y * (X - Y)) / K
sigyy(X, x, Y, y, m, pT, f, η, π, R0, K) = (2y * (1 - η/2 - y * (1 - η)) + R0 * Y * (x - y)) / K

sigXx(X, x, Y, y, m, pT, f, η, π, R0, K) = (2x * (1 - X)) / K
sigYy(X, x, Y, y, m, pT, f, η, π, R0, K) = (2y * (1 - η) * (1 - Y) + R0 * Y * (x - y)) / K

sigXY(X, x, Y, y, m, pT, f, η, π, R0, K) = (2 - η) * (1 - X) * Y / K
sigxy(X, x, Y, y, m, pT, f, η, π, R0, K) = (2 - η) * (1 - x) * y / K

sigXy(X, x, Y, y, m, pT, f, η, π, R0, K) = (2 - η) * (1 - X) * y / K
sigxY(X, x, Y, y, m, pT, f, η, π, R0, K) = (2 - η) * (y - x * Y) / K

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

    sXX = sqrt(abs(sigXX(X, x, Y, y, m, pT, f, η, π, R0, K)))
    sYY = sqrt(abs(sigYY(X, x, Y, y, m, pT, f, η, π, R0, K)))
    
    sxx = sqrt(abs(sigxx(X, x, Y, y, m, pT, f, η, π, R0, K)))
    syy = sqrt(abs(sigyy(X, x, Y, y, m, pT, f, η, π, R0, K)))
   
    sXx = sqrt(abs(sigXx(X, x, Y, y, m, pT, f, η, π, R0, K)))
    sYy = sqrt(abs(sigYy(X, x, Y, y, m, pT, f, η, π, R0, K)))

    sXY = sqrt(abs(sigXY(X, x, Y, y, m, pT, f, η, π, R0, K)))
    sxy = sqrt(abs(sigxy(X, x, Y, y, m, pT, f, η, π, R0, K)))

    sXy = sqrt(abs(sigXy(X, x, Y, y, m, pT, f, η, π, R0, K)))
    sxY = sqrt(abs(sigxY(X, x, Y, y, m, pT, f, η, π, R0, K)))

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