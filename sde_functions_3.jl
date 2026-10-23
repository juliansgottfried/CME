module helper

μx(x, y, K, m, xbar, R0, η, π) = m * (xbar - x)
μy(x, y, K, m, xbar, R0, η, π) = R0 * y * (x - y) + m * (1 - η) * (xbar * π - y) - η * y

σxx(x, y, K, m, xbar, R0, η, π) = 2x * (1 - x) / K
σyy(x, y, K, m, xbar, R0, η, π) = R0 * y * (x - y) / K + (2 * y * (1 - η / 2 - y * (1 - η))) / K
σxy(x, y, K, m, xbar, R0, η, π) = (2 - η) * y * (1 - x) / K

function f!(du, u, p, t)
    K, m, xbar, R0, η, π = values(p)    

    x = u[1]
    y = u[2]

    #=
    if x < 0 x = 0 end
    if x > 1 x = 1 end

    if y < 0 y = 0 end
    if y > x y = x end
    =#
    
    du[1] = μx(x, y, K, m, xbar, R0, η, π)
    du[2] = μy(x, y, K, m, xbar, R0, η, π)
end

function g!(du, u, p, t)
    K, m, xbar, R0, η, π = values(p)
    
    x = u[1]
    y = u[2]

    #=
    if x < 0 x = 0 end
    if x > 1 x = 1 end

    if y < 0 y = 0 end
    if y > x y = x end
    =#

    #=
    sxx = sqrt(max(σxx(x, y, K, m, xbar, R0, η, π), 0))
    syy = sqrt(max(σyy(x, y, K, m, xbar, R0, η, π), 0))
    sxy = sqrt(max(σxy(x, y, K, m, xbar, R0, η, π), 0))
    =#

    sxx = sqrt(abs(σxx(x, y, K, m, xbar, R0, η, π)))
    syy = sqrt(abs(σyy(x, y, K, m, xbar, R0, η, π)))
    sxy = sqrt(abs(σxy(x, y, K, m, xbar, R0, η, π)))

    #=
    sxx = sqrt(σxx(x, y, K, m, xbar, R0, η, π))
    syy = sqrt(σyy(x, y, K, m, xbar, R0, η, π))
    sxy = sqrt(σxy(x, y, K, m, xbar, R0, η, π))
    =#
    
    du[1, 1] = sxx
    du[1, 2] = sxy
    du[2, 1] = sxy
    du[2, 2] = syy
end

function getsol(sol)
    outu = (reduce(hcat, [sol.u[i] for i in eachindex(sol.t)])')
    (outu[:, 1], outu[:, 2])
end

function gettime(sol, J)
    T = size(sol.t)[1]
    stepper = div(T, J) == 0 ? 1 : div(T, J)
    min(div(T, 2), J):stepper:T
    # 1:stepper:T
end

end