module helper

import StochasticDiffEq as SDE
import DiffEqCallbacks as CB

μX(X, Xbar, κ) = κ * (Xbar - X)
μx(x, xbar, κ) = κ * (xbar - x)

μY(X, Y, Xbar, κ, π, R0, η) = κ * (Xbar * π - Y) - η * Y + R0 * Y * (X - Y)
μy(x, Y, y, xbar, κ, π, R0, η) = κ * (xbar * π - y) - η * y + R0 * Y * (x - y)

sigXX(X) = 2X * (1 - X)
sigxx(x) = 2x * (1 - x)

sigYY(Y) = 2Y * (1 - Y)
sigyy(y) = 2y * (1 - y)

sigXx(X, x) = 2x * (1 - X)
sigYy(Y, y) = 2y * (1 - Y)

sigXY(X, Y, Xbar) = -Y * (Xbar - X) + 2(1 - X) * Y
sigxy(x, y, xbar) = -y * (xbar - x) + 2(1 - x) * y

sigXy(X, y, Xbar) = -y * (Xbar - X) + 2(1 - X) * y
sigxY(x, Y, y, xbar) = -Y * (xbar - x) + 2(y - x * Y)

function f!(du, u, p, t)
    Xbar, xbar, κ, π, R0, η = values(p)

    X = u[1]
    x = u[2]
    Y = u[3]
    y = u[4]

    du[1] = μX(X, Xbar, κ)
    du[2] = μx(x, xbar, κ)
    du[3] = μY(X, Y, Xbar, κ, π, R0, η)
    du[4] = μy(x, Y, y, xbar, κ, π, R0, η)
end

function g!(du, u, p, t)
    Xbar, xbar, κ, π, R0, η = values(p)

    X = u[1]
    x = u[2]
    Y = u[3]
    y = u[4]

    sXX = sigXX(X)
    sYY = sigYY(Y)
    
    sxx = sigxx(x)
    syy = sigyy(y)
   
    sXx = sigXx(X, x)
    sYy = sigYy(Y, y)

    sXY = sigXY(X, Y, Xbar)
    sxy = sigxy(x, y, xbar)

    sXy = sigXy(X, y, Xbar)
    sxY = sigxY(x, Y, y, xbar)

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

    du .= real.(sqrt(du))
end

function getsol(sol)
    outu = reduce(hcat, [sol.u[i] for i in eachindex(sol.t)])'

    outX = outu[:, 1]
    outx = outu[:, 2]
    outY = outu[:, 3]
    outy = outu[:, 4]

    #=
    outX[outX .> 1] .= 1
    outX[outX .< 0] .= 0
    outY[outY .> outX] .= outX[outY .> outX]
    outY[outY .< 0] .= 0

    outx[outx .> 1] .= 1
    outx[outx .< 0] .= 0
    outy[outy .> outx] .= outx[outy .> outx]
    outy[outy .< 0] .= 0
    =#

    (outX, outx, outY, outy)
end

function gettime(sol, J)
    T = size(sol.t)[1]
    stepper = div(T, J) == 0 ? 1 : div(T, J)
    min(div(T, 2), J):stepper:T
end

function runner(i, j, Xbar, xbar, κ, π, R0, η, J)
    p = (Xbar, xbar[j], κ, π, R0[i], η)

    u0 = [Xbar; xbar[j]; Xbar * π; xbar[j] * π]

    prob = SDE.SDEProblem(helper.f!, helper.g!, 
        u0, (0.0, 10.0), p, noise_rate_prototype = ones(Float64, 4, 4))
    sol = SDE.solve(prob, SDE.ImplicitEM(), 
        callback = CB.PositiveDomain())

    outX, outx, outY, outy = helper.getsol(sol)
    timeidx = helper.gettime(sol, J)

    idx = sortperm(outY)

    (outX[idx][timeidx], outx[idx][timeidx], outY[idx][timeidx], outy[idx][timeidx])
end

function overbins(i, Xbar, xbar, κ, π, R0, η, J)
    print(1)
    outX, outx, outY, outy = runner(i, 1, Xbar, xbar, κ, π, R0, η, J)
    collector = [outX;; outY;; outx;; outy]
    for j in 2:size(xbar)[1]
        print(j)
        outX, outx, outY, outy = runner(i, j, Xbar, xbar, κ, π, R0, η, J)
        collector = [collector;; outx;; outy]
    end
    collector
end

end