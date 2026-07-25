module sim1

using Distributions

birth(pX, μ, Δt) = rand(Poisson(pX * μ * Δt))
death(X, Δt) = rand(Binomial(X, 1 - exp(-Δt)))

function FOI(β, S, I, Δt)
    if S == 0 && I == 0 return 0 end
    rand(Binomial(S, 1 - exp(-β * I / (S + I) * Δt)))
end

saver!(g, Is, I, inter) = Is[Int(g / inter)] = I

function iterate(β, μ, p, S, I, Δt)
    muSI = FOI(β, S, I, Δt)
    S -= muSI
    I += muSI
    S += birth(1 - p, μ, Δt) - death(S, Δt)
    I += birth(p, μ, Δt) - death(I, Δt)
    (S, I)
end

function loop!(Is, β, μ, p, Δt, G, inter)
    S = 100
    I = 0
    for g in 1:G
        if g % inter == 0 saver!(g, Is, I, inter) end
        S, I = iterate(β, μ, p, S, I, Δt)
    end
end    

function replication(β, μ, p, Δt, G, J, inter)
    nrun = Int(G / inter)
    all = zeros(Int, J * nrun)
    for j in 1:J
        idx = nrun * (j - 1) + 1
        @views loop!(all[idx:(idx + nrun - 1)], 
            β, μ, p, Δt, G, inter)
    end
    all
end

end
