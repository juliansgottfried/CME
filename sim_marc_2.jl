module Sim

function iterate(g, S, I, R0, π, μ)
    FOI = R0 * S * I
    λ = μ + FOI + S + I
    g -= 1 / λ * log(rand())
    dart = rand()
    if dart < (1 - π) * μ / λ S += 1
    elseif dart < μ / λ I += 1
    elseif dart < (μ + FOI) / λ
        S -= 1
        I += 1
    elseif dart < (μ + FOI + S) / λ S -= 1
    else I -= 1
    end
    (g, S, I)
end

function loop!(Is, nN, G, inter, t, R0, π, μ)
    g = 0
    S = 25
    I = 0
    t .= 0
    while g < G
        N = S + I
        if N > 0 && N <= nN
            if g ≥ inter * t[N]
                t[N] += 1
                Is[N, t[N]] = I
            end
        end
        g, S, I = iterate(g, S, I, R0, π, μ)
    end
end    

function replication!(all, J, G, inter, nN, len, t, R0, π, μ)
    all .= 0
    for j in 1:J
        idx = len * (j - 1) + 1
        @views loop!(all[:, idx:(idx + len - 1)],
            nN, G, inter, t, R0, π, μ)
    end
end

end
