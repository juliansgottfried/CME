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

function loop!(Is, G, inter, R0, π, μ)
    g = 0
    S = 25
    I = 0
    t = 0
    while g < G
        N = S + I
        if N == μ
            if g ≥ inter * t
                t += 1
                Is[t] = I
            end
        end
        g, S, I = iterate(g, S, I, R0, π, μ)
    end
end    

function replication!(all, J, G, inter, len, R0, π, μ)
    all .= 0
    for j in 1:J
        idx = len * (j - 1) + 1
        @views loop!(all[idx:(idx + len - 1)],
            G, inter, R0, π, μ)
    end
end

end
