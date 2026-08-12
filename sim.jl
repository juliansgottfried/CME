module Sim

function iterate(g, S, I, α, f, R0)
    FOI = R0 * S * I / (S + I)
    λ = α + FOI + S + I
    g -= 1 / λ * log(rand())
    dart = rand()
    if dart < (1 - f) * α / λ S += 1
    elseif dart < α / λ I += 1
    elseif dart < (α + FOI) / λ
        S -= 1
        I += 1
    elseif dart < (α + FOI + S) / λ S -= 1
    else I -= 1
    end
    (g, S, I)
end

function loop!(Is, nN, G, inter, t, α, f, R0)
    g = 0
    S = 25
    I = 0
    t .= 0
    while g < G
        N = S + I
        if N > 0 & N <= nN
            if g ≥ inter * t[N]
                t[N] += 1
                Is[N, t[N]] = I
            end
        end
        g, S, I = iterate(g, S, I, α, f, R0)
    end
end    

function replication!(all, J, G, inter, nN, len, t, α, f, R0)
    all .= 0
    for j in 1:J
        idx = len * (j - 1) + 1
        @views loop!(all[:, idx:(idx + len - 1)],
            nN, G, inter, t, α, f, R0)
    end
end

end
