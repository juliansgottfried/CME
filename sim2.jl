module sim2

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

function loop!(Is, G, inter, α, f, R0)
    S = 25
    I = 0
    g = 0
    t = 0
    while g < G
        if g ≥ inter * t 
            t += 1
            Is[t] = I
        end
        g, S, I = iterate(g, S, I, α, f, R0)
    end
end    

function replication(J, G, inter, α, f, R0)
    len = Int(G / inter)
    all = zeros(Int, J * len)
    for j in 1:J
        idx = len * (j - 1) + 1
        @views loop!(all[idx:(idx + len - 1)],
            G, inter, α, f, R0)
    end
    all
end

end
