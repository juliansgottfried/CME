module sim

function iterate(g, S, I, αS, αI, R0)
    FOI = R0 * S * I / (S + I)
    λ = αS + αI + FOI + S + I
    g -= 1 / λ * log(rand())
    dart = rand()
    if dart < αS / λ S += 1
    elseif dart < (αS + αI) / λ I += 1
    elseif dart < (αS + αI + FOI) / λ
        S -= 1
        I += 1
    elseif dart < (αS + αI + FOI + S) / λ S -= 1
    else I -= 1
    end
    (g, S, I)
end

function loop!(Is, G, inter, αS, αI, R0)
    S = 25
    I = 0
    g = 0
    t = 0
    while g < G
        if g ≥ inter * t 
            t += 1
            Is[t] = I
        end
        g, S, I = iterate(g, S, I, αS, αI, R0)
    end
end    

function replication(J, G, inter, αS, αI, R0)
    len = Int(G / inter)
    all = zeros(Int, J * len)
    for j in 1:J
        idx = len * (j - 1) + 1
        @views loop!(all[idx:(idx + len - 1)],
            G, inter, αS, αI, R0)
    end
    all
end

end
