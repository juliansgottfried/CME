module Sim

function iterate(g, S, I, β, π, ϕ, μ, λ)
    N = S + I
    Sup = (1 - π) * (μ + λ * N)
    Iup = π * (μ + λ * N)
    Sdown = (1 + λ) * S
    Idown = (1 + λ) * I
    SdownIup = (π * ϕ * λ + β * I) * S
    IdownSup = (1 - π) * ϕ * λ * I
    rate = Sup + Iup + Sdown + Idown + SdownIup + IdownSup
    g -= 1 / rate * log(rand())
    dart = rand()
    if dart < Sup / rate S += 1
    elseif dart < (Sup + Iup) / rate I += 1
    elseif dart < (Sup + Iup + Sdown) / rate S -= 1
    elseif dart < (Sup + Iup + Sdown + Idown) / rate I -= 1
    elseif dart < (Sup + Iup + Sdown + Idown + SdownIup) / rate
        S -= 1
        I += 1
    else
        I -= 1
        S += 1
    end
    (g, S, I)
end

function loop!(Is, conds, nN, G, inter, t, β, π, ϕ, μ, λ)
    g = 0
    S = 25
    I = 0
    t .= 0
    while g < G
        N = S + I
        Nidx = (1:nN)[N .== conds]
        if size(Nidx)[1] > 0
            if g ≥ inter * t[Nidx[1]]
                t[Nidx[1]] += 1
                Is[Nidx[1], t[Nidx[1]]] = I
            end
        end
        g, S, I = iterate(g, S, I, β, π, ϕ, μ, λ)
    end
end    

function replication!(all, J, G, inter, conds, nN, len, t, β, π, ϕ, μ, λ)
    all .= 0
    for j in 1:J
        idx = len * (j - 1) + 1
        @views loop!(all[:, idx:(idx + len - 1)],
            conds, nN, G, inter, t, β, π, ϕ, μ, λ)
    end
end

end
