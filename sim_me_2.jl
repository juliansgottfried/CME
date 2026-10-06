module Sim

import StatsBase, Distributions

function iterate(g, S, I, β, π, ϕ, μ, λ, γ)
    N = S + I
    Sup = (1 - π) * (μ + λ * N)
    Iup = π * (μ + λ * N)
    Sdown = (1 + λ) * S
    Idown = (1 + λ) * I
    SdownIup = (π * ϕ * λ + β * I) * S
    IdownSup = ((1 - π) * ϕ * λ + γ) * I
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

function loop!(Is, nN, G, inter, t, β, π, ϕ, μ, λ, γ)
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
        g, S, I = iterate(g, S, I, β, π, ϕ, μ, λ, γ)
    end
end    

function replication!(all, J, G, inter, nN, len, t, β, π, ϕ, μ, λ, γ)
    all .= 0
    for j in 1:J
        idx = len * (j - 1) + 1
        @views loop!(all[:, idx:(idx + len - 1)],
            nN, G, inter, t, β, π, ϕ, μ, λ, γ)
    end
end

function trans(θ, lim)
    log.(θ ./ (lim .- θ))
end

function invtrans(x, lim)
    lim ./ (1 .+ exp.(-x))
end

function proposal(d, cov)
   Distributions.MvNormal(zeros(Float64, d), cov)
end

function estimatecov(d, track, cov, i, lim)
    est = StatsBase.cov(eachcol(broadcast(trans, track[1:d, 1:(i - 1)], lim)))
    (1 - 0.05) ^ 2 .* 2.38 ^ 2 .* est ./ d .+ 0.05 ^ 2 .* cov
end

function freqs!(all, nbin, nN, counts, J, G, inter)
    all .= min.(all, nbin - 1) .+ 1
    for u in 1:nN
        @views counts[:, u] .= StatsBase.counts(all[u, :], 1:nbin)
    end
    counts .= log.((counts .+ 0.1) / (J * G / inter + 0.1 * nbin))
end

function mcmcbody!(i,
                θ, λ, π, d, cov,
                all, J, G, inter, nN, len, t, μ,
                nbin, conds, counts, slimmed,
                Idat, fac, lim, ϵ,
                loglik, track, rate, M)

    println(θ)
    if i <= 2d
        propdistr = proposal(d, cov)
    else
        propdistr = proposal(d, estimatecov(d, track, cov, i, lim))
    end

    if (i > 1)
        newθ = invtrans(trans(θ, lim) .+ Distributions.rand(propdistr), lim)
    else newθ = θ
    end

    # β = newθ[1] * (1 + newθ[2] + λ * (1 + newθ[3] * (1 - π)))
    β = newθ[1] * (1 + λ)

    # replication!(all, J, G, inter, nN, len, t, β, π, newθ[3], μ, λ, newθ[2])
    replication!(all, J, G, inter, nN, len, t, β, π, 0, μ, λ, 0)
    freqs!(all, nbin, nN, counts, J, G, inter)
    
    slimmed .= counts[:, conds]
    newloglik = sum(slimmed .* Idat) + fac
    if rand() < min(exp(newloglik - loglik), 1) 
        θ = newθ
        θ[θ .>= lim] .= (lim .- ϵ)[θ .>= lim]
        θ[θ .<= 0] .= ϵ
        loglik = newloglik
        rate += 1 / M
    end
    track[:, i] .= [θ; loglik]
    (θ, loglik, rate)
end

function cmedat(dat, idx, nbin)
    slimdat = dat[dat[:, 1] .== idx, 2:3]

    nN = maximum(slimdat[:, 1])

    conds = sort(unique(slimdat[:, 1]))
    ncond = size(conds)[1]
    condidx = [(1:ncond)[conds .== n][1] for n in slimdat[:, 1]]

    μ = StatsBase.mean(slimdat[:, 1])
    σ2 = StatsBase.var(slimdat[:, 1])
    λ = σ2 / μ - 1
    π = StatsBase.mean(slimdat[:, 2]) / μ

    p = μ / σ2
    r = μ * p / (1 - p)
    if λ > 0 law = Distributions.NegativeBinomial(r, p)
    else law = Distributions.Poisson(μ)
    end

    Idat = zeros(Int, nbin, ncond)
    for n in 1:ncond Idat[:, n] .= StatsBase.counts(slimdat[condidx .== n, 2], 0:(nbin - 1)) end

    fac = convert(Float64, sum(log.(factorial.(big.(sum(eachrow(Idat))))) .- sum(eachrow(log.(factorial.(big.(Idat)))))))
    fac = fac + sum(log.(Distributions.pdf(law, conds)))

    (nN, conds, ncond, μ, λ, π, Idat, fac)
end

function loopmcmc!(θ, J, G, inter, len,
                    nbin, d, cov, lim, ϵ,
                    loglik, track, rate, M,
                    idx, dat)

    nN, conds, ncond, μ, λ, π, Idat, fac = cmedat(dat, idx, nbin)

    counts = zeros(Float64, nbin, nN)
    slimmed = zeros(Float64, nbin, ncond)
    all = zeros(Int, nN, J * len)
    t = zeros(Int, nN)

    for i in 1:M
        println("$idx, $i")
        (θ, loglik, rate) = mcmcbody!(i,
                    θ, λ, π, d, cov,
                    all, J, G, inter, nN, len, t, μ,
                    nbin, conds, counts, slimmed,
                    Idat, fac, lim, ϵ,
                    loglik, track, rate, M)
    end
    rate
end

end
