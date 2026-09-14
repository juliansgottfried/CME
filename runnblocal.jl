include("simnb.jl")
import StatsBase, DelimitedFiles

R0 = 5 - 0.05
π = 0.01
ϕ = 0:0.1:(10 - 0.1)

μ = 97
σ2 = 516
λ = σ2 / μ - 1

J = 100
G = 50
inter = 1
nbin = 50

nN = 160

@time for i in eachindex(R0)
    counts = zeros(Int, nbin * nN, length(π) * length(ϕ))
    len = Int(G / inter)
    all = zeros(Int, nN, J * len)
    t = zeros(Int, nN)
    for j in eachindex(π)
        for k in eachindex(ϕ)
            println("i: $i, j: $j, k: $k")
            β = R0[i] * (1 + λ * (1 + ϕ[k] * (1 - π[j])))
            Sim.replication!(all, J, G, inter, nN, len, t, β, π[j], ϕ[k], μ, λ)
            all[all .> nbin - 1] .= nbin - 1
            all .+= 1
            for u in 1:nN
                counts[nbin * (u - 1) .+ (1:nbin), length(ϕ) * (j - 1) + k] .= StatsBase.counts(all[u, :], 1:nbin)
            end
        end
    end
    # DelimitedFiles.writedlm("/scratch/users/jgottf/CME/results/counts_$(i).csv", counts, ',')
end
