#!/usr/bin/env julia

using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere include("/scratch/users/jgottf/CME/sim_me_2.jl")
@everywhere import StatsBase, DelimitedFiles

@everywhere plugins = DelimitedFiles.readdlm("/scratch/users/jgottf/CME/plugins.csv", ',')
@everywhere μ = plugins[:, 1]
@everywhere λ = plugins[:, 2] ./ plugins[:, 1] .- 1

@everywhere R0 = (0:0.0004:(0.2 - 0.0004))
@everywhere π = 0:0.0002:(0.02 - 0.0002)
@everywhere ϕ = 0

@everywhere J = 100
@everywhere G = 50
@everywhere inter = 1
@everywhere nbin = 25

@everywhere conds = DelimitedFiles.readdlm("/scratch/users/jgottf/CME/conds.csv", Int)
@everywhere nN = size(conds)[1]

pmap(eachindex(μ)) do i
    counts = zeros(Int, nbin * nN, length(R0) * length(π))
    len = Int(G / inter)
    all = zeros(Int, nN, J * len)
    t = zeros(Int, nN)
    for j in eachindex(R0)
        for k in eachindex(π)
            println("i: $i, j: $j, k: $k")
            β = R0[j] * (1 + λ[i] * (1 + ϕ * (1 - π[k])))
            Sim.replication!(all, J, G, inter, conds, nN, len, t, β, π[k], ϕ, μ[i], λ[i])
            all[all .> nbin - 1] .= nbin - 1
            all .+= 1
            for u in 1:nN
                counts[nbin * (u - 1) .+ (1:nbin), length(π) * (j - 1) + k] .= StatsBase.counts(all[u, :], 1:nbin)
            end
        end
    end
    DelimitedFiles.writedlm("/scratch/users/jgottf/CME/results/counts_$(i).csv", counts, ',')
end
