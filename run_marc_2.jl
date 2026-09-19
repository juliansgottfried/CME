#!/usr/bin/env julia

using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere include("/scratch/users/jgottf/CME/sim_marc_2.jl")
@everywhere import StatsBase, DelimitedFiles

@everywhere R0 = (0:0.002:(0.2 - 0.002))
@everywhere π = 0:0.0005:(0.05 - 0.0005)
@everywhere μ = 1:150

@everywhere J = 100
@everywhere G = 50
@everywhere inter = 1
@everywhere nbin = 50

@everywhere nN = 160

pmap(eachindex(R0)) do i
    counts = zeros(Int, nbin * nN, length(π) * length(μ))
    len = Int(G / inter)
    all = zeros(Int, nN, J * len)
    t = zeros(Int, nN)
    for j in eachindex(π)
        for k in eachindex(μ)
            println("i: $i, j: $j, k: $k")
            Sim.replication!(all, J, G, inter, nN, len, t, R0[i], π[j], μ[k])
            all[all .> nbin - 1] .= nbin - 1
            all .+= 1
            for u in 1:nN
                counts[nbin * (u - 1) .+ (1:nbin), length(μ) * (j - 1) + k] .= StatsBase.counts(all[u, :], 1:nbin)
            end
        end
    end
    DelimitedFiles.writedlm("/scratch/users/jgottf/CME/results/counts_$(i).csv", counts, ',')
end
