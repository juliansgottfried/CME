#!/usr/bin/env julia

using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere include("/scratch/users/jgottf/CME/sim.jl")
@everywhere import StatsBase, DelimitedFiles

@everywhere R0 = (0:0.02:(2 - 0.02)) ./ 100
@everywhere π = 0:0.0001:(0.01 - 0.0001)
@everywhere μ = [54; 57; 63; 64; 66;
                74; 78; 79; 82; 85;
                86; 89; 90; 92; 95; 
                96; 98; 101; 105; 107; 
                108; 110; 111; 112; 118; 
                119; 121; 123; 125; 128; 
                136; 160]

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
