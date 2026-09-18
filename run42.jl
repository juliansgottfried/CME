#!/usr/bin/env julia

using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere include("/scratch/users/jgottf/CME/sim42.jl")
@everywhere import StatsBase, DelimitedFiles

@everywhere R0 = (0:0.001:(0.1 - 0.001))
@everywhere π = 0:0.0005:(0.05 - 0.0005)
@everywhere μ = [1:136; 138; 139; 141; 142; 
                144; 147; 150; 151; 153; 
                154; 156; 157; 160; 161; 
                168]

@everywhere J = 100
@everywhere G = 50
@everywhere inter = 1
@everywhere nbin = 50

pmap(eachindex(R0)) do i
    counts = zeros(Int, nbin, length(π) * length(μ))
    len = Int(G / inter)
    all = zeros(Int, J * len)
    for j in eachindex(π)
        for k in eachindex(μ)
            println("i: $i, j: $j, k: $k")
            Sim.replication!(all, J, G, inter, len, R0[i], π[j], μ[k])
            all[all .> nbin - 1] .= nbin - 1
            all .+= 1
            counts[:, length(μ) * (j - 1) + k] .= StatsBase.counts(all, 1:nbin)
        end
    end
    DelimitedFiles.writedlm("/scratch/users/jgottf/CME/results/counts_$(i).csv", counts, ',')
end
