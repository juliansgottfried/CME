#!/usr/bin/env julia

using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere include("/scratch/users/jgottf/CME/sim.jl")
@everywhere import StatsBase, DelimitedFiles

@everywhere α = [76:80; 84:88]
@everywhere f = 0:0.01:0.99
@everywhere R0 = 0:0.1:9.9

@everywhere J = 200
@everywhere G = 7000
@everywhere inter = 100
@everywhere nbin = 50
@everywhere nN = 168

pmap(eachindex(R0)) do i
    counts = zeros(Int, nbin * nN, length(f) * length(α))
    len = Int(G / inter)
    all = zeros(Int, nN, J * len)
    t = zeros(Int, nN)
    for j in eachindex(α)
        for k in eachindex(f)
            println("i: $i, j: $j, k: $k")
            Sim.replication!(all, J, G, inter, nN, len, t, α[j], f[k], R0[i])
            all[all .> nbin - 1] .= nbin - 1
            all .+= 1
            for u in 1:nN
                counts[nbin * (u - 1) .+ (1:nbin), length(f) * (j - 1) + k] .= StatsBase.counts(all[u, :], 1:nbin)
            end
        end
    end
    println(isa(counts, Array))
    DelimitedFiles.writedlm("/scratch/users/jgottf/CME/tmpdir4/counts_$i.csv", counts, ',')
end
