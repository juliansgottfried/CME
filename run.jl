#!/usr/bin/env julia

using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere include("/scratch/users/jgottf/CME/sim.jl")
@everywhere using JLD2, DelimitedFiles

@everywhere β = 0:0.05:4.95
@everywhere μ = 0:1:99
@everywhere p = 0:0.001:0.099

@everywhere Δt = 0.05
@everywhere G = 25000
@everywhere J = 400
@everywhere inter = 500

pmap(1:100) do i
    counts = zeros(Int, 100, 100 ^ 2)
    for j in 1:100
        for k in 1:100
            println("i: $i, j: $j, k: $k")
            all = sim.replication(β[i], μ[j], p[k], Δt, G, J, inter)
            all[all .> 99] .= 99
            all .+= 1
            [counts[l, 100 * (j - 1) + k] += 1 for l in all]
        end
    end
    writedlm("/scratch/users/jgottf/CME/JLresults4/counts_$i.csv", counts, ',')
end
