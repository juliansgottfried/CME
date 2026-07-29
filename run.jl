#!/usr/bin/env julia

using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere include("/scratch/users/jgottf/CME/sim.jl")
@everywhere using JLD2, DelimitedFiles

# @everywhere α = 0:1:99
# @everywhere f = 0:0.005:0.495
# @everywhere R0 = 0:0.1:9.9

@everywhere αS = 0:1:99
@everywhere αI = 0:0.1:9.9
@everywhere R0 = 0:0.1:9.9

@everywhere J = 200
@everywhere G = 7000
@everywhere inter = 100

pmap(1:100) do i
    counts = zeros(Int, 50, 100 ^ 2)
    for j in 1:100
        for k in 1:100
            println("i: $i, j: $j, k: $k")
            all = sim.replication(J, G, inter, αS[i], αI[j], R0[k])
            all[all .> 49] .= 49
            all .+= 1
            [counts[l, 100 * (j - 1) + k] += 1 for l in all]
        end
    end
    writedlm("/scratch/users/jgottf/CME/JLresults6/counts_$i.csv", counts, ',')
end
