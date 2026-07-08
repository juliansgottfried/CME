#!/usr/bin/env julia

using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere include("sim.jl")
@everywhere using JLD2, DelimitedFiles

@everywhere β = 0:0.01:0.99
@everywhere μ = 1:10:1000
@everywhere p = 0.0001:0.01:1

@everywhere Δt = 0.05
@everywhere G = 25000
@everywhere J = 400
@everywhere inter = 500

pmap(1:100) do i
    counts = zeros(Int, 100, 100 ^ 2)
    for j in 1:100
        for k in 1:100
            all = sim.replication(β[i], μ[j], p[k], Δt, G, J, inter)
            all[all .> 99] .= 99
            all .+= 1
            [counts[l, 100 * (j - 1) + k] += 1 for l in all]
        end
    end
    save_object("/scratch/users/jgottf/cme/results/counts_$i.csv", counts)
end
