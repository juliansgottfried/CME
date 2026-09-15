#!/usr/bin/env julia

using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere include("/scratch/users/jgottf/CME/simnb.jl")
@everywhere import StatsBase, DelimitedFiles

# @everywhere R0 = 0:0.05:(5 - 0.05)
# @everywhere π = 0:0.0025:(0.25 - 0.0025)
# @everywhere ϕ = 0:0.1:(10 - 0.1)

@everywhere R0 = 0:0.025:(2.5 - 0.025)
@everywhere π = 0:0.0001:(0.01 - 0.0001)
@everywhere ϕ = 0

@everywhere μ = 97
@everywhere σ2 = 516
@everywhere λ = σ2 / μ - 1

@everywhere J = 100
@everywhere G = 50
@everywhere inter = 1
@everywhere nbin = 50

@everywhere nN = 160

pmap(eachindex(R0)) do i
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
    DelimitedFiles.writedlm("/scratch/users/jgottf/CME/results/counts_$(i).csv", counts, ',')
end
