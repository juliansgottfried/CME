#!/usr/bin/env julia

using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere include("/scratch/users/jgottf/CME/convertfunctions.jl")
@everywhere import DelimitedFiles

@everywhere mixtures = DelimitedFiles.readdlm("/scratch/users/jgottf/CME/mixtures.csv", ',')[2:115, :]

@everywhere R0 = (0:0.002:(0.2 - 0.002))
@everywhere π = 0:0.0005:(0.05 - 0.0005)
@everywhere μ = 1:150
@everywhere nbin = 50
@everywhere nN = 168

pmap(1:size(mixtures)[1]) do idx
    shapeliest = Helper.makeshapely(mixtures[idx, ], R0, π, μ, nbin, nN)
    DelimitedFiles.writedlm("/scratch/users/jgottf/CME/cme_results/counts_$(idx).csv", shapeliest, ',')
end
