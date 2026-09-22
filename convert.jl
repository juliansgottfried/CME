#!/usr/bin/env julia

using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere include("/scratch/users/jgottf/CME/helpconvert.jl")
@everywhere import DelimitedFiles

@everywhere mixtures = DelimitedFiles.readdlm("/scratch/users/jgottf/CME/mixtures.csv", ',')

@everywhere R0 = (0:0.002:(0.2 - 0.002))
@everywhere π = 0:0.0005:(0.05 - 0.0005)
@everywhere μ = 1:150
@everywhere nbin = 50
@everywhere nN = 168

@everywhere l = 10
@everywhere chunk = 9

pmap(1:size(mixtures)[1]) do i
    shapeliest = zeros(Float64, nN * length(π) * l, nbin + 3)
    Helper.makeshapely!(mixtures[i, :], R0, π, μ, nbin, nN, shapeliest, l, chunk)
    DelimitedFiles.writedlm("/scratch/users/jgottf/CME/cme_results_marc_2/counts_$(i)_$(chunk).csv", shapeliest, ',')
end
