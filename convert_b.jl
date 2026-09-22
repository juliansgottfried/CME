#!/usr/bin/env julia

using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere include("/scratch/users/jgottf/CME/helpconvert.jl")
@everywhere import DelimitedFiles

@everywhere R0 = (0:0.002:(0.2 - 0.002))
@everywhere π = 0:0.0005:(0.05 - 0.0005)
@everywhere nbin = 50
@everywhere nN = 168

pmap(1:114) do i
    counts = zeros(Float64, nN * length(R0) * length(π), nbin + 3)
    for j in 1:10
        println("$i, $j")
        Helper.loop!(counts, i, j)
    end
    DelimitedFiles.writedlm("/scratch/users/jgottf/CME/results_marc_2_reshape/counts_$i.csv", counts, ',')
end
