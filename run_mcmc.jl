#!/usr/bin/env julia

using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere include("/scratch/users/jgottf/CME/sim_me_2.jl")
@everywhere import DelimitedFiles, LinearAlgebra

@everywhere dat = DelimitedFiles.readdlm("/scratch/users/jgottf/CME/NI.csv", ',')
@everywhere dat = Int.(dat)

@everywhere ncme = maximum(dat[:, 1])

@everywhere nbin = 25

@everywhere J = 100
@everywhere G = 50
@everywhere inter = 1
@everywhere len = Int(G / inter)

@everywhere M = 6000

@everywhere loglik = -Inf

@everywhere d = 1
@everywhere cov = 0.01 ^ 2 .* LinearAlgebra.I(d) ./ d
@everywhere θ = 0.01 .* ones(Float64, d)

@everywhere lim = 10 * ones(Float64, d)
@everywhere ϵ = 0.0001

pmap(1:ncme) do idx
    track = zeros(Float64, d + 1, M)
    rate = Sim.loopmcmc!(θ, J, G, inter, len,
            nbin, d, cov, lim, ϵ, loglik, track, 
            M, idx, dat)
    DelimitedFiles.writedlm("/scratch/users/jgottf/CME/results/fit_$(idx).csv", track, ',')
end
