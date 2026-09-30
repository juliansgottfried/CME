#!/usr/bin/env julia

using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere include("/scratch/users/jgottf/CME/sim_me_2.jl")
@everywhere import StatsBase, DelimitedFiles, Distributions

@everywhere dat = DelimitedFiles.readdlm("/scratch/users/jgottf/CME/NI.csv", ',')
@everywhere dat = Int.(dat)

@everywhere ncme = maximum(dat[:, 1])

@everywhere nbin = 25

@everywhere J = 100
@everywhere G = 50
@everywhere inter = 1
@everywhere len = Int(G / inter)

@everywhere M = 6000

@everywhere ϵ1 = 0.005
@everywhere ϵ2 = 0.0005
@everywhere proposal = Distributions.MvNormal([0; 0], [ϵ1; 0;; 0; ϵ2])

@everywhere ϕ = 0
@everywhere θ = [0.001; 0.001]

@everywhere loglik = -Inf
@everywhere rate = 0

pmap(1:ncme) do idx
    track = zeros(Float64, 3, M)
    Sim.loopmcmc!(proposal, θ, ϕ,
            J, G, inter, len,
            nbin, loglik, track, 
            rate, M, idx, dat)
    DelimitedFiles.writedlm("/scratch/users/jgottf/CME/results/fit_$(idx).csv", track, ',')
end
