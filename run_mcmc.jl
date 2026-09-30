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

@everywhere idcs = [2; 4; 8; 10; 17; 18; 20; 23; 24; 25; 26; 27; 29; 35; 40; 42; 67; 71; 74; 75; 78; 80; 81; 84; 85; 87; 90; 91; 92; 94; 97; 98; 99; 103; 113; 116; 123; 124; 125; 127; 129; 134; 135; 139; 140; 143; 146; 147; 151; 153; 154; 155; 157; 159; 162; 163; 165; 166; 171; 176; 179; 180; 181; 183; 186; 192; 193; 195; 196; 201; 204; 207; 208; 210; 213; 214; 216; 219; 224; 225; 227; 229; 233; 237; 240; 245; 248; 252; 254; 259; 261; 263; 264; 268; 269; 273; 282; 283; 285; 287; 291; 294; 295; 297; 304; 306; 308; 310; 313; 318; 319; 320; 324; 326; 329; 330; 331; 335; 337; 348; 353]

pmap(idcs) do idx
    track = zeros(Float64, 3, M)
    Sim.loopmcmc!(proposal, θ, ϕ,
            J, G, inter, len,
            nbin, loglik, track, 
            rate, M, idx, dat)
    DelimitedFiles.writedlm("/scratch/users/jgottf/CME/results/fit_$(idx).csv", track, ',')
end
