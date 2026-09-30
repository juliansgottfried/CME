include("sim_me_2.jl")
import StatsBase, DelimitedFiles, Distributions

dat = DelimitedFiles.readdlm("NI.csv", ',')
dat = Int.(dat)

ncme = maximum(dat[:, 1])

nbin = 25

J = 100
G = 50
inter = 1
len = Int(G / inter)

M = 6000

ϵ1 = 0.005
ϵ2 = 0.0005
proposal = Distributions.MvNormal([0; 0], [ϵ1; 0;; 0; ϵ2])

ϕ = 0
θ = [0.001; 0.001]

loglik = -Inf
rate = 0

@time for idx in 207:207
    track = zeros(Float64, 3, M)
    Sim.loopmcmc!(proposal, θ, ϕ,
            J, G, inter, len,
            nbin, loglik, track, 
            rate, M, idx, dat)
    # DelimitedFiles.writedlm("/scratch/users/jgottf/CME/results/fit_$(idx).csv", track, ',')
end
