include("sim_me_2.jl")

ncme = 353

M = 6000

track = zeros(Float64, 3 * ncme, M)

for idx in 1:ncme
    println(idx)
    @views Sim.rw!(track[(1:3) .+ 3 * (idx - 1), :], idx)
end

DelimitedFiles.writedlm("mcmc.csv", track, ',')
