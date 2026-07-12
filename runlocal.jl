include("sim.jl")
using JLD2, DelimitedFiles

β = 0:0.01:0.99
μ = 0:10:990
p = 0:0.005:0.495

Δt = 0.05
G = 25000
J = 400
inter = 500

for i in 1:100
    counts = zeros(Int, 100, 100 ^ 2)
    for j in 1:100
        for k in 1:100
            println("i: $i, j: $j, k: $k")
            all = sim.replication(β[i], μ[j], p[k], Δt, G, J, inter)
            all[all .> 99] .= 99
            all .+= 1
            [counts[l, 100 * (j - 1) + k] += 1 for l in all]
        end
    end
    writedlm("JLresults2/counts_$i.csv", counts, ',')
end
