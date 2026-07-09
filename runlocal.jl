include("sim.jl")
using JLD2, DelimitedFiles

β = 0:0.01:0.99
μ = 1:10:1000
p = 0.0001:0.01:1

Δt = 0.05
G = 25000
J = 400
inter = 500

for i in 1:2
    counts = zeros(Int, 2, 2 ^ 2)
    for j in 1:2
        for k in 1:2
            println("i: $i, j: $j, k: $k")
            all = sim.replication(β[i], μ[j], p[k], Δt, G, J, inter)
            all[all .> 1] .= 1
            all .+= 1
            [counts[l, 2 * (j - 1) + k] += 1 for l in all]
        end
    end
    writedlm("JLresults/counts_$i.csv", counts, ',')
end
