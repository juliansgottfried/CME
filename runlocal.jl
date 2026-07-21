include("sim2.jl")
using JLD2, DelimitedFiles

α = 0:1:99
f = 0:0.005:0.495
R0 = 0:0.1:9.9

J = 200
G = 7000
inter = 100

for i in 1:100
    counts = zeros(Int, 50, 100 ^ 2)
    for j in 1:100
        for k in 1:100
            println("i: $i, j: $j, k: $k")
            all = sim2.replication(J, G, inter, α[i], f[j], R0[k])
            all[all .> 49] .= 49
            all .+= 1
            [counts[l, 100 * (j - 1) + k] += 1 for l in all]
        end
    end
    writedlm("JLresults/counts_$i.csv", counts, ',')
end
