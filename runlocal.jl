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

for i in 1:100
    pars = zeros(Float64, 3, 100 ^ 2)
    for j in 1:100
        for k in 1:100
            pars[:, 100 * (j - 1) + k] = [β[i], μ[j], p[k]]
        end
    end
    writedlm("writepars/pars_$i.csv", pars, ',')
end

using Plots
Ss = zeros(Int, G)
Is = zeros(Int, G)
Ss[1] = 100
for g in 2:G
    Ss[g], Is[g] = sim.iterate(0.73, 11, 0.0001, Ss[g - 1], Is[g - 1], Δt)
end
plot(Is)
