include("sim.jl")
using JLD2, DelimitedFiles

α = [3; 9:16; 26; 28; 30; 33; 58:67; 69:80; 84:89; 91:94; 96:97]
f = 0:0.01:0.99
R0 = 0:0.1:9.9

α = 0:1
f = 0:1
R0 = 0:1

J = 200
G = 7000
inter = 100
nbin = 50
nN = 168

for i in eachindex(R0)
    counts = zeros(Int, nbin * nN, length(f) * length(α))
    len = Int(G / inter)
    all = zeros(Int, nN, J * len)
    t = zeros(Int, nN)
    for j in eachindex(α)
        for k in eachindex(f)
            println("i: $i, j: $j, k: $k")
            Sim.replication!(all, J, G, inter, nN, len, t, α[j], f[k], R0[i])
            all[all .> nbin - 1] .= nbin - 1
            all .+= 1
            for u in 1:nN
                [counts[nbin * (u - 1) + l, length(f) * (j - 1) + k] += 1 for l in all[u, :]]
            end
        end
    end
    writedlm("JLresultstmp/counts_$i.csv", counts, ',')
end
