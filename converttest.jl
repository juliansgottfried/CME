include("convertfunctionstest.jl")
import DelimitedFiles

mixtures = DelimitedFiles.readdlm("mixtures.csv", ',')[2:115, :]

R0 = (0:0.002:(0.2 - 0.002))
π = 0:0.0005:(0.05 - 0.0005)
μ = 1:150
nbin = 50
nN = 160

l = 10
chunk = 5

@time for i in 1:1
    shapeliest = zeros(Float64, nN * length(π) * l, nbin + 3)
    Helper.makeshapely!(mixtures[i, :], R0, π, μ, nbin, nN, shapeliest, l , chunk)
    # DelimitedFiles.writedlm("/scratch/users/jgottf/CME/cme_results/counts_$(i).csv", shapeliest, ',')
end
