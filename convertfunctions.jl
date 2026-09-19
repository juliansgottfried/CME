module Helper

import DelimitedFiles

function converter(i, mixture, R0, π, μ, nbin, nN)
    tmpcounts = DelimitedFiles.readdlm("/scratch/users/jgottf/CME/results/counts_$(i).csv", ',')'

    mixed = zeros(Float64, length(π), nbin * nN + 2)

    for j in eachindex(π)
        tmptmp = tmpcounts[(1:length(μ)) .+ length(μ) .* (j - 1), 1:(nbin * nN)]
        for k in 1:(nbin * nN)
            mixed[j, k] = sum(tmptmp[:, k] .* mixture ./ sum(mixture))
        end
    end
    mixed[:, nbin * nN + 1] .= R0[i]
    mixed[:, nbin * nN + 2] = π

    shapely = zeros(Float64, length(π) * nN, nbin + 3)
    for j in 1:nN
        shapely[length(π) * (j - 1) .+ (1:length(π)), nbin + 1] .= j
        shapely[length(π) * (j - 1) .+ (1:length(π)), [1:nbin; nbin .+ (2:3)]] .= 
            mixed[:, [nbin * (j - 1) .+ (1:nbin); nbin * nN .+ (1:2)]]
    end
    shapely
end

function makeshapely(mixture, R0, π, μ, nbin, nN)
    shapeliest = zeros(Float64, nN * length(π) * length(R0), nbin + 3)
    for i in eachindex(R0)
        shapeliest[(1:(nN * length(π))) .+ nN * length(π) * (i - 1), :] .= converter(i, mixture, R0, π, μ, nbin, nN)
    end
    shapeliest
end

end