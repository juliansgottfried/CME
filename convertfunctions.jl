module Helper

import DelimitedFiles

function converter!(idx, mixture, R0, π, μ, nbin, nN, tmpcounts, mixed, shapely, input)
    tmpcounts .= DelimitedFiles.readdlm("results/counts_$(idx).csv", ',')'

    for j in eachindex(π)
        @views tmptmp = tmpcounts[(1:length(μ)) .+ length(μ) .* (j - 1), 1:(nbin * nN)]
        mixed[j, 1:(nbin * nN)] .= tmptmp' * (mixture ./ sum(mixture))
    end
    mixed[:, nbin * nN + 1] .= R0[idx]
    mixed[:, nbin * nN + 2] .= π

    for j in 1:nN
        shapely[length(π) * (j - 1) .+ (1:length(π)), nbin + 1] .= j
        shapely[length(π) * (j - 1) .+ (1:length(π)), [1:nbin; nbin .+ (2:3)]] .= 
            mixed[:, [nbin * (j - 1) .+ (1:nbin); nbin * nN .+ (1:2)]]
    end
    
    input .= shapely
end

function makeshapely!(mixture, R0, π, μ, nbin, nN, shapeliest, l, chunk)
    tmpcounts = zeros(Float64, length(π) * length(μ), nbin * nN)
    mixed = zeros(Float64, length(π), nbin * nN + 2)
    shapely = zeros(Float64, length(π) * nN, nbin + 3)
    for i in 1:l
        idx = i + l * (chunk - 1)
        @views converter!(idx, mixture, R0, π, μ, nbin, nN, tmpcounts, mixed, shapely, shapeliest[(1:(nN * length(π))) .+ nN * length(π) * (i - 1), :])
        tmpcounts .= 0
        mixed .= 0
        shapely .= 0
    end
end

end