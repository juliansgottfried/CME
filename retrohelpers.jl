module helpers

using JLD2, DelimitedFiles
import Distributions

bootstrap! = function(b, m, preserve, tmpdata, tmpcounts, N, O, pars, pars0)
    if (b == 1) tmptmpdat = tmpdata[:, m]
    else 
        tmptmpcounts = Distributions.sample(tmpcounts, N, replace = true)
        tmptmpdat = [sum(tmptmpcounts .== o) for o in 1:O]
    end
    factor = Float64(log(factorial(big(N)) / prod(factorial.(big.(tmptmpdat)))))
    loglik = pars * tmptmpdat
    loglik0 = pars0 * tmptmpdat
    preserve[:, b] = [argmax(loglik); maximum(loglik) + factor; argmax(loglik0); maximum(loglik0) + factor]
end

replication! = function(m, k, tmpdata, tmpdistr, N, O, J, pars, pars0, path)
    tmpdata[:, m] = rand(Distributions.Multinomial(N, tmpdistr))
    tmpcounts = reduce(vcat, [repeat([o], inner = tmpdata[o, m]) for o in 1:O])
    preserve = zeros(Float64, 4, J + 1)
    for b in 1:(J + 1)
        bootstrap!(b, m, preserve, tmpdata, tmpcounts, N, O, pars, pars0)
    end
    writedlm("$path/retrodiction/stats_$(k)_$m.csv", preserve, ',')
end

processpar = function(i, j, N, O, M, J, parsper, pars, pars0, selected, path)
    k = parsper * (i - 1) + j
    tmpdistr = selected[k, :]
    tmpdistr ./= sum(tmpdistr)
    tmpdata = zeros(Int, O, M)
    for m in 1:M
        replication!(m, k, tmpdata, tmpdistr, N, O, J, pars, pars0, path)
    end
    writedlm("$path/retrodiction/dat_$(k).csv", tmpdata, ',')
end

end
