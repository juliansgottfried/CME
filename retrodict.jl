#!/usr/bin/env julia

using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere using JLD2, DelimitedFiles
@everywhere import Distributions

@everywhere N = 42
@everywhere O = 50
@everywhere M = 100
@everywhere J = 1000
@everywhere nthreads = 200
@everywhere parsper =  60

@everywhere pars = readdlm("/scratch/users/jgottf/CME/pars.csv", ',', Any, '\n')
@everywhere npars = size(pars)[1]
@everywhere logpars = log.(pars)
@everywhere logpars0 = logpars[pars[:, O + 3] .== 0, :]

@everywhere idx = Distributions.sample(1:npars, nthreads * parsper, replace = false)
@everywhere selected = pars[idx, :]
@everywhere writedlm("/scratch/users/jgottf/CME/retrodiction/selected.csv", selected[:, (O + 1):(O + 3)], ',')

pmap(1:nthreads) do i
    for j in 1:parsper
        k = parsper * (i - 1) + j
        tmpdistr = convert(Array{Float64, 1}, selected[k, 1:O])
        tmpdistr ./= sum(tmpdistr)
        tmpdata = zeros(Int, O, M)
        for m in 1:M
            tmpdata[:, m] = rand(Distributions.Multinomial(N, tmpdistr))
            tmpcounts = reduce(vcat, [repeat([o], inner = tmpdata[o, m]) for o in 1:O])
            preserve = zeros(Float64, 4, J + 1)
            for b in 1:(J + 1)
                if (b == 1) tmptmpdat = tmpdata[:, m]
                else 
                    tmptmpcounts = Distributions.sample(tmpcounts, N, replace = true)
                    tmptmpdat = [sum(tmptmpcounts .== o) for o in 1:O]
                end
                factor = Float64(log(factorial(big(N)) / prod(factorial.(big.(tmptmpdat)))))
                loglik = logpars[:, 1:O] * tmptmpdat
                loglik0 = logpars0[:, 1:O] * tmptmpdat
                preserve[:, b] = [argmax(loglik); maximum(loglik) + factor; argmax(loglik0); maximum(loglik0) + factor]
            end
            writedlm("/scratch/users/jgottf/CME/retrodiction/stats_$(k)_$m.csv", preserve, ',')
        end
        writedlm("/scratch/users/jgottf/CME/retrodiction/dat_$(k).csv", tmpdata, ',')
    end
end
