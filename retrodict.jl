#!/usr/bin/env julia

using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere path = "/scratch/users/jgottf/CME"

@everywhere include("$path/retrohelpers.jl")
@everywhere using JLD2, DelimitedFiles
@everywhere import Distributions

@everywhere N = 42
@everywhere O = 50
@everywhere M = 100
@everywhere J = 1000
@everywhere nthreads = 200
@everywhere parsper =  5
@everywhere npars = 1000000

@everywhere pars = readdlm("$path/pars.csv", ',', Float64, '\n')

@everywhere bool0 = pars[:, O + 3] .== 0
@everywhere idx0 = Distributions.sample((1:npars)[bool0], div(nthreads * parsper, 2), replace = false)
@everywhere idxnon0 = Distributions.sample((1:npars)[.!bool0], div(nthreads * parsper, 2), replace = false)
@everywhere idx = [idx0; idxnon0]
@everywhere selected = pars[idx, 1:O]

@everywhere writedlm("$path/retrodiction/selected.csv", [idx pars[idx, (O + 1):(O + 3)]], ',')

@everywhere pars = log.(pars[:, 1:O])
@everywhere pars0 = pars[bool0, :]

pmap(1:nthreads) do i 
    for j in 1:parsper 
        helpers.processpar(i, j, N, O, M, J, parsper, pars, pars0, selected, path) 
    end 
end
