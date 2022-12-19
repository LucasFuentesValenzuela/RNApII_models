#!/bin/sh
#
julia feasible_analysis.jl --type narrow --ntimes 10 --Omega 1
julia feasible_analysis.jl --type narrow --ntimes 10 --Omega 2
julia feasible_analysis.jl --type narrow --ntimes 10 --Omega 3
julia feasible_analysis.jl --type narrow --ntimes 10 --Omega 4

julia feasible_analysis.jl --type wide --ntimes 10 --Omega 1
julia feasible_analysis.jl --type wide --ntimes 10 --Omega 2
julia feasible_analysis.jl --type wide --ntimes 10 --Omega 3
julia feasible_analysis.jl --type wide --ntimes 10 --Omega 4


