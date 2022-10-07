#!/bin/sh
julia sweep_params.jl --param gamma
julia sweep_params.jl --param gamma --L 10
julia sweep_params.jl --param gamma --L 20
julia sweep_params.jl --param gamma --L 35

julia sweep_params.jl --param L
julia sweep_params.jl --param L --gamma 0.01
julia sweep_params.jl --param L --gamma 0.1
julia sweep_params.jl --param L --gamma 1.
julia sweep_params.jl --param L --gamma 10.
