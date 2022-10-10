#!/bin/sh
julia sweep_params.jl --param gamma --L 1 --beta2 3
julia sweep_params.jl --param gamma --L 10 --beta2 3
julia sweep_params.jl --param gamma --L 20 --beta2 3
julia sweep_params.jl --param gamma --L 35 --beta2 3
