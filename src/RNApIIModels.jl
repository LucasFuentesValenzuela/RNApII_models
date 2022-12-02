module RNApIIModels

# Module imports
using Distributions
using ProgressBars
using StatsBase
using Random
using ProgressBars
using Plots
using Statistics
using DataFrames
using CSV
using Interpolations
using JLD2
using ArgParse

# Exports
export J, get_regime, ρ, effective_α, check_γ
export Params, DEFAULT_PARAMS, OCCUPANCY_PARAMS, LITERATURE_PARAMS
export avg_cell_size
export k_on_vec_screen, α_vec_screen, n_times
export run_walker, set_Δt, get_total_occupancy
export run_occupancy_simulation, get_feasible_pts
export load_ChIP_data, RNA_free_frac, Rpb1_occupancy_haploid_interp
export CV_to_RNAfree_interp, RNAfree_to_CV_interp


# Files
include("parameters.jl")

include("RNAp_model.jl")
include("theory.jl")

include("experiments.jl")
include("data_utils.jl")
# include("plot_utils.jl")

end # module
