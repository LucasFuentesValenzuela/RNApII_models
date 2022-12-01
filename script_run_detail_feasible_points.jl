using JLD2
using CSV
using DataFrames

include("model.jl")
include("experiments.jl")
include("data_utils.jl")
include("parameters_feasibility.jl")

# 1. unpack feasible points
results_fs = JLD2.load("results/feasible_pts.jld2")
occupancy = results_fs["occupancy"]
promoter_occ = results_fs["promoter_occ"]
params_occ = results_fs["params_occ"]
params_iter = results_fs["params_iter"]
feasible = results_fs["feasible"]

# 2. Load the ChIP data
RNA_free_avgCell = CV_to_RNAfree_interp(avg_cell_size)

# 3. Run the feasible points
#   we run it n_times because we want to get some statistics

n_times = 10
n_kon_pts = 10

params_iter_feasible = []
kon_to_CV_interps = []
CV_interps = LinRange(30, 200, n_kon_pts)
RNA_free_interps = CV_to_RNAfree_interp.(CV_interps)

k_on_vec = params_iter[1][3]
α_vec = params_iter[1][1]

for (idx_k, idx_α) in Tuple.(findall(feasible .== 1))

    k_crt = params_iter[1][3][idx_k]
    α_crt = params_iter[idx_α][1]

    kon_C_crt = k_crt / RNA_free_avgCell

    k_on_vec_crt = RNA_free_interps .* kon_C_crt

    push!(params_iter_feasible, (α_crt, β, k_on_vec_crt))

    kon_to_CV_interp = linear_interpolation(
            RNA_free_frac * kon_C_crt, df[!, :cell_volume_fL]
    )
    push!(kon_to_CV_interps, kon_to_CV_interp)
end


occupancy_feasible = []
promoter_occ_feasible = []
params_occ_feasible = []

println("Computing the occupancy as a function of kon for feasible points.")
for iter in ProgressBar(1:n_times)
    occupancy_feasible_crt, promoter_occ_feasible_crt, params_occ_feasible_crt = run_occupancy_simulation(
        params_iter_feasible, Ω, δ, γ, L, Δt, nothing, n_sites, n_end_sites; n_events=n_events
    )
    push!(occupancy_feasible, occupancy_feasible_crt)
    push!(promoter_occ_feasible, promoter_occ_feasible_crt)
    push!(params_occ_feasible, params_occ_feasible_crt)
end

# 5. Save the data
JLD2.jldsave(
    "results/feasible_points_detailed_test.jld2"; 
    occupancy_feasible, promoter_occ_feasible, params_occ_feasible, params_iter_feasible,
    kon_to_CV_interps
)