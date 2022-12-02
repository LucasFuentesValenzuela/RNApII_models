# running the model for feasible points for a large range of kon
# values so that one can get to large and small values of occupancy
using JLD2
using CSV
using DataFrames
using ProgressBars

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

# 3. Run the feasible points

params_iter_feasible = []

for (idx_k, idx_α) in Tuple.(findall(feasible .== 1))

    k_crt = params_iter[1][3][idx_k]
    α_crt = params_iter[idx_α][1]

    k_on_vec_crt = k_crt .* (10 .^(LinRange(-1.5, 1.5, 10)))

    push!(params_iter_feasible, (α_crt, β, k_on_vec_crt))

end


n_times = 10

occupancy_feasible = []
promoter_occupancy_feasible = []
params_occ_feasible = []

println("Computing the occupancy as a function of kon for feasible points.")
for iter in ProgressBar(1:n_times)
    
    occupancy_feasible_crt, promoter_occ_feasible_crt, params_occ_feasible_crt = run_occupancy_simulation(
        params_iter_feasible, Ω, δ, γ, L, Δt, nothing, n_sites, n_end_sites; n_events=n_events
    )
    push!(occupancy_feasible, occupancy_feasible_crt)
    push!(promoter_occupancy_feasible, promoter_occ_feasible_crt)
    push!(params_occ_feasible, params_occ_feasible_crt)
end

# 5. Save the data
JLD2.jldsave(
    "results/feasible_points_wide_test.jld2"; 
    occupancy_feasible, promoter_occ_feasible, params_occ_feasible, params_iter_feasible,
)