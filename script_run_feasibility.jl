using JLD2
using CSV
using DataFrames

include("model.jl")
include("experiments.jl")
include("data_utils.jl")
include("parameters_feasibility.jl")


n_kon_pts = 10
k_on_vec = 10 .^(LinRange(
    log10(min_k_on/1.5), log10(max_k_on*1.5), n_kon_pts
));

n_α_values = 8
α_vec = LinRange(min_α, max_α, n_α_values)

# 1. Simulations to determine feasible points in a narrow range
params_iter = collect(Iterators.product(α_vec, [β], [k_on_vec]))

println("Computing the feasible points.")
occupancy, promoter_occ, params_occ = run_occupancy_simulation(
    params_iter, Ω, δ, γ, L, Δt, nsteps, n_sites, n_end_sites
)

# 2. Determine feasible points
occ_mat = reduce(hcat, occupancy)
occ_mat = (occ_mat .< max_ρ_g) .& (occ_mat .> min_ρ_g)
prom_occ_mat = reduce(hcat, promoter_occ)
prom_occ_mat = (prom_occ_mat .< max_ρ_p) .& (prom_occ_mat .> min_ρ_p)
feasible = (occ_mat .& prom_occ_mat)

# 3. Save the data
JLD2.jldsave(
    "results/feasible_pts.jld2"; 
    occupancy, promoter_occ, params_occ, params_iter, 
    feasible
)

