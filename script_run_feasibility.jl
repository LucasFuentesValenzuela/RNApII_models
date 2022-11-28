using JLD2
using CSV
using DataFrames

include("model.jl")
include("experiments.jl")
include("data_utils.jl")

# Define hyper-parameters for the simulations
Ω = min_Ω #(min_Ω + max_Ω)/2
ρ_p = (min_ρ_p + max_ρ_p)/2
Ψ = (min_Ψ + max_Ψ) / 2
β = 33

γ = 10
nsteps = 1e6
L = 1
δ = 35 / L # how much we "compress the phenomenon"
n_sites = Int(round(gL/δ))
n_end_sites = Int(round(300/δ))
Δt = 1e-2
length_gene = gL


n_kon_pts = 10
k_on_vec = 10 .^(LinRange(
    log10(min_k_on/1.5), log10(max_k_on*1.5), n_kon_pts
));

n_α_values = 8
α_vec = LinRange(min_α, max_α, n_α_values)

# 1. Simulations to determine feasible points in a narrow range
params_iter = collect(Iterators.product(α_vec, [β], [k_on_vec]))

occupancy, promoter_occ, params_occ = run_occupancy_simulation(
    params_iter, Ω, δ, γ, L, Δt, nsteps, n_sites, n_end_sites
)

# 2. Determine feasible points
occ_mat = reduce(hcat, occupancy)
occ_mat = (occ_mat .< max_ρ_g) .& (occ_mat .> min_ρ_g)
prom_occ_mat = reduce(hcat, promoter_occ)
prom_occ_mat = (prom_occ_mat .< max_ρ_p) .& (prom_occ_mat .> min_ρ_p)
feasible = (occ_mat .& prom_occ_mat)

# 3. Load the ChIP data
RNA_free_avgCell = CV_to_RNAfree_interp(avg_cell_size)

# 4. Run the feasible points
params_iter_feasible = []
kon_to_CV_interps = []
for (idx_k, idx_α) in Tuple.(findall(feasible .== 1))

    k_crt = k_on_vec[idx_k]
    α_crt = α_vec[idx_α]

    kon_C_crt = k_crt / RNA_free_avgCell

    kon_to_CV_interp = linear_interpolation(
            RNA_free_frac * kon_C_crt, df[!, :cell_volume_fL]
    )

    k_on_vec_crt = (LinRange(
        minimum(RNA_free_frac)/RNA_free_avgCell * 1.01, 
        maximum(RNA_free_frac)/RNA_free_avgCell * .99, 
        n_kon_pts
        )) .* k_crt

    push!(params_iter_feasible, (α_crt, β, k_on_vec_crt))
    push!(kon_to_CV_interps, kon_to_CV_interp)
end


occupancy_feasible, promoter_occ_feasible, params_occ_feasible = run_occupancy_simulation(
    params_iter_feasible, Ω, δ, γ, L, Δt, nsteps, n_sites, n_end_sites
)

# 5. Save the data
JLD2.jldsave(
    "results/feasible_points.jld2"; 
    occupancy, promoter_occ, params_occ, params_iter, 
    feasible, 
    occupancy_feasible, promoter_occ_feasible, params_occ_feasible, params_iter_feasible,
    kon_to_CV_interps
)
