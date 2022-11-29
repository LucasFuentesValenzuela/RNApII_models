include("model.jl")

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