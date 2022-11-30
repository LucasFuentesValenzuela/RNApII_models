include("model.jl")

# Define hyper-parameters for the simulations
Ω = min_Ω #(min_Ω + max_Ω)/2
ρ_p = (min_ρ_p + max_ρ_p)/2
Ψ = (min_Ψ + max_Ψ) / 2
β = 33

nsteps = 1e6
L = 1
δ = 35 / L # how much we "compress the phenomenon"
n_sites = Int(round(gL/δ))
n_end_sites = Int(round(10*δ/δ))
length_gene = gL

# we set them to nothing so that they are adaptively set
γ = 10
Δt = 1e-2

# how many binding events you expect on average
n_events = 5e3