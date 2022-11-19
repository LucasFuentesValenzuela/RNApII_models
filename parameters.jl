"""
Default parameters 
"""
LARGE_γ = 1000

α_default = 0.0033
β_default = 0.57
γ_default = 0.014
Δt_default = .01
ratio_β2_default = 8
β2_default = β_default/ratio_β2_default
L_default = 1
kon_default = 1
koff_default = 1

DEFAULT_nsteps = 5e5
DEFAULT_n_sites = 42
DEFAULT_n_end_sites = 10

"""
Ranges of acceptable metrics
"""

# gene length for average gene
gL = 1000

# residence time of mRNA molecules on promoter [s]
min_Ω = 2
max_Ω = 4

# number of mRNA molecules on promoter [-]
min_ρ_p = 0.02
max_ρ_p = 0.04


# number of mRNA molecules in gene body [-]
min_ρ_g = 0.15
max_ρ_g = .4

# total residence time on promoter + gene body [-]
min_Ψ = 23
max_Ψ = 26

# elongation rates [s-1]
min_β = 13
max_β = 40

# residence time in the gene body
min_Λ = min_Ψ - max_Ω
max_Λ = max_Ψ - min_Ω

# on rate or the promoter
max_k_on = max_ρ_p/(1-max_ρ_p)/min_Ω
min_k_on = min_ρ_p/(1-min_ρ_p)/max_Ω

# initiation rate
max_α = max_Λ * max_β/min_Ω/gL
min_α = min_Λ * min_β/max_Ω/gL

# off rate from the promoter
max_k_off = 1/min_Ω-min_α
min_k_off = max(0, 1/max_Ω - max_α)
