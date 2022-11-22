include("model.jl")
include("plot_utils.jl")

"""
Run simulations of sweeping kon and alpha
"""
function run_occupancy_simulation(params_iter, k_on_vec, Ω, δ, γ, L, Δt, nsteps, nsites, n_end_sites)
    occupancy = []
	promoter_occ = []
	params_occ = []

	for (α, β) in params_iter
			k_off = max(0, 1/Ω-α)
			
			occupancy_crt = []
			prom_occ_crt = []
	
			for k_on in k_on_vec
			
				params_crt = Params(
				α, β/δ, γ, L, k_on, k_off, 
				Δt, nsteps, 
				nsites, n_end_sites, 
				β/8/δ
				)
			
				_, density, _, _ = run_walker(params_crt);
			
				push!(
					occupancy_crt, 
					get_total_occupancy(density, params_crt; start_bp=2)
				)
				push!(
					prom_occ_crt, 
					get_total_occupancy(
						density, params_crt; start_bp=1, end_bp = 1
					)
				)
				
			end
	
			push!(occupancy, occupancy_crt)
			push!(promoter_occ, prom_occ_crt)
			push!(params_occ, "α=$(round(α; digits=3)), β=$(round(β; digits=3))")
			
	end

    return occupancy, promoter_occ, params_occ
end