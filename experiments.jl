include("model.jl")
include("plot_utils.jl")

"""
Run simulations of sweeping kon and alpha
"""
function run_occupancy_simulation(params_iter, Ω, δ, γ, L, Δt, nsteps, nsites, n_end_sites; n_events=nothing)

    occupancy = []
	promoter_occ = []
	params_occ = []

	for (α, β, k_on_vec) in params_iter
			k_off = max(0, 1/Ω-α)
			
			occupancy_crt = []
			prom_occ_crt = []
			params_list = []
	
			for k_on in k_on_vec

				# adaptively select the Δt and γ
				if Δt === nothing
					Δt_crt = set_Δt(α, β/δ, β/8/δ, k_on, k_off, γ)
				end
				if γ === nothing
					γ_crt = 1/Δt_crt
				end
				if n_events !== nothing
					nsteps = Int(round(n_events/k_on))
				end

			
				params_crt = Params(
					α, β/δ, γ_crt, L, k_on, k_off, 
					Δt_crt, nsteps, 
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
				push!(params_list, params_crt)
				
			end
	
			push!(occupancy, occupancy_crt)
			push!(promoter_occ, prom_occ_crt)
			push!(params_occ, params_list)
			
	end

    return occupancy, promoter_occ, params_occ
end