"""
	run_walker(
		α, β, γ, L, kon, koff, Δt, n_steps, n_sites, n_end_sites; 
		β2=nothing, ratio_β2=5
	)

Runs the simulation for a specified number of steps.

Parameters:
	`α`: initiation rate
	`β`: elongation rate
	`γ`: termination rate
	`L`: particle footprint
	`kon`: binding rate to the promoter
	`koff`: unbinding rate from the promoter
	`Δt`: timestep
	`n_steps`: number of timesteps to run the simulation
	`n_sites`: number of gene sites
	`n_end_sites`: number of sites past the termination point
	`β2`: elongation rate in the termination region
	`ratio_β2`: ratio of elongation rates between the gene body and the termination region
"""
function run_walker(
	α, β, γ, L, kon, koff, Δt, n_steps, n_sites, n_end_sites; 
	β2=nothing, ratio_β2=5
)

	if β2===nothing
		β2 = β/ratio_β2 #slower rate in the second part
	end

	exits = zeros(n_steps)
	gene = zeros(n_sites + n_end_sites); # initializing the sites
	density = zeros(n_sites + n_end_sites); 

	tracker_end = Dict(
		"current" => zeros(n_end_sites), 
		"terminated" => [], 
	) # tracks every RNA getting into the last strand
	
	n_inits = 0
	for k in 1:n_steps # number of timesteps
	
		gene, finish_flag, tracker_end, n_inits = step(
			α, β, β2, γ, L, kon, koff, gene, n_sites, Δt, tracker_end, n_inits
		)
		density += gene

		exits[k] = finish_flag

	end

	return exits, density, gene, tracker_end, n_inits
	
end

run_walker(params) = run_walker(
	params.α, params.β, params.γ, params.L, params.kon, params.koff,
	params.Δt, params.n_steps, params.n_sites, 
	params.n_end_sites; β2=params.β2
)

"""
	step(
		α, β, β2, γ, L, kon, koff, gene, n_sites, Δt, tracker_end, n_inits
	)

Takes a step forward, for a given model.

Parameters:
	`α`: initiation rate
	`β`: elongation rate
	`β2`: elongation rate in the termination region
	`γ`: termination rate
	`L`: particle footprint
	`kon`: binding rate to the promoter
	`koff`: unbinding rate from the promoter
	`gene`: vector representing the occupancy of the gene
	`n_sites`: number of gene sites
	`Δt`: timestep
	`tracker_end`: Dict() tracking the terminations
	`n_inits`: counter of the number of attempted initiations
"""
function step(
	α, β, β2, γ, L, kon, koff, gene, n_sites, Δt, tracker_end, n_inits
)

	# update tracker for RNAs present on the second strand
	tracker_end["current"][tracker_end["current"].>0] .+= 1

	occupied_sites = findall(gene .== 1)
	occupied_sites = occupied_sites[occupied_sites .> 1]
	process_order = occupied_sites[randperm(length(occupied_sites))]

	finish_flag = false

	for j in process_order
		# detachment can occur at any point in the termination zone
		if (j <= length(gene)) & (j >= n_sites+1)

			if (γ === nothing) || (rand(Bernoulli(γ*Δt))) # detach
				gene[j] = 0

				push!(tracker_end["terminated"], tracker_end["current"][j-n_sites])
				tracker_end["current"][j-n_sites] = 0
				
			elseif (j<=length(gene)-L) && (gene[j+L]==0) && (rand(Bernoulli(β2*Δt))) # advance

				@assert all(gene[j+1:j+L] .== 0)

				gene[j+1] = 1
				gene[j] = 0

				tracker_end["current"][j-n_sites+1] = tracker_end["current"][j-n_sites]
				tracker_end["current"][j-n_sites] = 0

			end
			
		end
		
		# in the gene body
		if (j<=n_sites) & (j>=2)
			if (gene[j+L]==0.) & (rand(Bernoulli(β*Δt)))

				@assert all(gene[j+1:j+L] .== 0)

				gene[j+1]=1.
				gene[j]=0.

				if j == n_sites
					finish_flag = true
					tracker_end["current"][1]=1
				end
				
			end
		end

	end

	# initiation
	if (gene[1]==1.)
		s = wsample(["off", "init", "nothing"], [koff*Δt, α*Δt, 1-Δt*(koff+α)])

		if s=="init"
			n_inits +=1
		end

		if s=="off"
			gene[1]=0.
		elseif (s=="init") & (gene[2]==0)
			gene[2]=1
			gene[1]=0
		end
	end


	# entrance of new RNAp
	if (rand(Bernoulli(kon*Δt))) & all(gene[1:L] .== 0.)
		gene[1]=1.
	end


	return gene, finish_flag, tracker_end, n_inits
	
end

"""
	get_t_d(γ, n_end_sites, β2)

Compute the characteristic time for detachment after reaching the termination point, from
the termination rate `γ`, the number of termination sites `n_end_sites` and the elongation rate
in the termination region `β2`.
"""
get_t_d(γ, n_end_sites, β2) = 1/γ - (n_end_sites)/β2

get_t_d(params) = get_t_d(params.γ, params.n_end_sites, params.β2)

"""
	get_trans_rate(exits, β, Δt)

Compute the transcription rate from the number of termination `exits`, the elongation rate `β`
and the timestep `Δt`.
"""
get_trans_rate(exits, β, Δt) = mean(exits)/β/Δt

get_trans_rate(exits, params) = get_trans_rate(exits, params.β, params.Δt)

"""
	set_Δt(α, β, β2, kon, koff, γ)

Set the Δt for the simulation, from the initiation rate `α`, the elongation rate `β`, 
the elongation rate in the termination region `β2`, the binding rate `kon`, the unbinding
rate `koff` and the termination rate `γ`.
"""
function set_Δt(α, β, β2, kon, koff, γ)

	if γ === nothing
		return .3 / max(α, β, β2, kon, koff)
	else
		return .3 / max(α, β, β2, kon, koff, γ)
	end
end

"""
	get_total_occupancy(density, params; start_bp=1, end_bp=nothing)

Compute the total occupancy on the gene, with `density` representing the occupied sites on the gene, 
`params` containing the parameters for the current simulation, `start_bp` representing the site after the promoter region
and `end_bp` the termination site.
"""
function get_total_occupancy(density, params; start_bp=1, end_bp=nothing)

	if end_bp === nothing
		end_bp = params.n_sites
	end
	n_steps = params.n_steps
	L = params.L

	total_occupancy = sum(density[start_bp:end_bp])/(n_steps)*L

	return total_occupancy
end

# get_total_occupancyδ(α_vec, p, densities, params) = [
# 	get_total_occupancy(densities[p][α], params[p][1]) for α in α_vec
# ]
