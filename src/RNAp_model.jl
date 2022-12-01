"""
Runs the walker for a specified number of steps.
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
	
	for k in 1:n_steps # number of timesteps
	
		gene, finish_flag, tracker_end = step(
			α, β, β2, γ, L, kon, koff, gene, n_sites, Δt, tracker_end
		)
		density += gene

		exits[k] = finish_flag

	end

	return exits, density, gene, tracker_end
	
end

run_walker(params) = run_walker(
	params.α, params.β, params.γ, params.L, params.kon, params.koff,
	params.Δt, params.n_steps, params.n_sites, 
	params.n_end_sites; β2=params.β2
)

"""
Takes a step forward, for a given model. The kwarg `model` defines the behavior of the RNAs on the second strand. 


Notes: 
* Right now the function is run linearly. But for a real MC approach you would have to pick the location at random, for a sufficient number of times that it actually represents one real time step across the entire strand. 
* I have to check the Poisson process document again to make sure the probabilities are properly encoded. 
* The RNAs are indicated by the left-most site they occupy
"""
function step(
	α, β, β2, γ, L, kon, koff, gene, n_sites, Δt, tracker_end
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

			if (rand(Bernoulli(γ*Δt))) # detach
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
		# fix the thing with L > 1
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

	# fix the thing with L > 1
	if (gene[1]==1.)
		s = wsample(["off", "init", "nothing"], [koff*Δt, α*Δt, 1-Δt*(koff+α)])
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


	return gene, finish_flag, tracker_end
	
end

"""
"""
get_t_d(γ, n_end_sites, β2) = 1/γ - (n_end_sites)/β2

get_t_d(params) = get_t_d(params.γ, params.n_end_sites, params.β2)

"""Compute the transcription rate.
"""
get_trans_rate(exits, β, Δt) = mean(exits)/β/Δt

get_trans_rate(exits, params) = get_trans_rate(exits, params.β, params.Δt)

"""Set the Δt for the simulation.
"""
function set_Δt(α, β, β2, kon, koff, γ)

	if γ === nothing
		return .3 / max(α, β, β2, kon, koff)
	else
		return .3 / max(α, β, β2, kon, koff, γ)
	end
end

"""Get the total occupancy of the gene
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

get_total_occupancyδ(α_vec, p, densities, params) = [
	get_total_occupancy(densities[p][α], params[p][1]) for α in α_vec
]
