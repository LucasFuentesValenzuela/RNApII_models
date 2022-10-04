using Distributions

"""
Define a parameter object
"""
struct Params
	α::Float32
	β::Float32
	γ::Float32
	L::Int
	Δt::Float32
	n_steps::Int
	n_sites::Int
	n_end_sites::Int
	β2::Float32
end

Params(
	α, β, γ, Δt, n_steps, n_sites, n_end_sites, β2
) = Params(
	α, β, γ, Δt, 1, n_steps, n_sites, n_end_sites, β2
)

"""
Runs the walker for a specified number of steps.
"""
function run_walker(
	α, β, γ, L, Δt, n_steps, n_sites, n_end_sites; 
	β2=nothing, model="continuous_detachment", ratio_β2=5
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
			α, β, β2, γ, L, gene, n_sites, Δt, tracker_end; 
			model=model
		)
		density += gene

		exits[k] = finish_flag

	end

	return exits, density, gene, tracker_end
	
end

run_walker(params, model) = run_walker(
	params.α, params.β, params.γ, params.L, params.Δt, params.n_steps, params.n_sites, 
	params.n_end_sites; β2=params.β2, model=model
)

"""
Takes a step forward, for a given model. The kwarg `model` defines the behavior of the RNAs on the second strand. 


Notes: 
* Right now the function is run linearly. But for a real MC approach you would have to pick the location at random, for a sufficient number of times that it actually represents one real time step across the entire strand. 
* I have to check the Poisson process document again to make sure the probabilities are properly encoded. 
* The RNAs are indicated by the left-most site they occupy
"""
function step(
	α, β, β2, γ, L, gene, n_sites, Δt, tracker_end; 
	model="continuous_detachment"
)

	# update tracker for RNAs present on the second strand
	tracker_end["current"][tracker_end["current"].>0] .+= 1
	
	# iterate on the termination strand
	if model == "continuous_detachment"
		
		# detachment can occur at any point
		for j in length(gene):-1:n_sites+1

			p_detach = γ*Δt
	
			if (gene[j]==1) & (rand(Bernoulli(p_detach))) # detach
				gene[j] = 0

				push!(tracker_end["terminated"], tracker_end["current"][j-n_sites])
				tracker_end["current"][j-n_sites] = 0
				
			elseif (j<=length(gene)-L) && (gene[j]==1) && (gene[j+L]==0) &&(rand(Bernoulli(β2*Δt))) # advance

				@assert all(gene[j+1:j+L] .== 0)

				gene[j+1] = 1
				gene[j] = 0

				tracker_end["current"][j-n_sites+1] = tracker_end["current"][j-n_sites]
				tracker_end["current"][j-n_sites] = 0

			end
			
		end
		
	end

	
	finish_flag = false
	
	for j in n_sites:-1:1
		if (gene[j]==1.) & (gene[j+L]==0.) & (rand(Bernoulli(β*Δt)))

			@assert all(gene[j+1:j+L] .== 0)

			gene[j+1]=1.
			gene[j]=0.

			if j == n_sites
				finish_flag = true
				tracker_end["current"][1]=1
			end
			
		end
	end

	# entrance of new RNAp
	if (rand(Bernoulli(α*Δt))) & all(gene[1:L]==0.)
		gene[1]=1.
	end

	return gene, finish_flag, tracker_end
	
end

get_t_d(δ, n_end_sites, β2) = 1/γ - (n_end_sites)/β2

get_t_d(params) = get_t_d(params.γ, params.n_end_sites, params.β2)

get_trans_rate(exits, β, Δt) = mean(exits)/β/Δt

get_trans_rate(exits, params) = get_trans_rate(exits, params.β, params.Δt)
