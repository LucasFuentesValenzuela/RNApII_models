using Distributions
using ProgressBars
using StatsBase

include("parameters.jl")

"""
Define a parameter object
"""
struct Params
	α::Float32 # initiation rate from promoter to gene body
	β::Float32 # elongation rate
	γ::Float32 # termination rate
	L::Int # footprint of the body
	kon::Float32 # rate of binding to the promoter
	koff::Float32 # rate of unbinding to the promoter
	Δt::Float32 # timestep
	n_steps::Int # number of steps in the simulation
	n_sites::Int # number of sites on the gene
	n_end_sites::Int # number of sites where termination is possible
	β2::Float32 # elongation rate on the termination sites
end

Params(
	α, β, γ, Δt, n_steps, n_sites, n_end_sites, β2
) = Params(
	α, β, γ, 1, Δt, n_steps, n_sites, n_end_sites, β2
)


DEFAULT_PARAMS = Params(
	α_default, β_default, γ_default, L_default, kon_default, koff_default, Δt_default, 
	DEFAULT_nsteps, DEFAULT_n_sites, DEFAULT_n_end_sites, β2_default
)

"""
Runs the walker for a specified number of steps.
"""
function run_walker(
	α, β, γ, L, kon, koff, Δt, n_steps, n_sites, n_end_sites; 
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
			α, β, β2, γ, L, kon, koff, gene, n_sites, Δt, tracker_end; 
			model=model
		)
		density += gene

		exits[k] = finish_flag

	end

	return exits, density, gene, tracker_end
	
end

run_walker(params, model) = run_walker(
	params.α, params.β, params.γ, params.L, params.kon, params.koff,
	params.Δt, params.n_steps, params.n_sites, 
	params.n_end_sites; β2=params.β2, model=model
)

run_walker(params) = run_walker(params, "continuous_detachment")

"""
Takes a step forward, for a given model. The kwarg `model` defines the behavior of the RNAs on the second strand. 


Notes: 
* Right now the function is run linearly. But for a real MC approach you would have to pick the location at random, for a sufficient number of times that it actually represents one real time step across the entire strand. 
* I have to check the Poisson process document again to make sure the probabilities are properly encoded. 
* The RNAs are indicated by the left-most site they occupy
"""
function step(
	α, β, β2, γ, L, kon, koff, gene, n_sites, Δt, tracker_end; 
	model="continuous_detachment"
)

	# update tracker for RNAs present on the second strand
	tracker_end["current"][tracker_end["current"].>0] .+= 1

	# iterate on the termination strand
	if model == "continuous_detachment"
		
		# detachment can occur at any point in the termination zone
		for j in length(gene):-1:n_sites+1

			if (gene[j]==1) & (rand(Bernoulli(γ*Δt))) # detach
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
	
	for j in n_sites:-1:2
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

	# fix the thing with L > 1
	if (gene[1]==1.)
		s = wsample(["off", "init", "nothing"], [koff*Δt, α*Δt, 1-Δt*(koff+α)])
		if s=="off"
			gene[1]=0.
		elseif s=="init"
			gene[2]=1
			gene[1]=0
		end
	end

	# entrance of new RNAp
	if (rand(Bernoulli(kon*Δt))) & all(gene[1:L].==0.)
		gene[1]=1.
	end

	return gene, finish_flag, tracker_end
	
end

get_t_d(γ, n_end_sites, β2) = 1/γ - (n_end_sites)/β2

get_t_d(params) = get_t_d(params.γ, params.n_end_sites, params.β2)

get_trans_rate(exits, β, Δt) = mean(exits)/β/Δt

get_trans_rate(exits, params) = get_trans_rate(exits, params.β, params.Δt)

function sweep_params(α_vec, p_vec, param_name; params=DEFAULT_PARAMS)

	params_dict = Dict()
	trans_rates = Dict()
	residence_times = Dict() # in number of steps
	densities = Dict()

	for p in ProgressBar(p_vec)

		params_dict[p] = []
		trans_rates[p] = []
		residence_times[p] = []
		densities[p] = Dict()
	
		for α in α_vec

			if param_name == "γ"
				β = params.β*params.L
				γ = p 
				L = params.L
				n_steps = params.n_steps
				n_sites = params.n_sites*params.L
				n_end_sites = params.n_end_sites*params.L
				β2 = params.β2*params.L
			elseif param_name == "L"
				β = params.β*p
				γ = params.γ
				L = p 
				n_steps = params.n_steps
				n_sites = params.n_sites*p
				n_end_sites = params.n_end_sites*p
				β2 = params.β2*p
			end
			
			Δt = 1/3 * minimum([1/α, 1/β, 1/γ])

			params_crt = Params(
				α, β, γ, L, Δt, n_steps, n_sites, n_end_sites, β2
			)

			exits_, density_, _, tracker_ = run_walker(params_crt);

			push!(params_dict[p], params_crt)
			push!(trans_rates[p], get_trans_rate(exits_, params_crt))
			push!(residence_times[p], mean(tracker_["terminated"]))
			densities[p][α] = density_
		
		end
		
	end

	return params_dict, trans_rates, residence_times, densities
end