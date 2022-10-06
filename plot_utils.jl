using Plots
using Statistics
using Interpolations

include("theory.jl")
include("model.jl")

LARGE_γ = 1000

function plot_tracker_end(tracker_end, params)
	histogram(tracker_end["terminated"], label="Data", normed=true)
	vline!([mean(tracker_end["terminated"])], linewidth=3., label="Mean")
	vline!([1/params.γ/params.Δt], linewidth=3., linestyle=:dash, label="Expected mean")

	xlabel!("Number of steps taken by RNAs")
	ylabel!("Count")
end

function plot_density(density, n_steps, n_sites)

	p = plot(density/n_steps, label="")
	xlabel!("basepair")
	ylabel!("Steady-state density")
	vline!([n_sites], linestyle=:dash, label="")
    vline!([1], linestyle=:dash, label="")
	ylims!(0, 1)
	# plot!(legend=:outertopleft)

	return p
end

plot_density(density, params) = plot_density(density, params.n_steps, params.n_sites)

# ======================= Plotting for param sweep ================================

"""
"""
function plot_transcription_rate_sweep(α_vec, p_vec, params, param_name, trans_rates, DEFAULT_PARAMS)

	color_palette = palette([:blue, :green], length(p_vec)+1)

	p1 = plot()

	# Analytical solution
	if param_name=="γ"
		plot!(
			α_vec, 
			J.(α_vec, DEFAULT_PARAMS.β, LARGE_γ, params[p_vec[1]][1].L), 
			label="Theory: γ >> 1, L=$(params[p_vec[1]][1].L)", linewidth=2,
			color=:firebrick
		)

	elseif param_name=="L"

		plot!(
			α_vec, 
			J.(α_vec, DEFAULT_PARAMS.β, params[DEFAULT_PARAMS.L][1].γ, DEFAULT_PARAMS.L), 
			label="Theory: L=$(DEFAULT_PARAMS.L), γ = $(params[1][1].γ)", linewidth=2,
			color=:firebrick
		)

		plot!(
			α_vec, 
			J.(α_vec, DEFAULT_PARAMS.β*p_vec[end], params[p_vec[end]][1].γ, p_vec[end]), 
			label="Theory:  L=$(p_vec[end]),  γ = $(params[1][end].γ)", linewidth=2,
			color=:orange
		)

	end

	for (k, p) in enumerate(p_vec)

		if k%trunc(Int, length(p_vec)/4)==0
			label = "$(param_name) = $(round(p; digits=3))"
		else
			label=""
		end

		# the below is assuming that all elements of list params[p] have the same β
		plot!(α_vec, trans_rates[p]*params[p][1].β, 
			label=label, linestyle=:dash, linewidth=2, color=color_palette[k]
		)

		# scatter!(
		# 	α_vec, trans_rates[p]*params[p][1].β, 
		# 	label="", markersize=5, color=color_palette[k]
		# )
	end

	xlabel!("Initiation rate α [1/s]")
	ylabel!("Transcription rate [1/s]")
	plot!(legend=:topleft)
	plot!(xscale=:log)

	return p1
end

plot_transcription_rate_sweep(results, param_name, DEFAULT_PARAMS) = plot_transcription_rate_sweep(
	results["α_vec"], results["p_vec"], results["params_dict"], param_name, results["trans_rates"], 
	DEFAULT_PARAMS
)

"""
test
"""
function plot_density_sweep(α_vec, p_plot, params, param_name, densities)

	color_palette = palette([:blue, :green], length(α_vec))

	p1 = plot()
	
	for (k, α) in enumerate(α_vec)

		if k%trunc(Int, length(α_vec)/5)==0
			label = "α/β=$(round(α/(params[p_plot][1].β); digits=3))"
		else
			label=""
		end

		plot!(
			LinRange(0, 1, params[p_plot][1].n_sites + params[p_plot][1].n_end_sites),
			densities[p_plot][α]/params[p_plot][1].n_steps, 
			label=label, 
			color=color_palette[k], linewidth=2.
		)
		xlabel!("Relative position")
		ylabel!("Steady-state density")
		
	end
	vline!([params[p_plot][1].n_sites/(params[p_plot][1].n_sites + params[p_plot][1].n_end_sites)], linestyle=:dash, label="", linewidth=2)

	plot!(legend=:outertopright)

	title!("Occupancy as a function of α, $(param_name) = $(round(p_plot; digits=3))")

	return p1
end

plot_density_sweep(results, p_plot, param_name) = plot_density_sweep(
	results["α_vec"], p_plot, results["params_dict"], param_name, results["densities"]
	)

"""
"""
function plot_residence_times_sweep(α_vec, p, params, param_name, residence_times)

	p1 = plot()

	if param_name=="γ"
		γ = p
	else
		γ = params[p][1].γ
	end
	expected = 1/γ/params[p][1].Δt

	plot!(α_vec, residence_times[p], label="", linestyle=:dash)
	scatter!(α_vec, residence_times[p], label="Data")
	hline!(
		[expected], 
		label="Expected", linestyle=:dash, linewidth=3, 
	)

	xlabel!("Initiation rate [1/s]")
	ylabel!("Number of steps on second strand")
	plot!(legend=:bottomright)
	ylims!(.6*expected, 1.4*expected)

	title!("Residence time for γ = $(round(γ; digits=3))")

	return p1
end

plot_residence_times_sweep(results, p, param_name) =  plot_residence_times_sweep(
	results["α_vec"], p, results["params_dict"], param_name, results["residence_times"]
	)

"""
"""

function get_total_occupancy(density, params)

	n_sites = params.n_sites
	n_steps = params.n_steps
	L = params.L

	total_occupancy = sum(density[1:n_sites])/(n_steps*n_sites)*L

	return total_occupancy
end

get_total_occupancy(α_vec, p, densities, params) = [
	get_total_occupancy(densities[p][α], params[p][1]) for α in α_vec
]

function plot_occupancy_sweep(α_vec, p, params, param_name, densities)

	p1 = plot()

	total_occupancy = get_total_occupancy(α_vec, p, densities, params)

	plot!(α_vec, total_occupancy, linestyle=:dash, label="")
	scatter!(α_vec, total_occupancy, linestyle=:dash, label="")

	if param_name == "γ"
		vline!([p], label="α = γ")
	else
		vline!([params[p][1].γ], label="α = γ")
	end

	vline!([params[p][1].β], label="α = β")

	xlabel!("Initiation rate α [1/s]")
	ylabel!("Total occupancy")

	plot!(xscale=:log)
	ylims!(-.1, 1.1)
	plot!(legend=:bottomright)

	title!("""Occupancy vs α, $(param_name) = $(round(p; digits=3))""")
	return p1
end

plot_occupancy_sweep(results, p, param_name) = plot_occupancy_sweep(
	results["α_vec"], p, results["params_dict"], param_name, results["densities"]
	)


function plot_occupancy_fold_change(
		results, param_name; DEFAULT_PARAMS=DEFAULT_PARAMS, max_fold_change_α=3
	)
	
		# values of α at which to evaluate the fold change
		α_interp = results["α_vec"][1] .* 10. .^(
			LinRange(
				0, 
				log10(results["α_vec"][end]/max_fold_change_α/results["α_vec"][1]), 
				100
			)
		)
	
		
		px = plot()
	
		color_palette = palette([:orange, :green], length(results["p_vec"]))
		
		for (k,p) in enumerate(results["p_vec"])
			
			occupancy = get_total_occupancy(
				 results["α_vec"], p, results["densities"], results["params_dict"]
			);
			occupancy_interp = linear_interpolation(results["α_vec"], occupancy)
	
			fold_changes = occupancy_interp.(max_fold_change_α*α_interp)./occupancy_interp.(α_interp)
	
			if k%trunc(Int, length(results["p_vec"])/4)==0
				label = "$param_name = $(round(p; digits=3))"
			else
				label=""
			end
			
			plot!(
				α_interp, fold_changes, linestyle=:dash, label=label, 
				color = color_palette[k], linewidth=2
			)
		end
		
		vline!([DEFAULT_PARAMS.α], label="α = α_default", color=:red)
		xlabel!("α")
		ylabel!("Occupancy fold change")
		plot!(xscale=:log)
		plot!(yscale=:log)
	
	end


function plot_occupancy_fold_change_scaling(
		results, p, param_name; max_fold_change_α=3, colors = [:orange, :green]
	)
	
		if param_name == "γ"
			γ_plot = p
		elseif param_name == "L"
			γ_plot = results["params_dict"][p][1].γ
		end
	
	
		lower_ = results["α_vec"][1]/γ_plot * 1.1
		upper_ = results["α_vec"][end]/γ_plot/(max_fold_change_α+.1)
		
		α_interp = γ_plot .* 10. .^(
			collect(
				LinRange(
					log10(lower_), 
					log10(upper_), 
					15
				)
			)
		)
		
		color_palette = palette(colors, length(α_interp))
	
		x = LinRange(1, max_fold_change_α, 30) # scaling factors
	
		occupancy = get_total_occupancy(
			 results["α_vec"], p, results["densities"], results["params_dict"]
		);
	
		occupancy_interp = linear_interpolation(results["α_vec"], occupancy)
	
		fold_changes = [
			occupancy_interp.(x*α)./occupancy_interp.(α) for α in α_interp
		]
	
		pocc_x = plot(size=(800, 400))
		
		for (k, α) in enumerate(α_interp)
			plot!(
				x, fold_changes[k], 
				linestyle=:dash, label="α/γ = $(round(α/γ_plot; digits=3))",
				color=color_palette[k]
			)
		end
		
		xlabel!("α fold change")
		ylabel!("occupancy fold change")
		plot!(legend=:outertopleft)
		ylims!(0, 10)
		title!("Occupancy fold change for $param_name = $(round(p; digits=3))")
	
		pocc_x
	end


"""
"""
function plot_comparison_J_th_sim(results, p, param_name)
	
	pe = plot()

	# load results
	α_vec = results["α_vec"]
	params = results["params_dict"][p][1]
	trans_rates = results["trans_rates"][p]

	# plot
	plot!(
		α_vec, 
		J.(α_vec, params.β, params.γ, params.L), 
		label="Theory: γ = $(params.γ), L=$(params.L)", linewidth=2,
		color=:firebrick
	)

	# the below is assuming that all elements of list params[p] have the same β
	plot!(α_vec, trans_rates*params.β, 
		label="", linestyle=:dash, linewidth=2
	)

	scatter!(
		α_vec, trans_rates*params.β, 
		label="$param_name = $(round(p; digits=3))", markersize=5
	)
	
	plot!(legend=:bottomright)
	xlabel!("α")
	ylabel!("J")

	return pe
end