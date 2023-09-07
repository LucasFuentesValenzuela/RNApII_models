"""
Series of plotting functions added over the course of the project. 
They were mostly used in the earlier parts of the project to debug and build intuition.
"""

"""
"""
function plot_tracker_end(tracker_end, params)
	histogram(tracker_end["terminated"], label="Data", normed=true)
	vline!([mean(tracker_end["terminated"])], linewidth=3., label="Mean")
	vline!([1/params.γ/params.Δt], linewidth=3., linestyle=:dash, label="Expected mean")

	xlabel!("Number of steps taken by RNAs")
	ylabel!("Count")
end

"""
"""
function plot_density(density, n_steps, n_sites, L; normalize=true)

	p = plot(density/n_steps*L, label="")
	xlabel!("basepair")
	ylabel!("Steady-state density")
	vline!([n_sites], linestyle=:dash, label="")
    vline!([1], linestyle=:dash, label="")
	if normalize
		ylims!(0, 1)
	end

	return p
end
plot_density(density, params) = plot_density(density, params.n_steps, params.n_sites, params.L)

# ======================= Plotting for param sweep ================================

"""
"""
function plot_transcription_rate_sweep(α_vec, p_vec, params, param_name, trans_rates)

	color_palette = palette([:blue, :green], length(p_vec)+1)

	p1 = plot()

	# Analytical solution
	if param_name=="γ"

		β = params[p_vec[1]][1].β
		γ = LARGE_γ
		L = params[p_vec[1]][1].L

		plot!(
			α_vec, 
			J.(α_vec, β, γ, L), 
			label="Theory: γ >> 1, L=$(L)", linewidth=2,
			color=:firebrick
		)

	elseif param_name=="L"

		# parameters for the small value of L
		β = params[p_vec[1]][1].β
		γ = params[p_vec[1]][1].γ
		L = p_vec[1]

		plot!(
			α_vec, 
			J.(α_vec, β, γ, L), 
			label="Theory: L=$(L), γ = $(γ)", linewidth=2,
			color=:firebrick
		)

		# parameters for the small value of L
		βL = params[p_vec[end]][1].β
		γL = params[p_vec[end]][1].γ
		LL = p_vec[end]

		plot!(
			α_vec, 
			J.(α_vec, βL, γL, LL), 
			label="Theory:  L=$(LL),  γ = $(γL)", linewidth=2,
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
		plot!(
			α_vec, trans_rates[p]*params[p][1].β, 
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
	plot!(legend=false)

	return p1
end

plot_transcription_rate_sweep(results, param_name) = plot_transcription_rate_sweep(
	results["α_vec"], results["p_vec"], results["params_dict"], param_name, results["trans_rates"]
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

		# tmp fix: densities need to be rescaled
		density_crt = copy(densities[p_plot][α])
		L_crt = params[p_plot][1].L
		for j in 1:length(densities[p_plot][α])-L_crt
			for l in 1:L_crt-1
				density_crt[j+l] += densities[p_plot][α][j]
			end
		end

		plot!(
			LinRange(0, 1, params[p_plot][1].n_sites + params[p_plot][1].n_end_sites),
			density_crt/params[p_plot][1].n_steps, 
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
	expected = [1/γ/params[p][k].Δt for k in 1:length(params[p])]

	plot!(α_vec, residence_times[p], label="", linestyle=:dash)
	scatter!(α_vec, residence_times[p], label="Data")
	plot!(α_vec, expected,
		label="Expected", linestyle=:dash, linewidth=3, 
	)

	xlabel!("Initiation rate [1/s]")
	ylabel!("Number of steps on second strand")
	plot!(legend=:topleft)
	# ylims!(.6*expected, 1.4*expected)
	plot!(xscale=:log)
	plot!(yscale=:log)

	title!("Residence time for γ = $(round(γ; digits=3))")

	return p1
end

plot_residence_times_sweep(results, p, param_name) =  plot_residence_times_sweep(
	results["α_vec"], p, results["params_dict"], param_name, results["residence_times"]
	)

"""
"""
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
		results, param_name; α_default=α_default, max_fold_change_α=3
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
		
		vline!([α_default], label="α = α_default", color=:red)
		xlabel!("α")
		ylabel!("Occupancy fold change")
		plot!(xscale=:log)
		plot!(yscale=:log)

		return px
	
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