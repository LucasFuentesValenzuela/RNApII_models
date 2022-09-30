using Plots
using Statistics

function plot_tracker_end(tracker_end, params)
	histogram(tracker_end["terminated"], label="Data", normed=true)
	vline!([mean(tracker_end["terminated"])], linewidth=3., label="Mean")
	vline!([1/params.δ/params.Δt], linewidth=3., linestyle=:dash, label="Expected mean")

	xlabel!("Number of steps taken by RNAs")
	ylabel!("Count")
end

function plot_density(density, n_steps, n_sites)

	p = plot(density/n_steps, label="")
	xlabel!("basepair")
	ylabel!("Steady-state density")
	vline!([n_sites], linestyle=:dash, label="")
    vline!([1], linestyle=:dash, label="")
	# plot!(legend=:outertopleft)

	return p
end

plot_density(density, params) = plot_density(density, params.n_steps, params.n_sites)