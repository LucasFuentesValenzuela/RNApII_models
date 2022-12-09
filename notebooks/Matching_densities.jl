### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 645482f3-eca3-4a0b-966b-640a502a9051
begin
	using Revise
	using Pkg
	Pkg.activate()
	using RNApIIModels
end

# ╔═╡ 3d756f7c-24d9-49fd-9e46-6bcee639ecf9
using Plots

# ╔═╡ 8bb557b4-d3f0-4bca-b96a-4f893bbe94e2
using Statistics

# ╔═╡ 5d691fb8-9dcc-48b1-8d98-cbaf85bc086c
using PlutoUI

# ╔═╡ 7fa20536-710b-11ed-1d6e-d5ea9ba16f39
md"""
# Description 

Here we validate that the theoretical model can be used for at least some parameters to make predictions quickly about the results of simulations. 
"""

# ╔═╡ 63606801-a5f5-472b-80b0-fccc2287c634
md"""
TODO: 
- Check why there is discrepancy for larger values of k_on: In simulations for the feasibility study we see that for large values of k_on there is a discrepancy between simulations and theory. Why is that? It seems that the discrepancy starts to appear around kon = 1e-1. Let's see. I think the plot we want is error as a function of kon and alpha.  I think it is because I implement the number of iterations based on the number of *binding* events. While I should implement it as a function of the number of *initiation* events. I think that makes more sense. Why? Right now the number of iterations is n_events/kon. But if we are saturating the promoter, we are decreasing the number of iterations without changing the effective initiation rate, really. So I think we want at least the same number of steps as before, not fewer. And right now, increasing kon changes the number of iterations by the exact same amount!

- I thought this was solved --- but apparently this is not the case in the results - why? 


- How does the result change with δ?
"""

# ╔═╡ bcb3969d-271a-4be7-aa4c-9d08b8d60ec8
TableOfContents()

# ╔═╡ 5b0f57dd-af58-446f-8aac-0fbec9ba3217
md"""# Load"""

# ╔═╡ 2bfcb6de-c716-4bee-8d6d-c054ddc37bb8
TableOfContents()

# ╔═╡ c27f81e6-eb47-4bb3-9248-962e22d278fe
md"""# Define parameters"""

# ╔═╡ 6aa4adfc-beb4-4615-8ae2-1d0eff4684ce
begin
	α = 1/2 # LITERATURE_PARAMS["min_α"]
	β = 1# LITERATURE_PARAMS["max_β"] / OCCUPANCY_PARAMS["δ"]
	Ω = 2 # OCCUPANCY_PARAMS["Ω"]
	γ = nothing
	L = OCCUPANCY_PARAMS["L"]

	SIM_PARAMS = copy(OCCUPANCY_PARAMS)
	SIM_PARAMS["n_events"] = 1e4

	ntimes = 5


	n_kon_pts = 8
	k_on_vec = 10 .^(LinRange(
		-2, 1.5, n_kon_pts
	));
	
	params_iter = collect(
		Iterators.product([α], [β], [k_on_vec])
	)
end;

# ╔═╡ 5e80ac50-88df-4a32-b0f4-a2fdc1a8c509
koff = 1/Ω - α

# ╔═╡ d749ba8e-8267-4545-9cc4-0c0cd8b544a0
# WHY is it taking so long when I decrease α? to me it should take less long? ... 
# oooh not if you define a number of initiation events!!!

# ╔═╡ 1ed4e2b1-4c8d-491c-b2fc-ee7005627aa8
begin
	occupancy = []
	promoter_occ = []
	params_occ = []
	α_eff_sim = []
	for k in 1:ntimes
		occupancy_crt, promoter_occ_crt, params_occ_crt, α_eff_crt = run_occupancy_simulation(
			params_iter, SIM_PARAMS
		)
		push!(occupancy, occupancy_crt)
		push!(promoter_occ, promoter_occ_crt)
		push!(params_occ, params_occ_crt)
		push!(α_eff_sim, α_eff_crt)
	end
end;

# ╔═╡ 2f71e54d-e676-4713-afd6-e49ed5f8736d
begin
	occ_sim = get_quantile(occupancy, .5);
	occ_q_up = get_quantile(occupancy, .9)
	occ_q_dwn = get_quantile(occupancy, .1)
	
	α_eff_sim_median = get_quantile(α_eff_sim, .5)
	α_eff_q_up = get_quantile(α_eff_sim, .9)
	α_eff_q_dwn = get_quantile(α_eff_sim, .1)

	prom_occ_sim = get_quantile(promoter_occ, .5)
	prom_occ_q_up = get_quantile(promoter_occ, .9)
	prom_occ_q_dwn = get_quantile(promoter_occ, .1)
end;

# ╔═╡ 6a6b01b9-dc76-46df-9d1d-2687125b4781
begin
	
	# @show koff
	α_eff = effective_α.(k_on_vec, koff, α)

	ρ_th = map(f -> (ρ.(f, β, γ, L))[2], α_eff)

	occ_th = ρ_th .* OCCUPANCY_PARAMS["n_sites"]

	prom_occ_th = ρp.(k_on_vec, koff, α)

	prom_occ_MC = ρp_MC.(α, β, k_on_vec, koff; order=5)

	occ_MC = map(f -> (ρ.(f, β, γ, L))[2], α * prom_occ_MC) .* OCCUPANCY_PARAMS["n_sites"]

end;

# ╔═╡ 6fe14578-5c4c-49e9-9020-97f8ac53c010
md"""
The above hints at the fact that the promoter saturation pattern might not quite match the theoretical one in some instances. However we have checked that it did. Why would that not be the case: 
- β is too low, creating a jam for very high values of the kon. 
- Equivalently, the residence time is too short
"""

# ╔═╡ 9b57f836-d6dc-4295-a75f-48b74d402ffc
begin
	α_eff_prom = α .* prom_occ_sim[1]
	occ_from_prom_occ = map(
		f -> (
			ρ.(f, β, OCCUPANCY_PARAMS["γ"], OCCUPANCY_PARAMS["L"]))[2] .* OCCUPANCY_PARAMS["n_sites"] , α_eff_prom
	)

	occ_from_α_eff_sim = map(
		f -> (
			ρ.(f, β, OCCUPANCY_PARAMS["γ"], OCCUPANCY_PARAMS["L"]))[2] .* OCCUPANCY_PARAMS["n_sites"] , α_eff_sim_median[1]
	)
end;

# ╔═╡ 3f858e05-5067-45fd-b5ea-77f726d5fe33
let
	p_gene = plot(
		k_on_vec, occ_th, 
		xscale=:log10, label="", 
		color=:firebrick
	)
	scatter!(
		k_on_vec, occ_th, 
		xscale=:log10, 
		label="Theory", 
		color=:firebrick
	)

	w = (occ_q_up .- occ_q_dwn)
	plot!(
		k_on_vec, (occ_q_up .+ occ_q_dwn)./2, 
		ribbon = w[1], label="", fillalpha=.3
	)
	scatter!(
		k_on_vec, occ_sim, 
		xscale=:log10, 
		label="Simulations", 
		color=:forestgreen
	)
	plot!(
		k_on_vec, occ_MC, 
		xscale=:log10, label="", 
		color=:darkorange
	)
	scatter!(
		k_on_vec, occ_MC, 
		xscale=:log10, 
		label="MC", 
		color=:darkorange
	)
	
	plot!(
		k_on_vec, occ_from_prom_occ, 
		linestyle=:dash, linewidth=2, label="α_eff from prom occupancy"
	)
	plot!(
		k_on_vec, occ_from_α_eff_sim,
		linestyle=:dash, linewidth=2, label="α_eff from simulations"
	)

	plot!(xlabel="kon", ylabel="occupancy", title="Gene body")
	# plot!(yscale=:log10)
	plot!(legend=:bottomright)
	vline!([α], label="α", linewidth=2, color=:red)


	p_prom = plot(
		k_on_vec, prom_occ_th, 
		xscale=:log10, label="", 
		color=:firebrick
	)
	scatter!(
		k_on_vec, prom_occ_th, 
		xscale=:log10, 
		label="Theory", 
		color=:firebrick
	)

	w = (prom_occ_q_up .- prom_occ_q_dwn)
	plot!(
		k_on_vec, (prom_occ_q_up .+ prom_occ_q_dwn)./2, 
		ribbon = w[1], label="", fillalpha=.3
	)
	scatter!(
		k_on_vec, prom_occ_sim, 
		xscale=:log10, 
		label="Simulations", 
		color=:forestgreen
	)

	plot!(
		k_on_vec, prom_occ_MC, 
		xscale=:log10, label="", 
		color=:darkorange
	)
	scatter!(
		k_on_vec, prom_occ_MC, 
		xscale=:log10, 
		label="MC", 
		color=:darkorange
	)

	plot!(xlabel="kon", ylabel="occupancy", title="promoter")
	# plot!(yscale=:log10)
	plot!(legend=:bottomright)
	vline!([α], label="α", linewidth=2, color=:red)

	p_α_eff = plot(
		k_on_vec, α_eff_sim_median, label=""
	)
	scatter!(k_on_vec, α_eff_sim_median, ribbon = (α_eff_q_up .- α_eff_q_dwn)[1], label="simulation")

	plot!(
		k_on_vec, effective_α.(k_on_vec, koff, α), 
		linestyle=:dash, label="theory"
	)

	plot!(
		k_on_vec, α_eff_prom, label="α_from_promoter", linewidth=2, color=:blue
	)
	plot!(
		k_on_vec, α * prom_occ_MC, label="MC", color=:orange
	)

	plot!(xscale=:log10)
	plot!(xlabel="kon", ylabel="α_eff")
	plot!(legend=:topleft)
	hline!([β/2], label="β/2", linestyle=:dash, linewidth=2, color=:firebrick)
	vline!([α], label="α", linewidth=2, color=:red)
	plot!(legend=:bottomright)


	l = @layout [
grid(3, 1)
	]
	plot([p_gene, p_prom, p_α_eff]..., layout=l, size=(500, 800))
end

# ╔═╡ 2f92f886-10ae-4e58-ab24-c4f5cb8eea9c
plot(
	k_on_vec, α_eff_sim_median[1] ./ α_eff_prom
)

# ╔═╡ 21617063-2d2b-4400-871e-552d500c5750
begin

end

# ╔═╡ 858d3a09-dfb2-41f6-b0ac-bab387f6c147
md"""
We see that the promoter occupancy is over-estimated from the simulations. Therefore, the effective α is also overestimated, which can only increase the occupancy. This could probably explain **part** of the discrepancy? 

- It does seem that part of it is explained by changing the rule of promoter occupancy. it might be because the law is not right, or breaks down in some regimes
- 
"""

# ╔═╡ c32fa3dd-0a49-4973-a32a-5b968eb9d5ca
begin
	plot(k_on_vec, α_eff ./ (β/ 1 + sqrt(L)), label="", xlabel="kon", ylabel="α_eff/β")
	plot!(xscale=:log10)
end

# ╔═╡ 2758e8ff-03fd-446f-a53c-e222c274be7e
md"""
# Impact of the n_event parameter

It seems to play a significant role (unsurprisingly). But do we get a saturation effect? let's investigate for a single value of kon
"""

# ╔═╡ 60631acc-623d-4324-89ea-68d567b6f93d
# ╠═╡ disabled = true
#=╠═╡
let
	
	kon = 50
	n_events_vec = [1e2, 5e2, 1e3, 5e3]
	ntimes = 5

	α_eff = effective_α(kon, koff, α)
	ρ_th = ρ.(α_eff, β, OCCUPANCY_PARAMS["γ"], OCCUPANCY_PARAMS["L"])[2]
	occ_th = ρ_th .* OCCUPANCY_PARAMS["n_sites"]

	occ_sim = []
	occ_std = []
	for nevents in n_events_vec
		
		SIM_PARAMS = copy(OCCUPANCY_PARAMS)
		SIM_PARAMS["n_events"] = nevents

		occupancy, promoter_occ, params_occ = run_occupancy_simulation(
			[(α, β, kon .* ones(ntimes))], SIM_PARAMS
		)

		push!(occ_sim, median(occupancy[1]))
		push!(occ_std, std(occupancy[1]))

	end
	

	# @show occ_sim, occ_th
	q = plot()
	plot!(n_events_vec, (occ_sim .- occ_th) ./ occ_th, ribbon=occ_std./occ_th, label="", fillalpha=.3)
	scatter!(n_events_vec, (occ_sim .- occ_th) ./ occ_th, label="")
	plot!(xlabel="nevents", ylabel="error")
	hline!([0], linestyle=:dash, label="", linewidth=2)
	plot!(xscale=:log10)

	q
	
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═7fa20536-710b-11ed-1d6e-d5ea9ba16f39
# ╠═63606801-a5f5-472b-80b0-fccc2287c634
# ╠═bcb3969d-271a-4be7-aa4c-9d08b8d60ec8
# ╟─5b0f57dd-af58-446f-8aac-0fbec9ba3217
# ╠═3d756f7c-24d9-49fd-9e46-6bcee639ecf9
# ╠═645482f3-eca3-4a0b-966b-640a502a9051
# ╠═8bb557b4-d3f0-4bca-b96a-4f893bbe94e2
# ╠═5d691fb8-9dcc-48b1-8d98-cbaf85bc086c
# ╠═2bfcb6de-c716-4bee-8d6d-c054ddc37bb8
# ╟─c27f81e6-eb47-4bb3-9248-962e22d278fe
# ╠═6aa4adfc-beb4-4615-8ae2-1d0eff4684ce
# ╠═5e80ac50-88df-4a32-b0f4-a2fdc1a8c509
# ╠═d749ba8e-8267-4545-9cc4-0c0cd8b544a0
# ╠═1ed4e2b1-4c8d-491c-b2fc-ee7005627aa8
# ╠═2f71e54d-e676-4713-afd6-e49ed5f8736d
# ╠═6a6b01b9-dc76-46df-9d1d-2687125b4781
# ╠═3f858e05-5067-45fd-b5ea-77f726d5fe33
# ╟─6fe14578-5c4c-49e9-9020-97f8ac53c010
# ╠═9b57f836-d6dc-4295-a75f-48b74d402ffc
# ╠═2f92f886-10ae-4e58-ab24-c4f5cb8eea9c
# ╠═21617063-2d2b-4400-871e-552d500c5750
# ╠═858d3a09-dfb2-41f6-b0ac-bab387f6c147
# ╠═c32fa3dd-0a49-4973-a32a-5b968eb9d5ca
# ╟─2758e8ff-03fd-446f-a53c-e222c274be7e
# ╠═60631acc-623d-4324-89ea-68d567b6f93d
