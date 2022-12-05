### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 645482f3-eca3-4a0b-966b-640a502a9051
begin
	using Revise
	using Pkg
	# Pkg.develop("RNApIIModels")
	Pkg.activate()
	using RNApIIModels
end

# ╔═╡ 3d756f7c-24d9-49fd-9e46-6bcee639ecf9
using Plots

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

# ╔═╡ 7fe13856-cf5c-4870-8fbc-aa1389f4a66d
xx = 2

# ╔═╡ 06fc6fe5-d724-4131-a651-cc0b2e948a2e
xx == 2 ? y = true : y = false

# ╔═╡ a3b6fe9a-e566-424e-8307-b0e658579a07
y

# ╔═╡ 5b0f57dd-af58-446f-8aac-0fbec9ba3217
md"""# Load"""

# ╔═╡ 6aa4adfc-beb4-4615-8ae2-1d0eff4684ce
begin
	α = LITERATURE_PARAMS["min_α"]
	β = LITERATURE_PARAMS["max_β"] / OCCUPANCY_PARAMS["δ"]
end;

# ╔═╡ e6611b5e-0c16-4e31-9724-fe726bca1943
begin
	n_kon_pts = 5
	k_on_vec = 10 .^(LinRange(
		-1, 1.3, n_kon_pts
	));
	
	params_iter = collect(
		Iterators.product([α], [β], [k_on_vec])
	)
end;

# ╔═╡ 1ce46953-b363-4d4c-a235-af379cf285d2
begin
	SIM_PARAMS = copy(OCCUPANCY_PARAMS)
	SIM_PARAMS["n_events"] = 1e5
end

# ╔═╡ 1ed4e2b1-4c8d-491c-b2fc-ee7005627aa8
begin
	occupancy, promoter_occ, params_occ = run_occupancy_simulation(
		params_iter, SIM_PARAMS
	)
end;

# ╔═╡ 2f71e54d-e676-4713-afd6-e49ed5f8736d
occ_sim = max.(occupancy[1], 1e-5)

# ╔═╡ 6a6b01b9-dc76-46df-9d1d-2687125b4781
begin
	
	koff = 1/OCCUPANCY_PARAMS["Ω"] - α
	# @show koff
	α_eff = effective_α.(k_on_vec, koff, α)

	ρ_th = map(f -> (ρ.(f, β, OCCUPANCY_PARAMS["γ"], OCCUPANCY_PARAMS["L"]))[2], α_eff)

	occ_th = ρ_th .* OCCUPANCY_PARAMS["n_sites"]

end;

# ╔═╡ 3f858e05-5067-45fd-b5ea-77f726d5fe33
begin
	plot(k_on_vec, occ_th, xscale=:log10, label="", color=:firebrick)
	scatter!(k_on_vec, occ_th, xscale=:log10, label="Theory", color=:firebrick)

	plot!(k_on_vec, occ_sim, xscale=:log10, label="", color=:forestgreen)
	scatter!(k_on_vec, occ_sim, xscale=:log10, label="Simulations", color=:forestgreen)

	plot!(xlabel="kon", ylabel="occupancy", title="absolute densities")
	plot!(yscale=:log10)
	plot!(legend=:bottomright)
end

# ╔═╡ Cell order:
# ╠═7fa20536-710b-11ed-1d6e-d5ea9ba16f39
# ╠═63606801-a5f5-472b-80b0-fccc2287c634
# ╠═7fe13856-cf5c-4870-8fbc-aa1389f4a66d
# ╠═06fc6fe5-d724-4131-a651-cc0b2e948a2e
# ╠═a3b6fe9a-e566-424e-8307-b0e658579a07
# ╟─5b0f57dd-af58-446f-8aac-0fbec9ba3217
# ╠═3d756f7c-24d9-49fd-9e46-6bcee639ecf9
# ╠═645482f3-eca3-4a0b-966b-640a502a9051
# ╠═6aa4adfc-beb4-4615-8ae2-1d0eff4684ce
# ╠═e6611b5e-0c16-4e31-9724-fe726bca1943
# ╠═1ce46953-b363-4d4c-a235-af379cf285d2
# ╠═1ed4e2b1-4c8d-491c-b2fc-ee7005627aa8
# ╠═2f71e54d-e676-4713-afd6-e49ed5f8736d
# ╠═6a6b01b9-dc76-46df-9d1d-2687125b4781
# ╠═3f858e05-5067-45fd-b5ea-77f726d5fe33
