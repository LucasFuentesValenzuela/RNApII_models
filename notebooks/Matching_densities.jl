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
Here we validate that the theoretical model can be used for at least some parameters to make predictions quickly about the results of simulations. 
"""

# ╔═╡ 6aa4adfc-beb4-4615-8ae2-1d0eff4684ce
begin
	α = LITERATURE_PARAMS["min_α"]
	β = LITERATURE_PARAMS["max_β"] / OCCUPANCY_PARAMS["δ"]
end;

# ╔═╡ e6611b5e-0c16-4e31-9724-fe726bca1943
begin
	n_kon_pts = 10
	k_on_vec = 10 .^(LinRange(
		-2, 0, n_kon_pts
	));
	
	params_iter = collect(
		Iterators.product([α], [β], [k_on_vec])
	)
end;

# ╔═╡ 1ed4e2b1-4c8d-491c-b2fc-ee7005627aa8
begin
	occupancy, promoter_occ, params_occ = run_occupancy_simulation(
		params_iter, OCCUPANCY_PARAMS
	)
end;

# ╔═╡ 4612767c-6317-4453-80f7-b4200bfdae03
begin
	plot(k_on_vec, occupancy[1], xscale=:log10, label="")
	scatter!(k_on_vec, occupancy[1], xscale=:log10, label="")

	plot!(xlabel="kon", ylabel="total occupancy")
end

# ╔═╡ 6a6b01b9-dc76-46df-9d1d-2687125b4781
begin
	
	koff = 1/OCCUPANCY_PARAMS["Ω"] - α

	α_eff = effective_α.(k_on_vec, koff, α)

	ρ_th = map(f -> (ρ.(f, β, OCCUPANCY_PARAMS["γ"], OCCUPANCY_PARAMS["L"]))[2], α_eff)

	occ_th = ρ_th .* OCCUPANCY_PARAMS["n_sites"]

end;

# ╔═╡ 3f858e05-5067-45fd-b5ea-77f726d5fe33
begin
	plot(k_on_vec, occ_th, xscale=:log10, label="", color=:firebrick)
	scatter!(k_on_vec, occ_th, xscale=:log10, label="Theory", color=:firebrick)

	plot!(k_on_vec, occupancy[1], xscale=:log10, label="", color=:forestgreen)
	scatter!(k_on_vec, occupancy[1], xscale=:log10, label="Simulations", color=:forestgreen)

	plot!(xlabel="kon", ylabel="occupancy", title="absolute densities")
	plot!(yscale=:log10)
	plot!(legend=:bottomright)
end

# ╔═╡ Cell order:
# ╟─7fa20536-710b-11ed-1d6e-d5ea9ba16f39
# ╠═3d756f7c-24d9-49fd-9e46-6bcee639ecf9
# ╠═645482f3-eca3-4a0b-966b-640a502a9051
# ╠═6aa4adfc-beb4-4615-8ae2-1d0eff4684ce
# ╠═e6611b5e-0c16-4e31-9724-fe726bca1943
# ╠═1ed4e2b1-4c8d-491c-b2fc-ee7005627aa8
# ╟─4612767c-6317-4453-80f7-b4200bfdae03
# ╠═6a6b01b9-dc76-46df-9d1d-2687125b4781
# ╠═3f858e05-5067-45fd-b5ea-77f726d5fe33
