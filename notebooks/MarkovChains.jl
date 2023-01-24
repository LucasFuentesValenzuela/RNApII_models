### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 5a1439f1-758f-40d5-a46a-6bf14bf81f5c
begin
	using Revise
	import Pkg
	Pkg.activate()
	using RNApIIModels
end

# ╔═╡ cd72c8d0-7747-11ed-0da4-a9933ab03aa2
begin
	using Plots
	using PlutoUI
	using Distributions
end

# ╔═╡ 11e4e2f4-2dfa-48e9-8b09-ba1d9d61075e
md"""
States are: 
- 1: 00
- 2: 01
- 3: 10
- 4: 11

**Promoter occupancy**: what proportion of the time is the first site occupied? that is, what is the probability to be in site 3 or 4? 

**Effective α**: what is the the reduced probability? (if we reduce the transition matrix)? 


"""

# ╔═╡ f95ab251-74b4-4015-9bd8-688cebf79a75
kon_vec = 10. .^(LinRange(-5, -.8, 20))

# ╔═╡ e889769e-e079-4583-b4ea-df10de284df2
begin
	α = 1e-2
end

# ╔═╡ 97a32a85-7542-4f01-98de-706c2f6e3294
koff = 1/2 - α

# ╔═╡ e00e1492-177f-4578-af99-4623e66fa483
@bind β Slider(10. .^(collect(LinRange(-4, 0, 10))))

# ╔═╡ 3272530f-5435-4d65-b3e1-f3a4f93e9605
β

# ╔═╡ 0ebf4aa7-1b9e-4926-b9df-507dd1bb4af8
begin
	plot(kon_vec, ρp.(kon_vec, koff, α), label="Theory", linewidth=2)

	for order in 2:5
		plot!(
			kon_vec, 
			ρp_MC.(α, β, kon_vec, koff; order=order), label="MC - $order", linestyle=:dash, linewidth=2
		)
	end
	# plot!(
	# 	kon_vec, 
	# 	ρp_MC.(α, β, kon_vec, koff; order=1), label="Markov Chain", linestyle=:dash, linewidth=2
	# )
	# plot!(
	# 	kon_vec, 
	# 	ρp_MC.(α, β, kon_vec, koff; order=2), 
	# 	label="Markov Chain-2", linestyle=:dash, linewidth=2
	# )
	plot!(xscale=:log10)
	plot!(xlabel="kon", ylabel="ρp")
	plot!(legend=:topleft)
end

# ╔═╡ 7e72f382-9e81-4066-b03d-0b1541945c5e
kon = kon_vec[1]

# ╔═╡ 7646279c-50b6-4b7a-943a-f247095725b1
Δt = 1e-2

# ╔═╡ 60c1852f-2809-4790-9d24-ab5f8151f73c
begin
	T(α, β, kon, koff, Δt) = transpose(hcat([
		[1 - kon*Δt, β*Δt, koff*Δt, 0],
		[0, 1 - (β + kon) * Δt, α*Δt, 0],
		[kon*Δt, 0, 1- (α+koff) * Δt, β*Δt],
		[0, kon*Δt, 0, 1-β*Δt]
	]...))
	
	T2(α, β, kon, koff, Δt) = transpose(hcat(
		[
			[1-kon * Δt, β * Δt, 0, 0, koff * Δt, 0, 0, 0], 
			[0, 1-(β+kon)* Δt, β* Δt, 0, 0, 0, 0, 0], 
			[0, 0, 1-(β+kon)* Δt, β* Δt, α* Δt, 0, 0, 0], 
			[0, 0, 0, 1- (β+kon) * Δt, 0, α * Δt, 0, 0], 
			[kon * Δt, 0, 0, 0, 1 - (α + koff) * Δt, β * Δt, 0, 0], 
			[0, kon * Δt, 0, 0, 0, 1 - (α + β) * Δt, β * Δt, 0], 
			[0, 0, kon * Δt, 0, 0, 0, 1 - β * Δt, β * Δt], 
			[0, 0, 0, kon * Δt, 0, 0, 0, 1 - β * Δt]
		]...
	))
end

# ╔═╡ 7a7f0863-fc4a-4770-a907-0a79f4280b0f
T(α, β, kon, koff, Δt)

# ╔═╡ bc3a4f7c-36c9-4f84-9e7b-575182d02f5c
T2(α, β, kon, koff, Δt)

# ╔═╡ Cell order:
# ╠═5a1439f1-758f-40d5-a46a-6bf14bf81f5c
# ╠═cd72c8d0-7747-11ed-0da4-a9933ab03aa2
# ╟─11e4e2f4-2dfa-48e9-8b09-ba1d9d61075e
# ╠═f95ab251-74b4-4015-9bd8-688cebf79a75
# ╠═e889769e-e079-4583-b4ea-df10de284df2
# ╠═97a32a85-7542-4f01-98de-706c2f6e3294
# ╠═3272530f-5435-4d65-b3e1-f3a4f93e9605
# ╠═e00e1492-177f-4578-af99-4623e66fa483
# ╠═0ebf4aa7-1b9e-4926-b9df-507dd1bb4af8
# ╠═7e72f382-9e81-4066-b03d-0b1541945c5e
# ╠═7646279c-50b6-4b7a-943a-f247095725b1
# ╠═7a7f0863-fc4a-4770-a907-0a79f4280b0f
# ╠═bc3a4f7c-36c9-4f84-9e7b-575182d02f5c
# ╠═60c1852f-2809-4790-9d24-ab5f8151f73c
