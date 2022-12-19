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

# ╔═╡ 6fa6e566-eb31-4803-924f-8e4e1a10b9f6
begin
	using Revise
	import Pkg
	Pkg.activate()
	using RNApIIModels
end

# ╔═╡ 3cb614e0-cbfa-40d6-b2fd-5e2a925a3216
begin
	using Plots
	using PlutoUI
end

# ╔═╡ 798003cc-abd9-42f2-8f28-a52713ac72f7
md"""# Description

We want to replicate the same type of analysis we do with the simulations purely from a theoretical basis. This allows faster iterations, results verification, and more fine-grained analysis

"""

# ╔═╡ 77bf5141-2715-4016-bd3d-e7aeaeee3a51
md"""# Load"""

# ╔═╡ 9b253a2e-baa0-4e23-838b-d437d89e0109
md"""
# Compression factor δ
"""

# ╔═╡ 50e76ee8-fdf5-4044-a233-3e063455c195
@bind δ Select([1, 35])

# ╔═╡ a13c3565-5dd5-4727-8d70-a38516f48804
md"""δ = $δ"""

# ╔═╡ 8e66639e-0b2f-4c0b-ac34-1c70ce7c3565
OCCUPANCY_PARAMS["δ"]

# ╔═╡ 089db03c-d59a-49ec-a4d9-08489347715d
L = Int(round(35/δ))

# ╔═╡ 50f02860-1673-4a38-b690-b9d8ecbd6a2c
md"""
# Feasible points for the average gene. 

Remember, those limits are only for the average gene in the average cell. 
"""

# ╔═╡ 693e7709-159e-4dce-ab68-4badf8f3f95f
md"""
The vertical 1/Ω line gives the maximum value of α beyond which we are not able to ensure that the RNAp actually remains the amount of time we tell it to (Ω). 
"""

# ╔═╡ 5a82ec54-7d95-4281-b84c-fd918bf29d70
Ω = 1

# ╔═╡ 7914de87-192d-4176-b4ca-b3ad13a5be95
# range of parameters
begin
	n_kon_pts = 100
	n_α_pts = 100
	
	k_on_vec = LinRange(
		LITERATURE_PARAMS["min_k_on"]/1.2, LITERATURE_PARAMS["max_k_on"]*1.2, n_kon_pts
	)
	α_vec = LinRange(LITERATURE_PARAMS["min_α"], min.(LITERATURE_PARAMS["max_α"] * 1.5, 1/Ω), n_α_pts)
	β_screen = LITERATURE_PARAMS["max_β"] / δ
end;

# ╔═╡ 99508ec5-1b8b-4402-92f1-8c86f001e5c6
koffs = max.(1/Ω .- α_vec, 0)

# ╔═╡ ce561af2-24e2-40fc-9664-dcb23ffbe861
α_eff = [effective_α.(k_on_vec, koffs[i], α_vec[i]) for i in 1:length(α_vec)]

# ╔═╡ 751e8d87-69cb-492b-ac48-d8ae6eb630a3
# dims are n_k_points x n_α_pts
begin
	
	ρ_g_th = reduce(hcat, 
		[map(
		f -> (ρ.(f, β_screen, OCCUPANCY_PARAMS["γ"], L))[2] * OCCUPANCY_PARAMS["n_sites"] * L, 
		α_eff[k]
	) for k in 1:length(α_vec)]
	);

	ρ_p_th = reduce(
	hcat, [
		map(f-> ρp.(f, koffs[i], α_vec[i]), k_on_vec) for i in 1:length(α_vec)
	]
	);

end;

# ╔═╡ c064f0e1-891b-4418-bb6e-400cce570957
feasible = is_feasible.(ρ_g_th, ρ_p_th);

# ╔═╡ afca47ff-3552-4b86-940b-770a6311f48b
begin
	heatmap(α_vec, k_on_vec, feasible)
	plot!(xlabel="α", ylabel="k_on")
	hline!([RNApIIModels.min_k_on, RNApIIModels.max_k_on], label="limits kon", linewidth=2)
	vline!([RNApIIModels.min_α, RNApIIModels.max_α], label="limits α", linewidth=2)
	vline!([1/Ω], label="1/Ω", linewidth=4, color=:red)
	plot!(title="Ω = $Ω")
	plot!(xlim=(0, 1))
end

# ╔═╡ c2e302d2-c02f-444f-bf92-134d6a7467be
md"""
Very minor differences in terms of which points are feasible in between the two represenations
"""

# ╔═╡ e6ae1115-a82a-4e95-a7f9-d9f9308c1b72
md"""
# Fold changes
"""

# ╔═╡ 753f983b-2612-49d9-8fa8-94bc60ab81ca
@bind Ωfc Slider(1:5)

# ╔═╡ 621e9282-c17a-422f-b291-01a0ae9cff8a
md"""Ωfc = $Ωfc"""

# ╔═╡ 2093cec3-c705-4f01-a38c-07f6d386deb2
md"""
Ωfc does have an impact because it scales the initiation rate (effective). BUT it also scales the feasible points (because it means higher occupancy, effectively). We see that if we rescale the kons with Ωfc it does not change anything to the distribution of fold changes. This is because in this regime everything is linear. 
"""

# ╔═╡ 7e4226a5-2223-4916-8482-2f1966c346fe
@bind αfc Slider(LinRange(LITERATURE_PARAMS["min_α"], min.(LITERATURE_PARAMS["max_α"] * 1.5, 1/Ωfc), n_α_pts))

# ╔═╡ e00ea34c-0954-44aa-be03-7085c8562908
md"""αfc = $(αfc)"""

# ╔═╡ 42aa81d2-270e-45b4-8434-dd1f2fe1fa81
md"""
We see that αfc has no impact on the fold changes. This is true because the density is linear in the initiation rate. Therefore, fold changes will be constnat. 
"""

# ╔═╡ 8fcddff7-97b8-4edd-a322-faec03a25e87
cs = LinRange(50, 150, 10)

# ╔═╡ 5e9bc538-7c92-4351-9fea-16da63ad275d
RNAs = RNApIIModels.CV_to_RNAfree_interp().(cs)

# ╔═╡ 23a735c9-1921-45ad-afe6-9909efe46beb
begin
	kons = 10. .^(LinRange(-3, -.5, 10))
	α_eff_fc(x) = effective_α(x, 1/Ωfc - αfc, αfc)
	
	ρs(x) = map(
		f -> ρ.(f, β_screen, OCCUPANCY_PARAMS["γ"], OCCUPANCY_PARAMS["L"])[2], α_eff_fc.(x)
	)
	
	xx = RNAs/RNAs[1]
	
	fcs = [ρs.(xx .* kon) ./ ρs(kon) for kon in kons]
end;

# ╔═╡ 47cdcca2-f7f9-429c-8146-3931c407cf06
begin
	colors = palette([:purple, :orange, :green], length(fcs))
	q = plot()
	for (k,fc) in enumerate(fcs)
		plot!(cs, fc, label="", color=colors[k])
	end
	plot!(xlabel="Cell volume [fl]", ylabel="Occ. FC")
	plot!(title="Ω = $Ωfc, α = $αfc")
	q
end

# ╔═╡ e6c51b0b-cf93-4934-9615-3ad3ae35d019
md"""
# Saturation away from the average gene

What we have so far analyses the feasibility for the average gene in the average cell. We can probably ask similar questions for higher expressed genes. 

But remember, we will always be constrained by how long the RNA has to remain on the promoter. However we can probably play with this. 
"""

# ╔═╡ 694cdb4b-2434-4a4c-aea9-d7d317ce6389
md"""
One way to ask this is in what regime you're in! From theory, we have three regimes: 
- Entry limited (the case for most of our analysis so far)
- Exit limited (never, as we assume the termination rate is very fast)
- Maximum current (limited by the elongation rate on the gene itself) 

In our setting, we have one constraint: that the RNA remains on the promoter for a duration Ω. For a given value of Ω, we have a range of acceptable value of α. This gives us a range of possible values for the effective `on` rate. And therefore, we can get a map of what regimes one system is in versus another. 
"""

# ╔═╡ f72a9a4b-8e48-4aeb-8fdb-00dc1ad10862
md"""
If we have to derive it formally, we have that the effective on rate is 

$α_{eff} = \frac{k_{on}α}{k_{on} + α + k_{off}}$, 

and that the condition at which we switch from a regime to another is 
$α_{eff}/β = 1/2$. 

That defines two regions of the plane, with different regimes. Let us see where we are. 

We can also say that at very high values of kon, the effective value converges toward $α$. The highest value that $α$ can ever have is `1/Ω`. 

Therefore, a criterion for the gene to be saturated is: 

$1/Ω \geq β/2 \Leftrightarrow Ω \leq 2/β$. If β is on the order of 10-50, then, it means that $Ω \leq .2$ ==> that is one order of magnitude in times shorter residence time. 
"""

# ╔═╡ 6adea95f-42b5-4703-b772-d3ce7d2d5752
@show δ, L

# ╔═╡ c40e2baa-2f7a-4d2e-b37a-561ddaf52f0f
@bind Ωtest Slider(LinRange(RNApIIModels.min_Ω/2, RNApIIModels.max_Ω*2, 100))

# ╔═╡ 31a759e6-5752-4850-8ef2-56e7cad02678
md"""Ωtest = $(Ωtest)"""

# ╔═╡ b4dc6631-83d7-4ee4-8a98-8fa566798ce5
@bind βtest Slider(
	LinRange(
		RNApIIModels.min_β/δ, RNApIIModels.max_β/δ, 
		100
	)
)

# ╔═╡ 9e44df78-4983-4d46-aec3-2c6a35ff1bfc
md"""
β = $βtest
"""

# ╔═╡ 2fc1d07b-e3de-47ea-bfc9-0b76d1c44f27
let

	α_vec = LinRange(RNApIIModels.min_α, 1/Ωtest, n_α_pts);

	# QUESTION; What is a reasonable range over which we can actually consider the kons
	k_on_vec = LinRange(0, 100 * RNApIIModels.max_k_on, n_kon_pts); 

	α_eff = reduce(hcat, [
		effective_α.(k_on_vec, 1/Ωtest - α_vec[i], α_vec[i]) for i in 1:length(α_vec)
	]
	)

	th = βtest/(1 + sqrt(L))

	regime = α_eff .< th

	p1 = plot()
	plot!(k_on_vec, α_eff[:, 1], label = "α = $(α_vec[1])")
	plot!(k_on_vec, α_eff[:, end], label="α = $(α_vec[end])")
	plot!(xlabel="kon", ylabel="α_eff")
	
	hline!([βtest/(1 + sqrt(L))], linewidth=2, linestyle=:dash, label="α_c")
	plot!(legend=:topleft)

	h1 = heatmap(α_vec, k_on_vec, regime)
	plot!(xlabel="α", ylabel="kon")

	plot([p1, h1]..., layout=(2, 1))
	
end

# ╔═╡ 179e1209-859d-46ce-bb70-3adfc208aaee
md"""
Ok, so this is interesting. It does indeed seem that in limit cases (very small residence time on the promoter, maximum value of α, very very large value of kon, you can get into a regime that is gene saturated). 

*However*, I am confused because the answer seems to change if you change the value of L, right? 
- β gets rescaled, as you need to account for it moving more slowly/rapidly. 
- The treshold for regime change also shifts, but not as fast as β changes I think? Indeed, it is $\frac{β}{1 + \sqrt{L}}$. So it seems that the threshold actually moves as $\sqrt{L}$, right? 

There are two things here: 
1. Is that what we see in simulations (we simulate with $L=1$, so we expect the theory should match). Do we observe something different for some values of kon? We have seen there is a discrepancy at high values of kon -- is that due to a potential close regime change? Why would that be? 
2. If the answer is significantly different depending on the way we model it, we probably need to at least investigate from the theoretical model what the answer is with `L = 35` and `β = 33-40`. 

"""

# ╔═╡ Cell order:
# ╠═798003cc-abd9-42f2-8f28-a52713ac72f7
# ╟─77bf5141-2715-4016-bd3d-e7aeaeee3a51
# ╟─6fa6e566-eb31-4803-924f-8e4e1a10b9f6
# ╠═3cb614e0-cbfa-40d6-b2fd-5e2a925a3216
# ╟─9b253a2e-baa0-4e23-838b-d437d89e0109
# ╟─a13c3565-5dd5-4727-8d70-a38516f48804
# ╠═50e76ee8-fdf5-4044-a233-3e063455c195
# ╠═8e66639e-0b2f-4c0b-ac34-1c70ce7c3565
# ╠═089db03c-d59a-49ec-a4d9-08489347715d
# ╟─50f02860-1673-4a38-b690-b9d8ecbd6a2c
# ╠═7914de87-192d-4176-b4ca-b3ad13a5be95
# ╠═99508ec5-1b8b-4402-92f1-8c86f001e5c6
# ╠═ce561af2-24e2-40fc-9664-dcb23ffbe861
# ╠═751e8d87-69cb-492b-ac48-d8ae6eb630a3
# ╠═c064f0e1-891b-4418-bb6e-400cce570957
# ╠═693e7709-159e-4dce-ab68-4badf8f3f95f
# ╠═5a82ec54-7d95-4281-b84c-fd918bf29d70
# ╠═afca47ff-3552-4b86-940b-770a6311f48b
# ╠═c2e302d2-c02f-444f-bf92-134d6a7467be
# ╟─e6ae1115-a82a-4e95-a7f9-d9f9308c1b72
# ╟─621e9282-c17a-422f-b291-01a0ae9cff8a
# ╠═753f983b-2612-49d9-8fa8-94bc60ab81ca
# ╠═2093cec3-c705-4f01-a38c-07f6d386deb2
# ╟─e00ea34c-0954-44aa-be03-7085c8562908
# ╟─7e4226a5-2223-4916-8482-2f1966c346fe
# ╠═42aa81d2-270e-45b4-8434-dd1f2fe1fa81
# ╠═23a735c9-1921-45ad-afe6-9909efe46beb
# ╠═47cdcca2-f7f9-429c-8146-3931c407cf06
# ╠═5e9bc538-7c92-4351-9fea-16da63ad275d
# ╠═8fcddff7-97b8-4edd-a322-faec03a25e87
# ╟─e6c51b0b-cf93-4934-9615-3ad3ae35d019
# ╟─694cdb4b-2434-4a4c-aea9-d7d317ce6389
# ╟─f72a9a4b-8e48-4aeb-8fdb-00dc1ad10862
# ╠═6adea95f-42b5-4703-b772-d3ce7d2d5752
# ╟─31a759e6-5752-4850-8ef2-56e7cad02678
# ╠═c40e2baa-2f7a-4d2e-b37a-561ddaf52f0f
# ╟─9e44df78-4983-4d46-aec3-2c6a35ff1bfc
# ╠═b4dc6631-83d7-4ee4-8a98-8fa566798ce5
# ╟─2fc1d07b-e3de-47ea-bfc9-0b76d1c44f27
# ╠═179e1209-859d-46ce-bb70-3adfc208aaee
