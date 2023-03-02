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

# ╔═╡ 9d253cf6-a8c1-4f01-8604-b05044e4a72f
begin
	using Revise
	import Pkg
	Pkg.activate()
	using RNApIIModels
end

# ╔═╡ 78b6d800-c9b7-4487-aaec-c804410278f3
begin
	using Plots
	using PlutoUI
	using Distributions
end

# ╔═╡ 9c79970a-f98c-40fb-bd47-65a6164699bf
md"""
# Questions

* What can we really extract from the analytical models: scaling of the current, critical parameters, scaling of the densities, ...? 

* What is the impact of the termination rate? Do we **really** care about the way the end is modeled? I don't think so. However, if we want to model the occupancy, we probably need it, right? Or we approximate densities from Eqs (25) in Lazaros, Chou and take ρ_(N/2) as an approximation for the steady-state density and ignore the boundary effects. 

* What is the real impact of δ and L? How do results actually change, in what regime? How to make sense of everything? 
"""

# ╔═╡ 94c158e1-4de0-4b20-89e8-b318b6ca3923
md"""# Models"""

# ╔═╡ bfe95384-e816-4ae8-afc2-631fc824bbef
md"""We can analyze the limit regimes from (Klumpp, Hwa - 2008). 

$J(α) = α (ε - α) / (ε + α (L-1)),$

$J_{max}= ε/(1+L^{1/2})^2.$ 

In the above: 
- L is the object size
- ε is the attempt rate of elongation
- α is the initiation rate
"""

# ╔═╡ 7359697b-cc44-47b6-8a8c-2230dbee2eeb
md"""
**Note**: the general model is actually derived in Lakatos, Chou, 2003. 
"""

# ╔═╡ ca87554a-0149-446a-99e9-443ca7dfd389
md"""
### Transcription rates for d = 1
"""

# ╔═╡ b5342ffd-0d37-44d1-9656-7929e0418156
#illustration of the theoretical transcription rate curves for large γ 

let

	plot()


	rates = collect(LinRange(0.01, 0.1, 20))
	for k in 1:length(rates)

		β = rates[k]
		color = palette([:blue, :green], length(rates))[k]
		plot!(rates/β, J.(vec(rates), β, nothing, 1)/β, label="", linestyle=:dash, color=color)	
		scatter!(rates/β, J.(vec(rates), β, nothing, 1)/β, label="", color=color)
	end
	
	# ylims!(0, .5)
	xlabel!("α")
	ylabel!("Current J")
	title!("Currents")
	
end

# ╔═╡ c4edf93a-b67e-415c-8d73-43d2c44050ee
md"""
The invariance can be understood because 
$J(α)/ε = α/ε (1 - α/ε) / (1 + α/ε (L-1)).$

This is the case regardless of $L$.
"""

# ╔═╡ 6ee4187d-f8fb-4617-83fe-b7b38a43404d
md"""
### Playing with transcription rates
"""

# ╔═╡ e8577416-39cc-4dac-a5b3-e04ae43d3d7b
rates = 10. .^[-3, -2, -1, 0, 1, 2, 3]

# ╔═╡ a2780f8f-20e9-49d0-a181-46ad13def02d
@bind β Slider(rates)

# ╔═╡ 5d36ad32-3e1d-44a0-a524-998731ef0f27
md"""β = $β"""

# ╔═╡ ff9ea155-5664-433c-8a4d-3cc8f7f857b6
@bind γ Slider(rates)

# ╔═╡ 5a8876aa-4b87-487a-8e39-427e48a8e10d
md"""γ = $γ"""

# ╔═╡ 6c541268-e624-4ea6-be16-830b2e3fccc9
@bind L Slider([1, 5, 10, 20, 30, 50])

# ╔═╡ 93984140-86e3-4422-a96e-7a4f313030dd
md"""L = $L"""

# ╔═╡ 52cc9dc6-b816-4d30-a9b2-abe8599e9616
let 
	α_vec = 10. .^(collect(LinRange(-3, 1, 100))) .* β

	J = J.(α_vec, β, γ, L)

	plot(α_vec ./β, J, linewidth=2, label="")

	# plot!(xscale=:log)
	vline!([1], label="α=β", linestyle=:dash)
	xlabel!("α/β")
	ylabel!("J")
	plot!(legend=:bottomright)
	
end

# ╔═╡ a2fc5a3f-d5c8-4a74-8d83-c9639cc08a32
md"""
### Parameter study for different L, same β
"""

# ╔═╡ 21f9f9a1-bc6a-4de2-99d6-b48738caca2c
md"""
The critical parameter is expressed as

$$α_c = \frac{β}{1 + \sqrt{L}}$$
"""

# ╔═╡ 4173cf58-f5fc-4758-b1cd-1919cd34265e
md"""
Therefore, there is a very strong dependence on the size of the particle for the critical rate at which we expect current to reach a maximum, for a same value of β. Namely, for a size $L' = δL$, we have that $α_c'/α_c = \frac{1 + \sqrt{L}}{1+ \sqrt{δL}} \sim 1/\sqrt{δ}$ for large L. 
"""

# ╔═╡ b4b30449-3d6f-438e-a398-57d799888a68
md"""
Here below we assume that the rate β does not change (i.e. the translation from one site to the next keeps the same rate/probability). We only change the size/footprint of the particle of interest. 
"""

# ╔═╡ f067ebdf-3ab8-4efd-ad93-e26ecadfe11a
# current in regimes across α (assume γ is large)
let
	L_vec = [1, 2, 5, 10, 20]

	β = .1

	α_vec = collect(LinRange(0.01, 1, 30)) .* β

	p = plot()

	for L in L_vec
		plot!(α_vec/β, J.(α_vec, β, nothing, L)/β, label="L=$L")
	end
	xlabel!("α/β")
	ylabel!("J/β")
	p
end

# ╔═╡ c463821d-8a79-4fca-b59f-b4a354ceb606
md"""We see that the transition point decreases with $L$ and that the current is smaller with increasing L"""

# ╔═╡ cc28a585-97b2-48f4-b133-b06d8ddf8927
md"""## functional forms with rescaling"""

# ╔═╡ cdbd7292-cdbe-4e40-bb68-f5115ede9c16
# transition point α_c
let
	β = .1 
	L̄ = 50
	L_vec = LinRange(1, L̄, 100)
	δ = L̄ ./ L_vec
	β_vec = β ./ δ
	α_c = β_vec ./ (1 .+ sqrt.(L_vec))

	plot(L_vec, α_c, label="", xlabel="L", ylabel = "α_c", title="α_c with β rescaling")

end

# ╔═╡ e4e9dd81-6380-402c-ba2e-6147b0f32493
md"""we see that the critical value for α with a rescaling of the elongation rate β changes as √L"""

# ╔═╡ f654bd6b-6850-4bb0-b4bb-d45311d21769
md"""**question** when does the transition appear if we make one assumption vs another? what is the scale of the difference/mistake/error?

For instance, let's say that we want to compare the point of transition between the two regimes for 
* the real parameters
* a rescaled version of the parameters

We have 
$$α_c = \frac{β̄}{1 + \sqrt{L̄}}$$ and $$α_c' = \frac{β̄/δ}{1 + \sqrt{L̄ /δ}}$$

The error is
$$α_c' / α_c = \frac{1 + \sqrt{L̄}}{δ + \sqrt{δL̄}}$$

"""

# ╔═╡ 3fff6894-bbe5-4589-93c8-aee7d06769a5
αc_err(δ, L̄) = (1 + sqrt(L̄)) ./ (δ + sqrt.(δ*L̄))

# ╔═╡ 889ff8cf-d1a7-4ade-adea-16c216694d2b
let
	p = plot()
	
	L̄_vec = [2, 5, 10, 20, 30, 50]
	
	for L̄ in L̄_vec
		δ = LinRange(1, L̄, 100)
		plot!(δ./L̄, αc_err(δ, L̄), label="L=$L̄")
	end
	plot!(xlabel="δ/L̄", ylabel="α_c'/α_c")

	p
end

# ╔═╡ 31ad0b42-7a71-493f-89ec-0d980e95f7ef
md"""We see that we make more and more error if we compress the system (i.e. δ is large). More particularly, let's look only at δ = L̄ (i.e. L = 1), as a function of L̄."""

# ╔═╡ b031a2aa-c3e2-4be3-a049-df7a303972d7
let
	L̄ = LinRange(1, 100, 1000)

	e = αc_err.(L̄, L̄)
	
	plot(
		L̄, e, 
		label="", 
		xlabel="L̄", ylabel="α_c'/α_c", 
		title="error for δ = L̄"
	)
end

# ╔═╡ 3eae7ee1-4f9a-4509-82c2-b41a0f4c3b9c
md"""
In our case, we have δ = 35, L̄ = 35. We see that α_c will be underestimated by one order of magnitude! What this means, therefore, is that you need to actually be able to go 10x larger in your initiation rate to actually cause a change in densities.
"""

# ╔═╡ 56668918-92ed-49a7-b3fe-38f6b9215cd7
@show αc_err(35, 35)

# ╔═╡ 5275d868-09ab-43ea-a336-7267e2084353
md"""
### Different L, rescaled β
"""

# ╔═╡ 2e9e5fd9-9e3e-44dd-8a47-61a647097325
md"""
Here we want to grasp the impact of introducing granularity into the model. 

Right now we are assuming that the particles take 35 bp jumps. That is, all spatial coordinates have been rescaled by 35 (x' = x/35). Such that β is actually 35 x too small (it moves slower because takes larger steps). 

Let us denote by δ the rescaling factor, i.e. δ = 35/L, and β = β̄ / δ. 

Therefore, we have to compare three results: 
* δ = 35 (i.e. L = 1), β = β̄/δ: everything scaled so that the RNAp takes only 1bp jumps (i.e. what we have mostly been modelling). 
* δ = 1 (i.e. L = 35), β = β̄: large particles moving at the given proper rate (what we should be modelling to be really accurate)
* δ = 35 (i.e. L=1), β = β̄: what is the impact of changing the size of the representation but keeping the actual rate

Note: exploring the fourth possibility (i.e. δ = 1, L = 35 with β = β̄) does not make a lot of sense, as there is no reason it would help us. 
"""

# ╔═╡ 2adc96f0-aeae-4c8a-80e6-1ed6c9ddffcc
@bind β_rescale Slider(rates)

# ╔═╡ ef2863e3-ad2f-45fc-9bee-3846581818fe
md""" β = $β_rescale"""

# ╔═╡ 825da303-dd0b-494c-bc46-eb7b8cfaf444
γ_rescale = nothing

# ╔═╡ 8dfbb51d-5d37-4783-8467-6faec8111653
md""" γ = $γ_rescale"""

# ╔═╡ fd8625d5-0065-43a8-9dcf-1f0dabf043a7
# currents, again assuming that we are not in γ limited regime

let
	β̄ = β_rescale
	L̄ = 35
	L_vec = [1, 10, 20, 25, 30, L̄]
	δ = L̄ ./ L_vec
	β_vec = β̄ ./ δ
		
	α_ratio_vec = 10. .^(collect(LinRange(-3, 3, 1000)))

	p = plot()

	for k in 1:length(β_vec)
		β = β_vec[k]
		L = L_vec[k]
		αs = α_ratio_vec ./ β
		
		plot!(
			αs, J.(αs, β, γ_rescale, L), 
			label="L=$L, β=$(round(β; digits=3))", linewidth=2
		)

	end
	
	xlabel!("α")
	ylabel!("J")
	plot!(legend=:outertopright)
	plot!(xscale=:log10)
	# xlims!(0, 5*β)
	# ylims!(0, β)
	# hline!([β], label="J=β", linestyle=:dash, linewidth=3)
	p
end

# ╔═╡ 425c6f37-936b-466a-acae-d5c2efdd5f3d
md"""
We see that the maximum achievable J is higher for larger $L$. The reason is that the maximum current scales as 

$J \propto β(L)/(\sqrt{L}+1)^2)$

If we impose that β(L) = β̄L, then we see that

$J \propto \bar{β}L/(\sqrt{L}+1)^2)$

Therefore the current really maxes out at β, and gets there as $L$ increases

"""

# ╔═╡ 54bf5fb0-c219-4385-8ee6-88411821ac32
md"""
# Densities
"""

# ╔═╡ 772e015a-fa71-4cab-83d9-ac24762b132b
md"""
there must be something useful/interesting to extract regarding densities, even if it is only the limit + medium densities. You can probably make some assumptions purely from the theoretical models. 
"""

# ╔═╡ c2062357-22d5-48e5-976a-8cb2b473d0f0
md"""### Entry limited"""

# ╔═╡ 4cfa3406-a838-4449-9ec6-ebe033d4fae6
# for γ very large

let
	
	β =.1
	γ = nothing

	x = [0, .5, 1]

	color_palette = palette([:blue, :green], 3)

	plots = []
	
	for (k,α) in enumerate([.1, 1, 10] .* β)
		p = plot()
		
		for (i,L) in enumerate([1, 10, 30])
			ρL, ρN, ρR = ρ(α, β, γ, L)
			hline!([L*ρN], linestyle=:dash, color=color_palette[i], label="α/β=$(round(α/β; digits=2)), L=$L")
			scatter!(x, [L*ρL, L*ρN, L*ρR], color=color_palette[i], label="", alpha=.3)
		end
		xlabel!("Relative length")
		ylabel!("Lρ")
		plot!(legend=:outertopright)
		xlims!(-.1, 1.1)
		ylims!(-.1, 1.1)

		push!(plots, p)
	
	end

	plot(plots..., layout = (length(plots), 1))
	
end

# ╔═╡ 42312f09-54d6-4ae9-a67b-d06f11c85c31
md"""
**Note** I don't understand why the densities are very low at the beginning if α is very large...? I would expect the opposite... 
"""

# ╔═╡ 4a9578ac-1580-40ce-a72d-26815fea6129
md"""### Exit limited"""

# ╔═╡ de58ca4e-7134-4910-9532-6defa8e7689c
# for α very large

let
	
	β =.1
	α = 1000

	x = [0, .5, 1]

	color_palette = palette([:blue, :green], 3)

	plots = []
	
	for (k,γ) in enumerate([.1, 1, 10] .* β)
		p = plot()
		
		for (i,L) in enumerate([1, 20, 30])
			ρL, ρN, ρR = ρ(α, β, γ, L)
			hline!([L*ρN], linestyle=:dash, color=color_palette[i], label="γ/β=$(round(γ/β; digits=2)), L=$L")
			scatter!(x, [L*ρL, L*ρN, L*ρR], color=color_palette[i], label="", alpha=.3)
		end
		xlabel!("Relative length")
		ylabel!("Lρ")
		plot!(legend=:outertopright)
		xlims!(-.1, 1.1)
		ylims!(-.1, 1.1)

		push!(plots, p)
	
	end

	plot(plots..., layout = (length(plots), 1))
	
end

# ╔═╡ 7e82dc32-10e0-4e55-84fd-80c6eca84514
md"""### Occupancy

Can we say anything about the occupancy over the strand, and whether or not we see this abrupt transition?"""

# ╔═╡ 2d2c0d45-5d3b-4256-a9d8-86a44675db01
# do we observe the transition as in both models I have analyzed last week? 

let
	β = 1
	γ = .1
	Ls = [1, 10, 20, 30]
	
	color_palette = palette([:blue, :green], length(Ls))

	αs = 10 .^collect(LinRange(-3, 1, 50)).* γ
	p = plot()
	
	for (i,L) in enumerate(Ls)
		ρs = ρ.(αs, β, γ, L)
		ρN = [ρ[2] for ρ in ρs]
		plot!(
			αs./γ, ρN, linestyle=:dash, label="L=$L", color=color_palette[i]
		)
		scatter!(
			αs./γ, ρN, linestyle=:dash, label="", color=color_palette[i]
		)
	end

	plot!(xscale=:log10, yscale=:log10)
	vline!([1], label="", color=:red)
	plot!(legend=:topleft)
	plot!(title="Occupancy as a function of particle size")

	xlabel!("α/γ")
	ylabel!("~ Occupancy")


end

# ╔═╡ 7e0462d8-e251-4b1c-9d22-3c905b89f95d
# do we observe the transition as in both models I have analyzed last week? 

let
	β = 1
	γ = nothing
	Ls = [1, 2, 10, 20, 30]
	
	color_palette = palette([:blue, :green], length(Ls))

	αs = 10 .^collect(LinRange(-3, 2, 100))
	p = plot()
	
	for (i,L) in enumerate(Ls)
		ρs = ρ.(αs, β*L, γ, L)
		ρN = [ρ[2] for ρ in ρs]
		plot!(
			αs, ρN * L, linestyle=:dash, label="L=$L, β = $(round(β*L; digits=3))", color=color_palette[i], linewidth=2
		)

		vline!([β*L / (1 + L^(1/2))], color=:gray, linewidth=.5, label="")
		# scatter!(
		# 	αs, ρN, linestyle=:dash, label="", color=color_palette[i]
		# )
	end

	plot!(xscale=:log10, yscale=:log10)
	plot!(legend=:topleft)
	plot!(title="Occupancy as a function of particle size")

	xlabel!("α")
	ylabel!("~ Occupancy")

end

# ╔═╡ 1c909f4b-7743-4ac8-bfc8-32f2de8fc63b
md"""
We see that as long as the gene is not saturated (as long as we are in the entry-limited regime) there is little variation between the different representations.

**But** there is about one order of magnitude between the different transitions -- and the high δ cases saturate much faster
"""

# ╔═╡ fa404ca1-049b-49fe-b894-0cba3080dc01
md"""
## Where are we on the phase space? 

We have values for all the parameters, can we make a statement about where we are in the phase space? 
"""

# ╔═╡ c7ca281b-e980-44f1-8963-f2c04290225d
md"""
The current parameters are: 
* β = 20
* γ = 1/70
* α = 0.0033
* L = 35
"""

# ╔═╡ d15931d9-2263-4b34-b972-6cb1c7a7ae1b
md"""
Therefore, the normalized parameters (wrt β) have the following value:

* α = 0.000165
* γ = 0.0007
* L = 35
"""

# ╔═╡ 32741bf0-ef1a-4ba1-9a48-2363bdb52e57
md"""
According to (Lakatos, Chou) the critical parameter values for α and γ are

$α_c = γ_c = \frac{1}{\sqrt{d} + 1} = 0.15$

Therefore, we see that $α \ll α_c$ and $\gamma \ll \gamma_c$. 
"""

# ╔═╡ a0fb1d2c-56a8-42e1-bfeb-aae5040df1cd
md"""To determine whether we are in the entry or exit limited regime, we have to determine whether $α \lesseqgtr \gamma$.

Basically, it really seems that $\gamma$ will act as the threshold rate at which we get a strong change in density.
""" 

# ╔═╡ 29246f54-4220-421c-ae41-62d9593686cb
md"""
# Dependence of transcription rate on γ

This can be useful to study the second strand problem
"""

# ╔═╡ e87db36e-70bf-4006-bd08-9eecb2f396d5
ranges = 10. .^(collect(LinRange(-3, 1, 15)));

# ╔═╡ ce96621b-bdf2-4a8c-b127-6306b395e5c9
@bind α_ Slider(ranges)

# ╔═╡ 5cdb2d51-257a-4f82-80c3-308a7d624d06
md"""α_ = $α_"""

# ╔═╡ e3bbb1f0-52f3-4279-b426-5c902305709b
@bind β_ Slider(ranges)

# ╔═╡ d22e249b-8bce-47d1-8dd9-dbb75fb15997
md"""β_ = $β_"""

# ╔═╡ 9560d375-ed9b-4ef6-8cf1-9f7b509f2af3
# currents, again assuming that we are not in γ limited regime

let

	γ_vec = 10. .^(collect(LinRange(-3, 3, 1000)))

	p = plot()

	plot!(
		γ_vec, J.(α_, β_, γ_vec, 1), 
		label="", linewidth=2
	)

	
	xlabel!("γ")
	ylabel!("J")
	plot!(legend=:topleft)
	plot!(xscale=:log)
	# xlims!(0, )
	ylims!(0, 1.)
	hline!([α_], label="J=α", linestyle=:dash, linewidth=3)
	p
end

# ╔═╡ 2e8e2d60-ec3a-440a-b178-33f7c541dea7
md"""
## What is a probable termination rate? 
"""

# ╔═╡ 89ecf69e-3d2c-4aac-abea-556fbadba1d5
θ = collect(LinRange(0, 1, 1000))

# ╔═╡ 7bc355fd-fb4d-4c18-ab2a-486b8ba21680
Tt(θ, Et) = (21 - 5θ)/(1-θ) - Et

# ╔═╡ 16db0c80-a9e4-439d-8c8b-fbab797e9690
let
	p = plot()
	for Et in LinRange(30, 80, 10)
		plot!(θ, Tt.(θ, Et), label="", linewidth=2, color=:gray)
	end
	hline!([0], linestyle=:dash, label="", linewidth=2)
	ylims!(0, 100)
	xlabel!("Relative fraction of promoter to gene-body bound RNApII")
	ylabel!("Tt (s)")
end

# ╔═╡ 7d66c043-19bf-4bbf-ab4f-8ab597dc4ee6
md"""
## Second strand models
"""

# ╔═╡ 0babe808-4ec9-424c-8ffa-9e686f0ca4d9
md"""
## Trying to understand how the occupancy should really scale with `α` or `kon`. 
"""

# ╔═╡ 1f759a4f-5dca-46d1-b8d1-482e7a610f79
let 
	
	xx = collect(10. .^(collect(LinRange(-5, 2, 100)))) # kon
	kin = 1/2
	kout = 1/4
	α(xx) = xx * kin ./ (xx .+ kin .+ kout) # effective alpha

	p0 = plot(xx, α.(xx), label="")
	plot!(xlabel="kon", ylabel="α_eff")
	plot!(title="theoretical model")
	plot!(xscale=:log10)
	
	ρ(α) = 1/2 .* (1 .- sqrt.(1 .- 4α .* (1 .- α)))

	p1 = plot(xx, ρ.(α.(xx)), label="")
	plot!(yscale=:log10, xscale=:log10)
	plot!(xlabel="kon", ylabel="occupancy at N/2")

	fc(x) = ρ(α(2*x))/ ρ(α(x))

	p2 = plot(xx, fc.(xx), label="", ylim=(0, 3.), xscale=:log10)
	plot!(xlabel="kon", ylabel="occ f.c.")

	plot([p0, p1, p2]..., layout=(3, 1), size=(800, 900))
end

# ╔═╡ fd8d3155-226d-43e1-b545-35e520673e9b


# ╔═╡ 3f2c4518-a44b-4d79-bba0-51187c84dad9
md"""
# Promoter model
"""

# ╔═╡ f08b27e9-5328-4361-90d9-9ac944b090c0
md"""
We model the system with a promoter, which has us lose a lot of effective computation time. It is important if you have promoter jamming before gene body jamming (we will see if that is the case separately, from simulations), but maybe there is a way to link an effective on-rate which would enable us to make links with theory. 

I believe that the effective initiation rate 

$\bar{α} = \frac{k_{on}α}{k_{on} + k_{out} + \alpha}$

Let us test that, and see when this breaks down (which I guess will be whenever the occupancy of the promoter region is too high. )
"""

# ╔═╡ 5036ec33-5c94-4b5e-bcdf-089365166343
function simulate_promoter(k_on_vec, α, β, k_out, n_steps, Δt)
	sites = [0, 0]
	n_inits = []
	n_bindings = []
	promoter_occ = []
	for k_on in k_on_vec

		n_inits_crt = 0
		n_bindings_crt = 0
		promoter_occ_crt = 0
		
		
		for _ in 1:n_steps

			if sites[1] == 1
				promoter_occ_crt += 1
			end
			
			if (rand(Bernoulli(k_on * Δt))) & (sites[1] == 0)
				sites[1] = 1
				n_bindings_crt +=1 
			end
			
			if (rand(Bernoulli(β * Δt))) & (sites[2] == 1)
				sites[2] = 0
			end

			if (sites[1]==1)
				s = wsample(
					["off", "init", "nothing"], [k_out*Δt, α*Δt, 1-Δt*(k_out+α)]
				)
				if s=="init"
					n_inits_crt += 1
				end
				
				if (s=="off")
					sites[1]=0
				elseif (s=="init") & (sites[2]==0)
					sites[2]=1
					sites[1]=0
				end
			end
			
		end

		push!(n_inits, n_inits_crt)
		push!(n_bindings, n_bindings_crt)
		push!(promoter_occ, promoter_occ_crt)
		
	end

	return n_inits, n_bindings, promoter_occ
end

# ╔═╡ a46cd130-d46a-4a3c-b99d-a8fe69640ea8
begin
	k_on_vec = 10. .^(LinRange(-1, 2, 20))
	Δtp = 1e-2
	βp = 1

	Ω = 2
	αp = 1e0/1 * 1/Ω
	k_out_p = 1/Ω - αp
	@show k_out_p
	n_steps = 1e6
	
	n_inits, n_bindings, promoter_occ = simulate_promoter(k_on_vec, αp, βp, k_out_p, n_steps, Δtp)
end

# ╔═╡ da26766c-f908-41da-bde4-ee2899237dac
let
	p1 = scatter(k_on_vec, n_inits ./ (n_steps * Δtp), label="simulations")
	plot!(k_on_vec, k_on_vec .* αp ./ (k_on_vec .+ αp .+ k_out_p), label="theory")
	plot!(legend=:bottomright)
	# plot!(xscale=:log10)
	plot!(xlabel="kon", ylabel="effective initiation rate")
	plot!(xscale=:log10)
	# hline!([βp], label="β/2")

	p2 = scatter(k_on_vec, promoter_occ ./ (n_steps), label="simulations")
	plot!(k_on_vec, k_on_vec ./ (k_on_vec .+ αp .+ k_out_p), label="theory")
	plot!(legend=:bottomright)
	plot!(xlabel="kon", ylabel="promoter occupancy")
	plot!(xscale=:log10)

	plot([p1, p2]..., layout=(1, 2))
end

# ╔═╡ 0b5f1ed0-51c0-41d0-9946-9e3f6c91b524
md"""
The things that can change the validity of the predictions are:
- value of β
- value of α
- value of Ω?
"""

# ╔═╡ 097625a8-43a1-44a7-86b3-04f0b718f505
md"""
What dictates the discrepancy between these values?  It seems that as I decrease β or increase α we have a difference between theory and prediction. 

Let us check that for a single value of kon, many values of β, α
"""

# ╔═╡ ec875061-2f1d-4d9b-8c0a-d6a6925dc4a4
function α_eff_map(Ω, α, β; kon_map = 1)
	
	α_eff = zeros(length(α), length(β))
	α_eff_th = copy(α_eff)

	k_out = 1/Ω .- α

	Δt = 1e-2
	n_steps = 5e6

	for (i,j) in collect(Iterators.product(1:length(α), 1:length(β)))
	
		n_inits, n_bindings, promoter_occ = simulate_promoter(
			kon_map, α[i], β[j], k_out[i], n_steps, Δt
		)
		α_eff[i,j] = n_inits[1] / (n_steps * Δt)

		α_eff_th[i,j] = effective_α(kon_map, k_out[i], α[i])
	end

	return α_eff ./ α_eff_th
	
end

# ╔═╡ e8703636-bbfe-46e2-9191-bd70f45d4a57
begin
	kon_map = 1

	Ω_map = 1
	α_map = (10. .^(collect(LinRange(-2, 0, 10)))) .* 1/Ω_map
	β_map = 10. .^(collect(LinRange(-2, 1, 5)))
end;

# ╔═╡ 06fc9bd4-4f61-4c6f-b096-4fd94c7d26f3
α_map

# ╔═╡ 7345c8db-73ce-4a90-a1b2-9478ab4d8a41
β_map

# ╔═╡ 757777f2-889e-4156-afd7-79baf084b20d
α_rat = α_eff_map(Ω_map, α_map, β_map; kon_map=kon_map);

# ╔═╡ 24d55bf6-aed0-4d87-9682-00eb828b19bd
heatmap(
	β_map, α_map, α_rat, 
	ylabel="α", xlabel="β", 
	xscale=:log10, yscale=:log10, 
	c=cgrad(:matter, 10, categorical = true)
)

# ╔═╡ 4960f6f2-7569-48f7-82f2-d5fd6af151ee
md"""
We see that the prediction fails for 
- low values of β
- high values of α

We are not able to grasp the real value of the effective initiation rate
"""

# ╔═╡ 5af99272-9b87-4b25-8eb4-e75e0a14cfc3
md"""# Influence of Ω"""

# ╔═╡ 24ae324c-8e00-4749-8440-1fb471ddf2b7
@bind αΩ Slider(
	LinRange(
		LITERATURE_PARAMS["min_α"], 
		min(LITERATURE_PARAMS["max_α"]), 
		10
	)
)

# ╔═╡ ca29e38b-171b-4dd7-b7ec-b3d31f66a19a
koffs(x) = 1/x - αΩ

# ╔═╡ 2585179f-14b1-4d28-81a3-f8e86b299f35
let 
	q = plot()
	for Ω in 1:5
		α_eff = effective_α.(RNApIIModels.k_on_vec_screen, koffs(Ω), αΩ)

		plot!(RNApIIModels.k_on_vec_screen, α_eff, label="Ω=$Ω")
	end
	plot!(ylim=(0, 0.01))
	plot!(xlabel="kon", ylabel="α_eff", title="α = $αΩ")
	q
end

# ╔═╡ fe6b9eff-0963-4008-abc5-2efb99748bb2
md"""
We see that the initiation rate (effective) increases with its residence time. The reason is that koff decreases, so it is more likely to actually initiate. 
"""

# ╔═╡ Cell order:
# ╠═4ddbce24-6f3d-4e54-83cb-efddb662546e
# ╠═78b6d800-c9b7-4487-aaec-c804410278f3
# ╠═9d253cf6-a8c1-4f01-8604-b05044e4a72f
# ╟─94c158e1-4de0-4b20-89e8-b318b6ca3923
# ╟─bfe95384-e816-4ae8-afc2-631fc824bbef
# ╟─7359697b-cc44-47b6-8a8c-2230dbee2eeb
# ╟─ca87554a-0149-446a-99e9-443ca7dfd389
# ╟─b5342ffd-0d37-44d1-9656-7929e0418156
# ╟─c4edf93a-b67e-415c-8d73-43d2c44050ee
# ╟─6ee4187d-f8fb-4617-83fe-b7b38a43404d
# ╠═e8577416-39cc-4dac-a5b3-e04ae43d3d7b
# ╟─5d36ad32-3e1d-44a0-a524-998731ef0f27
# ╟─a2780f8f-20e9-49d0-a181-46ad13def02d
# ╟─5a8876aa-4b87-487a-8e39-427e48a8e10d
# ╟─ff9ea155-5664-433c-8a4d-3cc8f7f857b6
# ╟─93984140-86e3-4422-a96e-7a4f313030dd
# ╟─6c541268-e624-4ea6-be16-830b2e3fccc9
# ╟─52cc9dc6-b816-4d30-a9b2-abe8599e9616
# ╟─a2fc5a3f-d5c8-4a74-8d83-c9639cc08a32
# ╠═21f9f9a1-bc6a-4de2-99d6-b48738caca2c
# ╟─4173cf58-f5fc-4758-b1cd-1919cd34265e
# ╟─b4b30449-3d6f-438e-a398-57d799888a68
# ╠═f067ebdf-3ab8-4efd-ad93-e26ecadfe11a
# ╟─c463821d-8a79-4fca-b59f-b4a354ceb606
# ╠═cc28a585-97b2-48f4-b133-b06d8ddf8927
# ╠═cdbd7292-cdbe-4e40-bb68-f5115ede9c16
# ╟─e4e9dd81-6380-402c-ba2e-6147b0f32493
# ╠═f654bd6b-6850-4bb0-b4bb-d45311d21769
# ╠═3fff6894-bbe5-4589-93c8-aee7d06769a5
# ╠═889ff8cf-d1a7-4ade-adea-16c216694d2b
# ╟─31ad0b42-7a71-493f-89ec-0d980e95f7ef
# ╠═b031a2aa-c3e2-4be3-a049-df7a303972d7
# ╟─3eae7ee1-4f9a-4509-82c2-b41a0f4c3b9c
# ╠═56668918-92ed-49a7-b3fe-38f6b9215cd7
# ╟─5275d868-09ab-43ea-a336-7267e2084353
# ╟─2e9e5fd9-9e3e-44dd-8a47-61a647097325
# ╟─ef2863e3-ad2f-45fc-9bee-3846581818fe
# ╠═2adc96f0-aeae-4c8a-80e6-1ed6c9ddffcc
# ╟─8dfbb51d-5d37-4783-8467-6faec8111653
# ╠═825da303-dd0b-494c-bc46-eb7b8cfaf444
# ╠═fd8625d5-0065-43a8-9dcf-1f0dabf043a7
# ╟─425c6f37-936b-466a-acae-d5c2efdd5f3d
# ╟─54bf5fb0-c219-4385-8ee6-88411821ac32
# ╟─772e015a-fa71-4cab-83d9-ac24762b132b
# ╟─c2062357-22d5-48e5-976a-8cb2b473d0f0
# ╟─4cfa3406-a838-4449-9ec6-ebe033d4fae6
# ╟─42312f09-54d6-4ae9-a67b-d06f11c85c31
# ╟─4a9578ac-1580-40ce-a72d-26815fea6129
# ╟─de58ca4e-7134-4910-9532-6defa8e7689c
# ╠═7e82dc32-10e0-4e55-84fd-80c6eca84514
# ╟─2d2c0d45-5d3b-4256-a9d8-86a44675db01
# ╠═7e0462d8-e251-4b1c-9d22-3c905b89f95d
# ╟─1c909f4b-7743-4ac8-bfc8-32f2de8fc63b
# ╟─fa404ca1-049b-49fe-b894-0cba3080dc01
# ╟─c7ca281b-e980-44f1-8963-f2c04290225d
# ╟─d15931d9-2263-4b34-b972-6cb1c7a7ae1b
# ╠═32741bf0-ef1a-4ba1-9a48-2363bdb52e57
# ╟─a0fb1d2c-56a8-42e1-bfeb-aae5040df1cd
# ╠═29246f54-4220-421c-ae41-62d9593686cb
# ╠═e87db36e-70bf-4006-bd08-9eecb2f396d5
# ╟─5cdb2d51-257a-4f82-80c3-308a7d624d06
# ╟─ce96621b-bdf2-4a8c-b127-6306b395e5c9
# ╟─d22e249b-8bce-47d1-8dd9-dbb75fb15997
# ╟─e3bbb1f0-52f3-4279-b426-5c902305709b
# ╠═9560d375-ed9b-4ef6-8cf1-9f7b509f2af3
# ╠═2e8e2d60-ec3a-440a-b178-33f7c541dea7
# ╟─89ecf69e-3d2c-4aac-abea-556fbadba1d5
# ╠═7bc355fd-fb4d-4c18-ab2a-486b8ba21680
# ╠═16db0c80-a9e4-439d-8c8b-fbab797e9690
# ╠═7d66c043-19bf-4bbf-ab4f-8ab597dc4ee6
# ╠═0babe808-4ec9-424c-8ffa-9e686f0ca4d9
# ╠═1f759a4f-5dca-46d1-b8d1-482e7a610f79
# ╠═fd8d3155-226d-43e1-b545-35e520673e9b
# ╟─3f2c4518-a44b-4d79-bba0-51187c84dad9
# ╟─f08b27e9-5328-4361-90d9-9ac944b090c0
# ╠═5036ec33-5c94-4b5e-bcdf-089365166343
# ╠═a46cd130-d46a-4a3c-b99d-a8fe69640ea8
# ╟─da26766c-f908-41da-bde4-ee2899237dac
# ╠═0b5f1ed0-51c0-41d0-9946-9e3f6c91b524
# ╠═097625a8-43a1-44a7-86b3-04f0b718f505
# ╠═ec875061-2f1d-4d9b-8c0a-d6a6925dc4a4
# ╠═e8703636-bbfe-46e2-9191-bd70f45d4a57
# ╠═06fc9bd4-4f61-4c6f-b096-4fd94c7d26f3
# ╠═7345c8db-73ce-4a90-a1b2-9478ab4d8a41
# ╠═757777f2-889e-4156-afd7-79baf084b20d
# ╠═24d55bf6-aed0-4d87-9682-00eb828b19bd
# ╠═4960f6f2-7569-48f7-82f2-d5fd6af151ee
# ╠═5af99272-9b87-4b25-8eb4-e75e0a14cfc3
# ╠═24ae324c-8e00-4749-8440-1fb471ddf2b7
# ╠═ca29e38b-171b-4dd7-b7ec-b3d31f66a19a
# ╠═2585179f-14b1-4d28-81a3-f8e86b299f35
# ╠═fe6b9eff-0963-4008-abc5-2efb99748bb2
