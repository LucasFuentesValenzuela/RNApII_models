### A Pluto.jl notebook ###
# v0.19.26

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

# ╔═╡ 87cfce1e-2624-46d7-8540-0c552508d9c4
begin
	import Pkg
	Pkg.activate("../")
	using Revise
	using RNApIIModels
end

# ╔═╡ fde65485-580c-4aab-b2be-104f35ea3e53
using PlutoUI

# ╔═╡ 3d6e8c09-dc5c-47a0-83c1-83fbcbd7bf58
using JLD2

# ╔═╡ c83b06d1-d6fa-4f15-9515-183999289582
using Plots

# ╔═╡ 38975f19-7e86-491b-8794-b2e7841b1967
using Interpolations

# ╔═╡ d0a7578b-39c9-44c7-af78-2123e763e4c5
using CubicSplines

# ╔═╡ 683257e7-3140-4bda-b88a-c25e0401f773
using DataFrames

# ╔═╡ 4feb63d1-7a61-47bb-97d6-d1ff7c5aad87
using CSV

# ╔═╡ eee93bb1-22ef-45e8-9bbd-dd69ebb1ddf9
using TOML

# ╔═╡ 3f79c487-6617-43e2-b2cb-48d31da2f4dc
using Tables

# ╔═╡ 9581b94f-0c7d-4c4d-bfc6-2c9824ba3ac7
TableOfContents()

# ╔═╡ e45a3666-e398-4ec3-94dc-79c66f1e4fd0
config = TOML.parsefile("../config.toml")

# ╔═╡ e5de0f99-a58d-4952-824e-f6ada230bef8
PATH = config["DATADIR"]

# ╔═╡ 94cec09f-d40e-4a80-8cb8-e6a1c2c3b6ba
md"""# Description

Main analysis notebook for the paper. We want to show that feasible points (i.e. points that give an acceptable/within range) do behave as measured. 

The process is: 
- Find feasible points
- Calibrate the kon so that it matches occupancy for a 50 fL cell
- Show that the model applies throughout the board

Note: the simulations are not run on this notebook --- they should be run beforehand.
"""

# ╔═╡ b19fd665-236a-458f-ab90-e223d0cc2cda
md"""# Feasible points screen

Goal: determine what points are feasible from simulations. 
"""

# ╔═╡ 77fe4330-421f-463b-b21f-542a3aaa68c4
md"""
Select the value of Ω for which to load the dataset.

**Note**: the simulations should be run beforehand, using `scripts/feasible_analysis.jl`
"""

# ╔═╡ 629b3e35-3331-4bed-8186-541016e9fa58
@bind Ω_val Select([1.7, 1, 2, 3, 4], default=1.7)

# ╔═╡ 4bd04289-7f6c-4fa5-9e10-6de3b2040c23
fnm(x) = joinpath(config["RESULTSDIR"], "feasible_pts_$(x)_Omega$(Ω_val).jld2")

# ╔═╡ 35a0da09-d22b-4e08-82e3-7edb88e003a8
begin 
	fnm_screen = fnm("screen")
	results = JLD2.load(
		fnm_screen; 
		typemap=Dict("Main.Params" => RNApIIModels.Params)
	);

	# taking the first element in }the array to facilitate analysis
	occupancy = get_quantile(results["occupancy"], .5)
	promoter_occ = get_quantile(results["promoter_occ"], .5)
	params_iter = results["params_iter"]
end;

# ╔═╡ 287994aa-7329-4a18-a7ac-9556a485306c
feasible, feasible_pts = get_feasible_pts(fnm_screen);

# ╔═╡ 5c24344c-ff17-4ec2-84c8-f981b36f55e3
feasible_pts

# ╔═╡ 67783a0e-82fb-46b2-b69b-c36a619a89bb
begin
	α_vec_screen = [params_iter[k][1] for k in 1:length(params_iter)]
	kon_vec_screen = params_iter[1][3]
	heatmap(α_vec_screen, kon_vec_screen, feasible)
	plot!(xlabel="α", ylabel="k_on")
	hline!([RNApIIModels.min_k_on, RNApIIModels.max_k_on], label="limits kon", linewidth=2)
	vline!([RNApIIModels.min_α, RNApIIModels.max_α], label="limits α", linewidth=2)
	vline!([1/Ω_val], label="1/Ω", linewidth=4, color=:red)
	plot!(title="Ω = $Ω_val")
end

# ╔═╡ 46e6137c-cecd-4763-8f45-b1076ae36e13
@bind iter_nb Slider(1:length(params_iter))

# ╔═╡ a4b4db5a-483a-4af2-90f9-5ed46c1892b8
let
	plot()

	colors = palette([:orange, :forestgreen, :firebrick], length(occupancy))

	for (k, occ) in enumerate(occupancy)
		k_on_vec = params_iter[k][3]
		plot!(k_on_vec, occ, label="", color=colors[k])
		scatter!(k_on_vec, occ, label="", color=colors[k])
	end
	xlabel!("k_on")
	ylabel!("occupancy")
	plot!(xscale=:log)
	plot!(title="gene body occupancy")
	vline!([RNApIIModels.min_k_on, RNApIIModels.max_k_on], label="kon, average gene")
	hline!([RNApIIModels.min_ρ_g, RNApIIModels.max_ρ_g], label = "occupancy, average gene", linewidth=2)


	params_ = params_iter[iter_nb]
	α_ = params_[1]
	β_ = params_[2]
	konvec = params_[3]
	koff_ = max(1/Ω_val - α_, 0)

	α_eff = effective_α.(konvec, koff_, α_)

	ρ_th = map(f -> (ρ.(f, β_, OCCUPANCY_PARAMS["γ"], OCCUPANCY_PARAMS["L"]))[2], α_eff)


	plot!(konvec, ρ_th * OCCUPANCY_PARAMS["n_sites"], xscale=:log10, linestyle=:dash, linewidth=2, label="theory", color=:blue)
	
	plot!(legend=:topleft)
end

# ╔═╡ 3d778cb5-17dc-46bf-8523-e7cd220f276d
let
	plot()

	colors = palette([:orange, :forestgreen, :firebrick], length(promoter_occ))

	for (k, occ) in enumerate(promoter_occ)
		k_on_vec = params_iter[k][3]
		plot!(k_on_vec, occ, label="", color=colors[k])
		scatter!(k_on_vec, occ, label="", color=colors[k])
	end
	
	xlabel!("k_on")
	ylabel!("occupancy")
	plot!(xscale=:log10)
	plot!(title="Promoter occupancy")

	params_ = params_iter[iter_nb]
	α_ = params_[1]
	β_ = params_[2]
	konvec = params_[3]
	koff_ = max(1/Ω_val - α_, 0)


	ρ_p_th = RNApIIModels.ρp.(konvec, koff_, α_)

	plot!(konvec, ρ_p_th, xscale=:log10, linestyle=:dash, linewidth=2, label="theory", color=:blue)
	
	
	plot!(legend=:topleft)
end

# ╔═╡ c2a6b9a6-bab2-45cf-93ab-8c6751b0254c
md"""
# Average gene

Here we test whether the points identified as feasible behave, for the average gene, in accordance with the measurements.

We therefore use data in the `DATADIR` and compare it with simulation results
"""

# ╔═╡ d8069c70-d820-49f8-add6-9e85f5ca88da
RNA_free_avgCell = CV_to_RNAfree_interp(config["DATADIR"])(avg_cell_size);

# ╔═╡ e22c410c-eeb0-4e5d-b435-8d54943936d1
# loading the "narrow" dataset
begin
	fnm_narrow = fnm("narrow")

	results_nw = JLD2.load(
		fnm_narrow;
		typemap=Dict("Main.Params" => RNApIIModels.Params)
	);

	occupancy_nw = get_quantile(results_nw["occupancy"], .5)
	promoter_occ_nw = get_quantile(results_nw["promoter_occ"], .5)
	params_iter_nw = results_nw["params_iter"]

	q_up_occ_nw = get_quantile(results_nw["occupancy"], .9)
	q_dwn_occ_nw = get_quantile(results_nw["occupancy"], .1)

end;

# ╔═╡ d23a3327-4e8c-4fc5-aeb0-ef1c166bb1b5
# plotting the results for the feasible points
# also saving the data in a DataFrame
begin
	
	# plot the result
	plots_feasible = []
	sfs = []
	plots_feasible_kon = []

	df = load_ChIP_data(config["DATADIR"])

	dfs_avg_gene = []
	
	for idx_feasible_point in 1:length(feasible_pts)
		
		kon_vec_ = params_iter_nw[idx_feasible_point][3]
		k_crt = feasible_pts[idx_feasible_point][1] # kon at 50fL
		CV_crt = RNAfree_to_CV_interp(config["DATADIR"]).(kon_vec_ / k_crt * RNA_free_avgCell)
		
		q_up = q_up_occ_nw[idx_feasible_point]
		q_dwn = q_dwn_occ_nw[idx_feasible_point]

		mid = (q_up+q_dwn)/2
		w = q_up - q_dwn
		
		occ_crt = occupancy_nw[idx_feasible_point]


		
		###########
		#
		# Plotting occupancy as a function of kon
		# 
		###########

		p_crt = plot(CV_crt, mid, ribbon=w, fillalpha=.3, label="")
		plot!(CV_crt, occ_crt, label="")
		scatter!(CV_crt, occ_crt, label="Simulations", markerstrokewidth=0.5)
	
		occupancy_crt_interp = linear_interpolation(CV_crt, occ_crt)
		sf_crt = occupancy_crt_interp(RNApIIModels.avg_cell_size)/Rpb1_occupancy_haploid_interp(config["DATADIR"]).(RNApIIModels.avg_cell_size)

		push!(sfs, sf_crt)
		
		plot!(
				df[!, :cell_volume_fL], df[!, :Rpb1_occupancy_haploid_fit] .* sf_crt, 
				label="haploid occupancy", linewidth=2
		)
	
		hline!([RNApIIModels.min_ρ_g, RNApIIModels.max_ρ_g], linestyle=:dash, linewidth=2, label="average gene occupancy at 50fL")
		
		if idx_feasible_point == 1
			plot!(legend=:topright)
		else
			plot!(legend = false)
		end
		plot!(xlabel="Cell Size [fL]", ylabel="Gene body occupancy")
		plot!(ylims=(0 * RNApIIModels.min_ρ_g, 2 * RNApIIModels.max_ρ_g))
	
		push!(plots_feasible, p_crt)

		###########
		#
		# Same plot as above but with kon on the x axis (not cell size)
		# 
		###########

		p_kon = plot(kon_vec_, mid, ribbon=w, fillalpha=.3, label="")
		plot!(kon_vec_, occ_crt, label="")
		scatter!(kon_vec_, occ_crt, label="Simulations", markerstrokewidth=0.5)
		
		plot!(
				kon_vec_, Rpb1_occupancy_haploid_interp(config["DATADIR"]).(CV_crt) .* sf_crt, 
				label="haploid occupancy", linewidth=2
		)
	
		hline!([RNApIIModels.min_ρ_g, RNApIIModels.max_ρ_g], linestyle=:dash, linewidth=2, label="average gene occupancy at 50fL")
		
		if idx_feasible_point == 1
			plot!(legend=:topright)
		else
			plot!(legend = false)
		end
		plot!(xlabel="kon", ylabel="Gene body occupancy")
		plot!(ylims=(0 * RNApIIModels.min_ρ_g, 2 * RNApIIModels.max_ρ_g))
	
		push!(plots_feasible_kon, p_kon)

		
		df_avg_gene_crt = DataFrame(
			id = idx_feasible_point, 
			kon = kon_vec_, 
			cell_size = CV_crt, 
			occ90 = q_up, 
			occ10 = q_dwn, 
			occ50 = occ_crt, 
			dyn_eq_haploid = Rpb1_occupancy_haploid_interp(config["DATADIR"]).(CV_crt) .* sf_crt,
		)

		push!(dfs_avg_gene, df_avg_gene_crt)
	end
	
end

# ╔═╡ 01462590-ba62-4796-9d5c-8f0127de58d0
let
	p_all = plot((plots_feasible[1:10])..., layout = (5, 2))
	plot!(size=(1000, 1000))

	p_all
end

# ╔═╡ 121e3d11-9a12-44a6-9403-0ab7dde865fc
let
	p_all = plot((plots_feasible_kon[1:8])..., layout = (4, 2))
	plot!(size=(1000, 1000))

	p_all
end

# ╔═╡ 510d3054-85f1-4942-8068-467cdc3a7e4b
md"""
# Gene bins

Now that we have studied the behavior of the system for the average gene, we want to study it for genes that are more or less expressed. Mostly, we are interested in the behavior at higher expression levels, as we see hints of promoter saturation in the data. 
"""

# ╔═╡ 77919038-a1a2-4baa-be85-c71f2d598773
df_bins = load_gene_bins(config["DATADIR"]);

# ╔═╡ ac788865-7558-4736-8920-cbcd93efe83b
# loading the "wide" dataset
begin
	fnm_wd = fnm("wide")
	
	results_wd = JLD2.load(
		fnm_wd;
		typemap=Dict("Main.Params" => RNApIIModels.Params)
	);

	occupancy_wd = get_quantile(results_wd["occupancy"], .5)
	promoter_occ_wd = get_quantile(results_wd["promoter_occ"], .5)
	params_iter_wd = results_wd["params_iter"]

	q_up_occ_wd = get_quantile(results_wd["occupancy"], .9)
	q_dwn_occ_wd = get_quantile(results_wd["occupancy"], .1)

	q_up_prom_occ_wd = get_quantile(results_wd["promoter_occ"], .9)
	q_dwn_prom_occ_wd = get_quantile(results_wd["promoter_occ"], .1)

end;

# ╔═╡ f91bb98a-7442-40b9-bd24-07e554fb65b5
# kon fold changes
# these apply across the board (for any bin/kon) because it is simply a change in 
# [RNA]_free
begin
	kon_fc = CV_to_RNAfree_interp(config["DATADIR"]).(df_bins.cell_volume_fL) ./ CV_to_RNAfree_interp(config["DATADIR"]).(df_bins.cell_volume_fL[3]);
end

# ╔═╡ d92ad4b6-f36b-4ffb-b55f-7ff76d5c0c03
colors = palette([:purple, :orange, :green], size(df_bins, 2)-1);

# ╔═╡ 286b0fab-7a1c-48b8-99db-e885629d5651
# a linear interpolation is not really possible because the function is not monotonous (therefore, the inverse is not defined)
# we will bracket the value of occupancy and linearly interpolate the kon

function invert_occupancy(x, occ, kon)

	idxs = findall(occ .> x)
	occ_up, ii = findmin(occ[idxs])
	kon_up = kon[idxs][ii]

	
	idxs_dwn = findall(occ .< x)
	occ_dwn, jj = findmax(occ[idxs_dwn])
	kon_dwn = kon[idxs_dwn][jj]

	return kon_dwn .+ (x .- occ_dwn)/(occ_up-occ_dwn) * (kon_up - kon_dwn)

end

# ╔═╡ 3142e3a8-5b7c-46e4-9de1-4568cb001dfa
# run the analysis and the comparative plotting for all the feasible points
begin
	# params for the plotting
	ref_idx = 3 # idx of the reference cell size in the df_bins DataFrame
	ylim = (0.8, 1.7)
	
	plots_fold_changes = []
	plots_fold_changes_th = []
	kons_lims = []

	dfs_occ_model = []
	dfs_occ_model_fc = [] # with fold change instead of absolute occupancies
	
	for idx_ in 1:length(feasible_pts)

		α_wide = params_iter_wd[idx_][1]
		β_wide = params_iter_wd[idx_][2]
		kon_wide = params_iter_wd[idx_][3]
		
		occ_model(x) = (CubicSpline(kon_wide, occupancy_wd[idx_]))[x]
		
		occ_avg_gene = occ_model(feasible_pts[idx_][1])
		
		# rescale the occupancy data for each gene bin
		df_bins_ = copy(df_bins)
		df_bins_[:, 2:end] .= df_bins_[:, 2:end] .* occ_avg_gene

		# kon at 54 fL for each bin under the selected feasible point
		kons_54fL = map(
			f -> invert_occupancy(
				f, occupancy_wd[idx_], kon_wide
				), 
				collect(df_bins_[3, 2:end])
		)
		
		push!(
			kons_lims, [
			(kons_54fL[1] .* kon_fc[1], kons_54fL[1] .* kon_fc[end]), 
			(kons_54fL[end] .* kon_fc[1], kons_54fL[end] .* kon_fc[end])
		]
		)
	
		# DataFrame with the occupancy from the model
		df_kons = copy(df_bins_)
		df_occ_model = copy(df_bins_)
		df_occ_model_fc = copy(df_bins_)
		for bin_nbr in 1:(length(kons_54fL))
			kon_values = kons_54fL[bin_nbr] .* kon_fc
			df_kons[:, bin_nbr+1] = kon_values
			occ_values = occ_model.(kon_values)
			df_occ_model[:, bin_nbr+1] = occ_values
			df_occ_model_fc[:, bin_nbr+1] = occ_values ./ occ_values[3]
		end

		df_occ_model[!, :id] .= idx_
		df_occ_model_fc[!, :id] .= idx_
		
		push!(dfs_occ_model, df_occ_model)
		push!(dfs_occ_model_fc, df_occ_model_fc)

		
	
		#######
		#
		# plotting 
		# 
		#######
		
	
		p_crt = plot()
		for (k, bin) in enumerate(names(df_bins)[2:end])
			plot!(
				df_occ_model[:, "cell_volume_fL"], 
				df_occ_model[:, bin] / df_occ_model[ref_idx, bin], 
				label="", color=colors[k]
			)
		end
		plot!(xlabel="Cell Volume [fL]", ylabel="Occupancy fold change")
		plot!(ylim=ylim)
		plot!(title="Model")
	
	
		push!(plots_fold_changes, p_crt)

		p_th_crt = plot()
		for (k, bin) in enumerate(names(df_bins)[2:end])
			kons = df_kons[:, bin]
			α_eff = effective_α.(kons, 1/Ω_val - α_wide, α_wide)

			ρs = [(ρ.(
				α, β_wide, OCCUPANCY_PARAMS["γ"], OCCUPANCY_PARAMS["L"])
			)[2] for α in α_eff]
			
			plot!(
				df_kons[:, "cell_volume_fL"], 
				ρs ./ ρs[3], 
				label="", 
				color = colors[k]
			)
		end
		plot!(xlabel="Cell Size [fL]", ylabel="Occupancy fold change", title="Theory")
		plot!(ylim=ylim)

		push!(plots_fold_changes_th, p_th_crt)
	
	
	end
end

# ╔═╡ abddf0af-85e8-44f3-8e68-5f0c011cb3b8
begin
	p_fc_data = plot()
	for (k, bin) in enumerate(names(df_bins)[2:end])
		plot!(df_bins[:, "cell_volume_fL"], df_bins[:, bin] / df_bins[ref_idx, bin], label="", color=colors[k])
	end
	plot!(xlabel="Cell Volume [fL]", ylabel="Occupancy fold change")
	plot!(ylim=ylim)
	plot!(title="Data")
end;

# ╔═╡ 9d6f40f3-cb34-42ad-b97f-f185ea45fada
@bind idx_p Slider(1:(length(plots_fold_changes)))

# ╔═╡ cd113d08-d1bd-4b7b-a2b9-3441c54dbdbb
md"""
**Note**: In the below plot: 
- sim: results of simulation
- theory: results of theory for a particle of size 1
- MC: correction of theory with the Markov Chain approach
- L: results of theory for a particle of size 35 (actual footprint)
"""

# ╔═╡ 402294f8-2ccc-46d8-86c8-1778d679d1bf
let

	normalize_occ = true
	n_sites = OCCUPANCY_PARAMS["n_sites"]
	
	ε = 1e-4

	kon_wide = params_iter_wd[idx_p][3]
	params_ = params_iter_wd[idx_p]
	α_ = params_[1]
	β_ = params_[2]
	konvec = params_[3]
	koff_ = max(1/Ω_val - α_, 0)
	α_eff = effective_α.(konvec, koff_, α_)
	ρ_th = map(f -> (ρ.(f, β_, OCCUPANCY_PARAMS["γ"], OCCUPANCY_PARAMS["L"]))[2], α_eff)
	ρp_th = ρp.(konvec, koff_, α_)

	ρp_MC_ = ρp_MC.(α_, β_, konvec, koff_; order=4)
	α_eff_MC = α_ * ρp_MC_
	ρ_MC = map(f -> (ρ.(f, β_, OCCUPANCY_PARAMS["γ"], OCCUPANCY_PARAMS["L"]))[2], α_eff_MC)

	δ = OCCUPANCY_PARAMS["δ"]
	ρ_L = map(f -> (ρ.(f, β_ * δ, OCCUPANCY_PARAMS["γ"], OCCUPANCY_PARAMS["L"] * δ))[2], α_eff) .* δ


	q_up = q_up_occ_wd[idx_p]
	q_dwn = q_dwn_occ_wd[idx_p]
	mid = (q_up+q_dwn)/2
	w = q_up - q_dwn

	q_prom_up = q_up_prom_occ_wd[idx_p]
	q_prom_dwn = q_dwn_prom_occ_wd[idx_p]
	mid_prom = (q_prom_up + q_prom_dwn)/2
	w_prom = (q_prom_up - q_prom_dwn)
	
	p0 = plots_fold_changes[idx_p]

	p2 = plot(kon_wide, occupancy_wd[idx_p] ./ n_sites, label="")
	plot!(kon_wide, mid ./ n_sites, ribbon=w ./ n_sites, fillalpha=.3, label="")
	scatter!(kon_wide, occupancy_wd[idx_p] ./ n_sites, label="gene body (sim)")
	plot!(konvec, ρ_th, linestyle=:dash, linewidth=3, label="gene body (theory)", color=:firebrick)
	plot!(konvec, ρ_MC, linestyle=:dot, linewidth=2, label="gene body (MC)")
	plot!(konvec, ρ_L, linewidth=2, label="gene body (L)")
	plot!(
		xlabel="kon", ylabel="normalized occupancy", 
		# xscale=:log10,
		# yscale=:log10
		legend=:bottomright
	)
	
	vline!([kons_lims[idx_p][1][1], kons_lims[idx_p][1][2]], label="")
	vline!([kons_lims[idx_p][2][1], kons_lims[idx_p][2][2]], label="")
	
	# plot!(legend=:topleft)
	
	p3 = plot(kon_wide, promoter_occ_wd[idx_p], label="")
	plot!(kon_wide, mid_prom, ribbon=w_prom, fillalpha=.3, label="")
	scatter!(kon_wide, promoter_occ_wd[idx_p], label="promoter (sim)")
	plot!(konvec, ρp_th, linestyle=:dash, linewidth=3, label="promoter (theory)")
	plot!(konvec, ρp_MC_, linestyle=:dot, linewidth=2, label="promoter (Monte Carlo correction)")
	
	plot!(
		xlabel="kon", ylabel="normalized occupancy", 
		xscale=:log10,
		# yscale=:log10
	)
	
	vline!([kons_lims[idx_p][1][1], kons_lims[idx_p][1][2]], label="")
	vline!([kons_lims[idx_p][2][1], kons_lims[idx_p][2][2]], label="")
	
	# plot!(legend=:topleft)
	
	err_gene = abs.((occupancy_wd[idx_p] ./ n_sites .- ρ_MC) ./ ρ_MC)
	err_prom = abs.((promoter_occ_wd[idx_p] .- ρp_MC_) ./ ρp_MC_)
	
	err_plot = plot()
	plot!(kon_wide, err_gene, label="")
	scatter!(kon_wide, err_gene, label="gene")
	plot!(kon_wide, err_prom, label="")
	scatter!(kon_wide, err_prom, label="promoter")
	plot!(xlabel="kon", ylabel="relative diff.")
	plot!(legend=:topleft)
	plot!(xscale=:log10)
	plot!(ylim=(0, .2))
	
	q=plot()
	for kk in 1:length(feasible_pts)
		if kk == idx_p
			c = :blue
			lw = 2
			
		else
			c = :gray
			lw = 1
		end
		
		plot!(kon_wide, max.(occupancy_wd[kk], ε), label="", color=c, linewidth=lw)
	end
	plot!(kon_wide, mid, ribbon=w, fillalpha=.3, label="")

	# plot the theory

	plot!(konvec, ρ_th * n_sites, xscale=:log10, linestyle=:dash, linewidth=3, label="theory", color=:firebrick)
	
	plot!(xscale=:log10, yscale=:log10)
	plot!(xlabel="kon", ylabel="occupancy")
	plot!(ylim=(1e-3, 1e1))
	plot!(legend=:topleft)
	vline!([kons_lims[idx_p][1][1], kons_lims[idx_p][1][2]], label="")
	vline!([kons_lims[idx_p][2][1], kons_lims[idx_p][2][2]], label="")
	q
	
	plot([
		p_fc_data, p0, 
		plots_fold_changes_th[idx_p], p2, 
		err_plot, q
	]..., layout=(3, 2), size=(700, 700))	
end

# ╔═╡ 4b82c6d9-a148-4004-8451-58a80f99840f
let
	p_fc_all = plot((plots_fold_changes[1:8])..., layout = (4, 2))
	plot!(size=(1000, 1000))

	p_fc_all
end

# ╔═╡ 362c3c71-7cc6-4d56-bf4f-d2ae3cdfd486
md"""
# File outputs

We save the results as files. 
"""

# ╔═╡ 7648928d-5707-4140-9630-e7655aee6ba0
feasible_mat = hcat(
	vcat(
		["kon|alpha"], reshape(kon_vec_screen, (10, 1))
	), vcat(
		reshape(α_vec_screen, (1, 8)), feasible
	)
);

# ╔═╡ ad148034-87f3-4235-b6e8-46326f9328fb
CSV.write(
	joinpath(PATH, "data_out", "feasible_grid_Ω$(Ω_val).csv"), Tables.table(feasible_mat), writeheader=false
)

# ╔═╡ 1959a44a-23bb-435a-9470-392437fb94ed
# outputting the feasible points and their statistics
begin
	# output cells
	df_feasible_stats = DataFrame(
		id = [k for k in 1:length(feasible_pts)],
		kon = [pt[1] for pt in feasible_pts], 
		α = [pt[2] for pt in feasible_pts], 
		koff = [max.(1/Ω_val - pt[2], 0) for pt in feasible_pts],
		gene_body_occupancy = [
			occupancy[idx_α][idx_k] for (idx_k, idx_α) in Tuple.(findall(feasible .== 1))
		], 
		promoter_occupancy = [
			promoter_occ[idx_α][idx_k] for (idx_k, idx_α) in Tuple.(findall(feasible .== 1))
		],
	)
	
	CSV.write(
		joinpath(PATH, "data_out", "feasible_pts_Ω$(Ω_val).csv"), df_feasible_stats
	)
end

# ╔═╡ 99e21b1a-47e0-4b77-b799-0b9e7cf993f1
begin
	df_avg_gene = vcat(
		dfs_avg_gene...
	)
	
	CSV.write(
			joinpath(PATH, "data_out", "avg_gene_scaling_Ω$(Ω_val).csv"), df_avg_gene
	)
end

# ╔═╡ 93824239-5784-4c5b-9233-c66efb52611a
begin
	# write the occupancy model dataset
	
	dfs_occ_models = vcat(dfs_occ_model...)
	dfs_occ_models = dfs_occ_models[!, vcat([52], collect(1:51))]
	CSV.write(
		joinpath(PATH, "data_out", "occupancy_bins_Ω$(Ω_val).csv"),
		dfs_occ_models
		
	)

	dfs_occ_models_fc = vcat(dfs_occ_model_fc...)
	dfs_occ_models_fc = dfs_occ_models_fc[!, vcat([52], collect(1:51))]
	CSV.write(
		joinpath(PATH, "data_out", "occupancy_bins_foldChanges_Ω$(Ω_val).csv"),
		dfs_occ_models_fc
		
	)
	
end

# ╔═╡ Cell order:
# ╠═87cfce1e-2624-46d7-8540-0c552508d9c4
# ╠═fde65485-580c-4aab-b2be-104f35ea3e53
# ╠═3d6e8c09-dc5c-47a0-83c1-83fbcbd7bf58
# ╠═c83b06d1-d6fa-4f15-9515-183999289582
# ╠═38975f19-7e86-491b-8794-b2e7841b1967
# ╠═d0a7578b-39c9-44c7-af78-2123e763e4c5
# ╠═683257e7-3140-4bda-b88a-c25e0401f773
# ╠═4feb63d1-7a61-47bb-97d6-d1ff7c5aad87
# ╠═eee93bb1-22ef-45e8-9bbd-dd69ebb1ddf9
# ╠═3f79c487-6617-43e2-b2cb-48d31da2f4dc
# ╠═9581b94f-0c7d-4c4d-bfc6-2c9824ba3ac7
# ╠═e45a3666-e398-4ec3-94dc-79c66f1e4fd0
# ╠═e5de0f99-a58d-4952-824e-f6ada230bef8
# ╟─94cec09f-d40e-4a80-8cb8-e6a1c2c3b6ba
# ╟─b19fd665-236a-458f-ab90-e223d0cc2cda
# ╟─77fe4330-421f-463b-b21f-542a3aaa68c4
# ╠═629b3e35-3331-4bed-8186-541016e9fa58
# ╠═4bd04289-7f6c-4fa5-9e10-6de3b2040c23
# ╠═35a0da09-d22b-4e08-82e3-7edb88e003a8
# ╠═287994aa-7329-4a18-a7ac-9556a485306c
# ╠═5c24344c-ff17-4ec2-84c8-f981b36f55e3
# ╠═67783a0e-82fb-46b2-b69b-c36a619a89bb
# ╠═46e6137c-cecd-4763-8f45-b1076ae36e13
# ╠═a4b4db5a-483a-4af2-90f9-5ed46c1892b8
# ╠═3d778cb5-17dc-46bf-8523-e7cd220f276d
# ╟─c2a6b9a6-bab2-45cf-93ab-8c6751b0254c
# ╠═d8069c70-d820-49f8-add6-9e85f5ca88da
# ╠═e22c410c-eeb0-4e5d-b435-8d54943936d1
# ╟─d23a3327-4e8c-4fc5-aeb0-ef1c166bb1b5
# ╠═01462590-ba62-4796-9d5c-8f0127de58d0
# ╠═121e3d11-9a12-44a6-9403-0ab7dde865fc
# ╟─510d3054-85f1-4942-8068-467cdc3a7e4b
# ╠═77919038-a1a2-4baa-be85-c71f2d598773
# ╠═ac788865-7558-4736-8920-cbcd93efe83b
# ╠═f91bb98a-7442-40b9-bd24-07e554fb65b5
# ╠═d92ad4b6-f36b-4ffb-b55f-7ff76d5c0c03
# ╠═abddf0af-85e8-44f3-8e68-5f0c011cb3b8
# ╠═3142e3a8-5b7c-46e4-9de1-4568cb001dfa
# ╠═286b0fab-7a1c-48b8-99db-e885629d5651
# ╟─9d6f40f3-cb34-42ad-b97f-f185ea45fada
# ╟─cd113d08-d1bd-4b7b-a2b9-3441c54dbdbb
# ╟─402294f8-2ccc-46d8-86c8-1778d679d1bf
# ╠═4b82c6d9-a148-4004-8451-58a80f99840f
# ╟─362c3c71-7cc6-4d56-bf4f-d2ae3cdfd486
# ╠═7648928d-5707-4140-9630-e7655aee6ba0
# ╠═ad148034-87f3-4235-b6e8-46326f9328fb
# ╠═1959a44a-23bb-435a-9470-392437fb94ed
# ╠═99e21b1a-47e0-4b77-b799-0b9e7cf993f1
# ╠═93824239-5784-4c5b-9233-c66efb52611a
