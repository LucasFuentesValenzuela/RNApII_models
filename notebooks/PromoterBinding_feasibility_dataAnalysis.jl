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

# ╔═╡ 87cfce1e-2624-46d7-8540-0c552508d9c4
begin
	using Revise
	import Pkg
	Pkg.activate()
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

# ╔═╡ ec58d9c0-ae9c-43cc-b2a8-bf9182204ebe
using Statistics

# ╔═╡ 9581b94f-0c7d-4c4d-bfc6-2c9824ba3ac7
TableOfContents()

# ╔═╡ 94cec09f-d40e-4a80-8cb8-e6a1c2c3b6ba
md"""# Description

This notebook is pretty much the same as `PromoterBinding.jl`. However, we flip it on its head. We want to show that feasible points (i.e. points that give an acceptable/within range) value of promoter and gene body occupancy. 

So the process is: 
- find feasible points
- calibrate the kon so that it matches occupancy for a 50 fL cell
- show that the model applies throughout the board


So it should contain two steps. 
Step 1. Calibrate for the "average gene in the average cell". 
Step 2. Extend it and use those parameters to run on large variations of kon. 
"""

# ╔═╡ d0a5c973-a446-4ce9-9151-d903295200a5
md"""
# TODO
- add percentiles (based on data from Matt)
- interpolate on the kon values? Now I am doing it on a narrow range, and therefore the discrepancies between two simulations seem to arise more quickly and notably. I could do it on a larger range and then select the values in between that correspond to desired cell sizes
- I should be able to run the simulations once and for all for all points that I started with and then select the plots based on those simulations. 

- visualize the promoter and gene occupancy as a function of kon and alpha

- careful with some parameter combinations where koff is negative because alpha is too large!

- build the same analysis from the theoretical perspective (feasible points, etc.)

- There is disagreement with the model when koff is negative ==> this means the constraint we set on the duration on the promoter cannot be guaranteed -- we need to see how to solve this. 
"""

# ╔═╡ e5de0f99-a58d-4952-824e-f6ada230bef8
PATH = "/Users/lucasfuentes/RNApII_models"

# ╔═╡ 4bd04289-7f6c-4fa5-9e10-6de3b2040c23
fnm(x) = joinpath(PATH, "results", "feasible_pts_$(x).jld2")

# ╔═╡ 35a0da09-d22b-4e08-82e3-7edb88e003a8
begin 
	fnm_screen = fnm("screen")
	results = JLD2.load(
		fnm_screen; 
		typemap=Dict("Main.Params" => RNApIIModels.Params)
	);

	# taking the first element in the array to facilitate analysis
	occupancy = get_quantile(results["occupancy"], .5)
	promoter_occ = get_quantile(results["promoter_occ"], .5)
	params_iter = results["params_iter"]
end;

# ╔═╡ 287994aa-7329-4a18-a7ac-9556a485306c
feasible, feasible_pts = get_feasible_pts(fnm_screen);

# ╔═╡ 67783a0e-82fb-46b2-b69b-c36a619a89bb
begin
	# heatmap(α_vec, k_on_vec, feasible)
	heatmap(feasible)
	# plot!(xlabel="α", ylabel="k_on")
	# hline!([ps.min_k_on, ps.max_k_on], label="limits kon", linewidth=2)
	# vline!([ps.min_α, ps.max_α], label="limits α", linewidth=2)
end

# ╔═╡ 46e6137c-cecd-4763-8f45-b1076ae36e13
@bind iter_nb Slider(1:length(params_iter))

# ╔═╡ b12bfc38-3f7f-40d2-ac84-a140cd8ddb13
md"""
*TODO* indicate which ones are the feasible points
"""

# ╔═╡ 09916819-dbf6-4b04-bf69-950a0c7e136b
md"""
At higher values of α the agreement with theory is not good... why? 
"""

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
	koff_ = max(1/OCCUPANCY_PARAMS["Ω"] - α_, 0)

	α_eff = effective_α.(konvec, koff_, α_)

	ρ_th = map(f -> (ρ.(f, β_, OCCUPANCY_PARAMS["γ"], OCCUPANCY_PARAMS["L"]))[2], α_eff)


	plot!(konvec, ρ_th * OCCUPANCY_PARAMS["n_sites"], xscale=:log10, linestyle=:dash, linewidth=2, label="theory", color=:blue)
	
	plot!(legend=:topleft)
end

# ╔═╡ 3d778cb5-17dc-46bf-8523-e7cd220f276d
let
	plot()

	colors = palette([:orange, :forestgreen, :firebrick], length(promoter_occ))
	# colors = [:orange]

	for (k, occ) in enumerate(promoter_occ)
		k_on_vec = params_iter[k][3]
		plot!(k_on_vec, occ, label="", color=colors[k])
		scatter!(k_on_vec, occ, label="", color=colors[k])
	end

	# plot!(
	# 	k_on_vec, occupancy_interp(k_on_vec), 
	# 	linestyle=:dash, linewidth=3, 
	# 	color=:blue
	# )
	xlabel!("k_on")
	ylabel!("occupancy")
	plot!(xscale=:log10, yscale=:log10)
	plot!(title="Promoter occupancy")
	# hline!([ps.min_ρ_p, ps.max_ρ_p], label = "occupancy, average gene", linewidth=2)
	# vline!([min_α, max_α], label="α, average gene")
	# hline!([min_ρ_g, max_ρ_g], label = "occupancy, average gene")
	# vline!([α], label="α")

	params_ = params_iter[iter_nb]
	α_ = params_[1]
	β_ = params_[2]
	konvec = params_[3]
	koff_ = max(1/OCCUPANCY_PARAMS["Ω"] - α_, 0)


	ρ_p_th = RNApIIModels.ρp.(konvec, koff_, α_)

	plot!(konvec, ρ_p_th, xscale=:log10, linestyle=:dash, linewidth=2, label="theory", color=:blue)
	
	
	plot!(legend=:topleft)
end

# ╔═╡ c2a6b9a6-bab2-45cf-93ab-8c6751b0254c
md"""
# Average gene

Remember, we compute all this for the average gene, average cell. 
- Average cell: 50fL. This is where we extract the ChIP validation data. 
- Average gene: matches the occupancy metrics that we have in the document. 
"""

# ╔═╡ d8069c70-d820-49f8-add6-9e85f5ca88da
RNA_free_avgCell = CV_to_RNAfree_interp()(avg_cell_size);

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
begin
	# plot the result
	plots_feasible = []
	sfs = []
	plots_feasible_kon = []

	df = load_ChIP_data()
	
	for idx_feasible_point in 1:length(params_iter_nw)
		kon_vec_ = params_iter_nw[idx_feasible_point][3]
		k_crt = feasible_pts[idx_feasible_point][1] # kon at 50fL
		CV_crt = RNAfree_to_CV_interp().(kon_vec_ / k_crt * RNA_free_avgCell)
		
		q_up = q_up_occ_nw[idx_feasible_point]
		q_dwn = q_dwn_occ_nw[idx_feasible_point]

		mid = (q_up+q_dwn)/2
		w = q_up - q_dwn
		
		occ_crt = occupancy_nw[idx_feasible_point]

		p_crt = plot(CV_crt, mid, ribbon=w, fillalpha=.3, label="")
		plot!(CV_crt, occ_crt, label="")
		scatter!(CV_crt, occ_crt, label="Simulations", markerstrokewidth=0.5)
	
		occupancy_crt_interp = linear_interpolation(CV_crt, occ_crt)
		sf_crt = occupancy_crt_interp(RNApIIModels.avg_cell_size)/Rpb1_occupancy_haploid_interp().(RNApIIModels.avg_cell_size)

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
				kon_vec_, Rpb1_occupancy_haploid_interp().(CV_crt) .* sf_crt, 
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
	end
	
end

# ╔═╡ 01462590-ba62-4796-9d5c-8f0127de58d0
let
	p_all = plot((plots_feasible[1:10])..., layout = (5, 2))
	plot!(size=(1000, 1000))
	savefig(joinpath(PATH, "figs", "feasible_points.png"))

	p_all
end

# ╔═╡ 121e3d11-9a12-44a6-9403-0ab7dde865fc
let
	p_all = plot((plots_feasible_kon[1:10])..., layout = (5, 2))
	plot!(size=(1000, 1000))
	savefig(joinpath(PATH, "figs", "feasible_points_kon.png"))

	p_all
end

# ╔═╡ 510d3054-85f1-4942-8068-467cdc3a7e4b
md"""
# Gene bins
"""

# ╔═╡ b7dcd8cf-f775-4e9a-ae65-d244623d1213
md"""
We have occupancy per gene bin. Basically, the genes are filtered by values of kon, and we have the scaling of kon as a function of size.
"""

# ╔═╡ 77919038-a1a2-4baa-be85-c71f2d598773
df_bins = load_gene_bins();

# ╔═╡ de2f8fde-98a1-4c9b-8317-86638fe4485b
md"""
Normalization: average gene at 50 fL = 1 (occupancy). 

What we can do: 
- scaling of occupancy at 50fL (how the different bins scale --> representative of how the kon changes with genes)
- you are expecting the change to follow the change from the model, 
"""

# ╔═╡ aef445ad-08ec-4b8e-8ac5-9a6f1d217cf2
md"""
The model gives us how the average gene in the average cell (50fL) scales in occupancy with size. We want something similar for all other genes. 
"""

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

end;

# ╔═╡ 884c8cb2-8863-4d4b-a9ab-d23e0be85814
md"""## rescaling the occupancies
we have computed the occupancy profile for a series of combinations of the parameters that are feasible. Let us try to use that for now"""

# ╔═╡ 680d1245-c19e-4ac5-ab47-c57abd2d2f83
# we load the new simulation results, that span a larger range of kon values
results_wide = JLD2.load("results/feasible_points_wide_test.jld2");

# ╔═╡ 4bdfc70d-01d9-4ebf-8001-1e971fd644dd
n_feas = sum(feasible .== 1)

# ╔═╡ 280e88ac-fe5b-46b6-915f-b4a30cd527c0
begin
	# occ_wide = results_wide["occupancy_feasible"]
	params_iter_wide = results_wide["params_iter_feasible"]
end;

# ╔═╡ ae1d9063-91ac-4205-84eb-fa5e282dfcb4
occupancy_wide = reshape(hcat(hcat(results_wide["occupancy_feasible"]...)...), 10, n_feas, :);

# ╔═╡ 4982d64c-f6ad-4029-b25d-7dec59b67fd1
occ_wide = [vec(median(occupancy_wide[:, ii, :], dims=2)) for ii in 1:n_feas];

# ╔═╡ f91bb98a-7442-40b9-bd24-07e554fb65b5
# kon fold changes
# these apply across the board (for any bin/kon) because it is simply a change in 
# [RNA]_free
kon_fc = data.CV_to_RNAfree_interp(df_bins.cell_volume_fL) ./ data.CV_to_RNAfree_interp(df_bins.cell_volume_fL[3])

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
begin
	# run the analysis and the comparative plotting for all the feasible points
	
	# params for the plotting
	ref_idx = 3
	ylim = (0.8, 1.7)
	
	colors = palette([:purple, :orange, :green], size(df_bins, 2)-1)
	
	plots_fold_changes = []

	p1 = plot()
	for (k, bin) in enumerate(names(df_bins)[2:end])
		plot!(df_bins[:, "cell_volume_fL"], df_bins[:, bin] / df_bins[ref_idx, bin], label="", color=colors[k])
	end
	plot!(xlabel="Cell Volume [fL]", ylabel="Occupancy fold change")
	plot!(ylim=ylim)
	plot!(title="Data")

	push!(plots_fold_changes, p1)
	
	
	for idx_ in 1:n_feas
		# @show idx_

		kon_wide = params_iter_wide[idx_][3]
		occ_model = linear_interpolation(kon_wide, occ_wide[idx_])
		
		occ_avg_gene = occ_model(kon_wide[1] / 10^-1.5)
		# occ_avg_gene = data.Rpb1_occupancy_haploid_interp(ps.avg_cell_size) * sfs[idx_]
		
		# rescale the occupancy data for each gene bin
		df_bins_ = copy(df_bins)
		df_bins_[:, 2:end] .= df_bins_[:, 2:end] .* occ_avg_gene
	
		

		# @show minimum(df_bins_[3, 2:end]), maximum(df_bins_[3, 2:end])
		# @show minimum(occ_wide[idx_]), maximum(occ_wide[idx_])
		# kon at 54 fL for each bin under the selected feasible point
		kons_54fL = map(
			f -> invert_occupancy(
				f, occ_wide[idx_], kon_wide
				), 
				collect(df_bins_[3, 2:end])
		)
	
		# DataFrame with the occupancy from the model
		df_occ_model = copy(df_bins_)
		for bin_nbr in 1:(length(kons_54fL)-1)
			# @show bin_nbr
			kon_values = kons_54fL[bin_nbr] .* kon_fc
			# @show minimum(kon_values), maximum(kon_values)
			# @show minimum(kon_wide), maximum(kon_wide)
			occ_values = occ_model.(kon_values)
			df_occ_model[:, bin_nbr+1] = occ_values
		end
	
	
		#######
		#
		# plotting 
		# 
		#######
		
	
		p_crt = plot()
		for (k, bin) in enumerate(names(df_bins)[2:end])
			plot!(df_occ_model[:, "cell_volume_fL"], df_occ_model[:, bin] / df_occ_model[ref_idx, bin], label="", color=colors[k])
		# scatter!(df_bins[:, "cell_volume_fL"], df_bins[:, bin] / df_bins[3, bin], label="")
		end
		plot!(xlabel="Cell Volume [fL]", ylabel="Occupancy fold change")
		plot!(ylim=ylim)
		plot!(title="Model")
	
	
		push!(plots_fold_changes, p_crt)
	
	
	end
end

# ╔═╡ 9d6f40f3-cb34-42ad-b97f-f185ea45fada
@bind idx_p Slider(1:(length(plots_fold_changes)-1))

# ╔═╡ e59cc0b2-2b5e-405b-9fc7-cc6e69952817
# idx_p = 1

# ╔═╡ 402294f8-2ccc-46d8-86c8-1778d679d1bf
let

	ε = 1e-4

	kon_wide = params_iter_wide[idx_p][3]
	
	p0 = plots_fold_changes[idx_p+1]
	
	p2 = plot(kon_wide, occ_wide[idx_p], label="", xscale=:log10)
	scatter!(kon_wide, occ_wide[idx_p], label="", xscale=:log10)
	plot!(xlabel="kon", ylabel="occupancy")
	vline!([kon_wide[1] * kon_fc[1], kon_wide[1] * kon_fc[end]], label="")
	vline!([kon_wide[end] * kon_fc[1], kon_wide[end] * kon_fc[end]], label="")


	q=plot()
	for kk in 1:13
		if kk == idx_p
			c = :blue
			lw = 2

			q_up = [
				quantile(occupancy_wide[k, kk, :], .7) for k in 1:size(occupancy_wide, 1)
			]
			q_down = [
				quantile(occupancy_wide[k, kk, :], .3) for k in 1:size(occupancy_wide, 1)
			]


			mid = (q_up+q_down)/2
			w = q_up - q_down

			plot!(params_iter_wide[kk][3], mid, ribbon=w, fillalpha=.3, label="")
		else
			c = :gray
			lw = 1
		end
		
		plot!(params_iter_wide[kk][3], max.(occ_wide[kk], ε), label="", color=c, linewidth=lw)
	end
	plot!(xscale=:log10, yscale=:log10)
	plot!(xlabel="kon", ylabel="occupancy")
	plot!(ylim=(1e-3, 1e1))
	q
	
	plot([p1, p0, p2, q]..., layout=(2, 2), size=(1000, 1000))	
end

# ╔═╡ 1cc3f01f-e036-4b8b-b726-734ddd3f217d
md"""Question: why are the vertical lines sometimes out of the range? how is it even eable to interpolate in that range then?"""

# ╔═╡ 4b82c6d9-a148-4004-8451-58a80f99840f
let
	p_fc_all = plot((plots_fold_changes[1:10])..., layout = (5, 2))
	plot!(size=(1000, 1000))
	savefig("figs/fold_changes.png")

	p_fc_all
end

# ╔═╡ 5a765471-87f2-406b-ace1-ee2b51cf35b6
md"""TODO: promoter occupancy plots for the `wide` dataset"""

# ╔═╡ Cell order:
# ╠═fde65485-580c-4aab-b2be-104f35ea3e53
# ╠═3d6e8c09-dc5c-47a0-83c1-83fbcbd7bf58
# ╠═c83b06d1-d6fa-4f15-9515-183999289582
# ╠═38975f19-7e86-491b-8794-b2e7841b1967
# ╠═9581b94f-0c7d-4c4d-bfc6-2c9824ba3ac7
# ╠═87cfce1e-2624-46d7-8540-0c552508d9c4
# ╟─94cec09f-d40e-4a80-8cb8-e6a1c2c3b6ba
# ╠═d0a5c973-a446-4ce9-9151-d903295200a5
# ╠═e5de0f99-a58d-4952-824e-f6ada230bef8
# ╠═4bd04289-7f6c-4fa5-9e10-6de3b2040c23
# ╠═35a0da09-d22b-4e08-82e3-7edb88e003a8
# ╠═ec58d9c0-ae9c-43cc-b2a8-bf9182204ebe
# ╠═67783a0e-82fb-46b2-b69b-c36a619a89bb
# ╠═287994aa-7329-4a18-a7ac-9556a485306c
# ╟─46e6137c-cecd-4763-8f45-b1076ae36e13
# ╠═b12bfc38-3f7f-40d2-ac84-a140cd8ddb13
# ╠═09916819-dbf6-4b04-bf69-950a0c7e136b
# ╟─a4b4db5a-483a-4af2-90f9-5ed46c1892b8
# ╠═3d778cb5-17dc-46bf-8523-e7cd220f276d
# ╟─c2a6b9a6-bab2-45cf-93ab-8c6751b0254c
# ╠═d8069c70-d820-49f8-add6-9e85f5ca88da
# ╠═e22c410c-eeb0-4e5d-b435-8d54943936d1
# ╠═01462590-ba62-4796-9d5c-8f0127de58d0
# ╟─d23a3327-4e8c-4fc5-aeb0-ef1c166bb1b5
# ╠═121e3d11-9a12-44a6-9403-0ab7dde865fc
# ╟─510d3054-85f1-4942-8068-467cdc3a7e4b
# ╟─b7dcd8cf-f775-4e9a-ae65-d244623d1213
# ╠═77919038-a1a2-4baa-be85-c71f2d598773
# ╠═de2f8fde-98a1-4c9b-8317-86638fe4485b
# ╠═aef445ad-08ec-4b8e-8ac5-9a6f1d217cf2
# ╠═ac788865-7558-4736-8920-cbcd93efe83b
# ╟─884c8cb2-8863-4d4b-a9ab-d23e0be85814
# ╠═680d1245-c19e-4ac5-ab47-c57abd2d2f83
# ╠═4bdfc70d-01d9-4ebf-8001-1e971fd644dd
# ╠═280e88ac-fe5b-46b6-915f-b4a30cd527c0
# ╠═ae1d9063-91ac-4205-84eb-fa5e282dfcb4
# ╠═4982d64c-f6ad-4029-b25d-7dec59b67fd1
# ╠═f91bb98a-7442-40b9-bd24-07e554fb65b5
# ╠═3142e3a8-5b7c-46e4-9de1-4568cb001dfa
# ╠═286b0fab-7a1c-48b8-99db-e885629d5651
# ╠═9d6f40f3-cb34-42ad-b97f-f185ea45fada
# ╠═e59cc0b2-2b5e-405b-9fc7-cc6e69952817
# ╠═402294f8-2ccc-46d8-86c8-1778d679d1bf
# ╟─1cc3f01f-e036-4b8b-b726-734ddd3f217d
# ╠═4b82c6d9-a148-4004-8451-58a80f99840f
# ╠═5a765471-87f2-406b-ace1-ee2b51cf35b6
