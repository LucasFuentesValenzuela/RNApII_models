"""Functions for use with the ChIP data"""

PATH = "/Users/lucasfuentes/Library/CloudStorage/GoogleDrive-lucasfv@stanford.edu/My Drive/Work/Skotheim/Data/RNApIIModels_data"

"""
Load the ChIP dataset.
"""
function load_ChIP_data()
	return DataFrame(
		CSV.File(
			joinpath(PATH, "data_out", "data_exp_nuclear_fraction.csv")
			)
	);
end

"""
Load the gene bins dataset.
"""
function load_gene_bins()
	return DataFrame(
		XLSX.readtable(
		joinpath(PATH, "data_in", "global_norm_filtered-noCC-noESR.xlsx"), "Sheet1"
		)
	);
end

# Series of functions used to compare data with simulations

RNA_free_frac(df) = (1 .-df[!, :Rpb1_bound_fraction_haploid_prediction]) ./ (df[!, :nuclear_volume] ./ df[!, :cell_volume_fL])

RNA_free_frac() = RNA_free_frac(load_ChIP_data())

Rpb1_occupancy_haploid_interp(df) = linear_interpolation(
	df[!, :cell_volume_fL], df[!, :Rpb1_occupancy_haploid_fit]
)
Rpb1_occupancy_haploid_interp() = Rpb1_occupancy_haploid_interp(load_ChIP_data())

CV_to_RNAfree_interp(df) = linear_interpolation(
		df[!, :cell_volume_fL], RNA_free_frac(df)
)
CV_to_RNAfree_interp() = CV_to_RNAfree_interp(load_ChIP_data())

RNAfree_to_CV_interp(df) = linear_interpolation(
		RNA_free_frac(df), df[!, :cell_volume_fL]
)
RNAfree_to_CV_interp() = RNAfree_to_CV_interp(load_ChIP_data())