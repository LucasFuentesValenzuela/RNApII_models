"""Functions for use with the ChIP data"""

"""
Load the ChIP dataset.
"""
function load_ChIP_data(PATH::String)
	return DataFrame(
		CSV.File(
			joinpath(PATH, "data_out", "data_exp_nuclear_fraction.csv")
			)
	);
end

"""
Load the gene bins dataset.
"""
function load_gene_bins(PATH::String)
	return DataFrame(
		XLSX.readtable(
		joinpath(PATH, "data_in", "global_norm_filtered-noCC-noESR.xlsx"), "Sheet1"
		)
	);
end

# Series of functions used to compare data with simulations

RNA_free_frac(df) = (1 .-df[!, :Rpb1_bound_fraction_haploid_prediction]) ./ (df[!, :nuclear_volume] ./ df[!, :cell_volume_fL])

RNA_free_frac() = RNA_free_frac(load_ChIP_data(PATH))

Rpb1_occupancy_haploid_interp(df::DataFrames.DataFrame) = linear_interpolation(
	df[!, :cell_volume_fL], df[!, :Rpb1_occupancy_haploid_fit]
)
Rpb1_occupancy_haploid_interp(PATH::AbstractString) = Rpb1_occupancy_haploid_interp(load_ChIP_data(PATH))

CV_to_RNAfree_interp(df::DataFrames.DataFrame) = linear_interpolation(
		df[!, :cell_volume_fL], RNA_free_frac(df)
)
CV_to_RNAfree_interp(PATH::AbstractString) = CV_to_RNAfree_interp(load_ChIP_data(PATH))

RNAfree_to_CV_interp(df::DataFrames.DataFrame) = linear_interpolation(
		RNA_free_frac(df), df[!, :cell_volume_fL]
)
RNAfree_to_CV_interp(PATH::AbstractString) = RNAfree_to_CV_interp(load_ChIP_data(PATH))