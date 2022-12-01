"""Functions for use with the chip data"""

"""
"""
function load_ChIP_data()
	return DataFrame(
		CSV.File(
			"/Users/lucasfuentes/RNApII_models/data/data_exp_nuclear_fraction.csv")
	);
end

# fraction * cell volume = amount
# amount/nuc_vol = concentration
RNA_free_frac = (1 .-df[!, :Rpb1_bound_fraction_haploid_prediction]) ./ (df[!, :nuclear_volume] ./ df[!, :cell_volume_fL])

Rpb1_occupancy_haploid_interp = linear_interpolation(
	df[!, :cell_volume_fL], df[!, :Rpb1_occupancy_haploid_fit]
)

CV_to_RNAfree_interp = linear_interpolation(
		df[!, :cell_volume_fL], RNA_free_frac
)

RNAfree_to_CV_interp = linear_interpolation(
		RNA_free_frac, df[!, :cell_volume_fL]
)