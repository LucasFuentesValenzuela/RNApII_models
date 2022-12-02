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