using DataFrames
using CSV
using Interpolations

df = DataFrame(CSV.File("Python/data_exp_nuclear_fraction.csv"));
# fraction * cell volume = amount
# amount/nuc_vol = concentration
RNA_free_frac = (1 .-df[!, :Rpb1_bound_fraction_haploid_prediction]) ./ (df[!, :nuclear_volume] ./ df[!, :cell_volume_fL])
Rpb1_occupancy_haploid_interp = linear_interpolation(
	df[!, :cell_volume_fL], df[!, :Rpb1_occupancy_haploid_fit]
)
CV_to_RNAfree_interp = linear_interpolation(
		df[!, :cell_volume_fL], RNA_free_frac
)