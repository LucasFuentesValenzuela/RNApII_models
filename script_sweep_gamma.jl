using JLD

include("model.jl")

n_points_ = 5
α_vec = 10. .^(LinRange(-3, 0, n_points_));
γ_vec = γ_default * 10 .^([-1, -.5, 0, .5, 1]);

params_dict_γ, trans_rates_γ, residence_times_γ, densities_γ = sweep_params(
	α_vec, γ_vec, DEFAULT_PARAMS, "γ"
);

