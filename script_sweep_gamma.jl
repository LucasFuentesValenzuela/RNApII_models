using JLD2

include("model.jl")

n_points_ = 25
α_vec = 10. .^(LinRange(-3, 1, n_points_));
param_name = "γ"
p_vec = DEFAULT_PARAMS.γ * 10. .^(collect(LinRange(-1.5, 1.5, 15)));

params_dict, trans_rates, residence_times, densities = sweep_params(
	α_vec, p_vec, DEFAULT_PARAMS, "γ"
);

JLD2.jldsave(
    "results/sweep_gamma_L35.jld2";
    α_vec, 
    p_vec,
    params_dict, 
    trans_rates, 
    residence_times, 
    densities, 
    )