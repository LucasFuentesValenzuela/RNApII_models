using JLD2

include("model.jl")

n_points_ = 5
α_vec = 10. .^(LinRange(-3, 0, n_points_));
param_name = "L"
p_vec = [1, 2, 5, 10]

params_dict, trans_rates, residence_times, densities = sweep_params(
	α_vec, p_vec, DEFAULT_PARAMS, param_name
);

JLD2.jldsave(
    "results/sweep_L.jld2";
    α_vec, 
    p_vec,
    params_dict, 
    trans_rates, 
    residence_times, 
    densities
    )
