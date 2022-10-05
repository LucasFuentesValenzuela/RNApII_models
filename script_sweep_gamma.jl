using JLD2

include("model.jl")

n_points_ = 5
α_vec = 10. .^(LinRange(-3, 0, n_points_));
param_name = "γ"
p_vec = DEFAULT_PARAMS.γ * 10. .^([-1, -.5, 0, .5, 1]);

params_dict, trans_rates, residence_times, densities = sweep_params(
	α_vec, p_vec, DEFAULT_PARAMS, "γ"
);

JLD2.jldsave(
    "results/sweep_gamma_test.jld2";
    α_vec, 
    p_vec,
    params_dict, 
    trans_rates, 
    residence_times, 
    densities, 
    test_var = Dict(1 => [2, 3])
    )