using JLD2
using ArgParse

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--param"
            help = "parameter to iterate upon"
            arg_type = String
            default = "gamma"
        "--L"
            help = "value of L"
            arg_type = Int
            default = L_default
        "--gamma"
            help = "value of γ"
            arg_type = Float64
            default = γ_default
        "--beta2"
            help = "value of β2 ratio with β"
            arg_type = Int
            default = 8
    end
    return parse_args(s)
end

function edit_default_params(DEFAULT_PARAMS, parsed_args)

    L = parsed_args["L"]
    γ = parsed_args["gamma"]
    ratio_β2 = parsed_args["beta2"]

    params = Params(
        DEFAULT_PARAMS.α, 
        DEFAULT_PARAMS.β, 
        γ,
        L, 
        DEFAULT_PARAMS.Δt, 
        DEFAULT_PARAMS.n_steps, 
        DEFAULT_PARAMS.n_sites, 
        DEFAULT_PARAMS.n_end_sites,
        DEFAULT_PARAMS.β/ratio_β2
    )

    return params
end

# parameters
function main()

    @show parsed_args = parse_commandline()

    n_points_ = 25
    α_vec = 10. .^(LinRange(-3, 1, n_points_));

    # create default values
    if parsed_args["param"] == "gamma"
        p_vec = DEFAULT_PARAMS.γ * 10. .^(collect(LinRange(-1.5, 1.5, 15)));
        param_name = "γ"
        label = "$(parsed_args["param"])_L$(parsed_args["L"])_rat$(parsed_args["beta2"])"
    elseif parsed_args["param"] == "L"
        p_vec = [1, 2, 10, 20, 35]
        param_name = "L"
        label = "$(parsed_args["param"])_1000gamma$(trunc(Int, 1000*parsed_args["gamma"]))_rat$(parsed_args["beta2"])"
    end

    # create params based on kwargs
    params = edit_default_params(DEFAULT_PARAMS, parsed_args)

    @show params

    # run the function
    params_dict, trans_rates, residence_times, densities = sweep_params(
        α_vec, p_vec, param_name; params=params
    );

    JLD2.jldsave(
        "results/sweep_$(label).jld2";
        α_vec, 
        p_vec,
        params_dict, 
        trans_rates, 
        residence_times, 
        densities, 
        parsed_args
        )

end

main()