using RNApIIModels
using ArgParse
using JLD2
using ProgressBars
using Statistics

PATH = "/Users/lucasfuentes/RNApII_models"
fnm_screen = joinpath(PATH, "results", "feasible_pts_screen.jld2")


"""
"""
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--type"
            help = "type of analysis to run"
            arg_type = String
            default = "screen"
    end
    return parse_args(s)
end

"""
"""
function build_iteration_params(type)

    if type=="screen"

        # k_on_vec = k_on_vec_screen
        # α_vec = α_vec_screen
        params_iter = collect(
            Iterators.product(
                RNApIIModels.α_vec_screen, RNApIIModels.β_screen, [RNApIIModels.k_on_vec_screen]
                )
        )

    elseif type=="narrow"

        # unpack the results from the screen
        _, feasible_points = get_feasible_pts(fnm_screen)

        params_iter = []
        CV_interps = LinRange(30, 200, RNApIIModels.n_kon_pts_screen)
        RNA_free_interps = CV_to_RNAfree_interp().(CV_interps)
        RNA_free_avgCell = CV_to_RNAfree_interp()(avg_cell_size)

        for (k_crt, α_crt) in feasible_points


            kon_C_crt = k_crt / RNA_free_avgCell

            k_on_vec_crt = RNA_free_interps .* kon_C_crt

            push!(params_iter, (α_crt, RNApIIModels.β_screen, k_on_vec_crt))

        end

    elseif type=="wide"

        # unpack the results from the screen

        _, feasible_points = get_feasible_pts(fnm_screen)

        params_iter = []

        for (k_crt, α_crt) in feasible_points

            k_on_vec_crt = k_crt .* (10 .^(LinRange(-1.5, 1.5, 10)))

            push!(params_iter, (α_crt, RNApIIModels.β_screen, k_on_vec_crt))
        end

    end

    return params_iter

end



"""
"""
function main()

    @show parsed_args = parse_commandline()

    type = parsed_args["type"]

    params_iter = build_iteration_params(type)

    occupancy = []
    promoter_occ = []
    params_occ = []

    for _ in ProgressBar(1:n_times[type])

        occupancy_crt, promoter_occ_crt, params_occ_crt = run_occupancy_simulation(
            params_iter, OCCUPANCY_PARAMS
        )

        push!(occupancy, occupancy_crt)
        push!(promoter_occ, promoter_occ_crt)
        push!(params_occ, params_occ_crt)
    end
    
    # 2. Save the data
    JLD2.jldsave(
        joinpath(PATH, "results", "feasible_pts_$(type).jld2"); 
        occupancy, promoter_occ, params_occ, params_iter
    )

end

main()