using RNApIIModels
using ArgParse
using JLD2
using ProgressBars
using Statistics

PATH = "/Users/lucasfuentes/RNApII_models"

"""
"""
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--type"
            help = "type of analysis to run"
            arg_type = String
            default = "screen"
        "--fnm" 
            help = "output filename"
            arg_type = String
            default = nothing
        "--ntimes"
            help = "number of times each simulation is run"
            arg_type = Int
            default = 1
        "--Omega"
            help = "residence time of RNAp on the promoter"
            arg_type = Float64
            default = 1.7
    end
    return parse_args(s)
end

"""
"""
function build_iteration_params(type, fnm_screen)

    if type=="screen"

        params_iter = collect(
            Iterators.product(
                RNApIIModels.α_vec_screen, RNApIIModels.β_screen, [RNApIIModels.k_on_vec_screen]
                )
        )

    elseif type=="narrow"

        # unpack the results from the screen
        _, feasible_points = get_feasible_pts(fnm_screen)
        println("There are $(length(feasible_points)) feasible points.")

        params_iter = []
        CV_interps = LinRange(30, 200, RNApIIModels.n_kon_pts_screen)
        RNA_free_interps = CV_to_RNAfree_interp().(CV_interps)
        RNA_free_avgCell = CV_to_RNAfree_interp()(avg_cell_size)

        for (k_crt, α_crt) in feasible_points

            k_on_vec_crt = RNA_free_interps .* k_crt ./ RNA_free_avgCell

            push!(params_iter, (α_crt, RNApIIModels.β_screen, k_on_vec_crt))

        end

    elseif type=="wide"

        # unpack the results from the screen

        _, feasible_points = get_feasible_pts(fnm_screen)
        println("There are $(length(feasible_points)) feasible points.")

        params_iter = []

        for (k_crt, α_crt) in feasible_points

            k_on_vec_crt = k_crt .* (LinRange(10^(-1.5), 10^(1.5), 10))

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
    Ω = parsed_args["Omega"]
    ntimes = parsed_args["ntimes"]

    parsed_args["fnm"] === nothing ? fnm = "feasible_pts_$(type)_Omega$(Ω).jld2" : fnm = parsed_args["fnm"]

    OCCUPANCY_PARAMS_crt = copy(OCCUPANCY_PARAMS)
    OCCUPANCY_PARAMS_crt["Ω"] = Ω

    println("Running simulations with the current parameters:")
    @show OCCUPANCY_PARAMS_crt

    println(
        "Running analyis $(type) with $(ntimes) repetitions, Ω=$(Ω) and saving at $(fnm)"
    )

    fnm_screen = joinpath(PATH, "results", "feasible_pts_screen_Omega$(Ω).jld2")
    params_iter = build_iteration_params(type, fnm_screen)

    occupancy = []
    promoter_occ = []
    params_occ = []

    for _ in ProgressBar(1:ntimes)

        occupancy_crt, promoter_occ_crt, params_occ_crt = run_occupancy_simulation(
            params_iter, OCCUPANCY_PARAMS_crt
        )

        push!(occupancy, occupancy_crt)
        push!(promoter_occ, promoter_occ_crt)
        push!(params_occ, params_occ_crt)
    end
    
    JLD2.jldsave(
        joinpath(PATH, "results", fnm); 
        occupancy, promoter_occ, params_occ, params_iter, OCCUPANCY_PARAMS_crt
    )

end

main()