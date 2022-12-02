using RNApIIModels
PATH_TO_REPO = "/Users/lucasfuentes/RNApII_models"

n_kon_pts = 10
n_α_values = 8

# Parameters for a narrow screen with feasible and infeasible points
k_on_vec_screen = 10 .^(LinRange(
    log10(LITERATURE_PARAMS["min_k_on"]/1.5), log10(LITERATURE_PARAMS["max_k_on"]*1.5), n_kon_pts
))
α_vec_screen = LinRange(LITERATURE_PARAMS["min_α"], LITERATURE_PARAMS["max_α"], n_α_values)
β_screen = LITERATURE_PARAMS["max_β"] 

n_times = Dict(
    "screen" => 1, 
    "wide" => 1,
    "narrow" => 1
)

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

        k_on_vec = k_on_vec_screen
        α_vec = α_vec_screen
        params_iter = collect(
        Iterators.product(α_vec, β_screen, [k_on_vec])
        )

    elseif type=="narrow"

        # unpack the results from the screen
        results_fs = JLD2.load(joinpath(PATH, "results", "feasible_pts_screen.jld2"))
        params_iter = results_fs["params_iter"]
        feasible = results_fs["feasible"]

        params_iter = []
        # kon_to_CV_interps = [] # I think you can go around that
        CV_interps = LinRange(30, 200, n_kon_pts)
        RNA_free_interps = CV_to_RNAfree_interp().(CV_interps)
        RNA_free_avgCell = CV_to_RNAfree_interp()(avg_cell_size)

        for (idx_k, idx_α) in Tuple.(findall(feasible .== 1))

            k_crt = k_on_vec_screen[idx_k]
            α_crt = α_vec_screen[idx_α]

            kon_C_crt = k_crt / RNA_free_avgCell

            k_on_vec_crt = RNA_free_interps .* kon_C_crt

            push!(params_iter, (α_crt, β_screen, k_on_vec_crt))

            # kon_to_CV_interp = linear_interpolation(
            #         RNA_free_frac * kon_C_crt, df[!, :cell_volume_fL]
            # )
            # push!(kon_to_CV_interps, kon_to_CV_interp)
        end

    elseif type=="wide"

        # unpack the results from the screen
        results_fs = JLD2.load(joinpath(PATH, "results", "feasible_pts_screen.jld2"))
        params_iter = results_fs["params_iter"]
        feasible = results_fs["feasible"]

        params_iter = []

        for (idx_k, idx_α) in Tuple.(findall(feasible .== 1))

            k_crt = k_on_vec_screen[idx_k]
            α_crt = α_vec_screen[idx_α]

            k_on_vec_crt = k_crt .* (10 .^(LinRange(-1.5, 1.5, 10)))

            push!(params_iter, (α_crt, β_screen, k_on_vec_crt))
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

        occupancy_crt, promoter_occ_crt, params_occ = run_occupancy_simulation(
            params_iter, OCCUPANCY_PARAMS
        )

        push!(occupancy, occupancy_crt)
        push!(promoter_occ, promoter_occ_crt)
        push!(params_occ, params_occ_crt)
    end

    # 2. Determine feasible points
    occ_mat = reduce(hcat, occupancy)
    occ_mat = (occ_mat .< max_ρ_g) .& (occ_mat .> min_ρ_g)
    prom_occ_mat = reduce(hcat, promoter_occ)
    prom_occ_mat = (prom_occ_mat .< max_ρ_p) .& (prom_occ_mat .> min_ρ_p)
    feasible = (occ_mat .& prom_occ_mat)
    
    # 3. Save the data
    JLD2.jldsave(
        joinpath(PATH, "results", "feasible_pts_$(type).jld2"); 
        occupancy, promoter_occ, params_occ, params_iter, 
        feasible
    )



end

main()