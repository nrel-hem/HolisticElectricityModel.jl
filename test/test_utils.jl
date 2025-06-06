
test_dir = dirname(@__FILE__)

@testset "test read_param with county level index" begin
    index_z = Dimension("index_z", [:p129, :p130, :p131])
    index_h = Dimension("index_h",
        [:Com_50001, :Com_70001, :Com_90001, :Ind_70001, :Ind_90001, :Res_10001, :Res_50001, :Res_90001])
    index_m = Dimension("index_m", [:BTMPV])
    index_d = Dimension("index_d", [:d004])
    index_t = Dimension("index_t", [:h004, :h008])
    index_z_h_map = DimensionSet{2}("index_z_h_map", "index_z_h_map", "mapping from index_z to index_h",
        (index_z, index_h),
        [
            (:p129, :Res_10001),
            (:p129, :Com_50001),
            (:p129, :Ind_90001),
            (:p130, :Res_50001),
            (:p130, :Com_90001),
            (:p130, :Ind_70001),
            (:p131, :Res_90001),
            (:p131, :Com_70001),
            (:p131, :Ind_90001),
        ]
    )

    valid_param = read_param(
        "test_param",
        joinpath(test_dir, "data"),
        "test_input_param_county",
        index_t,
        [index_h, index_m, index_z, index_d]
    )

    for (z, h) in index_z_h_map, m in index_m, d in index_d, t in index_t
        @test !(ismissing(valid_param(h, m, z, d, t)))
    end

    @test_throws ArgumentError begin
        read_param(
            "test_param",
            joinpath(test_dir, "data"),
            "test_input_param_county_invalid",
            index_t,
            [index_h, index_m, index_z, index_d]
        )
    end

    param_df = CSV.read(
        joinpath(test_dir, "data", "test_input_param_county.csv"),
        DataFrame,
    )

    for (z, h) in index_z_h_map, m in index_m, d in index_d, t in index_t
        @test valid_param(h, m, z, d, t) == param_df[
            (param_df.index_h.==String(h)).&(param_df.index_z.==String(z)).&(param_df.index_m.==String(m)).&(param_df.index_d.==String(d)), String(t)
            ][1]
    end

end

@testset "test_save_config" begin

    hem_options = HEMOptions(
        WM(),
        DERAdoption(),
        NullUseCase(),
        DERAggregation()
        )

    basic_config = Dict(
        "regulator_options" => Dict(
            "rate_design" => "time_of_use",
            "net_metering_policy" => "excess_marginal_cost",
            "tou_suffix" => "NE2025",
        ),
        "ipp_options" => Dict(
            "ipp_algorithm" => "mppdcmer"
        ),
        "customer_options" => Dict(
            "pv_adoption_type" => "solar_plus_storage"
        ),
        "der_aggregator_options" => Dict(
            "incentive_curve" => 1,
            "frac_viu_cost_savings_as_revenue" => 0.5,
        ),
    )

    ipp_solvers = Dict()
    addsolvers_ipp!(ipp_solvers, :Ipopt)

    agent_options = AgentOptionsStore(
        Dict(
            Regulator => RegulatorOptions(
                TOU(),
                ExcessMarginalCost(),
            ),
            IPPGroup => IPPOptions(MPPDCMER(), ipp_solvers),
            CustomerGroup => CustomerOptions(
                SolarPlusStorageOnly(),
                JuMP.optimizer_with_attributes(() -> get_optimizer_for_solver(:Ipopt)),
            ),
            GreenDeveloper => GreenDeveloperOptions(
                JuMP.optimizer_with_attributes(() -> get_optimizer_for_solver(:Ipopt)),
            ),
            DERAggregator => DERAggregatorOptions(
                JuMP.optimizer_with_attributes(() -> get_optimizer_for_solver(:Ipopt));
            ),
        )
    )


    output_file = joinpath(test_dir, tempdir())
    save_config(
        hem_options,
        agent_options,
        output_file,
    )
    @test isfile(joinpath(output_file, "config.yaml"))
    config = YAML.load_file(joinpath(output_file, "config.yaml"))
    @test Set(keys(config)) == Set(["HEMOptions", "Regulator", "IPPGroup", "CustomerGroup", "GreenDeveloper", "DERAggregator"])
    @test Set(keys(config["HEMOptions"])) == Set(["market_structure", "der_use_case", "supply_choice_use_case", "der_aggregation_use_case"])
    @test Set(keys(config["Regulator"])) == Set(["rate_design", "net_metering_policy", "tou_suffix", "planning_reserve_margin", "allowed_return_on_investment"])
    @test Set(keys(config["IPPGroup"])) == Set(["ipp_algorithm"])
    @test Set(keys(config["CustomerGroup"])) == Set(["pv_adoption_type"])
end

    