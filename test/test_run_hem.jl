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


    output_file = joinpath(TEST_DIR, tempdir())
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