using Test
using YAML
using JuMP
using HolisticElectricityModel


include(joinpath("..", "script", "config.jl"))
include(joinpath("..", "script", "config_data.jl"))
include(joinpath("..", "script", "parse_options.jl"))

TEST_DIR = dirname(@__FILE__)

@testset "test_save_config" begin

    config = YAML.load_file(joinpath(TEST_DIR, "data", "test_input_config.yaml"))

    input_dir, = parse(config, "DataSelection", validators)
    solver, = parse(config, "SimulationParameters", validators)
    market_structure, der_use_case, supply_choice_use_case,
    der_aggregation_use_case = parse(config, "HEMOptions", validators)

    hem_opts = HEMOptions(
    market_structure,
    der_use_case,
    supply_choice_use_case,
    der_aggregation_use_case
    )
    agent_options = get_agent_options(config, hem_opts, solver)

    output_file = joinpath(TEST_DIR, tempdir())
    save_config(
        solver,
        input_dir,
        hem_opts,
        agent_options,
        output_file
    )
    @test isfile(joinpath(output_file, "config.yaml"))
    config = YAML.load_file(joinpath(output_file, "config.yaml"))
    @test Set(keys(config)) == Set(["HEMOptions", "Regulator", "IPPGroup", "CustomerGroup", "DERAggregator", "DataSelection", "RunOptions", "SimulationParameters"])
    @test Set(keys(config["HEMOptions"])) == Set(["market_structure", "der_use_case", "supply_choice_use_case", "der_aggregation_use_case"])
    @test Set(keys(config["Regulator"])) == Set(["rate_design", "net_metering_policy", "tou_suffix", "planning_reserve_margin", "allowed_return_on_investment"])
    @test Set(keys(config["IPPGroup"])) == Set(["ipp_algorithm"])
    @test Set(keys(config["CustomerGroup"])) == Set(["pv_adoption_type"])
    @test Set(keys(config["DERAggregator"])) == Set(["incentive_curve", "frac_viu_cost_savings_as_revenue"])
    @test Set(keys(config["DataSelection"])) == Set(["input_path"])
    @test Set(keys(config["RunOptions"])) == Set(["output_dir"])
    @test Set(keys(config["SimulationParameters"])) == Set(["solver"])
end