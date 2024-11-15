using Revise
# using Combinatorics
# using Ipopt
using JuMP
using YAML

using HolisticElectricityModel
import HolisticElectricityModelData

const HEMDataRepo = HolisticElectricityModelData

# Load config file
if length(ARGS) > 0 # otherwise, define this variable in the REPL
    config_fp = ARGS[0]
else
    config_fp = joinpath(@__DIR__, "configs", "hem_config.yaml")
end
config = YAML.load_file(config_fp)

# File locations
base_dir = abspath(joinpath(dirname(Base.find_package("HolisticElectricityModel")), ".."))
hem_data_dir = dirname(dirname(Base.find_package("HolisticElectricityModelData")))
input_path = joinpath(hem_data_dir, "inputs")

include(joinpath(base_dir, "script", "config.jl"))
include(joinpath(base_dir, "script", "config_data.jl"))

# ------------------------------------------------------------------------------
# Input Folder and Data
# ------------------------------------------------------------------------------

balancing_areas, base_year, num_future_years, num_ipps, folder_suffix, der_aggregator, set_nuclear_varcost_negative, sectors_with_county_level_load, skip_parse =
    parse(config, "data_selection", validators)

# Define scenario
future_years = [base_year + i for i in range(1, num_future_years)]
scenario = HEMDataRepo.DataSelection(balancing_areas, base_year, future_years, num_ipps)

# Define parse options
options = HEMDataRepo.ParseOptions(
    der_aggregator,
    set_nuclear_varcost_negative,
    sectors_with_county_level_load;
    load_profiles_subdir="NewEngland"
)

bal_areas_len = length(balancing_areas)
input_dir_name =
    "ba_" *
    "$bal_areas_len" *
    "_base_" *
    "$base_year" *
    "_future_" *
    "$num_future_years" *
    "_ipps_" *
    "$num_ipps" *
    "_" *
    "$folder_suffix"
input_dir = joinpath(hem_data_dir, "runs", input_dir_name)

if !skip_parse
    mkpath(input_dir)
    # For the following to work, someone needs to have done the following from HolisticElectricityModelData:
    #     - In julia: run(output_dir = PROFILES_DIRECTORY, user = "nguo", hostname = HOSTNAME, dbname = DATABASE, port = PORT, pca_ids = nothing) to get residential and commercial profiles
    #     - In command prompt: python inputs/write_industrial_profiles.py #ba to get industrial profiles
    HEMDataRepo.parse_inputs(input_path, input_dir, scenario, options)
end

# ------------------------------------------------------------------------------
# Model configuration
# ------------------------------------------------------------------------------

market_structure, der_use_case, supply_choice_use_case, der_aggregation_use_case = parse(config, "hem_options", validators)
ipp_algorithm, = parse(config, "ipp_options", validators)
rate_design, net_metering_policy = parse(config, "regulator_options", validators)
pv_adoption_type, = parse(config, "customer_options", validators)

# Define HEM run options
hem_opts = HEMOptions(market_structure, der_use_case, supply_choice_use_case, der_aggregation_use_case)

regulator_opts = RegulatorOptions(rate_design, net_metering_policy)

# Get the optimizer depending on the solver defined the config
mip_solver = get_optimizer_for_solver("Ipopt")
lp_solver = get_optimizer_for_solver(solver)

ipp_opts = IPPOptions(
    ipp_algorithm,
    Dict(
        "Lagrange_Sub_Investment_Retirement_Cap" => JuMP.optimizer_with_attributes(
            mip_solver,
            "print_level" => 0,
            # "tol" => 1e-6,
            # "max_iter" => 500,
        ),
        "Lagrange_Sub_Dispatch_Cap" => JuMP.optimizer_with_attributes(
            lp_solver,
            # "OUTPUTLOG" => 0,
        ),
        "Lagrange_Feasible_Cap" => JuMP.optimizer_with_attributes(
            lp_solver,
            "Presolve" => 0,
            # "OUTPUTLOG" => 0,
        ),
        "solve_agent_problem_ipp_cap" => JuMP.optimizer_with_attributes(
            lp_solver,
            "Presolve" => 1,
            # "OUTPUTLOG" => 0,
        ),
        "solve_agent_problem_ipp_mppdc" => JuMP.optimizer_with_attributes(
            lp_solver,
            "Presolve" => 1,
            "BarHomogeneous" => 1,
            # "NumericFocus" => 3,
            # "OUTPUTLOG" => 0,
        ),
        "solve_agent_problem_ipp_mppdc_mccormic_lower" =>
            JuMP.optimizer_with_attributes(
                lp_solver,
                "Presolve" => 1,
                "BarHomogeneous" => 1,
                # "OUTPUTLOG" => 0,
            ),
    ),
)

utility_opts = UtilityOptions(JuMP.optimizer_with_attributes(
    lp_solver,
    # "OUTPUTLOG" => 0,
))

green_developer_opts = GreenDeveloperOptions(
    JuMP.optimizer_with_attributes(
        lp_solver,
        # "OUTPUTLOG" => 0,
    ),
)

customer_opts = CustomerOptions(
    pv_adoption_type,
    JuMP.optimizer_with_attributes(
        lp_solver,
        # "OUTPUTLOG" => 0,
    ))

dera_opts = DERAggregatorOptions(
    JuMP.optimizer_with_attributes(
        () -> Xpress.Optimizer(),
        # "OUTPUTLOG" => 0,
    ),
)

jump_model = []

# Run HEM
output_dir = run_hem(
    input_dir,
    hem_opts;
    regulator_options=regulator_opts,
    ipp_options=ipp_opts,
    utility_options=utility_opts,
    green_developer_options=green_developer_opts,
    customer_options=customer_opts,
    dera_options=dera_opts,
    force=true,
    jump_model=jump_model,
)
# ------------------------------------------------------------------------------
