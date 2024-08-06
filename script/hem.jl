using Revise
# using Combinatorics
# using Ipopt
# using JuMP
using YAML

using HolisticElectricityModel
import HolisticElectricityModelData

const HEMDataRepo = HolisticElectricityModelData

# Load config file
if length(ARGS) > 0 # otherwise, define this variable in the REPL
    config_fp = ARGS[0]
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

balancing_areas, base_year, num_future_years, num_ipps, folder_suffix, skip_parse =
    parse(config, "data_selection", validators)

# Define scenario
future_years = [base_year + i for i in range(1, num_future_years)]
scenario = HEMDataRepo.DataSelection(balancing_areas, base_year, future_years, num_ipps)

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
    HEMDataRepo.parse_inputs(input_path, input_dir, scenario, false)
end

error("Continue implementation here")

# Acceptable options config_fields; maps config_fields to lists of acceptable options
option_dict = Dict{String,Vector{Any}}()
option_dict["balancing_areas"] = ["p129", "p130", "p131", "p132", "p133", "p134"]
option_dict["base_year"] = collect(2018:2024)  # unsure if these base years are reasonable
option_dict["num_future_years"] = collect(1:30)  # unsure if these horizons are reasonable
option_dict["market_structure"] = ["wholesale_market", "vertically_integrated_utility"]
option_dict["der_use_case"] = ["der_use_case", "null_use_case"]
option_dict["supply_choice_use_case"] = ["supply_choice_use_case", "null_use_case"]
option_dict["num_ipps"] = [1]  # unsure if more IPPs are desired
option_dict["ipp_algorithm"] =
    ["lagrange_decomposition", "mppdcmer_transportation_storage", "mpppcmer", "miqp"]
option_dict["rate_design"] = ["flat_rate", "time_of_use"]
option_dict["net_metering_policy"] =
    ["excess_retail_rate", "excess_marginal_cost", "excess_zero"]
option_dict["solver"] = ["Gurobi", "Xpress"]

function unpack_config_struct(config, config_field, option_dict)
    # if the config_field isn't supported
    if !(config_field in keys(config))
        error("Invalid config field: $config_field.")
        # elseif the option specified in the yaml for the config_field isn't an acceptable option
    elseif !(config[config_field] in option_dict[config_field])
        error("Invalid option $config[config_field] for config_field $config_field.")
        # all good -- unpack and return
    else
        config[config_field]
    end
end

bal_areas = option_dict["balancing_areas"]
base_year = unpack_config_struct(config, "base_year", option_dict)
num_future_years = unpack_config_struct(config, "num_future_years", option_dict)
market_structure = unpack_config_struct(config, "market_structure", option_dict)
der_use_case = unpack_config_struct(config, "der_use_case", option_dict)
supply_choice_use_case = unpack_config_struct(config, "supply_choice_use_case", option_dict)
num_ipps = unpack_config_struct(config, "num_ipps", option_dict)
ipp_algorithm = unpack_config_struct(config, "ipp_algorithm", option_dict)
rate_design = unpack_config_struct(config, "rate_design", option_dict)
net_metering_policy = unpack_config_struct(config, "net_metering_policy", option_dict)
solver = unpack_config_struct(config, "solver", option_dict)




# Restructure HEM options from user/config file input
# market structure
if market_structure == "wholesale_market"
    market_structure = WholesaleMarket()
elseif market_structure == "vertically_integrated_utility"
    market_structure = VerticallyIntegratedUtility()
end
# DER use case
if der_use_case == "der_use_case"
    der_use_case = DERUseCase()
elseif der_use_case == "null_use_case"
    der_use_case = NullUseCase()
end
# Supply Choice use case
if supply_choice_use_case == "supply_choice_use_case"
    supply_choice_use_case = SupplyChoiceUseCase()
elseif supply_choice_use_case == "null_use_case"
    supply_choice_use_case = NullUseCase()
end
# IPP algorithm
if ipp_algorithm == "lagrange_decomposition"
    ipp_algorithm = LagrangeDecomposition()
elseif ipp_algorithm == "mppdcmer_transportation_storage"
    ipp_algorithm = MPPDCMERTransStorage()
elseif ipp_algorithm == "mppdcmer"
    ipp_algorithm = MPPDCMER()
elseif ipp_algorithm == "miqp"
    ipp_algorithm = MIQP()
end
# Rate Design
if rate_design == "flat_rate"
    rate_design = FlatRate()
elseif rate_design == "time_of_use"
    rate_design = TOU()
end
# Net Metering Policy
if net_metering_policy == "excess_retail_rate"
    net_metering_policy = ExcessRetailRate()
elseif net_metering_policy == "excess_marginal_cost"
    net_metering_policy = ExcessMarginalCost()
elseif net_metering_policy == "excess_zero"
    net_metering_policy = ExcessZero()
end


# Define HEM run options
hem_opts = HEMOptions(market_structure, der_use_case, supply_choice_use_case)

regulator_opts = RegulatorOptions(rate_design, net_metering_policy)

ipp_opts = IPPOptions(
    ipp_algorithm,
    Dict(
        "Lagrange_Sub_Investment_Retirement_Cap" => JuMP.optimizer_with_attributes(
            Ipopt.Optimizer,
            "print_level" => 0,
            # "tol" => 1e-6,
            # "max_iter" => 500,
        ),
        "Lagrange_Sub_Dispatch_Cap" => JuMP.optimizer_with_attributes(
            () -> Gurobi.Optimizer(GUROBI_ENV),
            # "OUTPUTLOG" => 0,
        ),
        "Lagrange_Feasible_Cap" => JuMP.optimizer_with_attributes(
            () -> Gurobi.Optimizer(GUROBI_ENV),
            "Presolve" => 0,
            # "OUTPUTLOG" => 0,
        ),
        "solve_agent_problem_ipp_cap" => JuMP.optimizer_with_attributes(
            () -> Gurobi.Optimizer(GUROBI_ENV),
            "Presolve" => 1,
            # "OUTPUTLOG" => 0,
        ),
        "solve_agent_problem_ipp_mppdc" => JuMP.optimizer_with_attributes(
            () -> Gurobi.Optimizer(GUROBI_ENV),
            "Presolve" => 1,
            "BarHomogeneous" => 1,
            # "NumericFocus" => 3,
            # "OUTPUTLOG" => 0,
        ),
        "solve_agent_problem_ipp_mppdc_mccormic_lower" =>
            JuMP.optimizer_with_attributes(
                () -> Gurobi.Optimizer(GUROBI_ENV),
                "Presolve" => 1,
                "BarHomogeneous" => 1,
                # "OUTPUTLOG" => 0,
            ),
    ),
)

utility_opts = UtilityOptions(JuMP.optimizer_with_attributes(
    () -> Gurobi.Optimizer(GUROBI_ENV),
    # "OUTPUTLOG" => 0,
))

green_developer_opts = GreenDeveloperOptions(
    JuMP.optimizer_with_attributes(
        () -> Gurobi.Optimizer(GUROBI_ENV),
        # "OUTPUTLOG" => 0,
    ),
)

customers_opts = CustomersOptions(JuMP.optimizer_with_attributes(
    () -> Gurobi.Optimizer(GUROBI_ENV),
    # "OUTPUTLOG" => 0,
))

jump_model = []

# Run HEM
output_dir = run_hem(
    input_dir,
    hem_opts,
    regulator_options = regulator_opts,
    ipp_options = ipp_opts,
    utility_options = utility_opts,
    green_developer_options = green_developer_opts,
    customers_options = customers_opts,
    force = true,
    jump_model = jump_model,
)
# ------------------------------------------------------------------------------
