using Revise
using Gurobi
using Ipopt
using JuMP
using YAML

using HolisticElectricityModel
import HolisticElectricityModelData

const GUROBI_ENV = Gurobi.Env()
const HEMDataRepo = HolisticElectricityModelData

# File locations
base_dir = abspath(joinpath(dirname(Base.find_package("HolisticElectricityModel")), ".."))
hem_data_dir = dirname(dirname(Base.find_package("HolisticElectricityModelData")))
input_path = joinpath(hem_data_dir, "inputs")
test_data_dir = joinpath(base_dir, "test", "driver_outputs")

# Parse inputs
config = YAML.load_file("hem_config.yaml")
bal_areas = config["balancing_areas"]
base_year = config["base_year"]
num_future_years = config["num_future_years"]
market_structure = config["market_structure"]
der_use_case = config["der_use_case"]
supply_choice_use_case = config["supply_choice_use_case"]
num_ipps = config["num_ipps"]
ipp_algorithm = config["ipp_algorithm"]
rate_design = config["rate_design"]
net_metering_policy = config["net_metering_policy"]
solver = config["solver"]

# Define scenario
bal_areas_len = length(bal_areas)
future_years = [base_year + i for i in range(1, num_future_years)]
scenario = HEMDataRepo.DataSelection(bal_areas, base_year, future_years, num_ipps)

# need to run in julia: run(output_dir = PROFILES_DIRECTORY, user = "nguo", hostname = HOSTNAME, dbname = DATABASE, port = PORT, pca_ids = nothing) to get residential and commercial profiles
# also need to run in command prompt: python inputs/write_industrial_profiles.py #ba to get industrial profiles

input_dir_name =
    "ba_" *
    "$bal_areas_len" *
    "_base_" *
    "$base_year" *
    "_future_" *
    "$num_future_years" *
    "_ipps_" *
    "$num_ipps" *
    "_enhanced_test_full"
input_dir = joinpath(hem_data_dir, "runs", input_dir_name)
# mkpath(input_dir)

# HEMDataRepo.parse_inputs(input_path, input_dir, scenario)

# Restructure HEM options from user/config file input
# market structure
if market_structure == "wholesale_market"
    market_structure = WholesaleMarket()
elseif market_structure == "vertically_integrated_utility"
    market_structure = VerticallyIntegratedUtility()
else
    error("Invalid market structure: $market_structure.")
end
# DER use case
if der_use_case == "der_use_case"
    der_use_case = DERUseCase()
elseif der_use_case == "null_use_case"
    der_use_case = NullUseCase()
else
    error("Invalid DER use case: $der_use_case.")
end
# Supply Choice use case
if supply_choice_use_case == "supply_choice_use_case"
    supply_choice_use_case = SupplyChoiceUseCase()
elseif supply_choice_use_case == "null_use_case"
    supply_choice_use_case = NullUseCase()
else
    error("Invalid Supply Choice use case: $supply_choice_use_case.")
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
else
    error("Invalid IPP algorithm: $ipp_algorithm.")
end
# Rate Design
if rate_design == "flat_rate"
    rate_design = FlatRate()
elseif rate_design == "time_of_use"
    rate_design = TOU()
else
    error("Invalid Rate Design: $rate_design.")
end
# Net Metering Policy
if net_metering_policy == "excess_retail_rate"
    net_metering_policy = ExcessRetailRate()
elseif net_metering_policy == "excess_marginal_cost"
    net_metering_policy = ExcessMarginalCost()
elseif net_metering_policy == "excess_zero"
    net_metering_policy = ExcessZero()
else
    error("Invalid Net Metering Policy: $net_metering_policy.")
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
    regulator_options=regulator_opts,
    ipp_options=ipp_opts,
    utility_options=utility_opts,
    green_developer_options=green_developer_opts,
    customers_options=customers_opts,
    force=true,
    jump_model=jump_model,
)
# ------------------------------------------------------------------------------
