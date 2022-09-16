using Logging
using HolisticElectricityModel

# This is the driver script

# Define the solver ------------------------------------------------------------
# using Xpress
# MIP_solver = XpressSolver(Xpress)
using Ipopt
NLP_solver = Ipopt_Solver(Ipopt)

using Gurobi
const GRB_ENV = Gurobi.Env()
MIP_solver = Gurobi_Solver(Gurobi, GRB_ENV)
# ------------------------------------------------------------------------------

# Define the model run ---------------------------------------------------------

# File locations
base_dir = abspath(joinpath(dirname(Base.find_package("HolisticElectricityModel")), ".."))
hem_data_dir = joinpath(base_dir, "..", "HolisticElectricityModel-Data")

# input_filename =
#     joinpath(hem_data_dir, "inputs", "HEM_Parameters_ipp1_single_year_final.xlsx")     # HEM_Parameters_ipp1_single_year_final, HEM_Parameters_ipp1_two_year_test
include(joinpath(hem_data_dir, "inputs", "input_data_parsing.jl"))

input_path = joinpath(hem_data_dir, "inputs")
ba = ["p13"]
ba_len = length(ba)
base_year = 2018
future_years = [2019, 2020]
future_years_len = length(future_years)
ipp_number = 1
scenario = DataSelection(ba, base_year, future_years, ipp_number)

# need to run in julia: run(#ba, PROFILES_DIRECTORY, "nguo", HOSTNAME, DATABASE, PORT) to get residential and commercial profiles
# also need to run in command prompt: python inputs/write_industrial_profiles.py #ba to get industrial profiles

input_filename = joinpath(hem_data_dir, "inputs", "ba_"*"$ba_len"*"_base_"*"$base_year"*"_future_"*"$future_years_len"*"_ipps_"*"$ipp_number")
mkpath(input_filename)

main(input_path, input_filename, scenario)

export_file_path = joinpath(hem_data_dir, "outputs", "ba_"*"$ba_len"*"_base_"*"$base_year"*"_future_"*"$future_years_len"*"_ipps_"*"$ipp_number")
mkpath(export_file_path)

logger = configure_logging(
    console_level = Logging.Info,
    file_level = Logging.Info,
    filename = "driver.log",
)

hem_opts = HEMOptions(
    MIP_solver,                       # HEMSolver
    NLP_solver,
    WholesaleMarket(),    # MarketStructure    # VerticallyIntegratedUtility(), WholesaleMarket()
    DERUseCase(),                     # DERUseCase          
    NullUseCase(),                    # SupplyChoiceUseCase
)

regulator_opts = RegulatorOptions(
    TOU(),               # RateDesign       # FlatRate, #TOU
    ExcessZero(),  # NetMeteringPolicy    # ExcessRetailRate, ExcessMarginalCost, ExcessZero
)

ipp_opts = IPPOptions(
    LagrangeDecomposition(),              # LagrangeDecomposition, MIQP
)

# Load sets and parameters, define functions -----------------------------------
@info "Loading data"
model_data = HEMData(input_filename)
regulator = Regulator(input_filename, model_data)
utility = Utility(input_filename, model_data, regulator)
customers = CustomerGroup(input_filename, model_data)
ipp = IPPGroup(input_filename, model_data)
green_developer = GreenDeveloper(input_filename, model_data)
distribution_utility = DistributionUtility(input_filename, model_data)

max_iter = 100
window_length = 1

agents_and_opts = [
    AgentAndOptions(utility, NullAgentOptions()),
    AgentAndOptions(ipp, ipp_opts),
    AgentAndOptions(regulator, regulator_opts),
    AgentAndOptions(customers, NullAgentOptions()),
    AgentAndOptions(green_developer, NullAgentOptions()),
    AgentAndOptions(distribution_utility, NullAgentOptions()),
]

file_prefix = get_file_prefix(hem_opts, agents_and_opts)
@info "file_prefix: $(file_prefix)"

solve_equilibrium_problem!(
    hem_opts,
    model_data,
    agents_and_opts,
    export_file_path,
    file_prefix,
    max_iter,
    window_length
)

close(logger)
