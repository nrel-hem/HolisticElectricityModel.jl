using Logging
using HolisticElectricityModel

# This is the driver script

# Define the solver ------------------------------------------------------------
using Xpress
MIP_solver = XpressSolver(Xpress)

using Ipopt
NLP_solver = Ipopt_Solver(Ipopt)

#using Gurobi
#const GRB_ENV = Gurobi.Env()
#MIP_solver = Gurobi_Solver(Gurobi, GRB_ENV)
# ------------------------------------------------------------------------------

# Define the model run ---------------------------------------------------------

# File locations
base_dir = abspath(joinpath(dirname(Base.find_package("HolisticElectricityModel")), ".."))
hem_data_dir = joinpath(base_dir, "..", "HolisticElectricityModel-Data")
# 20220328@ETH - For now input_filename can be a folder of csvs or an xlsx file
input_filename =
    joinpath(hem_data_dir, "inputs", "HEM_Parameters_ipp1_single_year_final.xlsx")     # HEM_Parameters_ipp1_single_year_final, HEM_Parameters_ipp1_two_year_test
export_file_path = joinpath(hem_data_dir, "outputs")
mkpath(export_file_path)

logger = configure_logging(
    console_level = Logging.Info,
    file_level = Logging.Info,
    filename = "driver.log",
)

hem_opts = HEMOptions(
    MIP_solver,                       # HEMSolver
    NLP_solver,
    VerticallyIntegratedUtility(),    # MarketStructure    # VerticallyIntegratedUtility(), WholesaleMarket()
    DERUseCase(),                     # UseCase            # NullUseCase, DERUseCase
    SupplyChoiceUseCase(),            # UseCase            # NullUseCase, SupplyChoiceUseCase
)

regulator_opts = RegulatorOptions(
    FlatRate(),                       # RateDesign         # FlatRate(), TOU()
    ExcessRetailRate(),               # NetMeteringPolicy  # ExcessRetailRate(), ExcessZero()
)

ipp_opts = IPPOptions(
    LagrangeDecomposition(),          # LagrangeDecomposition, MIQP
)

# Load sets and parameters, define functions -----------------------------------
@info "Loading data"
model_data = HEMData(input_filename)
regulator = Regulator(input_filename, model_data)
utility = Utility(input_filename, model_data, regulator)
customers = CustomerGroup(input_filename, model_data)
ipp = IPPGroup(input_filename, model_data)
green_developer = GreenDeveloper(input_filename, model_data)

agents_and_opts = [
    AgentAndOptions(utility, NullAgentOptions()),
    AgentAndOptions(ipp, ipp_opts),
    AgentAndOptions(regulator, regulator_opts),
    AgentAndOptions(customers, NullAgentOptions()),
    AgentAndOptions(green_developer, NullAgentOptions()),
]

file_prefix = get_file_prefix(hem_opts, agents_and_opts)
@info "file_prefix: $(file_prefix)"

solve_equilibrium_problem!(
    hem_opts,
    model_data,
    agents_and_opts,
    export_file_path,
    file_prefix,
)

close(logger)
