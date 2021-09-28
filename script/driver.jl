using Logging
using HolisticElectricityModel

# This is the driver script

# Define the solver ------------------------------------------------------------
using Xpress
solver = XpressSolver(Xpress)

# using Gurobi
# const GRB_ENV = Gurobi.Env()
# solver = GurobiSolver(Gurobi, GRB_ENV)
# ------------------------------------------------------------------------------

# Define the model run ---------------------------------------------------------

# File locations
base_dir = abspath(joinpath(dirname(Base.find_package("HolisticElectricityModel")), ".."))
hem_data_dir = joinpath(base_dir, "..", "HolisticElectricityModel-Data")
input_filename = joinpath(hem_data_dir, "inputs", "HEM_Parameters_ReEDS_17_dGen_julia.xlsx")
export_file_path = joinpath(hem_data_dir, "outputs")
mkpath(export_file_path)

logger = configure_logging(
    console_level = Logging.Info,
    file_level = Logging.Info,
    filename = "driver.log",
)

hem_opts = HEMOptions(
    solver,                       # HEMSolver    
    VerticallyIntegratedUtility(), # MarketStructure    
)

regulator_opts = RegulatorOptions(
    TOU(),              # RateDesign
    ExcessRetailRate(),  # NetMeteringPolicy
)

# Load sets and parameters, define functions -----------------------------------
@info "Loading data"
model_data = HEMData(input_filename)
regulator = Regulator(input_filename, model_data)
utility = Utility(input_filename, model_data)
customers = CustomerGroup(input_filename, model_data)
ipps = IPPGroup(input_filename, model_data)

file_prefix = "Results_$(hem_opts.market_structure)_$(regulator_opts.rate_design)_$(regulator_opts.net_metering_policy)"

solve_equilibrium_problem!(
    hem_opts,
    model_data,
    [
        AgentAndOptions(regulator, regulator_opts),
        AgentAndOptions(utility, NullAgentOptions()),
        AgentAndOptions(customers, NullAgentOptions()),
        AgentAndOptions(ipps, NullAgentOptions()),
    ],
    export_file_path,
    file_prefix,
)

close(logger)
