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
hem_data_dir = joinpath(@__DIR__, "..", "..", "HolisticElectricityModel-Data")
input_filename = joinpath(hem_data_dir, "inputs", "HEM_Parameters_ReEDS_17_dGen_julia.xlsx")
exportfilepath = joinpath(hem_data_dir, "outputs")
if !isdir(exportfilepath)
    mkdir(exportfilepath)
end

hem_opts = HEMOptions(
    solver,                       # HEMSolver    
    VerticallyIntegratedUtility() # MarketStructure    
) 

regulator_opts = RegulatorOptions(
    TOU(),              # RateDesign
    ExcessRetailRate()  # NetMeteringPolicy
)

configure_logging(console_level = Logging.Info, file_level = Logging.Info, filename = "driver.log")

# Load sets and parameters, define functions -----------------------------------
@info "Loading data"
model_data = HEMData(input_filename)
regulator = Regulator(input_filename, model_data)
utility = Utility(input_filename, model_data)
customers = Customers(input_filename, model_data)
ipp = IPP(input_filename, model_data)

fileprefix = "Results_$(hem_opts.market_structure)_$(regulator_opts.rate_design)_$(regulator_opts.net_metering_policy)"

solve_equilibrium_problem(hem_opts, model_data, [
    AgentAndOptions(regulator, regulator_opts),
    AgentAndOptions(utility, NullAgentOptions()),
    AgentAndOptions(customers, NullAgentOptions()),
    AgentAndOptions(ipp, NullAgentOptions())], 
    exportfilepath, fileprefix)
