using HolisticElectricityModel

# This is the driver script

# Define the model run ---------------------------------------------------------

# File locations
hem_data_dir = joinpath(@__DIR__, "..", "..", "HolisticElectricityModel-Data")
input_filename = joinpath(hem_data_dir, "inputs", "HEM_Parameters_ReEDS_17_dGen_julia.xlsx")
exportfilepath = joinpath(hem_data_dir, "outputs")
if !isdir(exportfilepath)
    mkdir(exportfilepath)
end

# User can change retail rate making, DER net metering policy, and market structure
retailrate = TOU()     # FlatRate(), TOU()
dernetmetering = ExcessRetailRate()     # ExcessRetailRate(), ExcessMarginalCost(), ExcessZero()
marketstructure = WholesaleMarket()     # VerticallyIntegratedUtility(), WholesaleMarket()

# Logging options
#set_log_level(Logging.Debug) # if commented, will revert to Info

# Load sets and parameters, define functions -----------------------------------
@info "Loading data"
model_data = HEMData(input_filename)
regulator = Regulator(input_filename, model_data)
utility = Utility(input_filename, model_data)
customers = Customers(input_filename, model_data)
ipp = IPP(input_filename, model_data)

solve_equilibrium_problem(marketstructure, retailrate, dernetmetering,
    model_data, regulator, utility, customers, ipp, exportfilepath)
