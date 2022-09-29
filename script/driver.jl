using Logging
using HolisticElectricityModel
using JuMP

# This is the driver script

# Define the solver ------------------------------------------------------------
# using Xpress
using Ipopt
using Gurobi
const GUROBI_ENV = Gurobi.Env()
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

input_dir = joinpath(hem_data_dir, "inputs", "ba_"*"$ba_len"*"_base_"*"$base_year"*"_future_"*"$future_years_len"*"_ipps_"*"$ipp_number")
mkpath(input_dir)

main(input_path, input_dir, scenario)

# export_file_path = joinpath(hem_data_dir, "outputs", "ba_"*"$ba_len"*"_base_"*"$base_year"*"_future_"*"$future_years_len"*"_ipps_"*"$ipp_number")
# mkpath(export_file_path)

hem_opts = HEMOptions(
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
        )
    )
)

utility_opts = UtilityOptions(
    JuMP.optimizer_with_attributes(
        () -> Gurobi.Optimizer(GUROBI_ENV),
        # "OUTPUTLOG" => 0,
    ),
)

# Load sets and parameters, define functions -----------------------------------
run_hem(
    input_dir,
    hem_opts,
    regulator_options=regulator_opts,
    ipp_options=ipp_opts,
    utility_options=utility_opts,
)
