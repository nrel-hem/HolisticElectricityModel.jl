using Logging
using HolisticElectricityModel
using JuMP

# This is the driver script

# Define the solver ------------------------------------------------------------
using Ipopt
using Gurobi
const GUROBI_ENV = Gurobi.Env()
# ------------------------------------------------------------------------------

# Define the model run ---------------------------------------------------------

# File locations
base_dir = abspath(joinpath(dirname(Base.find_package("HolisticElectricityModel")), ".."))
hem_data_dir = joinpath(base_dir, "..", "HolisticElectricityModel-Data")
test_data_dir = joinpath(base_dir, "test", "driver_outputs")

# Create input data
include(joinpath(hem_data_dir, "inputs", "input_data_parsing.jl"))

#                                                # Test inputs
input_path = joinpath(hem_data_dir, "inputs")
ba = ["p13"]                                     # p13
ba_len = length(ba)
base_year = 2018                                 # 2018
future_years = [2019, 2020]                      # [2019, 2020]
future_years_len = length(future_years)
ipp_number = 1                                   # 1
scenario = DataSelection(ba, base_year, future_years, ipp_number)

# need to run in julia: run(#ba, PROFILES_DIRECTORY, "nguo", HOSTNAME, DATABASE, PORT) to get residential and commercial profiles
# also need to run in command prompt: python inputs/write_industrial_profiles.py #ba to get industrial profiles

input_dir_name = "ba_"*"$ba_len"*"_base_"*"$base_year"*"_future_"*"$future_years_len"*"_ipps_"*"$ipp_number"
input_dir = joinpath(hem_data_dir, "runs", input_dir_name)
mkpath(input_dir)

main(input_path, input_dir, scenario)

# Define the scenario and other run options
hem_opts = HEMOptions(
    VerticallyIntegratedUtility(),    # VerticallyIntegratedUtility(), WholesaleMarket()
    DERUseCase(),                     # DERUseCase(), NullUseCase()
    NullUseCase(),                    # SupplyChoiceUseCase(), NullUseCase()
)

regulator_opts = RegulatorOptions(
    FlatRate(),                       # FlatRate(), TOU()
    ExcessRetailRate(),               # ExcessRetailRate(), ExcessMarginalCost(), ExcessZero()
)

ipp_opts = IPPOptions(
    LagrangeDecomposition(),          # LagrangeDecomposition(), MIQP()
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

green_developer_opts = UtilityOptions(
    JuMP.optimizer_with_attributes(
        () -> Gurobi.Optimizer(GUROBI_ENV),
        # "OUTPUTLOG" => 0,
    ),
)
# ------------------------------------------------------------------------------

# Run HEM ----------------------------------------------------------------------
output_dir = run_hem(
    input_dir,
    hem_opts,
    regulator_options=regulator_opts,
    ipp_options=ipp_opts,
    utility_options=utility_opts,
    green_developer_options=green_developer_opts
)
# ------------------------------------------------------------------------------
