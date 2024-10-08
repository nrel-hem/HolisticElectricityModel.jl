using HolisticElectricityModel
import HolisticElectricityModelData
using JuMP

# This is the driver script

# Define the solver ------------------------------------------------------------
using Xpress
using Ipopt
# ------------------------------------------------------------------------------

const HEMDataRepo = HolisticElectricityModelData

# Define the model run ---------------------------------------------------------

# File locations
base_dir = abspath(joinpath(dirname(Base.find_package("HolisticElectricityModel")), ".."))
hem_data_dir = dirname(dirname(Base.find_package("HolisticElectricityModelData")))
test_data_dir = joinpath(base_dir, "test", "driver_outputs")

# Create input data
input_path = joinpath(hem_data_dir, "inputs")
ba = ["p13"]
ba_len = length(ba)
base_year = 2018
future_years = [2019, 2020]
future_years_len = length(future_years)
ipp_number = 1
scenario = HEMDataRepo.DataSelection(ba, base_year, future_years, ipp_number)

input_dir_name = "ba_"*"$ba_len"*"_base_"*"$base_year"*"_future_"*"$future_years_len"*"_ipps_"*"$ipp_number"
input_dir = joinpath(hem_data_dir, "runs", input_dir_name)
mkpath(input_dir)

HEMDataRepo.parse_inputs(input_path, input_dir, scenario)

# Define the scenario and other run options
hem_opts = HEMOptions(
    WholesaleMarket(),    # VerticallyIntegratedUtility(), WholesaleMarket()
    DERAdoption(),                    # DERAdoption(), NullUseCase()
    NullUseCase(),                    # SupplyChoice(), NullUseCase()
    DERAggregation(),                 # DERAggregation(), NullUseCase()
)

regulator_opts = RegulatorOptions(
    TOU(),                            # FlatRate(), TOU()
    ExcessZero(),                     # ExcessRetailRate(), ExcessMarginalCost(), ExcessZero()
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
            () -> Xpress.Optimizer(),
            # "OUTPUTLOG" => 0,
        ),
        "Lagrange_Feasible_Cap" => JuMP.optimizer_with_attributes(
            () -> Xpress.Optimizer()
        ),
        "solve_agent_problem_ipp_cap" => JuMP.optimizer_with_attributes(
            () -> Gurobi.Optimizer(GUROBI_ENV),
            "Presolve" => 0,
            # "OUTPUTLOG" => 0,
        ),
        "solve_agent_problem_ipp_mppdc" => JuMP.optimizer_with_attributes(
            () -> Gurobi.Optimizer(GUROBI_ENV),
            "Presolve" => 1,
            # "OUTPUTLOG" => 0,
        ),
        "solve_agent_problem_ipp_mppdc_mccormic_lower" => JuMP.optimizer_with_attributes(
            () -> Gurobi.Optimizer(GUROBI_ENV),
            "Presolve" => 1,
            # "OUTPUTLOG" => 0,
        )
    )
)

utility_opts = UtilityOptions(
    JuMP.optimizer_with_attributes(
        () -> Xpress.Optimizer(),
        # "OUTPUTLOG" => 0,
    ),
)

green_developer_opts = GreenDeveloperOptions(
    JuMP.optimizer_with_attributes(
        () -> Xpress.Optimizer(),
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
    green_developer_options=green_developer_opts,
    force=true,
)
# ------------------------------------------------------------------------------
