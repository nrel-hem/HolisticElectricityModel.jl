using Revise
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
ba = ["p129", "p130", "p131", "p132", "p133", "p134"]
ba_len = length(ba)
base_year = 2020
future_years = [2021, 2022, 2023]
future_years_len = length(future_years)
ipp_number = 1
scenario = HEMDataRepo.DataSelection(ba, base_year, future_years, ipp_number)

inputs_date = "20241119-ba"
input_dir_name = "$inputs_date"*"_ba_"*"$ba_len"*"_base_"*"$base_year"*"_future_"*"$future_years_len"*"_ipps_"*"$ipp_number"
input_dir = joinpath(hem_data_dir, "runs", input_dir_name)
# input_dir = joinpath(test_data_dir, input_dir_name)
# mkpath(input_dir)

# HEMDataRepo.parse_inputs(input_path, input_dir, scenario)

# Define the scenario and other run options
hem_opts = HEMOptions(
    WM(),                             # VIU(), WM()
    DERAdoption(),                    # DERAdoption(), NullUseCase()
    NullUseCase(),                    # SupplyChoice(), NullUseCase()
    DERAggregation(),                 # DERAggregation(), NullUseCase()
)

regulator_opts = RegulatorOptions(
    TOU(),                            # FlatRate(), TOU()
    ExcessRetailRate();               # ExcessRetailRate(), ExcessMarginalCost(), ExcessZero()
    tou_suffix="NE2025",
    planning_reserve_margin=0.129     # Value for New England from ReEDS-2.0/inputs/reserves/prm_annual.csv
)

ipp_opts = IPPOptions(
    MPPDCMERTransStorage(),          # LagrangeDecomposition(), MIQP()
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
            () -> Xpress.Optimizer(),
            # "Presolve" => 1,
            # "OUTPUTLOG" => 0,
        ),
        "solve_agent_problem_ipp_mppdc" => JuMP.optimizer_with_attributes(
            () -> Xpress.Optimizer(),
            # "Aggregate" => 0,
            # "Presolve" => 0,
            # "BarHomogeneous" => 1,
            # "FeasibilityTol" => 1e-3,
            # "Method" => 1
            # "NumericFocus" => 3,
            # "ScaleFlag" => 2,
            # "OUTPUTLOG" => 0,
        ),
        "solve_agent_problem_ipp_mppdc_mccormic_lower" => JuMP.optimizer_with_attributes(
            () -> Xpress.Optimizer(),
            # "Presolve" => 1,
            # "BarHomogeneous" => 1,
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

customer_opts = CustomerOptions(
    CompeteDERConfigs(),
    JuMP.optimizer_with_attributes(
        () -> Xpress.Optimizer(),
        # "OUTPUTLOG" => 0,
    ),
)

dera_opts = DERAggregatorOptions(
    JuMP.optimizer_with_attributes(
        () -> Xpress.Optimizer(),
        # "OUTPUTLOG" => 0,
    );
    incentive_curve=1,
    frac_viu_cost_savings_as_revenue=0.1
)

# ------------------------------------------------------------------------------

jump_model = []

# Run HEM ----------------------------------------------------------------------
output_dir = run_hem(
    input_dir,
    hem_opts,
    regulator_options=regulator_opts,
    ipp_options=ipp_opts,
    utility_options=utility_opts,
    green_developer_options=green_developer_opts,
    customer_options=customer_opts,
    dera_options=dera_opts,
    force=true,
    jump_model=jump_model
)
# ------------------------------------------------------------------------------
