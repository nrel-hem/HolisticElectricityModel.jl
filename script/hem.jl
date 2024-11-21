using Revise
# using Combinatorics
# using Ipopt
using JuMP
using YAML

using HolisticElectricityModel

# Load config file
if length(ARGS) > 0 # otherwise, define this variable in the REPL
    config_fp = ARGS[0]
else
    config_fp = joinpath(@__DIR__, "configs", "hem_config.yaml")
end
config = YAML.load_file(config_fp)

# File locations
base_dir = abspath(joinpath(dirname(Base.find_package("HolisticElectricityModel")), ".."))

include(joinpath(base_dir, "script", "config.jl"))
include(joinpath(base_dir, "script", "config_data.jl"))

input_dir, = parse(config, "data_selection", validators)

# configure solver
solver, = parse(config, "simulation_parameters", validators)
@info "Running on environment $(splitpath(Base.active_project())[end-1]) with solver $(solver)"
import_solver_package(solver)

# ------------------------------------------------------------------------------
# Model configuration
# ------------------------------------------------------------------------------

market_structure, der_use_case, supply_choice_use_case, der_aggregation_use_case = parse(config, "hem_options", validators)
ipp_algorithm, = parse(config, "ipp_options", validators)
rate_design, net_metering_policy = parse(config, "regulator_options", validators)
pv_adoption_type, = parse(config, "customer_options", validators)

# Define HEM run options
hem_opts = HEMOptions(market_structure, der_use_case, supply_choice_use_case, der_aggregation_use_case)

regulator_opts = RegulatorOptions(rate_design, net_metering_policy)

# Get the optimizer depending on the solver defined the config
ipp_opts = IPPOptions(
    ipp_algorithm,
    Dict(
        "Lagrange_Sub_Investment_Retirement_Cap" => JuMP.optimizer_with_attributes(
            () -> get_optimizer_for_solver("Ipopt"),
            "print_level" => 0,
            # "tol" => 1e-6,
            # "max_iter" => 500,
        ),
        "Lagrange_Sub_Dispatch_Cap" => JuMP.optimizer_with_attributes(
            () -> get_optimizer_for_solver(solver),
            # "OUTPUTLOG" => 0,
        ),
        "Lagrange_Feasible_Cap" => JuMP.optimizer_with_attributes(
            () -> get_optimizer_for_solver(solver),
            "Presolve" => 0,
            # "OUTPUTLOG" => 0,
        ),
        "solve_agent_problem_ipp_cap" => JuMP.optimizer_with_attributes(
            () -> get_optimizer_for_solver(solver),
            "Presolve" => 1,
            # "OUTPUTLOG" => 0,
        ),
        "solve_agent_problem_ipp_mppdc" => JuMP.optimizer_with_attributes(
            () -> get_optimizer_for_solver(solver),
            "Presolve" => 1,
            "BarHomogeneous" => 1,
            # "NumericFocus" => 3,
            # "OUTPUTLOG" => 0,
        ),
        "solve_agent_problem_ipp_mppdc_mccormic_lower" =>
            JuMP.optimizer_with_attributes(
                () -> get_optimizer_for_solver(solver),
                "Presolve" => 1,
                "BarHomogeneous" => 1,
                # "OUTPUTLOG" => 0,
            ),
    ),
)

utility_opts = UtilityOptions(JuMP.optimizer_with_attributes(
    () -> get_optimizer_for_solver(solver),
    # "OUTPUTLOG" => 0,
))

green_developer_opts = GreenDeveloperOptions(
    JuMP.optimizer_with_attributes(
        () -> get_optimizer_for_solver(solver)
        # "OUTPUTLOG" => 0,
    ),
)

customer_opts = CustomerOptions(
    pv_adoption_type,
    JuMP.optimizer_with_attributes(
        () -> get_optimizer_for_solver(solver),
        # "OUTPUTLOG" => 0,
    ),
)

dera_opts = DERAggregatorOptions(
    JuMP.optimizer_with_attributes(
        () -> get_optimizer_for_solver(solver),
        # "OUTPUTLOG" => 0,
    ),
)

jump_model = []

# Run HEM
output_dir = run_hem(
    input_dir,
    hem_opts;
    regulator_options=regulator_opts,
    ipp_options=ipp_opts,
    utility_options=utility_opts,
    green_developer_options=green_developer_opts,
    customer_options=customer_opts,
    dera_options=dera_opts,
    force=true,
    jump_model=jump_model,
)
# ------------------------------------------------------------------------------
