using Revise
# using Combinatorics
# using Ipopt
using JuMP
using YAML

using HolisticElectricityModel

# Load config file
if length(ARGS) > 0 # otherwise, define this variable in the REPL
    config_fp = ARGS[1]
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

# Define HEM run options

market_structure, der_use_case, supply_choice_use_case, 
    der_aggregation_use_case = parse(config, "hem_options", validators)

hem_opts = HEMOptions(
    market_structure, 
    der_use_case, 
    supply_choice_use_case, 
    der_aggregation_use_case
)

# Define agent options

rate_design, net_metering_policy, tou_suffix, planning_reserve_margin, 
    allowed_return_on_investment = parse(config, "regulator_options", validators)

regulator_opts = RegulatorOptions(
    rate_design, 
    net_metering_policy; 
    tou_suffix=tou_suffix,
    planning_reserve_margin=planning_reserve_margin,
    allowed_return_on_investment=allowed_return_on_investment,
)

# Get the optimizer depending on the solver defined the config
ipp_algorithm, = parse(config, "ipp_options", validators)

ipp_solvers = Dict()
addsolvers_ipp!(ipp_solvers,:Ipopt)
addsolvers_ipp!(ipp_solvers,solver)
ipp_opts = IPPOptions(ipp_algorithm, ipp_solvers)

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

pv_adoption_type, = parse(config, "customer_options", validators)

customer_opts = CustomerOptions(
    pv_adoption_type,
    JuMP.optimizer_with_attributes(
        () -> get_optimizer_for_solver(solver),
        # "OUTPUTLOG" => 0,
    ),
)

incentive_curve, frac_viu_cost_savings_as_revenue = 
    parse(config, "der_aggregator_options", validators)

dera_opts = DERAggregatorOptions(
    JuMP.optimizer_with_attributes(
        () -> get_optimizer_for_solver(solver),
        # "OUTPUTLOG" => 0,
    );
    incentive_curve=incentive_curve,
    frac_viu_cost_savings_as_revenue=frac_viu_cost_savings_as_revenue
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
