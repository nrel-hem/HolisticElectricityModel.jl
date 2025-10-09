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
include(joinpath(base_dir, "script", "parse_options.jl"))

input_dir, year_start, delta_t, stage_1_results_dir = parse(config, "DataSelection", validators)

if !(isnothing(stage_1_results_dir))
    @info "Using stage 1 results from path: $stage_1_results_dir"
end


# configure solver
solver, = parse(config, "SimulationParameters", validators)
@info "Running on environment $(splitpath(Base.active_project())[end-1]) with solver $(solver)"
import_solver_package(solver)

# ------------------------------------------------------------------------------
# Model configuration
# ------------------------------------------------------------------------------

# Define HEM run options

market_structure, der_use_case, supply_choice_use_case,
der_aggregation_use_case = parse(config, "HEMOptions", validators)

output_dir, = "RunOptions" in keys(config) ? parse(config, "RunOptions", validators) : (nothing,)

hem_opts = HEMOptions(
    market_structure,
    der_use_case,
    supply_choice_use_case,
    der_aggregation_use_case
)

# Define agent options

agent_options = get_agent_options(config, hem_opts, solver)

jump_model = []

# Run HEM
resolved_output_dir = run_hem(
    input_dir,
    stage_1_results_dir,
    hem_opts;
    agent_options,
    force=true,
    jump_model=jump_model,
    output_dir=output_dir
)

save_config(
    solver,
    input_dir,
    hem_opts,
    agent_options,
    resolved_output_dir
)
# ------------------------------------------------------------------------------
