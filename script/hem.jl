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

agent_options = get_agent_options(config, hem_opts)

jump_model = []

# Run HEM
output_dir = run_hem(
    input_dir,
    hem_opts;
    agent_options,
    force=true,
    jump_model=jump_model,
)
# ------------------------------------------------------------------------------
