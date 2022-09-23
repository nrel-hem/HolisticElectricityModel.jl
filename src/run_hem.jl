"""
Solve the problem with the given inputs.

# Arguments
- `input_dir::AbstractString`: Directory containing input data. Outputs will be recorded in
  a subdirectory.
- `options::HEMOptions`:
- `ipp_options::IPPOptions`: Independent Power Producer options. Defaults to NullAgentOptions.
- `regulator_options::RegulatorOptions`: Regulator options. Defaults to NullAgentOptions.
- `utility_options::RegulatorOptions`: Regulator options. Defaults to NullAgentOptions.
- `max_iterations::Int`: Max number of iterations to attempt a solution. Defaults to 100.
- `window_length::Int`:
- `force::Bool`: If true, overwrite results if a directory already exists.
"""
function run_hem(
    input_dir::AbstractString,
    options::HEMOptions;
    ipp_options=NullAgentOptions(),
    regulator_options=NullAgentOptions(),
    utility_options=NullAgentOptions(),
    max_iterations=100,
    window_length=1,
    force=false,
)
    model_data = HEMData(input_dir)
    regulator = Regulator(input_dir, model_data)
    utility = Utility(input_dir, model_data, regulator)
    customers = CustomerGroup(input_dir, model_data)
    ipp = IPPGroup(input_dir, model_data)
    green_developer = GreenDeveloper(input_dir, model_data)
    distribution_utility = DistributionUtility(input_dir, model_data)
    agents_and_opts = [
        AgentAndOptions(utility, utility_options),
        AgentAndOptions(ipp, ipp_options),
        AgentAndOptions(regulator, regulator_options),
        AgentAndOptions(customers, NullAgentOptions()),
        AgentAndOptions(green_developer, NullAgentOptions()),
        AgentAndOptions(distribution_utility, NullAgentOptions()),
    ]

    output_dir = joinpath(input_dir, get_file_prefix(options, agents_and_opts))
    if isdir(output_dir)
        if force
            rm(output_dir, recursive=true)
        else
            error("$output_dir already exists. Move the existing directory or set force=true to overwrite.")
        end
    end
    mkdir(output_dir)
    logger = configure_logging(
        console_level = Logging.Info,
        file_level = Logging.Info,
        filename = joinpath(output_dir, "run_hem.log"),
    )
    try
        @info "Output directory: $(output_dir)"
        solve_equilibrium_problem!(
            options,
            model_data,
            agents_and_opts,
            output_dir,
            max_iterations,
            window_length,
        )
    finally
        close(logger)
    end
end
