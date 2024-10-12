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
    ipp_options::Union{IPPOptions, NullAgentOptions}=NullAgentOptions(),
    regulator_options::Union{RegulatorOptions, NullAgentOptions}=NullAgentOptions(),
    utility_options::Union{UtilityOptions, NullAgentOptions}=NullAgentOptions(),
    green_developer_options::Union{GreenDeveloperOptions, NullAgentOptions}=NullAgentOptions(),
    customer_options::Union{CustomerOptions, NullAgentOptions}=NullAgentOptions(),
    dera_options::Union{DERAggregatorOptions, NullAgentOptions}=NullAgentOptions(),
    max_iterations=1,
    window_length=1,
    force=false,
    jump_model::Any
)
    model_data = HEMData(input_dir)

    regulator = Regulator(input_dir, model_data)
    utility = Utility(input_dir, model_data, regulator)
    customers = CustomerGroup(input_dir, model_data)
    ipp = IPPGroup(input_dir, model_data)
    green_developer = GreenDeveloper(input_dir, model_data)
    dera = DERAggregator(input_dir, model_data)
    # distribution_utility = DistributionUtility(input_dir, model_data)

    # the sequence of simulation matters a lot! (e.g., the year DER aggregation is picked is dependent on this)
    agents_and_opts = [
        AgentAndOptions(utility, utility_options),
        AgentAndOptions(ipp, ipp_options),
        AgentAndOptions(regulator, regulator_options),
        AgentAndOptions(customers, customer_options),
        AgentAndOptions(green_developer, green_developer_options),
        AgentAndOptions(dera, dera_options),
        # AgentAndOptions(distribution_utility, NullAgentOptions()),
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
            jump_model
        )
    finally
        close(logger)
    end
    
    return output_dir
end
