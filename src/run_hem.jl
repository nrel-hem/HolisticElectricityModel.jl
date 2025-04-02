function create_agents_and_opts(input_dir::AbstractString, model_data::HEMData, agent_options::Dict{<:AbstractString,AgentOptions}, ::HEMOptions{VIU})
    regulator = Regulator(input_dir, model_data, agent_options["regulator_options"])
    utility = Utility(input_dir, model_data)
    customers = CustomerGroup(input_dir, model_data)
    green_developer = GreenDeveloper(input_dir, model_data)
    dera = DERAggregator(input_dir, model_data, agent_options["der_aggregator_options"])
    # distribution_utility = DistributionUtility(input_dir, model_data)

    # the sequence of simulation matters a lot! (e.g., the year DER aggregation is picked is dependent on this)
    agents_and_opts = [
        AgentAndOptions(utility, agent_options["utility_options"]),
        AgentAndOptions(regulator, agent_options["regulator_options"]),
        AgentAndOptions(customers, agent_options["customer_options"]),
        AgentAndOptions(green_developer, agent_options["green_developer_options"]),
        AgentAndOptions(dera, agent_options["der_aggregator_options"]),
        # AgentAndOptions(distribution_utility, NullAgentOptions()),
    ]

    return agents_and_opts
end

function create_agents_and_opts(input_dir::AbstractString, model_data::HEMData, agent_options::Dict{<:AbstractString,AgentOptions}, ::HEMOptions{WM})
    regulator = Regulator(input_dir, model_data, agent_options["regulator_options"])
    ipp = IPPGroup(input_dir, model_data)
    customers = CustomerGroup(input_dir, model_data)
    green_developer = GreenDeveloper(input_dir, model_data)
    dera = DERAggregator(input_dir, model_data, agent_options["der_aggregator_options"])
    # distribution_utility = DistributionUtility(input_dir, model_data)

    # the sequence of simulation matters a lot! (e.g., the year DER aggregation is picked is dependent on this)
    agents_and_opts = [
        AgentAndOptions(ipp, agent_options["ipp_options"]),
        AgentAndOptions(regulator, agent_options["regulator_options"]),
        AgentAndOptions(customers, agent_options["customer_options"]),
        AgentAndOptions(green_developer, agent_options["green_developer_options"]),
        AgentAndOptions(dera, agent_options["der_aggregator_options"]),
        # AgentAndOptions(distribution_utility, NullAgentOptions()),
    ]

    return agents_and_opts
end

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
    agent_options::Dict{<:AbstractString,AgentOptions},
    max_iterations=1,
    window_length=1,
    force=false,
    jump_model::Any
)
    model_data = HEMData(input_dir)

    agents_and_opts = create_agents_and_opts(input_dir, model_data, agent_options, options)

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
        console_level=Logging.Info,
        file_level=Logging.Info,
        filename=joinpath(output_dir, "run_hem.log"),
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
