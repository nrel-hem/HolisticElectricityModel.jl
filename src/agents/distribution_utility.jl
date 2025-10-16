abstract type AbstractDistributionUtilityOptions <: AgentOptions end

struct DistributionUtilityOptions <: AbstractDistributionUtilityOptions
    solvers::HEMSolver
    stage_1_results_dir::String
    
    # Add other options as needed
end

function DistributionUtilityOptions(attributes::MOI.OptimizerWithAttributes, stage_1_results_dir::String)
    return DistributionUtilityOptions(AnySolver(attributes), stage_1_results_dir)
end

abstract type AbstractDistributionUtility <: AbstractAgent end

mutable struct DistributionUtility <: AbstractDistributionUtility
    id::String
    current_year::Symbol

    # other fields as needed
end

function DistributionUtility(input_dir::AbstractString, model_data::HEMData, distribution_utility_options::DistributionUtilityOptions; id = DEFAULT_ID)
    return DistributionUtility(id, first(model_data.index_y))
end

get_id(x::DistributionUtility) = x.id

function solve_agent_problem!(
    distribution_utility::DistributionUtility,
    distribution_utility_opts::DistributionUtilityOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{LocalDistributionAndDER},
    agent_store::AgentStore,
    w_iter,
    window_length,
    jump_model,
    export_file_path,
    update_results::Bool,
    output_intermediate_results::Bool
)
    # Implement the logic here
    @info("Solving problem for Distribution Utility: $(distribution_utility.id)")

    # get load profiles from the customer agent

    # compute distribution system (upgrade) costs

    # make sure costs are in the form expected by the LocalDistributionRateMaker
    
    return 0.0
end

function save_results(
    distribution_utility::DistributionUtility, 
    distribution_utility_opts::DistributionUtilityOptions,
    hem_opts::HEMOptions,
    export_file_path::AbstractString)

    @info("Saving results for Distribution Utility: $(distribution_utility.id)")
    # Implement the logic to save results here
end