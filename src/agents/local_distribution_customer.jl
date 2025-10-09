struct LocalDistributionCustomerOptions <: AbstractCustomerOptions
    solvers::HEMSolver
    adoption_rate_file_index::Int
    # Add other options as needed
end

function LocalDistributionCustomerOptions(
    attributes::MOI.OptimizerWithAttributes,
    adoption_rate_file_index::Int = 1,
)
    return LocalDistributionCustomerOptions(AnySolver(attributes), adoption_rate_file_index)
end

mutable struct LocalDistributionCustomer <: AbstractCustomerGroup
    id::String
    current_year::Symbol
    previous_year::Symbol

    # Sets
    "PV, EV, Heat pump technologies"
    index_m::Dimension

    # Parameters
    "number of customers for each customer type"
    gamma::ParamArray

    # other fields as needed
end

function LocalDistributionCustomer(input_dir::AbstractString, stage_1_results_dir::AbstractString, model_data::HEMData; id = DEFAULT_ID)
    # change this to read from input files
    index_m = Dimension("m_local", [:BTMPV, :EV, :HeatPump])
    gamma = ParamArray("gamma_local", (model_data.index_z, model_data.index_h), zeros((length(model_data.index_z), length(model_data.index_h))) )
    return LocalDistributionCustomer(id, first(model_data.index_y), first(model_data.index_y), index_m, gamma)
end

get_id(x::LocalDistributionCustomer) = x.id

function solve_agent_problem!(
    customers::LocalDistributionCustomer,
    customer_opts::LocalDistributionCustomerOptions,
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
    @info("Solving problem for Local Distribution Customer: $(customers.id)")

    # use exogenous DER adoption levels and year to calculate load profiles and get DER information

    # finalize DER profiles--could use utility rate information

    # calculate retail bill based on rates reported by the rate maker, load and DER profiles

    return 0.0
end


function save_results(
    customers::LocalDistributionCustomer,
    customer_opts::LocalDistributionCustomerOptions,
    hem_opts::HEMOptions,
    export_file_path::AbstractString)

    @info("Saving results for Local Distribution Customer: $(customers.id)")
    # Implement the logic to save results here
end