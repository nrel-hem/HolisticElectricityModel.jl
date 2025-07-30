# This module defines inputs that are held in common across all agents

const DEFAULT_ID = "default"
const HEM_TIMER = TimerOutputs.TimerOutput()

abstract type Options end

get_file_prefix(::Options) = String("")

# Struct with no fields used to dispatch -- this is the traits pattern
abstract type MarketStructure end
struct VIU <: MarketStructure end
struct WM <: MarketStructure end
struct RetailMarket <: MarketStructure end

abstract type UseCase end
struct NullUseCase <: UseCase end
struct DERAdoption <: UseCase end
struct DERAggregation <: UseCase end

struct HEMOptions{
    T<:MarketStructure,
    U<:Union{NullUseCase,DERAdoption},
    V<:NullUseCase,
    W<:Union{NullUseCase,DERAggregation}
} <: Options
    # market structure switch
    market_structure::T

    # use case switches
    der_use_case::U
    supply_choice_use_case::V
    der_aggregation_use_case::W
end

# TODO: Rethink file prefixes to create shorter directory names
#       If do this, might need to be able to provide mapping function when doing integration testing
function get_file_prefix(options::HEMOptions)
    return join([# "$(typeof(options.der_use_case))", 
            # "$(typeof(options.supply_choice_use_case))",
            "$(typeof(options.der_aggregation_use_case))",
            "$(typeof(options.market_structure))"], "_")
end


"""
    $(TYPEDEF)
Struct to store the data required for the HEM model.
$(TYPEDFIELDS)
"""
mutable struct HEMData
    # Configuration
    "iteration tolerance"
    epsilon::ParamScalar

    # Sets
    "Simulation years in a particular window (gets updated in solve_equilibrium_problem!)"
    index_y::Dimension
    """
    Represents the full simulation horizon (does not change).
    E.g., when we simulate years 2021-2030, "index_y_fix" will be [2021, ..., 2030]
    if the planning window is 5-year for utility or IPPs, so the first index_y will be
    [2021, ..., 2025], after solving the first window, index_y will be updated to [2022, ..., 2026] etc.
    """
    index_y_fix::Dimension
    "Representative day (from ReEDS)"
    index_d::Dimension
    "Representative hour index within each representative day"
    index_t::Dimension
    "Customer groups"
    index_h::Dimension
    "Zones (ReEDS BA modeled)"
    index_z::Dimension
    "Customer high level groups"
    index_sector::Dimension
    "Map from customer group to sector"
    index_h_sector_map::DimensionSet
    "Map from zone to customer group"
    index_z_h_map::DimensionSet

    # Parameters
    "Number of days per representative day"
    omega::ParamArray
    "Year"
    year::ParamArray
    "Time"
    time::ParamArray
    "Simulation start year"
    year_start::ParamScalar
    "Number of hours per representative hour"
    delta_t::ParamScalar
end

"""
Create a new `HEMData` instance.
"""
function HEMData(input_dir::String; year_start::Int=2020, delta_t::Int=4, epsilon::AbstractFloat=1.0E-3)

    index_y = read_set(
        input_dir,
        "index_y",
        "index_y",
        prose_name="simulation year index y",
        description="simulation years",
    )

    index_y_fix = read_set(
        input_dir,
        "index_y",
        "index_y_fix",
        prose_name="simulation year index y",
        description="simulation years",
    )

    index_d = read_set(
        input_dir,
        "index_d",
        "index_d",
        prose_name="representative day index d",
        description="ReEDS representative days",
    )

    index_t = read_set(
        input_dir,
        "index_t",
        "index_t",
        prose_name="time index t",
        description="ReEDS representative hour within each representative day",
    )

    index_h = read_set(
        input_dir,
        "index_h",
        "index_h",
        prose_name="customer group index h",
        description="customer groups",
    )

    index_z = read_set(
        input_dir,
        "index_z",
        "index_z",
        prose_name="zones index z",
        description="ReEDS BA modeled",
    )

    index_sector = read_set(
        input_dir,
        "index_sector",
        "index_sector";
        prose_name="customer group index sector",
        description="customer high level groups",
    )

    index_h_sector_map = read_set(
        "index_h_sector_map",
        input_dir,
        "index_h_sector_mapping",
        [index_h, index_sector];
        prose_name="Map from index_h to index_sector",
        description="Defines which customer groups (load and DER adoption) are in each sector (for ratemaking)"
    )

    index_z_h_map = read_set(
        "index_z_h_map",
        input_dir,
        "index_z_h_mapping",
        [index_z, index_h];
        prose_name="Map from index_z to index_h",
        description="Defines which customer groups (load and DER participation) are present in each zone (bulk power BA)"
    )

    omega = read_param(
        "omega",
        input_dir,
        "Omega",
        index_d,
        description="number of days per representative day"
    )
    year = read_param(
        "year",
        input_dir,
        "Year",
        index_y,
        description="Year"
    )
    time = read_param(
        "time",
        input_dir,
        "Time",
        index_t,
        description="Time"
    )

    year_start = ParamScalar("year_start", year_start, description="simulation start year")

    delta_t = ParamScalar("delta_t", delta_t, description="number of hours per representative hour")

    return HEMData(
        ParamScalar("epsilon", epsilon, description="iteration tolerance"),
        index_y,
        index_y_fix,
        index_d,
        index_t,
        index_h,
        index_z,
        index_sector,
        index_h_sector_map,
        index_z_h_map,
        omega,
        year,
        time,
        year_start,
        delta_t
    )
end

"""
Returns the current year and its Symbol representation from the model data.
"""
function get_reg_year(model_data::HEMData)
    reg_year = model_data.year(first(model_data.index_y))
    return reg_year, Symbol(Int(reg_year))
end

"""
Returns the current year and its Symbol representation from the model data,
taking into account the window iteration.
"""
function get_prev_reg_year(model_data::HEMData, w_iter::Integer)
    if w_iter >= 2
        prev_reg_year = model_data.year(first(model_data.index_y)) - 1
    else
        prev_reg_year = model_data.year(first(model_data.index_y))
    end
    return prev_reg_year, Symbol(Int(prev_reg_year))
end

"""
Returns the year two iterations before the current one and its Symbol representation,
taking into account the window iteration. If `w_iter` is 2, it returns the previous year.
If `w_iter` is 1, it returns the current year. 
"""
function get_prev_two_reg_year(model_data::HEMData, w_iter::Integer)
    if w_iter >= 3
        prev_reg_year = model_data.year(first(model_data.index_y)) - 2
    elseif w_iter == 2
        prev_reg_year = model_data.year(first(model_data.index_y)) - 1
    else
        prev_reg_year = model_data.year(first(model_data.index_y))
    end
    return prev_reg_year, Symbol(Int(prev_reg_year))
end

# TODO: Check that save_results argument names make sense

"""
    $(TYPEDEF)
Abstract type for agents.

Required interfaces:
- get_id(agent::Agent)::String
- get_current_year(agent::Agent, model_data::HEMData)::Tuple{Float64,Symbol}
- solve_agent_problem!(
      agents::AgentGroup,
      agent_opts::AgentOptions,
      model_data::HEMData,
      hem_opts::HEMOptions,
      agent_store::AgentStore,
  )
- save_results(
    agent::AbstractAgent,
    agent_opts::AgentOptions,
    hem_opts::HEMOptions,
    export_file_path::AbstractString,
    file_prefix::AbstractString,
)
"""
abstract type AbstractAgent end

"""
Default no-op method for updating cumulative parameters in AbstractAgents after
each model year.
"""
function update_cumulative!(model_data::HEMData, agent::AbstractAgent)
    return
end

"""
Get the last year for which agent's data have been updated. Provided so other 
agents can access the right data for agent.
"""
function get_current_year(agent::AbstractAgent, model_data::HEMData)
    yr = model_data.year(agent.current_year)
    return yr, Symbol(Int(yr))
end

# There is currently no behavioral difference between the structs AgentGroup and Agent, but
# there may be differences in the future.

"""
Abstract type for a group of individual agents.
"""
abstract type AgentGroup <: AbstractAgent end

"""
Abstract type for all individual agents.
"""
abstract type Agent <: AbstractAgent end

get_file_prefix(::AbstractAgent) = String("")

abstract type AgentOptions <: Options end
struct NullAgentOptions <: AgentOptions end

"""
Struct to store parsed agent options.
"""
struct AgentOptionsStore
    data::Dict{DataType,AgentOptions}
end

"""
Returns the options for the agent type `T` from the `options` store.
If no options are found for the agent type, an error is raised.
"""
function get_agent_option(::Type{T}, options::AgentOptionsStore) where T<:AbstractAgent
    if haskey(options.data, T)
        return options.data[T]
    else
        error("No agent options found for agent type: $(T).")
    end
end

"""
$(TYPEDEF)
Struct to store an agent and its options.
$(TYPEDFIELDS)
"""
struct AgentAndOptions{T<:AbstractAgent,U<:AgentOptions}
    agent::T
    options::U
end

AgentOrOptions = Union{AbstractAgent,Options}

"""
$(TYPEDEF)
Struct to store a collection of agents and their options.
$(TYPEDFIELDS)
This is passed to the agents to allow them to access other agents and their options.
"""
struct AgentStore
    data::OrderedDict{DataType,OrderedDict{String,AgentAndOptions}}
end

"""
Create an `AgentStore` from a vector of `AgentAndOptions`.
"""
function AgentStore(agents_and_opts::Vector{AgentAndOptions})
    data = OrderedDict{DataType,OrderedDict{String,AgentAndOptions}}()
    for item in agents_and_opts
        type = typeof(item.agent)
        id = get_id(item.agent)
        if haskey(data, type)
            sub_dict = data[type]
            haskey(sub_dict, id) && error("$type agent with ID = $id is already stored")
            sub_dict[id] = item
        else
            data[type] = OrderedDict{String,AgentAndOptions}()
            data[type][id] = item
        end
    end

    return AgentStore(data)
end

"""
Return the agent of the given type and ID from the store.

If there is only one agent of the given type then `id` is optional.
"""
function get_agent(::Type{T}, store::AgentStore, id=nothing) where {T<:AbstractAgent}
    !haskey(store.data, T) && error("No agents of type $T are stored.")
    agents_and_opts = store.data[T]

    if id === nothing
        if length(agents_and_opts) > 1
            error("Passing 'id' is required if more than one agent is stored.")
        end
        return first(values(agents_and_opts)).agent
    end

    !haskey(agents_and_opts, id) && error("No agent of type $T id = $id is stored")
    return agents_and_opts[id].agent
end

"""
Return the options for the agent of the given type and ID from the store.
If there is only one agent of the given type then `id` is optional.
"""
function get_option(::Type{T}, store::AgentStore, id=nothing) where {T<:AbstractAgent}
    !haskey(store.data, T) && error("No agents of type $T are stored.")
    agents_and_opts = store.data[T]

    if id === nothing
        if length(agents_and_opts) > 1
            error("Passing 'id' is required if more than one agent is stored.")
        end
        return first(values(agents_and_opts)).options
    end

    !haskey(agents_and_opts, id) && error("No agent of type $T id = $id is stored")
    return agents_and_opts[id].options
end

"""
Iterate over all agents and their options in the store and return a Tuple
    of the agent and its options.
"""
function iter_agents_and_options(store::AgentStore)
    return ((x.agent, x.options) for agents in values(store.data) for x in values(agents))
end

"""
Get the corresponding bulk system agent from the store based on the market structure within HEMOptions.
If the market structure is `VIU`, it returns the `Utility` agent.
If the market structure is `WM`, it returns the `IPPGroup` agent.
"""
function get_bulk_system_agent(store::AgentStore, ::HEMOptions{VIU})
    return get_agent(Utility, store)
end

function get_bulk_system_agent(store::AgentStore, ::HEMOptions{WM})
    return get_agent(IPPGroup, store)
end

function get_file_prefix(hem_opts::HEMOptions, agents_and_opts::Vector{AgentAndOptions})
    # create vector of items that may contribute information
    items = Vector{AgentOrOptions}()
    push!(items, hem_opts)
    for item in agents_and_opts
        push!(items, item.options, item.agent)
    end

    # call get_file_prefix on each item
    file_prefix = Vector{String}()
    for item in items
        val = get_file_prefix(item)
        if !isempty(val)
            push!(file_prefix, val)
        end
    end
    file_prefix = string("Results_", join(file_prefix, "_"))
    return file_prefix
end

"""
This method needs to be implemented by each agent type to save its results.
"""
function save_results(
    agent::AbstractAgent,
    agent_opts::AgentOptions,
    hem_opts::HEMOptions,
    export_file_path::AbstractString,
)
    @info "No results defined for $(typeof(agent)) agents when $hem_opts and $agent_opts"
    return
end

"""
Main function to run the HEM model.
It loops over the simulation years and calls the `solve_agent_problem!` method and the 
`save_results` method for each agent.
Arguments:
- `input_dir::AbstractString`: Directory containing input data. Outputs will be recorded in
  a subdirectory.
- `model_data::HEMData`: Data required for the model.
- `hem_opts::HEMOptions`: Options for the model.
- `agents_and_opts::Vector{AgentAndOptions}`: Vector of agents and their options.
- `export_file_path::AbstractString`: Path to export the results.
- `max_iter::Int64`: Maximum number of iterations to solve the problem. Set in `run_hem`.
- `window_length::Int64`: Length of the window for the simulation. Set in `run_hem`.
- `jump_model::Any`: Jump model to use for the simulation. Iniitialized in `run_hem`.
"""
function solve_equilibrium_problem!(
    hem_opts::HEMOptions,
    model_data::HEMData,
    agents_and_opts::Vector{AgentAndOptions},
    export_file_path::AbstractString,
    max_iter::Int64,
    window_length::Int64,
    jump_model::Any
)
    store = AgentStore(agents_and_opts)
    TimerOutputs.reset_timer!(HEM_TIMER)

    TimerOutputs.@timeit HEM_TIMER "solve_equilibrium_problem!" begin
        for w in 1:(length(model_data.index_y_fix)-window_length+1) # loop over windows
            model_data.index_y.elements =
                model_data.index_y_fix.elements[w:(w+window_length-1)]
            i = 0
            diff_iter = []
            for i in 1:max_iter
                diff_vec = []

                for (agent, options) in iter_agents_and_options(store)
                    TimerOutputs.@timeit HEM_TIMER "solve_agent_problem!" begin
                        @info "$(typeof(agent)), iteration $i"
                        diff_one = solve_agent_problem!(
                            agent,
                            options,
                            model_data,
                            hem_opts,
                            store,
                            w,
                            window_length,
                            jump_model,
                            export_file_path,
                            true,
                            false,
                        )
                    end
                    @assert !isnothing(diff_one) "Nothing returned by solve_agent_problem!($(typeof(agent))): $(diff_one)"
                    @assert !(diff_one isa Tuple) "Tuple returned by solve_agent_problem!($(typeof(agent))): $(diff_one)"
                    @info "$(diff_one)"

                    if diff_one isa HolisticElectricityModel.ParamArray
                        push!(diff_vec, maximum(diff_one.values))
                    else
                        push!(diff_vec, diff_one)
                    end
                end
                diff = maximum(diff_vec)
                @info "Iteration $i value: $diff"
                push!(diff_iter, diff)
                @info "Iteration $i value vector: $diff_iter"

                if diff < model_data.epsilon
                    break
                end
            end

            i >= max_iter && error("Reached max iterations $max_iter with no solution")
            @info "Problem solved!"

            update_cumulative!(model_data, agents_and_opts)
        end

    end

    for (agent, options) in iter_agents_and_options(store)
        save_results(agent, options, hem_opts, export_file_path)
    end

    @info "\n$(HEM_TIMER)\n"
end

function update_cumulative!(
    model_data::HEMData,
    agents_and_opts::Vector{AgentAndOptions},
)
    for item in agents_and_opts
        update_cumulative!(model_data, item.agent)
    end
end
