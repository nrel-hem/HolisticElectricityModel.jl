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
    T <: MarketStructure, 
    U <: Union{NullUseCase,DERAdoption},
    V <: NullUseCase,
    W <: Union{NullUseCase,DERAggregation}
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
   "$(typeof(options.market_structure))"],"_")
end


mutable struct HEMData
    # Configuration
    epsilon::ParamScalar # iteration tolerance

    # Sets
    index_y::Dimension # year index
    index_y_fix::Dimension # year index
    index_d::Dimension # representative day index
    index_t::Dimension # time index (within each representative day)
    index_h::Dimension # customer types
    index_z::Dimension # zone index
    index_sector::Dimension # sector index (for rate-making)

    index_h_sector_map::DimensionSet # map from customer group to sector
    index_z_h_map::DimensionSet      # map from zone to cutomer group

    # Parameters
    omega::ParamArray # number of hours per timeslice
    year::ParamArray
    time::ParamArray
    # TODO: Define with kwarg to constructor
    year_start::ParamScalar
end

# TODO: Change input_filename to input_dir
function HEMData(input_filename::String; epsilon::AbstractFloat = 1.0E-3)
    # simulation year index
    index_y = read_set(
        input_filename,
        "index_y",
        "index_y",
        prose_name = "simulation year index y",
        description = "simulation years",
    )
    # TODO: Move information like the below to places where it will get captured in documentation
    # "index_y_fix" represents the full simulation horizon (does not change)
    # "index_y" represents the simulation years in a particular window (gets updated in solve_equilibrium_problem!)
    # e.g., when we simulate years 2021-2030, "index_y_fix" will be [2021, ..., 2030]
    # if the planning window is 5-year for utility or IPPs, so the first index_y will be
    # [2021, ..., 2025], after solving the first window, index_y will be updated to [2022, ..., 2026] etc.
    index_y_fix = read_set(
        input_filename,
        "index_y",
        "index_y_fix",
        prose_name = "simulation year index y",
        description = "simulation years",
    )

    # TODO: Generalize descriptions
    # representative day and hour (from ReEDS)
    index_d = read_set(
        input_filename,
        "index_d",
        "index_d",
        prose_name = "representative day index d",
        description = "ReEDS representative days",
    )

    index_t = read_set(
        input_filename,
        "index_t",
        "index_t",
        prose_name = "time index t",
        description = "ReEDS representative hour within each representative day",
    )

    # customer group types
    index_h = read_set(
        input_filename,
        "index_h",
        "index_h",
        prose_name = "customer group index h",
        description = "customer groups",
    )

    # zones
    index_z = read_set(
        input_filename,
        "index_z",
        "index_z",
        prose_name = "zones index z",
        description = "ReEDS BA modeled",
    )

    # customer group types
    index_sector = read_set(
        input_filename,
        "index_sector",
        "index_sector";
        prose_name = "customer group index sector",
        description = "customer high level groups",
    )

    index_h_sector_map = read_set(
        "index_h_sector_map",
        input_filename,
        "index_h_sector_mapping",
        [index_h, index_sector];
        prose_name = "Map from index_h to index_sector",
        description = "Defines which customer groups (load and DER adoption) are in each sector (for ratemaking)"
    )

    index_z_h_map = read_set(
        "index_z_h_map",
        input_filename,
        "index_z_h_mapping",
        [index_z, index_h];
        prose_name = "Map from index_z to index_h",
        description = "Defines which customer groups (load and DER participation) are present in each zone (bulk power BA)"
    )

    omega = read_param(
        "omega",
        input_filename,
        "Omega",
        index_d,
        description = "number of days per representative day"
    )
    year = read_param(
        "year",
        input_filename,
        "Year",
        index_y,
        description = "Year"
    )
    time = read_param(
        "time",
        input_filename,
        "Time",
        index_t,
        description = "Time"
    )
    # TODO: Remove hard-coding. (This start year is also specified in HEMData.jl)
    # Perhaps requires loading the HEMData config, which currently isn't stored in
    # the input directory.
    year_start = ParamScalar("year_start", 2020, description = "simulation start year")

    # Return HEMData, passing the constructed h_to_sector
    return HEMData(
        ParamScalar("epsilon", epsilon, description = "iteration tolerance"),
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
        year_start
    )
end

# TODO: Maybe convert delta_t to a parameter, since index_t doesn't have to have the
#       structure implied by this function.
function get_delta_t(model_data::HEMData)
    return (
        parse(Int64, chop(string(model_data.index_t.elements[2]), head = 1, tail = 0)) - 
        parse(Int64, chop(string(model_data.index_t.elements[1]), head = 1, tail = 0))
    )
end

# TODO: Document the functions that follow.

function get_reg_year(model_data::HEMData)
    reg_year = model_data.year(first(model_data.index_y))
    return reg_year, Symbol(Int(reg_year))
end

function get_prev_reg_year(model_data::HEMData, w_iter::Integer)
    if w_iter >= 2
        prev_reg_year = model_data.year(first(model_data.index_y)) - 1
    else
        prev_reg_year = model_data.year(first(model_data.index_y))
    end
    return prev_reg_year, Symbol(Int(prev_reg_year))
end

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
    data::Dict{DataType, AgentOptions}
end

function get_agent_option(::Type{T}, options::AgentOptionsStore) where T <: AbstractAgent
    if haskey(options.data, T)
        return options.data[T]
    else
        error("No agent options found for agent type: $(T).")
    end
end

struct AgentAndOptions{T <: AbstractAgent, U <: AgentOptions}
    agent::T
    options::U
end

AgentOrOptions = Union{AbstractAgent, Options}

struct AgentStore
    data::OrderedDict{DataType, OrderedDict{String, AgentAndOptions}}
end

function AgentStore(agents_and_opts::Vector{AgentAndOptions})
    data = OrderedDict{DataType, OrderedDict{String, AgentAndOptions}}()
    for item in agents_and_opts
        type = typeof(item.agent)
        id = get_id(item.agent)
        if haskey(data, type)
            sub_dict = data[type]
            haskey(sub_dict, id) && error("$type agent with ID = $id is already stored")
            sub_dict[id] = item
        else
            data[type] = OrderedDict{String, AgentAndOptions}()
            data[type][id] = item
        end
    end

    return AgentStore(data)
end

"""
Return the agent of the given type and ID from the store.

If there is only one agent of the given type then `id` is optional.
"""
function get_agent(::Type{T}, store::AgentStore, id = nothing) where {T <: AbstractAgent}
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

function get_option(::Type{T}, store::AgentStore, id = nothing) where {T <: AbstractAgent}
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

function iter_agents_and_options(store::AgentStore)
    return ((x.agent, x.options) for agents in values(store.data) for x in values(agents))
end

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
    file_prefix = string("Results_",join(file_prefix, "_"))
    return file_prefix
end

function save_results(
    agent::AbstractAgent,
    agent_opts::AgentOptions,
    hem_opts::HEMOptions,
    export_file_path::AbstractString,
)
    @info "No results defined for $(typeof(agent)) agents when $hem_opts and $agent_opts"
    return
end

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
        for w in 1:(length(model_data.index_y_fix) - window_length + 1) # loop over windows
            model_data.index_y.elements =
                model_data.index_y_fix.elements[w:(w + window_length - 1)]
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
