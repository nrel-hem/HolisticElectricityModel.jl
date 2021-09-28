# This module defines inputs that are held in common across all agents

const DEFAULT_ID = "default"

struct HEMData
    # Configuration
    epsilon::ParamScalar # iteration tolerance

    # Sets
    index_t::Dimension # time index, currently 17 ReEDS timeslices
    index_h::Dimension # customer types

    # Parameters
    omega::ParamVector # number of hours per timeslice
end

function HEMData(input_filename::String; epsilon::AbstractFloat = 1.0E-3)
    # 17 timeslices (from ReEDS)
    index_t = read_set(
        input_filename,
        "index_t",
        "index_t",
        prose_name = "time index t",
        description = "ReEDS 17 timeslices representation",
    )

    return HEMData(
        ParamScalar("epsilon", epsilon, description = "iteration tolerance"),
        index_t,
        read_set(input_filename, "index_h", "index_h", prose_name = "customer types h"),
        read_param(
            "omega",
            input_filename,
            "Omega",
            index_t,
            description = "number of hours per timeslice",
        ),
    )
end

# Struct with no fields used to dispatch -- this is the traits pattern
abstract type MarketStructure end
struct VerticallyIntegratedUtility <: MarketStructure end
struct WholesaleMarket <: MarketStructure end

struct HEMOptions{T <: MarketStructure}
    solver::HEMSolver
    market_structure::T
end

"""
Abstract type for agents.

Required interfaces:
- get_id(agent::Agent)::String
- solve_agent_problem!(
      agents::Agents,
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

# There is currently no behavioral difference between the structs Agents and Agent, but
# there may be differences in the future.

"""
Abstract type for a group of individual agents.
"""
abstract type Agents <: AbstractAgent end

"""
Abstract type for all individual agents.
"""
abstract type Agent <: AbstractAgent end

abstract type AgentOptions end
struct NullAgentOptions <: AgentOptions end

struct AgentAndOptions{T <: AbstractAgent, U <: AgentOptions}
    agent::T
    options::U
end

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

function iter_agents_and_options(store::AgentStore)
    return ((x.agent, x.options) for agents in values(store.data) for x in values(agents))
end

function save_results(
    agent::AbstractAgent,
    agent_opts::AgentOptions,
    hem_opts::HEMOptions,
    export_file_path::AbstractString,
    file_prefix::AbstractString,
)
    @info "No results defined for $(typeof(agent)) agents when $hem_opts and $agent_opts"
    return
end

function solve_equilibrium_problem!(
    hem_opts::HEMOptions,
    model_data::HEMData,
    agents_and_opts::Vector{AgentAndOptions},
    export_file_path::AbstractString,
    file_prefix::AbstractString,
)
    store = AgentStore(agents_and_opts)
    max_iter = 10
    diff = 100.0

    i = 0
    for i in 1:max_iter
        diff = 0.0

        for (agent, options) in iter_agents_and_options(store)
            @info "$(typeof(agent)), iteration $i"
            diff += solve_agent_problem!(agent, options, model_data, hem_opts, store)
        end
        @info "Iteration $i value: $diff"

        if diff < model_data.epsilon
            break
        end
    end

    welfare = []
    for (agent, options) in iter_agents_and_options(store)
        # push!(welfare, welfare_calculation(agent, options, model_data, hem_opts, other_agents))
        save_results(agent, options, hem_opts, export_file_path, file_prefix)
    end

    # save_welfare(welfare, export_file_path, file_prefix)
    i >= max_iter && error("Reached max iterations $max_iter with no solution")
    @info "Problem solved!"
end

# TODO: Write the welfare calculation and saving more generally
function save_welfare(
    Supply::Any,
    Demand::Any,
    export_file_path::AbstractString,
    file_prefix::AbstractString,
)
    save_param(
        Demand[1],
        [:CustomerType, :DERTech],
        :PVNetCS_dollar,
        joinpath(
            export_file_path,
            "$(file_prefix)_$(marketstructure)_$(retailrate)_$(dernetmetering)_PVNetCS.csv",
        ),
    )
    save_param(
        Demand[2],
        [:CustomerType, :DERTech],
        :PVEnergySaving_dollar,
        joinpath(
            export_file_path,
            "$(file_prefix)_$(marketstructure)_$(retailrate)_$(dernetmetering)_PVSaving.csv",
        ),
    )
    save_param(
        Demand[3],
        [:CustomerType, :DERTech],
        :EnergyCost_dollar,
        joinpath(
            export_file_path,
            "$(file_prefix)_$(marketstructure)_$(retailrate)_$(dernetmetering)_EnergyCost.csv",
        ),
    )
    save_param(
        Demand[4],
        [:CustomerType, :DERTech],
        :NetCS_dollar,
        joinpath(
            export_file_path,
            "$(file_prefix)_$(marketstructure)_$(retailrate)_$(dernetmetering)_NetCS.csv",
        ),
    )
    SocialWelfare = DataFrame(
        SupplierRevenue = Supply[1],
        SupplierCost = Supply[2],
        TotalNetCS = Demand[5],
    )
    CSV.write(
        joinpath(
            export_file_path,
            "$(file_prefix)_$(marketstructure)_$(retailrate)_$(dernetmetering)_socialwelfare.csv",
        ),
        SocialWelfare,
    )
end
