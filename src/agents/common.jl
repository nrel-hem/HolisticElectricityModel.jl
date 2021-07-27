# This module defines inputs that are held in common across all agents

struct HEMData
    # Configuration
    epsilon::ParamScalar # iteration tolerance

    # Sets
    index_t::Set1D # time index, currently 17 ReEDS timeslices
    index_h::Set1D # customer types

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

abstract type Agent end
abstract type AgentOptions end
struct NullAgentOptions <: AgentOptions end

struct AgentAndOptions{T <: Agent, U <: AgentOptions}
    agent::T
    options::U
end

# make this be the data field of every agent type
# still need abstract type (no data) for Agent, then derive from it
# mutable struct AgentData
#     sets::Vector{MySetType}
#     params::Vector{MyParameterType}
#     vars::Vector{MyVariableType}
# end

SolveAgentCallInfo = @NamedTuple {agent::Agent, options::AgentOptions}

struct AgentStore
    data::Dict{DataType, <:AgentAndOptions}
end

function AgentStore(agents_and_opts::Vector{<:AgentAndOptions})
    data = Dict{DataType, AgentAndOptions}()
    for item in agents_and_opts
        type = typeof(item.agent)
        haskey(data, type) && error("$type cannot be stored multiple times")
        data[type] = item
    end

    return AgentStore(data)
end

get_agent(::Type{T}, store::AgentStore) where {T <: Agent} = store.data[T].agent
# TODO DT: could implement Base collection interfaces

function get_call_info(store::AgentStore)
    return [SolveAgentCallInfo((x.agent, x.options)) for x in values(store.data)]
end

function save_results(
    agent::Agent,
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
    return solve_equilibrium_problem!(
        hem_opts,
        model_data,
        store,
        export_file_path,
        file_prefix,
    )
end

function solve_equilibrium_problem!(
    hem_opts::HEMOptions,
    model_data::HEMData,
    agent_store::AgentStore,
    export_file_path::AbstractString,
    file_prefix::AbstractString,
)
    iter = 1
    max_iter = 10
    diff = 100.0

    call_info = get_call_info(agent_store)

    for i in 1:max_iter
        diff = 0.0

        for (agent, options) in call_info
            @info "$(typeof(agent)), iteration $iter"
            diff += solve_agent_problem!(agent, options, model_data, hem_opts, agent_store)
        end
        @info "Iteration $i value: $diff"

        if diff < model_data.epsilon
            break
        end
    end
    # TODO DT: does it matter if we reached max_iter?

    welfare = []
    for (agent, options) in call_info
        # push!(welfare, welfare_calculation(agent, options, model_data, hem_opts, other_agents))
        save_results(agent, options, hem_opts, export_file_path, file_prefix)
    end

    # save_welfare(welfare, export_file_path, file_prefix)
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
