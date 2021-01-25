# This module defines inputs that are held in common across all agents

struct HEMData
    # Configuration
    epsilon::AbstractFloat # iteration tolerance

    # Sets
    index_t::Set1D # time index, currently 17 ReEDS timeslices
    index_h::Set1D # customer types

    # Parameters
    omega::Dict{Symbol, AbstractFloat} # number of hours per timeslice
end

function HEMData(input_filename::String)
    # 17 timeslices (from ReEDS)
    index_t = read_set(input_filename, "index_t", "index_t",
                       prose_name = "time index t", 
                       description = "ReEDS 17 timeslices representation")
    
    return HEMData(
        1.0E-3,
        index_t, 
        read_set(input_filename, "index_h", "index_h", 
                 prose_name = "customer types h"),
        read_param(input_filename, "Omega", index_t) # number of hours per timeslice
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


function get_agent(agents::Vector{Agent}, agent_type::Type{T}) where {T <: Agent}
    result = filter(x -> x isa agent_type, agents)
    @assert length(result) == 1 "$agents does not have exactly one $T, it has $(length(result))"
    return first(result)
end


function save_results(
        agent::Agent,
        agent_opts::AgentOptions,
        hem_opts::HEMOptions,
        exportfilepath::AbstractString,
        fileprefix::AbstractString)
    @info "No results defined for $(typeof(agent)) agents when $hem_opts and $agent_opts"
    return
end


SolveAgentCallInfo = @NamedTuple{agent::Agent, options::AgentOptions, other_agents::Vector{Agent}}


function get_call_info(agents_and_opts::Vector{AgentAndOptions})
    result = Vector{SolveAgentCallInfo}()
    for (i, agent_and_opts) in enumerate(agents_and_opts)
        other_agents = Agent[val.agent for (j, val) in enumerate(agents_and_opts) if j != i]
        push!(result, SolveAgentCallInfo((agent_and_opts.agent, agent_and_opts.options, other_agents)))
    end
    return result
end


function solve_equilibrium_problem(hem_opts::HEMOptions, model_data::HEMData, 
        agents_and_opts::Vector{AgentAndOptions}, exportfilepath::AbstractString, 
        fileprefix::AbstractString)
    iter = 1
    max_iter = 10
    diff = 100.0

    call_info = get_call_info(agents_and_opts)

    while diff >= model_data.epsilon
        diff = 0.0

        for (agent, options, other_agents) in call_info
            @info "$(typeof(agent)), iteration $iter"
            diff += solve_agent_problem(agent, options, model_data, hem_opts, other_agents)
        end

        iter += 1
        @info "Iteration $iter value: $diff"
        if iter == max_iter
            stop()
        end
    end

    welfare = []
    for (agent, options, other_agents) in call_info
        # push!(welfare, welfare_calculation(agent, options, model_data, hem_opts, other_agents))
        save_results(agent, options, hem_opts, exportfilepath, fileprefix)
    end
    
    # save_welfare(welfare, exportfilepath, fileprefix)
    @info "Problem solved!"
end


# TODO: Write the welfare calculation and saving more generally
function save_welfare(Supply::Any, Demand::Any, exportfilepath::AbstractString, fileprefix::AbstractString)
    save_param(Demand[1], [:CustomerType, :DERTech], :PVNetCS_dollar, joinpath(exportfilepath, "$(fileprefix)_$(marketstructure)_$(retailrate)_$(dernetmetering)_PVNetCS.csv"))
    save_param(Demand[2], [:CustomerType, :DERTech], :PVEnergySaving_dollar, joinpath(exportfilepath, "$(fileprefix)_$(marketstructure)_$(retailrate)_$(dernetmetering)_PVSaving.csv"))
    save_param(Demand[3], [:CustomerType, :DERTech], :EnergyCost_dollar, joinpath(exportfilepath, "$(fileprefix)_$(marketstructure)_$(retailrate)_$(dernetmetering)_EnergyCost.csv"))
    save_param(Demand[4], [:CustomerType, :DERTech], :NetCS_dollar, joinpath(exportfilepath, "$(fileprefix)_$(marketstructure)_$(retailrate)_$(dernetmetering)_NetCS.csv"))
    SocialWelfare = DataFrame(SupplierRevenue = Supply[1], SupplierCost = Supply[2], TotalNetCS = Demand[5])
    CSV.write(joinpath(exportfilepath, "$(fileprefix)_$(marketstructure)_$(retailrate)_$(dernetmetering)_socialwelfare.csv"), SocialWelfare)
end
