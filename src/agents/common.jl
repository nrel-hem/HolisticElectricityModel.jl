# This module defines inputs that are held in common across all agents
using Xpress

Solver = Xpress

struct HEMData
    # Configuration
    epsilon::AbstractFloat # iteration tolerance

    # Sets
    index_t::Array{Symbol,1} # time index, currently 17 ReEDS timeslices
    index_h::Array{Symbol,1} # customer types

    # Parameters
    omega::Dict{Symbol, AbstractFloat} # number of hours per timeslice
end

function HEMData(input_filename::String)
    index_t = read_set(input_filename, "index_t") # 17 timeslices (from ReEDS)

    return HEMData(
        1.0E-3,
        index_t, 
        read_set(input_filename, "index_h"),
        read_param(input_filename, "Omega", index_t) # number of hours per timeslice
    )
end


# Struct with no fields used to dispatch -- this is the traits pattern
abstract type MarketStructure end
struct VerticallyIntegratedUtility <: MarketStructure end
struct WholesaleMarket <: MarketStructure end


struct HEMOptions{T <: MarketStructure}
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


function save_welfare(Supply::Any, Demand::Any, exportfilepath::AbstractString, fileprefix::AbstractString)
    save_param(Demand[1], [:CustomerType, :DERTech], :PVNetCS_dollar, joinpath(exportfilepath, "$(fileprefix)_$(marketstructure)_$(retailrate)_$(dernetmetering)_PVNetCS.csv"))
    save_param(Demand[2], [:CustomerType, :DERTech], :PVEnergySaving_dollar, joinpath(exportfilepath, "$(fileprefix)_$(marketstructure)_$(retailrate)_$(dernetmetering)_PVSaving.csv"))
    save_param(Demand[3], [:CustomerType, :DERTech], :EnergyCost_dollar, joinpath(exportfilepath, "$(fileprefix)_$(marketstructure)_$(retailrate)_$(dernetmetering)_EnergyCost.csv"))
    save_param(Demand[4], [:CustomerType, :DERTech], :NetCS_dollar, joinpath(exportfilepath, "$(fileprefix)_$(marketstructure)_$(retailrate)_$(dernetmetering)_NetCS.csv"))
    SocialWelfare = DataFrame(SupplierRevenue = Supply[1], SupplierCost = Supply[2], TotalNetCS = Demand[5])
    CSV.write(joinpath(exportfilepath, "$(fileprefix)_$(marketstructure)_$(retailrate)_$(dernetmetering)_socialwelfare.csv"), SocialWelfare)
end


SolveAgentCallInfo = @NamedTuple{agent::Agent, options::AgentOptions, other_agents::Vector{Agent}}


function get_call_info(agents_and_opts::Vector{AgentAndOptions})
    result = Vector{SolveAgentCallInfo}()
    for (i, agent_and_opts) in enumerate(agents_and_opts)
        other_agents = Agent[val.agent for (j, val) in enumerate(agents_and_opts) if not j == i]
        push!(result, SolveAgentCallInfo(agent_and_opts.agent, agent_and_opts.options, other_agents))
    end
end


function solve_equilibrium_problem(hem_opts::HEMOptions{T}, model_data::HEMData, 
        agents_and_opts::Vector{AgentAndOptions}, exportfilepath::AbstractString)
    iter = 1
    max_iter = 10
    diff = 100.0

    call_info = get_call_info(agents_and_opts)

    while diff >= model_data.epsilon
        diff = 0.0

        for (agent, options, other_agents) in call_info
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
        push!(welfare, welfare_calculation(agent, options, model_data, hem_opts, other_agents))
        save_results(agent, exportfilepath, "Results")
    end
    
    save_welfare(welfare, exportfilepath, "Results")
end


function solve_equilibrium_problem(marketstructure::VerticallyIntegratedUtility, retailrate, dernetmetering,
        model_data, regulator, utility, customers, ipp, exportfilepath)

    # ETH@20200818 - Dheepak thinks if we remove all of these global keywords 
    # the code will still work
    global iter = 1
    global max_iter = 10
    global diff = 100.0

    while diff >= model_data.epsilon
        global iter
        global max_iter
        global diff
    
        diff = 0.0
        diff += solve_agent_problem(utility, model_data, regulator, customers)
        diff += solve_agent_problem(regulator, marketstructure, model_data, utility, customers, retailrate, dernetmetering)
        diff += solve_agent_problem(customers, model_data, regulator)
    
        iter += 1
        @info "Iteration $iter value: $diff"
        if iter == max_iter
            stop()
        end
    end

    Welfare_Utility = welfare_calculation(utility, model_data, regulator, customers)
    Welfare_Consumer = welfare_calculation(customers, model_data, regulator)
    
    save_results(regulator, exportfilepath, "Results")
    save_results(utility, exportfilepath, "Results")
    save_results(customers, exportfilepath, "Results")
    save_welfare(Welfare_Utility, Welfare_Consumer, exportfilepath, "Results")
end


function solve_equilibrium_problem(marketstructure::WholesaleMarket, retailrate, dernetmetering, 
        model_data, regulator, utility, customers, ipp, exportfilepath)
    # Julia while-scope ridiculousness https://discourse.julialang.org/t/variable-scope-in-while-loops/15886/3
    global iter = 1
    global max_iter = 20
    global diff = 100.0
    @info "Iteration $iter value: $diff"
        
    while diff >= model_data.epsilon
        global iter
        global max_iter
        global diff
    
        diff = 0.0
        @info "Solving IPP problem"
        diff += solve_agent_problem(ipp, model_data, regulator, customers)
        @info "Solving Regulator problem"
        diff += solve_agent_problem(regulator, model_data, ipp, customers, retailrate, dernetmetering, marketstructure)
        @info "Solving Customers problem"
        diff += solve_agent_problem(customers, model_data, regulator)
    
        iter += 1
        @info "Iteration $iter value: $diff"
        if iter == max_iter
            stop()
        end
    end

    Welfare_IPP = welfare_calculation(ipp, model_data, regulator, customers)
    Welfare_Consumer = welfare_calculation(customers, model_data, regulator)
    
    save_results(regulator, exportfilepath, "Results")
    save_results(ipp, exportfilepath, "Results")
    save_results(customers, exportfilepath, "Results")
    save_welfare(Welfare_IPP, Welfare_Consumer, exportfilepath, "Results")
end

