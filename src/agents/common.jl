# This module defines inputs that are held in common across all agents

const DEFAULT_ID = "default"
const HEM_TIMER = TimerOutputs.TimerOutput()

struct HEMData
    # Configuration
    epsilon::ParamScalar # iteration tolerance

    # Sets
    index_y::Dimension # year index
    index_y_fix::Dimension # year index
    index_s::Dimension # year index (for new resources depreciation schedule)
    index_t::Dimension # time index, currently 17 ReEDS timeslices
    index_h::Dimension # customer types
    index_j::Dimension # green tariff technologies

    # Parameters
    omega::ParamAxisArray # number of hours per timeslice
    year::ParamAxisArray
    hour::ParamAxisArray
    year_start::ParamScalar
end

function HEMData(input_filename::String; epsilon::AbstractFloat = 1.0E-3)
    # simulation year index
    index_y = read_set(
        input_filename,
        "index_y",
        "index_y",
        prose_name = "simulation year index y",
        description = "simulation years",
    )

    # new resource depreciation year index
    index_s = read_set(
        input_filename,
        "index_s",
        "index_s",
        prose_name = "new resource depreciation year index s",
        description = "new resource depreciation years",
    )

    # 17 timeslices (from ReEDS)
    index_t = read_set(
        input_filename,
        "index_t",
        "index_t",
        prose_name = "time index t",
        description = "ReEDS 17 timeslices representation",
    )

    # customer group types
    index_h = read_set(
        input_filename,
        "index_h",
        "index_h",
        prose_name = "customer group index h",
        description = "customer groups",
    )

    # green technology types
    index_j = read_set(
        input_filename,
        "index_j",
        "index_j",
        prose_name = "green technologies index j",
        description = "green tariff technologies",
    )

    return HEMData(
        ParamScalar("epsilon", epsilon, description = "iteration tolerance"),
        index_y,
        index_y,
        index_s,
        index_t,
        index_h,
        index_j,
        read_param(
            "omega",
            input_filename,
            "Omega",
            index_t,
            description = "number of hours per timeslice",
        ),
        read_param("year", input_filename, "Year", index_y, description = "Year"),
        read_param("hour", input_filename, "Hour", index_t, description = "Hour"),
        ParamScalar("year_start", 2018, description = "simulation start year"),
    )
end

# Struct with no fields used to dispatch -- this is the traits pattern
abstract type MarketStructure end
struct VerticallyIntegratedUtility <: MarketStructure end
struct WholesaleMarket <: MarketStructure end

abstract type UseCase end
struct NullUseCase <: UseCase end
struct DERUseCase <: UseCase end
struct SupplyChoiceUseCase <: UseCase end

abstract type Options end

get_file_prefix(::Options) = String("")

struct HEMOptions{T <: MarketStructure, 
                  U <: Union{NullUseCase,DERUseCase}, 
                  V <: Union{NullUseCase,SupplyChoiceUseCase}} <: Options
    MIP_solver::HEMSolver
    NLP_solver::HEMSolver
    market_structure::T

    # use case switches
    der_use_case::U
    supply_choice_use_case::V
end

function get_file_prefix(options::HEMOptions)
    return join(["$(typeof(options.der_use_case))", 
                 "$(typeof(options.supply_choice_use_case))",
                 "$(typeof(options.market_structure))"],"_")
end

"""
Abstract type for agents.

Required interfaces:
- get_id(agent::Agent)::String
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

function iter_agents_and_options(store::AgentStore)
    return ((x.agent, x.options) for agents in values(store.data) for x in values(agents))
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
    max_iter = 100
    window_length = 1
    TimerOutputs.reset_timer!(HEM_TIMER)


    TimerOutputs.@timeit HEM_TIMER "solve_equilibrium_problem!" begin
        for w in 1:(length(model_data.index_y_fix) - window_length + 1)  # loop over windows
            model_data.index_y.elements =
                model_data.index_y_fix.elements[w:(w + window_length - 1)]
            i = 0
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
                        )
                    end
                    push!(diff_vec, diff_one)
                end
                diff = maximum(diff_vec)
                @info "Iteration $i value: $diff"

                if diff < model_data.epsilon
                    break
                end
            end

            # save_welfare(welfare, export_file_path, file_prefix)
            i >= max_iter && error("Reached max iterations $max_iter with no solution")
            @info "Problem solved!"

            update_cumulative!(model_data, agents_and_opts)
        end

    end

    for (agent, options) in iter_agents_and_options(store)
        # push!(welfare, welfare_calculation(agent, options, model_data, hem_opts, other_agents))
        save_results(agent, options, hem_opts, export_file_path, file_prefix)
    end

    if hem_opts.market_structure isa VerticallyIntegratedUtility
        x = store.data[Utility]["default"]
        Welfare_supply =
            welfare_calculation!(x.agent, x.options, model_data, hem_opts, store)
        y = store.data[CustomerGroup]["default"]
        Welfare_demand =
            welfare_calculation!(y.agent, y.options, model_data, hem_opts, store)
    elseif hem_opts.market_structure isa WholesaleMarket
        x = store.data[IPPGroup]["default"]
        Welfare_supply =
            welfare_calculation!(x.agent, x.options, model_data, hem_opts, store)
        y = store.data[CustomerGroup]["default"]
        Welfare_demand =
            welfare_calculation!(y.agent, y.options, model_data, hem_opts, store)
    end

    if hem_opts.supply_choice_use_case isa NullUseCase
        Welfare_green_developer = [AxisArray(
            [
                0.0
                for y in model_data.index_y_fix
            ],
            model_data.index_y_fix.elements,
        ), AxisArray(
            [
                0.0
                for y in model_data.index_y_fix
            ],
            model_data.index_y_fix.elements,
        ), AxisArray(
            [
                0.0
                for y in model_data.index_y_fix
            ],
            model_data.index_y_fix.elements,
        ), AxisArray(
            [
                0.0
                for y in model_data.index_y_fix
            ],
            model_data.index_y_fix.elements,
        ), AxisArray(
            [
                0.0
                for y in model_data.index_y_fix
            ],
            model_data.index_y_fix.elements,
        ), AxisArray(
            [
                0.0
                for y in model_data.index_y_fix
            ],
            model_data.index_y_fix.elements,
        ), AxisArray(
            [
                0.0
                for y in model_data.index_y_fix
            ],
            model_data.index_y_fix.elements,
        )]
    else
        z = store.data[GreenDeveloper]["default"]
        Welfare_green_developer =
            welfare_calculation!(z.agent, z.options, model_data, hem_opts, store)
    end

    save_welfare!(Welfare_supply, Welfare_demand, Welfare_green_developer, export_file_path, file_prefix)

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

# TODO: Write the welfare calculation and saving more generally
function save_welfare!(
    Supply::Any,
    Demand::Any,
    GreenDeveloper::Any,
    exportfilepath::AbstractString,
    fileprefix::AbstractString,
)
    save_param(
        Demand[1].values,
        [:Year, :CustomerType, :DERTech],
        :PVNetCS_dollar,
        joinpath(exportfilepath, "$(fileprefix)_PVNetCS.csv"),
    )
    save_param(
        Demand[2].values,
        [:Year, :CustomerType],
        :ConGreenPowerNetSurplus_dollar,
        joinpath(exportfilepath, "$(fileprefix)_GreenPowerNetCS.csv"),
    )
    # save_param(
    #     Demand[2],
    #     [:Year, :CustomerType, :DERTech],
    #     :PVEnergySaving_dollar,
    #     joinpath(exportfilepath, "$(fileprefix)_PVSaving.csv"),
    # )
    # save_param(
    #     Demand[3],
    #     [:Year, :CustomerType, :DERTech],
    #     :EnergyCost_dollar,
    #     joinpath(exportfilepath, "$(fileprefix)_EnergyCost.csv"),
    # )
    # save_param(
    #     Demand[4],
    #     [:Year, :CustomerType, :DERTech],
    #     :NetCS_dollar,
    #     joinpath(exportfilepath, "$(fileprefix)_NetCS.csv"),
    # )
    save_param(
        Demand[3],
        [:Year],
        :TotalNetCS_dollar,
        joinpath(exportfilepath, "$(fileprefix)_TotalNetCS.csv"),
    )
    # save_param(
    #     Demand[6],
    #     [:Year, :CustomerType, :DERTech],
    #     :NetCS_dollar,
    #     joinpath(exportfilepath, "$(fileprefix)_NetCS_per_customer.csv"),
    # )
    # save_param(
    #     Demand[7],
    #     [:Year, :CustomerType],
    #     :dollar,
    #     joinpath(exportfilepath, "$(fileprefix)_annual_bill_per_customer.csv"),
    # )
    # save_param(
    #     Demand[8],
    #     [:Year, :CustomerType],
    #     :dollar,
    #     joinpath(exportfilepath, "$(fileprefix)_average_bill_per_customer.csv"),
    # )
    save_param(
        Supply[1],
        [:Year],
        :SupplierRevenue_dollar,
        joinpath(exportfilepath, "$(fileprefix)_SupplierRevenue.csv"),
    )
    save_param(
        Supply[2],
        [:Year],
        :SupplierCost_dollar,
        joinpath(exportfilepath, "$(fileprefix)_SupplierCost.csv"),
    )
    save_param(
        Supply[3],
        [:Year],
        :DebtInterest_dollar,
        joinpath(exportfilepath, "$(fileprefix)_DebtInterest.csv"),
    )
    save_param(
        Supply[4],
        [:Year],
        :IncomeTax_dollar,
        joinpath(exportfilepath, "$(fileprefix)_IncomeTax.csv"),
    )
    save_param(
        Supply[5],
        [:Year],
        :OperationalCost_dollar,
        joinpath(exportfilepath, "$(fileprefix)_OperationalCost.csv"),
    )
    save_param(
        Supply[6],
        [:Year],
        :Depreciation_dollar,
        joinpath(exportfilepath, "$(fileprefix)_Depreciation.csv"),
    )
    save_param(
        Supply[7],
        [:Year],
        :Depreciation_dollar,
        joinpath(exportfilepath, "$(fileprefix)_Tax_Depreciation.csv"),
    )
    save_param(
        Supply[8],
        [:Year],
        :Metric_ton,
        joinpath(exportfilepath, "$(fileprefix)_Total_Emission.csv"),
    )
    save_param(
        GreenDeveloper[1],
        [:Year],
        :GreenDeveloperRevenue_dollar,
        joinpath(exportfilepath, "$(fileprefix)_GreenDeveloperRevenue.csv"),
    )
    save_param(
        GreenDeveloper[2],
        [:Year],
        :GreenDeveloperCost_dollar,
        joinpath(exportfilepath, "$(fileprefix)_GreenDeveloperCost.csv"),
    )
    save_param(
        GreenDeveloper[3],
        [:Year],
        :DebtInterest_dollar,
        joinpath(exportfilepath, "$(fileprefix)_GreenDeveloperDebtInterest.csv"),
    )
    save_param(
        GreenDeveloper[4],
        [:Year],
        :IncomeTax_dollar,
        joinpath(exportfilepath, "$(fileprefix)_GreenDeveloperIncomeTax.csv"),
    )
    save_param(
        GreenDeveloper[5],
        [:Year],
        :OperationalCost_dollar,
        joinpath(exportfilepath, "$(fileprefix)_GreenDeveloperOperationalCost.csv"),
    )
    save_param(
        GreenDeveloper[6],
        [:Year],
        :Depreciation_dollar,
        joinpath(exportfilepath, "$(fileprefix)_GreenDeveloperDepreciation.csv"),
    )
    save_param(
        GreenDeveloper[7],
        [:Year],
        :Depreciation_dollar,
        joinpath(exportfilepath, "$(fileprefix)_GreenDeveloper_Tax_Depreciation.csv"),
    )
end
