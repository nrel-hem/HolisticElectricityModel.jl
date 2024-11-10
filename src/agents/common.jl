# This module defines inputs that are held in common across all agents

const DEFAULT_ID = "default"
const HEM_TIMER = TimerOutputs.TimerOutput()

mutable struct HEMData
    # Configuration
    epsilon::ParamScalar # iteration tolerance

    # Sets
    index_y::Dimension # year index
    index_y_fix::Dimension # year index
    index_s::Dimension # year index (for new resources depreciation schedule)
    index_d::Dimension # representative day index
    index_t::Dimension # time index (within each representative day)
    index_h::Dimension # customer types
    index_j::Dimension # green tariff technologies
    index_z::Dimension # zone index
    index_sector::Dimension # zone index

    # Parameters
    omega::ParamArray # number of hours per timeslice
    year::ParamArray
    time::ParamArray
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
    # "index_y_fix" represents the full simulation horizon (does not change)
    # "index_y" represents the simulation years in a particular window (gets updated on line 296)
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

    # new resource depreciation year index
    index_s = read_set(
        input_filename,
        "index_s",
        "index_s",
        prose_name = "new resource depreciation year index s",
        description = "new resource depreciation years",
    )

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

    # green technology types
    index_j = read_set(
        input_filename,
        "index_j",
        "index_j",
        prose_name = "green technologies index j",
        description = "green tariff technologies",
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
        "index_sector",
        prose_name = "customer group index sector",
        description = "customer high level groups",
        
    )

    return HEMData(
        ParamScalar("epsilon", epsilon, description = "iteration tolerance"),
        index_y,
        index_y_fix,
        index_s,
        index_d,
        index_t,
        index_h,
        index_j,
        index_z,
        index_sector,
        read_param(
            "omega",
            input_filename,
            "Omega",
            index_d,
            description = "number of days per representative day",
        ),
        read_param("year", input_filename, "Year", index_y, description = "Year"),
        read_param("time", input_filename, "Time", index_t, description = "Time"),
        ParamScalar("year_start", 2020, description = "simulation start year"),
    )
end

function get_delta_t(model_data::HEMData)
    return (
        parse(Int64, chop(string(model_data.index_t.elements[2]), head = 1, tail = 0)) - 
        parse(Int64, chop(string(model_data.index_t.elements[1]), head = 1, tail = 0))
    )
end

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


# Struct with no fields used to dispatch -- this is the traits pattern
abstract type MarketStructure end
struct VerticallyIntegratedUtility <: MarketStructure end
struct WholesaleMarket <: MarketStructure end

abstract type UseCase end
struct NullUseCase <: UseCase end
struct DERAdoption <: UseCase end
struct SupplyChoice <: UseCase end
struct DERAggregation <: UseCase end

abstract type Options end

get_file_prefix(::Options) = String("")

struct HEMOptions{T <: MarketStructure, 
                  U <: Union{NullUseCase,DERAdoption},
                  V <: Union{NullUseCase,SupplyChoice},
                  W <: Union{NullUseCase,DERAggregation}} <: Options
    market_structure::T

    # use case switches
    der_use_case::U
    supply_choice_use_case::V
    der_aggregation_use_case::W
end

function get_file_prefix(options::HEMOptions)
    return join(["$(typeof(options.der_use_case))", 
                 "$(typeof(options.supply_choice_use_case))",
                 "$(typeof(options.der_aggregation_use_case))",
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
    jump_model::Any,
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
                            jump_model,
                            export_file_path,
                            true
                        )
                    end
                    @assert !isnothing(diff_one) "Nothing returned by solve_agent_problem!($(typeof(agent))): $(diff_one)"
                    @assert !(diff_one isa Tuple) "Tuple returned by solve_agent_problem!($(typeof(agent))): $(diff_one)"
                    @info "$(diff_one)"
                    push!(diff_vec, diff_one)
                end
                diff = maximum(diff_vec)
                @info "Iteration $i value: $diff"
                push!(diff_iter, diff)
                @info "Iteration $i value vector: $diff_iter"

                if diff < model_data.epsilon
                    break
                end
            end

            # save_welfare(welfare, export_file_path)
            i >= max_iter && error("Reached max iterations $max_iter with no solution")
            @info "Problem solved!"

            update_cumulative!(model_data, agents_and_opts)
        end

    end

    for (agent, options) in iter_agents_and_options(store)
        # push!(welfare, welfare_calculation(agent, options, model_data, hem_opts, other_agents))
        save_results(agent, options, hem_opts, export_file_path)
    end

    # if hem_opts.market_structure isa VerticallyIntegratedUtility
    #     x = store.data[Utility]["default"]
    #     Welfare_supply =
    #         welfare_calculation!(x.agent, x.options, model_data, hem_opts, store)
    #     y = store.data[CustomerGroup]["default"]
    #     Welfare_demand =
    #         welfare_calculation!(y.agent, y.options, model_data, hem_opts, store)
    # elseif hem_opts.market_structure isa WholesaleMarket
    #     x = store.data[IPPGroup]["default"]
    #     Welfare_supply =
    #         welfare_calculation!(x.agent, x.options, model_data, hem_opts, store)
    #     y = store.data[CustomerGroup]["default"]
    #     Welfare_demand =
    #         welfare_calculation!(y.agent, y.options, model_data, hem_opts, store)
    # end

    # if hem_opts.supply_choice_use_case isa NullUseCase
    #     Welfare_green_developer = [
    #         initialize_keyed_array(model_data.index_y_fix), 
    #         initialize_keyed_array(model_data.index_y_fix), 
    #         initialize_keyed_array(model_data.index_y_fix), 
    #         initialize_keyed_array(model_data.index_y_fix), 
    #         initialize_keyed_array(model_data.index_y_fix), 
    #         initialize_keyed_array(model_data.index_y_fix), 
    #         initialize_keyed_array(model_data.index_y_fix)]
    # else
    #     z = store.data[GreenDeveloper]["default"]
    #     Welfare_green_developer =
    #         welfare_calculation!(z.agent, z.options, model_data, hem_opts, store)
    # end

    # save_welfare!(Welfare_supply, Welfare_demand, Welfare_green_developer, export_file_path)

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
)
    save_param(
        Demand[1].values,
        [:Year, :CustomerType, :DERTech],
        :PVNetCS_dollar,
        joinpath(exportfilepath, "PVNetCS.csv"),
    )
    save_param(
        Demand[2].values,
        [:Year, :CustomerType],
        :ConGreenPowerNetSurplus_dollar,
        joinpath(exportfilepath, "GreenPowerNetCS.csv"),
    )
    # save_param(
    #     Demand[2],
    #     [:Year, :CustomerType, :DERTech],
    #     :PVEnergySaving_dollar,
    #     joinpath(exportfilepath, "PVSaving.csv"),
    # )
    # save_param(
    #     Demand[3],
    #     [:Year, :CustomerType, :DERTech],
    #     :EnergyCost_dollar,
    #     joinpath(exportfilepath, "EnergyCost.csv"),
    # )
    # save_param(
    #     Demand[4],
    #     [:Year, :CustomerType, :DERTech],
    #     :NetCS_dollar,
    #     joinpath(exportfilepath, "NetCS.csv"),
    # )
    save_param(
        Demand[3],
        [:Year],
        :TotalNetCS_dollar,
        joinpath(exportfilepath, "TotalNetCS.csv"),
    )
    # save_param(
    #     Demand[6],
    #     [:Year, :CustomerType, :DERTech],
    #     :NetCS_dollar,
    #     joinpath(exportfilepath, "NetCS_per_customer.csv"),
    # )
    # save_param(
    #     Demand[7],
    #     [:Year, :CustomerType],
    #     :dollar,
    #     joinpath(exportfilepath, "annual_bill_per_customer.csv"),
    # )
    # save_param(
    #     Demand[8],
    #     [:Year, :CustomerType],
    #     :dollar,
    #     joinpath(exportfilepath, "average_bill_per_customer.csv"),
    # )
    save_param(
        Supply[1],
        [:Year],
        :SupplierRevenue_dollar,
        joinpath(exportfilepath, "SupplierRevenue.csv"),
    )
    save_param(
        Supply[2],
        [:Year],
        :SupplierCost_dollar,
        joinpath(exportfilepath, "SupplierCost.csv"),
    )
    save_param(
        Supply[3],
        [:Year],
        :DebtInterest_dollar,
        joinpath(exportfilepath, "DebtInterest.csv"),
    )
    save_param(
        Supply[4],
        [:Year],
        :IncomeTax_dollar,
        joinpath(exportfilepath, "IncomeTax.csv"),
    )
    save_param(
        Supply[5],
        [:Year],
        :OperationalCost_dollar,
        joinpath(exportfilepath, "OperationalCost.csv"),
    )
    save_param(
        Supply[6],
        [:Year],
        :Depreciation_dollar,
        joinpath(exportfilepath, "Depreciation.csv"),
    )
    save_param(
        Supply[7],
        [:Year],
        :Depreciation_dollar,
        joinpath(exportfilepath, "Tax_Depreciation.csv"),
    )
    save_param(
        Supply[8],
        [:Year],
        :Metric_ton,
        joinpath(exportfilepath, "Total_Emission.csv"),
    )
    save_param(
        GreenDeveloper[1],
        [:Year],
        :GreenDeveloperRevenue_dollar,
        joinpath(exportfilepath, "GreenDeveloperRevenue.csv"),
    )
    save_param(
        GreenDeveloper[2],
        [:Year],
        :GreenDeveloperCost_dollar,
        joinpath(exportfilepath, "GreenDeveloperCost.csv"),
    )
    save_param(
        GreenDeveloper[3],
        [:Year],
        :DebtInterest_dollar,
        joinpath(exportfilepath, "GreenDeveloperDebtInterest.csv"),
    )
    save_param(
        GreenDeveloper[4],
        [:Year],
        :IncomeTax_dollar,
        joinpath(exportfilepath, "GreenDeveloperIncomeTax.csv"),
    )
    save_param(
        GreenDeveloper[5],
        [:Year],
        :OperationalCost_dollar,
        joinpath(exportfilepath, "GreenDeveloperOperationalCost.csv"),
    )
    save_param(
        GreenDeveloper[6],
        [:Year],
        :Depreciation_dollar,
        joinpath(exportfilepath, "GreenDeveloperDepreciation.csv"),
    )
    save_param(
        GreenDeveloper[7],
        [:Year],
        :Depreciation_dollar,
        joinpath(exportfilepath, "GreenDeveloper_Tax_Depreciation.csv"),
    )
end
