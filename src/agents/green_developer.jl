# This file defines the data and functions associated with the green developer

abstract type AbstractGreenDeveloper <: Agent end

struct GreenDeveloperOptions <: AgentOptions
    solvers::HEMSolver
    # solvers::Union{HEMSolver, Dict{String, <:HEMSolver}}
end

"""
Construct GreenDeveloperOptions with an MOI.OptimizerWithAttributes instance.
"""
function GreenDeveloperOptions(attributes::MOI.OptimizerWithAttributes)
    return GreenDeveloperOptions(AnySolver(attributes))
end

mutable struct GreenDeveloper <: AbstractGreenDeveloper
    id::String
    current_year::Symbol
    previous_year::Symbol
    index_j::Dimension # green tariff technologies

    "internal rate of return"
    irr::ParamScalar
    "ppa price"
    ppa_my::ParamArray
    "annual green tech buildout (under PPA)"
    green_tech_buildout_my::ParamArray
end

function GreenDeveloper(input_dir::AbstractString, model_data::HEMData; id = DEFAULT_ID)

    index_j = read_set(
        input_dir,
        "index_j",
        "index_j",
        prose_name = "green technologies index j",
        description = "green tariff technologies",
    )
    return GreenDeveloper(
        id,
        first(model_data.index_y),
        first(model_data.index_y),
        index_j,
        ParamScalar("irr", 0.12, description = "internal rate of return"),
        initialize_param(
            "ppa_my",
            model_data.index_y,
            model_data.index_h,
            value = 30.0,
            description = "multi-year ppa price",
        ),
        initialize_param(
            "green_tech_buildout_my",
            model_data.index_y,
            index_j,
            model_data.index_z,
            model_data.index_h,
            value = 0.0,
            description = "annual green tech buildout (under PPA)",
        ),
    )
end

get_id(x::GreenDeveloper) = x.id

function solve_agent_problem!(
    green_developer::GreenDeveloper,
    green_developer_opts::GreenDeveloperOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{WM, <:UseCase, SupplyChoice, <:UseCase},
    agent_store::AgentStore,
    w_iter,
    window_length,
    jump_model,
    export_file_path,
    update_results::Bool,
    output_intermediate_results::Bool
)   
    diff = 0.0
    ipp = get_agent(IPPGroup, agent_store)
    for p in ipp.index_p
        diff += solve_green_developer_problem(
            green_developer,
            green_developer_opts,
            model_data,
            hem_opts,
            agent_store,
            w_iter,
            p,
        )
    end

    return diff
end

function solve_agent_problem!(
    green_developer::GreenDeveloper,
    green_developer_opts::GreenDeveloperOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{VIU, <:UseCase, SupplyChoice, <:UseCase},
    agent_store::AgentStore,
    w_iter,
    window_length,
    jump_model,
    export_file_path,
    update_results::Bool,
    output_intermediate_results::Bool
)   
    diff = solve_green_developer_problem(
        green_developer,
        green_developer_opts,
        model_data,
        hem_opts,
        agent_store,
        w_iter,
        nothing,
    )

    return diff
end

function solve_agent_problem!(
    green_developer::GreenDeveloper,
    green_developer_opts::GreenDeveloperOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{<:MarketStructure, DERAdoption, NullUseCase, <:UseCase},
    agent_store::AgentStore,
    w_iter,
    window_length,
    jump_model,
    export_file_path,
    update_results::Bool,
    output_intermediate_results::Bool
)

    reg_year, reg_year_index = get_reg_year(model_data)
    reg_year_pre, reg_year_index_pre = get_prev_reg_year(model_data, w_iter)

    green_developer.current_year = reg_year_index
    green_developer.previous_year = reg_year_index_pre
    
    return 0.0
end

function get_objective_function(
    x_green::JuMP.Containers.DenseAxisArray,
    green_developer::GreenDeveloper,
    ipp::IPPGroup,
    model_data::HEMData,
    p::Symbol,
)

    objective_function = begin
        sum(
            ipp.fom_C_my(green_developer.current_year, p, z, j) * x_green[j, z, h] / (green_developer.irr * (1 + green_developer.irr)^20 / ((1 + green_developer.irr)^20 - 1)) +
            ipp.CapEx_my(green_developer.current_year, p, z, j) * x_green[j, z, h] for j in green_developer.index_j, (z, h) in model_data.index_z_h_map
        )
    end

    return objective_function
end


function get_objective_function(
    x_green::JuMP.Containers.DenseAxisArray,
    green_developer::GreenDeveloper,
    utility::Utility,
    model_data::HEMData,
    ::Nothing,
)

    objective_function = begin
        sum(
            utility.fom_C_my(green_developer.current_year, z, j) * x_green[j, z, h] / (green_developer.irr * (1 + green_developer.irr)^20 / ((1 + green_developer.irr)^20 - 1)) +
            utility.CapEx_my(green_developer.current_year, z, j) * x_green[j, z, h] for j in green_developer.index_j, (z, h) in model_data.index_z_h_map
        )
    end

    return objective_function
end

function get_green_tech_subscription(
    x_green::JuMP.Containers.DenseAxisArray,
    x_green_cumu::KeyedArray,
    z::Symbol,
    h::Symbol,
    green_developer::GreenDeveloper,
    ipp::IPPGroup,
    model_data::HEMData,
    p::Symbol,
)
    delta_t = model_data.delta_t.value

    green_tech_subscription = begin
        sum(
            model_data.omega(d) * delta_t * ipp.rho_C_my(p, j, z, d, t) * (x_green[j, z, h] + x_green_cumu(j, z, h)) for
            j in green_developer.index_j, t in model_data.index_t, d in model_data.index_d
        )
    end

    return green_tech_subscription

end

function get_green_tech_subscription(
    x_green::JuMP.Containers.DenseAxisArray,
    x_green_cumu::KeyedArray,
    z::Symbol,
    h::Symbol,
    green_developer::GreenDeveloper,
    utility::Utility,
    model_data::HEMData,
    ::Nothing,
)
    delta_t = model_data.delta_t.value

    green_tech_subscription = begin
        sum(
            model_data.omega(d) * delta_t * utility.rho_C_my(j, z, d, t) * (x_green[j, z, h] + x_green_cumu(j, z, h)) for
            j in green_developer.index_j, t in model_data.index_t, d in model_data.index_d
        )
    end

    return green_tech_subscription

end

function get_ppa_price(
    z::Symbol,
    h::Symbol,
    green_developer::GreenDeveloper,
    customers::CustomerGroup,
    ipp::IPPGroup,
    model_data::HEMData,
    p::Symbol,
)

    reg_year, reg_year_index = get_reg_year(model_data)
    ppa_price = begin
        (sum(sum(ipp.fom_C_my(reg_year_index, p, z, j) * green_developer.green_tech_buildout_my(reg_year_index, j, z, h) for j in green_developer.index_j) / (1+green_developer.irr)^n for n in 1:20) +
            sum(ipp.CapEx_my(reg_year_index, p, z, j) * green_developer.green_tech_buildout_my(reg_year_index, j, z, h) * (1 - ipp.ITC_new_my(reg_year_index, j)) for j in green_developer.index_j)) /
            (sum(customers.x_green_sub_incremental_my(reg_year_index, h, z) / ((1+green_developer.irr)^n) for n in 1:20))
    end

    return ppa_price

end

function get_ppa_price(
    z::Symbol,
    h::Symbol,
    green_developer::GreenDeveloper,
    customers::CustomerGroup,
    utility::Utility,
    model_data::HEMData,
    ::Nothing
)

    reg_year, reg_year_index = get_reg_year(model_data)
    ppa_price = begin
        (sum(sum(utility.fom_C_my(reg_year_index, z, j) * green_developer.green_tech_buildout_my(reg_year_index, j, z, h) for j in green_developer.index_j) / (1+green_developer.irr)^n for n in 1:20) +
            sum(utility.CapEx_my(reg_year_index, z, j) * green_developer.green_tech_buildout_my(reg_year_index, j, z, h) * (1 - utility.ITC_new_my(reg_year_index, j)) for j in green_developer.index_j)) /
            (sum(customers.x_green_sub_incremental_my(reg_year_index, h, z) / ((1+green_developer.irr)^n) for n in 1:20))
    end

    return ppa_price
end

"""
    Create the green developer optimization problem depending on the market structure
    and solve for the optimal green tech buildout. Calculate the ppa price from the optimization results.
"""
function solve_green_developer_problem(
    green_developer::GreenDeveloper,
    green_developer_opts::GreenDeveloperOptions,
    model_data::HEMData,
    hem_opts::HEMOptions,
    agent_store::AgentStore,
    w_iter::Int,
    p::Union{Symbol, Nothing} = nothing
)
    utility_or_ipp = get_bulk_system_agent(agent_store, hem_opts)
    customers = get_agent(CustomerGroup, agent_store)

    # the year green developer is solving PPA investment problem
    reg_year, reg_year_index = get_reg_year(model_data)
    reg_year_pre, reg_year_index_pre = get_prev_reg_year(model_data, w_iter)

    Green_Developer_model = get_new_jump_model(green_developer_opts.solvers)

    # x_green is the annual PPA buildout (x_green is indexed by j, z, h)
    @variable(Green_Developer_model, x_green[green_developer.index_j, model_data.index_z, model_data.index_h] >= 0)

    x_green_cumu = make_keyed_array(green_developer.index_j, model_data.index_z, model_data.index_h)
    for j in green_developer.index_j, (z, h) in model_data.index_z_h_map
        if reg_year == model_data.year(first(model_data.index_y_fix))
            x_green_cumu(j, z, h, :) .= 0.0
        else
            x_green_cumu(j, z, h, :) .= sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, z, h) 
                for y_symbol in model_data.year(first(model_data.index_y_fix)):(reg_year - 1))
        end
    end

    objective_function = get_objective_function(x_green, green_developer, utility_or_ipp, model_data, p)

    @objective(Green_Developer_model, Min, objective_function)

    @constraint(
        Green_Developer_model,
        Eq_ppa[h in model_data.index_h, z in model_data.index_z; (z, h) in model_data.index_z_h_map],
        get_green_tech_subscription(x_green, x_green_cumu, z, h, green_developer, utility_or_ipp, model_data, p) -
        customers.x_green_sub_my(reg_year_index, h, z) / (1 - utility_or_ipp.loss_dist) >=
        0
    )

    optimize!(Green_Developer_model)

    green_tech_buildout_before = ParamArray(green_developer.green_tech_buildout_my, "green_tech_buildout_before")

    for j in green_developer.index_j, (z, h) in model_data.index_z_h_map
        green_developer.green_tech_buildout_my(reg_year_index, j, z, h, :) .= value.(x_green[j, z, h])
    end

    for (z, h) in model_data.index_z_h_map
        # There is only one combination of (z, h) for each customer type h
        # so we can use the same ppa price for all j
        if sum(green_developer.green_tech_buildout_my(reg_year_index, j, z, h) for j in green_developer.index_j) > 0.0
            green_developer.ppa_my(reg_year_index, h, :) .= get_ppa_price(z, h, green_developer, customers, utility_or_ipp, model_data, p)
        else
            green_developer.ppa_my(reg_year_index, h, :) .= 0.0
        end
    end

    green_developer.current_year = reg_year_index
    green_developer.previous_year = reg_year_index_pre

    return compute_difference_percentage_one_norm([(green_tech_buildout_before, green_developer.green_tech_buildout_my)])
end

function save_results(
    green_developer::GreenDeveloper,
    green_developer_opts::GreenDeveloperOptions,
    hem_opts::HEMOptions{<:MarketStructure, U, SupplyChoice},
    export_file_path::AbstractString,
) where U <: Union{DERAdoption, NullUseCase}

    # Primal Variables
    save_param(
        green_developer.green_tech_buildout_my.values,
        [:Year, :GreenTech, :CustomerType],
        :Capacity_MW,
        joinpath(export_file_path, "green_tech_buildout.csv"),
    )
    
    save_param(
        green_developer.ppa_my.values,
        [:Year, :CustomerType],
        :PPAPrice,
        joinpath(export_file_path, "ppa.csv"),
    )
end
