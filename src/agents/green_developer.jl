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
            value = 10.0,
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
    hem_opts::HEMOptions{<:MarketStructure, <:UseCase, SupplyChoice, <:UseCase},
    agent_store::AgentStore,
    w_iter,
    window_length,
    jump_model,
    export_file_path,
    update_results::Bool,
    output_intermediate_results::Bool
)

    utility = get_agent(Utility, agent_store)
    customers = get_agent(CustomerGroup, agent_store)

    # the year green developer is solving PPA investment problem
    reg_year, reg_year_index = get_reg_year(model_data)
    reg_year_pre, reg_year_index_pre = get_prev_reg_year(model_data, w_iter)

    Green_Developer_model = get_new_jump_model(green_developer_opts.solvers)

    # x_green is the annual PPA buildout (x_green is indexed by h for rate-making purpose)
    @variable(Green_Developer_model, x_green[green_developer.index_j, model_data.index_h] >= 0)

    x_green_cumu = make_keyed_array(green_developer.index_j, model_data.index_h)
    for j in green_developer.index_j, h in model_data.index_h
        if reg_year == model_data.year(first(model_data.index_y_fix))
            x_green_cumu[j, h] = 0.0
        else
            x_green_cumu[j, h] = sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) 
                for y_symbol in model_data.year(first(model_data.index_y_fix)):(reg_year - 1))
        end
    end

    objective_function = begin
        sum(
            # fixed o&m
            utility.fom_C_my(reg_year_index, j) * x_green[j, h] / (green_developer.irr * (1 + green_developer.irr)^20 / ((1 + green_developer.irr)^20 - 1)) +
            # capital costs
            utility.CapEx_my(reg_year_index, j) * x_green[j, h] for j in green_developer.index_j, h in model_data.index_h
        )
    end

    @objective(Green_Developer_model, Min, objective_function)

    @constraint(
        Green_Developer_model,
        Eq_ppa[h in model_data.index_h],
        sum(
            model_data.omega(t) * utility.rho_C_my(j, t) * (x_green[j, h] + x_green_cumu[j, h]) for
            j in green_developer.index_j, t in model_data.index_t
        ) -
        customers.x_green_sub_my(reg_year_index, h) / (1 - utility.loss_dist) >=
        0
    )

    optimize!(Green_Developer_model)

    green_tech_buildout_before = ParamArray(green_developer.green_tech_buildout_my, "green_tech_buildout_before")

    for j in green_developer.index_j, h in model_data.index_h
        green_developer.green_tech_buildout_my(reg_year_index, j, h, :) .= value.(x_green[j, h])
    end

    for h in model_data.index_h
        if sum(green_developer.green_tech_buildout_my(reg_year_index, j, h) for j in green_developer.index_j) > 0.0
            green_developer.ppa_my(reg_year_index, h, :) .= (sum(sum(utility.fom_C_my(reg_year_index, j) * green_developer.green_tech_buildout_my(reg_year_index, j, h) for j in green_developer.index_j) / (1+green_developer.irr)^n for n in 1:20) +
            sum(utility.CapEx_my(reg_year_index, j) * green_developer.green_tech_buildout_my(reg_year_index, j, h) * (1 - utility.ITC_new_my(reg_year_index, j)) for j in green_developer.index_j)) /
            (sum(customers.x_green_sub_incremental_my(reg_year_index, h) / ((1+green_developer.irr)^n) for n in 1:20))
        else
            green_developer.ppa_my(reg_year_index, h, :) .= 0.0
        end
    end

    green_developer.current_year = reg_year_index
    green_developer.previous_year = reg_year_index_pre

    return compute_difference_percentage_one_norm([(green_tech_buildout_before, green_developer.green_tech_buildout_my)])

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

function save_results(
    green_developer::GreenDeveloper,
    green_developer_opts::GreenDeveloperOptions,
    hem_opts::HEMOptions{<:MarketStructure, <:UseCase, SupplyChoice, <:UseCase},
    export_file_path::AbstractString,
)

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
