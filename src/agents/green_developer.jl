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
    "internal rate of return"
    irr::ParamScalar
    "ppa price"
    ppa_my::ParamArray
    "annual green tech buildout (under PPA)"
    green_tech_buildout_my::ParamArray
end

function GreenDeveloper(input_filename::AbstractString, model_data::HEMData; id = DEFAULT_ID)
    return GreenDeveloper(
        id,
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
            model_data.index_j,
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
    hem_opts::HEMOptions{<:MarketStructure, <:Union{NullUseCase,DERUseCase}, SupplyChoiceUseCase},
    agent_store::AgentStore,
    w_iter,
    jump_model,
    export_file_path,
)

    utility = get_agent(Utility, agent_store)
    customers = get_agent(CustomerGroup, agent_store)

    # the year green developer is solving PPA investment problem
    reg_year, reg_year_index = get_reg_year(model_data)

    Green_Developer_model = get_new_jump_model(green_developer_opts.solvers)

    # x_green is the annual PPA buildout (x_green is indexed by h for rate-making purpose)
    @variable(Green_Developer_model, x_green[model_data.index_j, model_data.index_h] >= 0)

    x_green_cumu = make_keyed_array(model_data.index_j, model_data.index_h)
    for j in model_data.index_j, h in model_data.index_h
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
            utility.CapEx_my(reg_year_index, j) * x_green[j, h] for j in model_data.index_j, h in model_data.index_h
        )
    end

    @objective(Green_Developer_model, Min, objective_function)

    @constraint(
        Green_Developer_model,
        Eq_ppa[h in model_data.index_h],
        sum(
            model_data.omega(t) * utility.rho_C_my(j, t) * (x_green[j, h] + x_green_cumu[j, h]) for
            j in model_data.index_j, t in model_data.index_t
        ) -
        customers.x_green_sub_my(reg_year_index, h) / (1 - utility.loss_dist) >=
        0
    )

    optimize!(Green_Developer_model)

    green_tech_buildout_before = ParamArray(green_developer.green_tech_buildout_my, "green_tech_buildout_before")

    for j in model_data.index_j, h in model_data.index_h
        green_developer.green_tech_buildout_my(reg_year_index, j, h, :) .= value.(x_green[j, h])
    end

    for h in model_data.index_h
        if sum(green_developer.green_tech_buildout_my(reg_year_index, j, h) for j in model_data.index_j) > 0.0
            green_developer.ppa_my(reg_year_index, h, :) .= (sum(sum(utility.fom_C_my(reg_year_index, j) * green_developer.green_tech_buildout_my(reg_year_index, j, h) for j in model_data.index_j) / (1+green_developer.irr)^n for n in 1:20) +
            sum(utility.CapEx_my(reg_year_index, j) * green_developer.green_tech_buildout_my(reg_year_index, j, h) * (1 - utility.ITC_new_my(reg_year_index, j)) for j in model_data.index_j)) /
            (sum(customers.x_green_sub_incremental_my(reg_year_index, h) / ((1+green_developer.irr)^n) for n in 1:20))
        else
            green_developer.ppa_my(reg_year_index, h, :) .= 0.0
        end
    end

    return compute_difference_percentage_one_norm([(green_tech_buildout_before, green_developer.green_tech_buildout_my)])

end

function solve_agent_problem!(
    green_developer::GreenDeveloper,
    green_developer_opts::GreenDeveloperOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{<:MarketStructure, DERUseCase, NullUseCase},
    agent_store::AgentStore,
    w_iter,
    jump_model,
    export_file_path
)

    return 0.0, nothing

end

function save_results(
    green_developer::GreenDeveloper,
    green_developer_opts::AgentOptions,
    hem_opts::HEMOptions{<:MarketStructure, <:Union{NullUseCase,DERUseCase}, SupplyChoiceUseCase},
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

# TODO: welfare for green developer
function welfare_calculation!(
    green_developer::GreenDeveloper,
    green_developer_opts::AgentOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{<:MarketStructure, <:Union{NullUseCase,DERUseCase}, SupplyChoiceUseCase},
    agent_store::AgentStore,
)
    utility = get_agent(Utility, agent_store)
    customers = get_agent(CustomerGroup, agent_store)

    green_developer_Revenue = KeyedArray(
        [
            sum(green_developer.ppa_my(Symbol(Int(y_symbol)), h) * customers.x_green_sub_incremental_my(Symbol(Int(y_symbol)), h)
            for y_symbol in model_data.year(first(model_data.index_y_fix)):model_data.year(y), h in model_data.index_h)
            for y in model_data.index_y_fix
        ];
        [get_pair(model_data.index_y_fix)]...,
    )

    energy_cost = KeyedArray(
        [ 
            0.0
            for y in model_data.index_y_fix
        ];
        [get_pair(model_data.index_y_fix)]...,
    )
    fixed_om = KeyedArray(
        [
            sum(utility.fom_C_my(Symbol(Int(y_symbol)), j) * green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h)
            for y_symbol in model_data.year(first(model_data.index_y_fix)):model_data.year(y), h in model_data.index_h, j in model_data.index_j)
            for y in model_data.index_y_fix
        ];
        [get_pair(model_data.index_y_fix)]...
    )
    operational_cost = KeyedArray(
        [ 
            energy_cost(y) + fixed_om(y) for y in model_data.index_y_fix
        ];
        [get_pair(model_data.index_y_fix)]...
    )
    working_capital = KeyedArray(
        [ 
            utility.DaysofWC / 365 * operational_cost(y) 
            for y in model_data.index_y_fix
        ];
        [get_pair(model_data.index_y_fix)]...
    )

    # assume the green developer's new depreciation schedule is the same as utility's
    ADITNew = make_keyed_array(model_data.index_y_fix, model_data.index_j)
    for y in model_data.index_y_fix, j in model_data.index_j
        ADITNew[y, j] = sum(
            utility.CapEx_my(Symbol(Int(y_symbol)), j) *
            sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for h in model_data.index_h) *
            (
                utility.CumuTaxDepre_new_my(
                    Symbol(Int(model_data.year(y) - y_symbol + 1)),
                    j,
                ) - utility.CumuAccoutDepre_new_my(
                    Symbol(Int(model_data.year(y) - y_symbol + 1)),
                    j,
                )
            ) *
            utility.Tax +
            utility.ITC_new_my(Symbol(Int(y_symbol)), j) *
            utility.CapEx_my(Symbol(Int(y_symbol)), j) *
            sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for h in model_data.index_h) *
            (
                1 - utility.CumuITCAmort_new_my(
                    Symbol(Int(model_data.year(y) - y_symbol + 1)),
                    j,
                )
            ) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):model_data.year(y)
        )
    end
    RateBaseNoWC_new = make_keyed_array(model_data.index_y_fix, model_data.index_j)
    for y in model_data.index_y_fix, j in model_data.index_j
        RateBaseNoWC_new[y, j] =
            sum(
                utility.CapEx_my(Symbol(Int(y_symbol)), j) *
                sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for h in model_data.index_h) *
                (
                    1 - utility.CumuAccoutDepre_new_my(
                        Symbol(Int(model_data.year(y) - y_symbol + 1)),
                        j,
                    )
                ) for y_symbol in
                    model_data.year(first(model_data.index_y_fix)):model_data.year(y)
            ) - ADITNew[y, j]
    end

    rate_base = KeyedArray(
        [
            sum(RateBaseNoWC_new[y, j] for j in model_data.index_j) +
            working_capital(y) for y in model_data.index_y_fix
        ];
        [get_pair(model_data.index_y_fix)]...,
    )
    debt_interest = KeyedArray(
        [ 
            rate_base(y) * utility.DebtRatio * utility.COD 
            for y in model_data.index_y_fix
        ];
        [get_pair(model_data.index_y_fix)]...
    )

    depreciation = KeyedArray(
        [
            sum(
                utility.CapEx_my(Symbol(Int(y_symbol)), j) *
                sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for h in model_data.index_h) *
                utility.AnnualAccoutDepre_new_my(
                    Symbol(Int(model_data.year(y) - y_symbol + 1)),
                    j,
                 ) for y_symbol in
                    model_data.year(first(model_data.index_y_fix)):model_data.year(y),
                    j in model_data.index_j
            ) for y in model_data.index_y_fix
        ];
        [get_pair(model_data.index_y_fix)]...
    )

    depreciation_tax = KeyedArray(
        [
            sum(
                utility.CapEx_my(Symbol(Int(y_symbol)), j) *
                sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for h in model_data.index_h) *
                utility.AnnualTaxDepre_new_my(
                    Symbol(Int(model_data.year(y) - y_symbol + 1)),
                    j,
                ) for y_symbol in
                    model_data.year(first(model_data.index_y_fix)):model_data.year(y),
                    j in model_data.index_j
            ) for y in model_data.index_y_fix
        ];
        [get_pair(model_data.index_y_fix)]...
    )

    income_tax = KeyedArray(
        [
            (
                green_developer_Revenue(y) - debt_interest(y) - operational_cost(y) -
                depreciation_tax(y)
            ) * utility.Tax - 
            sum(
                utility.ITC_new_my(Symbol(Int(y_symbol)), j) *
                utility.CapEx_my(Symbol(Int(y_symbol)), j) *
                sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for h in model_data.index_h) *
                utility.AnnualITCAmort_new_my(Symbol(Int(model_data.year(y) - y_symbol + 1)), j) for
                y_symbol in model_data.year(first(model_data.index_y_fix)):model_data.year(y), j in model_data.index_j
            )
            for y in model_data.index_y_fix
        ];
        [get_pair(model_data.index_y_fix)]...
    )

    green_developer_Cost = KeyedArray(
        [
            debt_interest(y) +
            income_tax(y) +
            operational_cost(y) +
            depreciation(y) for y in model_data.index_y_fix
        ];
        [get_pair(model_data.index_y_fix)]...
    )

    return green_developer_Revenue,
    green_developer_Cost,
    debt_interest,
    income_tax,
    operational_cost,
    depreciation,
    depreciation_tax
end
