abstract type AbstractGreenDeveloper <: Agent end

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
            model_data.index_h,
            value = 0.0,
            description = "annual green tech buildout (under PPA)",
        ),
    )
end

get_id(x::GreenDeveloper) = x.id

function solve_agent_problem!(
    green_developer::GreenDeveloper,
    green_developer_opts::AgentOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{<:MarketStructure, <:Union{DERSupplyChoiceUseCase, SupplyChoiceUseCase}},
    agent_store::AgentStore,
    w_iter,
)

    utility = get_agent(Utility, agent_store)
    customers = get_agent(CustomerGroup, agent_store)

    # the year green developer is solving PPA investment problem
    reg_year = model_data.year[first(model_data.index_y)]
    reg_year_index = Symbol(Int(reg_year))

    Green_Developer_model = get_new_jump_model(hem_opts.MIP_solver)

    # x_green is the annual PPA buildout (x_green is indexed by h for rate-making purpose)
    @variable(Green_Developer_model, x_green[model_data.index_j, model_data.index_h] >= 0)

    if reg_year == model_data.year[first(model_data.index_y_fix)]
        x_green_cumu = Dict((j, h) => 0.0 for j in model_data.index_j, h in model_data.index_h)
    else
        x_green_cumu = Dict((j, h) => sum(green_developer.green_tech_buildout_my[y, j, h] 
            for y in model_data.year[first(model_data.index_y_fix)]:(reg_year - 1))
            for j in model_data.index_j, h in model_data.index_h)
    end

    objective_function = begin
        sum(
            # fixed o&m
            utility.fom_C_my[reg_year_index, j] * x_green[j, h] / (green_developer.irr * (1 + green_developer.irr)^20 / ((1 + green_developer.irr)^20 - 1)) +
            # capital costs
            utility.CapEx_my[reg_year_index, j] * x_green[j, h] for j in model_data.index_j, h in model_data.index_h
        )
    end

    @objective(Green_Developer_model, Min, objective_function)

    @constraint(
        Green_Developer_model,
        Eq_ppa[h in model_data.index_h],
        sum(
            model_data.omega[t] * utility.rho_C_my[j, t] * (x_green[j, h] + x_green_cumu[j, h]) for
            j in model_data.index_j, t in model_data.index_t
        ) -
        customers.x_green_sub_my[reg_year_index, h] >=
        0
    )

    optimize!(Green_Developer_model)

    green_tech_buildout_before = ParamArray(green_developer.green_tech_buildout_my, "green_tech_buildout_before")

    for j in model_data.index_j, h in model_data.index_h
        green_developer.green_tech_buildout_my[reg_year_index, j, h] = value.(x_green[j, h])
    end

    for h in model_data.index_h
        if sum(green_developer.green_tech_buildout_my[reg_year_index, j, h] for j in model_data.index_j) > 0.0
            green_developer.ppa_my[reg_year_index, h] = (sum(sum(utility.fom_C_my[reg_year_index, j] * green_developer.green_tech_buildout_my[reg_year_index, j, h] for j in model_data.index_j) / (1+green_developer.irr)^n for n in 1:20) +
            sum(utility.CapEx_my[reg_year_index, j] * green_developer.green_tech_buildout_my[reg_year_index, j, h] for j in model_data.index_j)) /
            (sum(sum(model_data.omega[t] * utility.rho_C_my[j, t] * green_developer.green_tech_buildout_my[reg_year_index, j, h] for j in model_data.index_j, t in model_data.index_t) / ((1+green_developer.irr)^n) for n in 1:20))
        else
            green_developer.ppa_my[reg_year_index, h] = 9999.0
        end
    end

    return compute_difference_one_norm([(green_tech_buildout_before, green_developer.green_tech_buildout_my)])

end

function solve_agent_problem!(
    green_developer::GreenDeveloper,
    green_developer_opts::AgentOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{<:MarketStructure, DERUseCase},
    agent_store::AgentStore,
    w_iter,
)

    return 0.0

end

function save_results(
    green_developer::GreenDeveloper,
    green_developer_opts::AgentOptions,
    hem_opts::HEMOptions{<:MarketStructure, <:Union{DERSupplyChoiceUseCase, SupplyChoiceUseCase}},
    export_file_path::AbstractString,
    fileprefix::AbstractString,
)

    # Primal Variables
    save_param(
        green_developer.green_tech_buildout_my.values,
        [:Year, :GreenTech, :CustomerType],
        :Capacity_MW,
        joinpath(export_file_path, "$(fileprefix)_green_tech_buildout.csv"),
    )
    
    save_param(
        green_developer.ppa_my.values,
        [:Year, :CustomerType],
        :PPAPrice,
        joinpath(export_file_path, "$(fileprefix)_ppa.csv"),
    )
end

# TODO: welfare for green developer

