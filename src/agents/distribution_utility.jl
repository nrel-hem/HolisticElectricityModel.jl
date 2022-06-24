using TableTransforms

mutable struct DistributionCostModel
    constant::ParamScalar
    saidi_coefficient::ParamScalar
    dpv_coefficient::ParamScalar
    total_sales_coefficient::ParamScalar
    residential_customer_coefficient::ParamScalar
    commercial_customer_coefficient::ParamScalar
    industrial_customer_coefficient::ParamScalar
end

abstract type AbstractDistributionUtility <: AgentGroup end

mutable struct DistributionUtility <: AbstractDistributionUtility
    id::String
    SAIDI::ParamAxisArray
    distribution_cost_model::DistributionCostModel
end

function DistributionUtility(
    input_filename::AbstractString,
    model_data::HEMData;
    id = DEFAULT_ID,
)

    distribution_cost_model = DistributionCostModel(
        ParamScalar("constant", -0.01257239),
        ParamScalar("saidi_coefficient", 0.04208811),
        ParamScalar("dpv_coefficient", 0.11860937),
        ParamScalar("total_sales_coefficient", -0.03801144),
        ParamScalar("residential_customer_coefficient", 0.54768282),
        ParamScalar("commercial_customer_coefficient", 0.0954236),
        ParamScalar("industrial_customer_coefficient", -0.05217491),
    )

    return DistributionUtility(
        id,
        read_param(
            "SAIDI",
            input_filename,
            "SAIDI",
            model_data.index_y,
            description = "SAIDI",
        ),
        distribution_cost_model,
    )

end

get_id(x::DistributionUtility) = x.id

function solve_agent_problem!(
    distribution_utility::DistributionUtility,
    distribution_utility_opts::AgentOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{<:MarketStructure,<:UseCase},
    agent_store::AgentStore,
    w_iter,
)

    utility = get_agent(Utility, agent_store)
    regulator = get_agent(Regulator, agent_store)
    customers = get_agent(CustomerGroup, agent_store)

    reg_year = model_data.year[first(model_data.index_y)]
    reg_year_index = Symbol(Int(reg_year))

    distribution_cost_model = distribution_utility.distribution_cost_model

    total_sale =
        sum(
            customers.gamma[h] * model_data.omega[t] * customers.d[h, t] for
            h in model_data.index_h, t in model_data.index_t
        ) +
        # export
        sum(
            model_data.omega[t] * utility.eximport_my[reg_year_index, t] for
            t in model_data.index_t
        ) -
        # DG
        sum(
            model_data.omega[t] * (
                customers.rho_DG[h, m, t] *
                customers.x_DG_E_my[first(model_data.index_y), h, m] + sum(
                    customers.rho_DG[h, m, t] * customers.x_DG_new_my[Symbol(Int(y)), h, m]
                    for y = model_data.year[first(model_data.index_y_fix)]:reg_year
                )
            ) for t in model_data.index_t, h in model_data.index_h, m in customers.index_m
        )

    distribution_cost_before = regulator.distribution_cost

    # Use static normalization numbers from the national database
    norm_inputs = DataFrame(
        id = [
            "saidi",
            "dpv",
            "total_sales",
            "residential",
            "commercial",
            "industrial",
            "distribution_cost",
        ],
        min = [0, 0.02, 33321, 3664, 820, 0, 1851740],
        max = [968.6, 3018.362, 131738016, 4489462, 619752, 33562, 3236569353],
    )

    function norm(data, norm_variable)
        data =
            (data - norm_inputs[norm_inputs.id.==norm_variable, "min"][1]) / (
                norm_inputs[norm_inputs.id.==norm_variable, "max"][1] -
                norm_inputs[norm_inputs.id.==norm_variable, "min"][1]
            )
        return data
    end

    function reverse_norm(data, norm_variable)
        data =
            data * (
                norm_inputs[norm_inputs.id.==norm_variable, "max"][1] -
                norm_inputs[norm_inputs.id.==norm_variable, "min"][1]
            ) + norm_inputs[norm_inputs.id.==norm_variable, "min"][1]
        return data
    end
    # customer module updates customers.x_DG_E to be the DPV at the beginning of the year,
    # to use end of the year, add customers.x_DG_new_my[reg_year_index, h, m]    
    dpv_pca =
        sum(customers.x_DG_E[h, m] for h in model_data.index_h, m in customers.index_m)

    distribution_cost_norm =
        distribution_cost_model.constant +
        distribution_cost_model.saidi_coefficient *
        norm(distribution_utility.SAIDI[reg_year_index], "saidi") +
        distribution_cost_model.dpv_coefficient * norm(dpv_pca, "dpv") +
        distribution_cost_model.total_sales_coefficient * norm(total_sale, "total_sales") +
        distribution_cost_model.residential_customer_coefficient *
        norm(customers.gamma[:Residential], "residential") +
        distribution_cost_model.commercial_customer_coefficient *
        norm(customers.gamma[:Commercial], "commercial") +
        distribution_cost_model.industrial_customer_coefficient *
        norm(customers.gamma[:Industrial], "industrial")

    regulator.distribution_cost[reg_year_index] =
        reverse_norm(distribution_cost_norm, "distribution_cost")

    return compute_difference_percentage_one_norm([
        (distribution_cost_before, regulator.distribution_cost),
    ])

end
