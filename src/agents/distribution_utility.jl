using TableTransforms

mutable struct DistributionCapexBalanceModel
    constant::ParamScalar
    total_sales_coefficient::ParamScalar
    residential_customer_coefficient::ParamScalar
    commercial_customer_coefficient::ParamScalar
    industrial_customer_coefficient::ParamScalar
end

mutable struct DistributionCapexAdditionModel
    constant::ParamScalar
    saidi_coefficient::ParamScalar
    dpv_coefficient::ParamScalar
    total_sales_coefficient::ParamScalar
    residential_customer_coefficient::ParamScalar
    commercial_customer_coefficient::ParamScalar
    industrial_customer_coefficient::ParamScalar
end

mutable struct DistributionOMCostModel
    constant::ParamScalar
    saidi_coefficient::ParamScalar
    total_sales_coefficient::ParamScalar
    residential_customer_coefficient::ParamScalar
    commercial_customer_coefficient::ParamScalar
    industrial_customer_coefficient::ParamScalar
end

abstract type AbstractDistributionUtility <: AgentGroup end

mutable struct DistributionUtility <: AbstractDistributionUtility
    id::String
    SAIDI::ParamAxisArray
    distribution_capex_balance_model::DistributionCapexBalanceModel
    distribution_capex_addition_model::DistributionCapexAdditionModel
    distribution_om_cost_model::DistributionOMCostModel
    norm_inputs::DataFrame
    beginning_balance_lifetime::Float64
    DistCumuTaxDepre_new_my::ParamAxisArray # cumulative tax depreciation of new capacity (%)
    DistCumuAccoutDepre_new_my::ParamAxisArray # cumulative accounting depreciation of new capacity (%)
    DistITC_new_my::ParamAxisArray # ITC of new capacity (%)
    DistCumuITCAmort_new_my::ParamAxisArray # ITC amortization of new capacity (%)
    DistAnnualAccoutDepre_new_my::ParamAxisArray # annual accounting depreciation of new capacity (%)
    DistAnnualTaxDepre_new_my::ParamAxisArray # annual tax depreciation of new capacity (%) 
    DistCapExAddition_new_my::ParamAxisArray
    DistOMCost_new_my::ParamAxisArray
    Tax::ParamScalar
    DaysofWC::ParamScalar
    DebtRatio::ParamScalar
    COD::ParamScalar
    COE::ParamScalar
end

function DistributionUtility(
    input_filename::AbstractString,
    model_data::HEMData;
    id = DEFAULT_ID,
)

    distribution_capex_balance_model = DistributionCapexBalanceModel(
        ParamScalar("constant", -0.00712743),
        ParamScalar("total_sales_coefficient", -0.09971823),
        ParamScalar("residential_customer_coefficient", 0.50846654),
        ParamScalar("commercial_customer_coefficient", 0.40584111),
        ParamScalar("industrial_customer_coefficient", 0.11309393),
    )

    distribution_capex_addition_model = DistributionCapexAdditionModel(
        ParamScalar("constant", -0.0193135),
        ParamScalar("saidi_coefficient", 0.03607853),
        ParamScalar("dpv_coefficient", 0.28150498),
        ParamScalar("total_sales_coefficient", -0.0284804),
        ParamScalar("residential_customer_coefficient", 0.53321818),
        ParamScalar("commercial_customer_coefficient", 0.07271089),
        ParamScalar("industrial_customer_coefficient", -0.02300505),
    )

    distribution_om_cost_model = DistributionOMCostModel(
        ParamScalar("constant", -0.00549679),
        ParamScalar("saidi_coefficient", 0.02759952),
        ParamScalar("total_sales_coefficient", -0.0783589),
        ParamScalar("residential_customer_coefficient", 0.4508443),
        ParamScalar("commercial_customer_coefficient", 0.04459026),
        ParamScalar("industrial_customer_coefficient", 0.16479862),
    )

    # Use static normalization numbers from the national database
    # Ella to update min and max for distribution_capex_balance, distribution_capex_addition, distribution_om_cost
    norm_inputs = DataFrame(
        id = [
            "saidi",
            "dpv",
            "total_sales",
            "residential",
            "commercial",
            "industrial",
            "distribution_capex_balance",
            "distribution_capex_addition",
            "distribution_om_cost",
        ],
        min = [0, 0.02, 33321, 3664, 820, 0, 1851740, 1851740, 1851740],
        max = [968.6, 3018.362, 131738016, 4489462, 619752, 33562, 3236569353, 3236569353, 3236569353],
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
        distribution_capex_balance_model,
        distribution_capex_addition_model,
        distribution_om_cost_model,
        norm_inputs,
        10,
        read_param(
            "DistCumuTaxDepre_new_my",
            input_filename,
            "DistCumuTaxDepreNewmy",
            model_data.index_s,
        ),
        read_param(
            "DistCumuAccoutDepre_new_my",
            input_filename,
            "DistCumuAccoutDepreNewmy",
            model_data.index_s,
        ),
        read_param(
            "DistITC_new_my",
            input_filename,
            "DistITCNewmy",
            model_data.index_y,
        ),
        read_param(
            "DistCumuITCAmort_new_my",
            input_filename,
            "DistCumuITCAmortNewmy",
            model_data.index_s,
        ),
        read_param(
            "DistAnnualAccoutDepre_new_my",
            input_filename,
            "DistAnnualAccoutDepreNewmy",
            model_data.index_s,
        ),
        read_param(
            "DistAnnualTaxDepre_new_my",
            input_filename,
            "DistAnnualTaxDepreNewmy",
            model_data.index_s,
        ),
        initialize_param(
            "DistCapExAddition_new_my",
            model_data.index_y,
            description = "distribution capex addition per year",
        ),
        initialize_param(
            "DistOMCost_new_my",
            model_data.index_y,
            description = "distribution O&M cost per year",
        ),
        ParamScalar("Tax", 0.26, description = "tax rate"),
        ParamScalar("DaysofWC", 45.0, description = "number of days of working capital"),
        ParamScalar("DebtRatio", 0.6, description = "debt ratio"),
        ParamScalar("COD", 0.06, description = "cost of debt"),
        ParamScalar("COE", 0.112, description = "cost of equity"),
    )

end

get_id(x::DistributionUtility) = x.id

function norm(data::Float64, norm_variable::String, norm_inputs::DataFrame)
    data =
        (data - norm_inputs[norm_inputs.id.==norm_variable, "min"][1]) / (
            norm_inputs[norm_inputs.id.==norm_variable, "max"][1] -
            norm_inputs[norm_inputs.id.==norm_variable, "min"][1]
        )
    return data
end

function reverse_norm(data::Float64, norm_variable::String, norm_inputs::DataFrame)
    data =
        data * (
            norm_inputs[norm_inputs.id.==norm_variable, "max"][1] -
            norm_inputs[norm_inputs.id.==norm_variable, "min"][1]
        ) + norm_inputs[norm_inputs.id.==norm_variable, "min"][1]
    return data
end

function existing_distribution_account(
    distribution_utility::DistributionUtility,
    model_data::HEMData,
    agent_store::AgentStore,
    reg_year::Int64,
)
    utility = get_agent(Utility, agent_store)
    customers = get_agent(CustomerGroup, agent_store)

    distribution_capex_balance_model = distribution_utility.distribution_capex_balance_model
    norm_inputs = distribution_utility.norm_inputs

    total_sale_initial =
        sum(
            customers.gamma[h] * model_data.omega[t] * customers.d_my[first(model_data.index_y_fix), h, t] for
            h in model_data.index_h, t in model_data.index_t
        ) +
        # export
        sum(
            model_data.omega[t] * utility.eximport_my[first(model_data.index_y_fix), t] for
            t in model_data.index_t
        ) -
        # DG
        sum(
            model_data.omega[t] * (
                customers.rho_DG[h, m, t] *
                customers.x_DG_E_my[first(model_data.index_y_fix), h, m]
            ) for t in model_data.index_t, h in model_data.index_h, m in customers.index_m
        )
    
    distribution_capex_balance_norm =
        distribution_capex_balance_model.constant +
        distribution_capex_balance_model.total_sales_coefficient * norm(total_sale_initial, "total_sales", norm_inputs) +
        distribution_capex_balance_model.residential_customer_coefficient *
        norm(customers.gamma[:Residential], "residential", norm_inputs) +
        distribution_capex_balance_model.commercial_customer_coefficient *
        norm(customers.gamma[:Commercial], "commercial", norm_inputs) +
        distribution_capex_balance_model.industrial_customer_coefficient *
        norm(customers.gamma[:Industrial], "industrial", norm_inputs)

    distribution_capex_balance_reverse = reverse_norm(distribution_capex_balance_norm, "distribution_capex_balance", norm_inputs)

    if reg_year - model_data.year_start > distribution_utility.beginning_balance_lifetime
        distribution_existing_balance = 0.0
        distribution_existing_annual_depreciation = 0.0
    else
        distribution_existing_balance = distribution_capex_balance_reverse * (1 - (reg_year - model_data.year_start) / distribution_utility.beginning_balance_lifetime)
        distribution_existing_annual_depreciation = distribution_capex_balance_reverse / distribution_utility.beginning_balance_lifetime
    end

    return distribution_existing_balance, distribution_existing_annual_depreciation

end


function new_distribution_account(
    distribution_utility::DistributionUtility,
    model_data::HEMData,
    agent_store::AgentStore,
    reg_year::Int64,
)
    utility = get_agent(Utility, agent_store)
    customers = get_agent(CustomerGroup, agent_store)
    reg_year_index = Symbol(Int(reg_year))

    distribution_capex_addition_model = distribution_utility.distribution_capex_addition_model
    distribution_om_cost_model = distribution_utility.distribution_om_cost_model
    norm_inputs = distribution_utility.norm_inputs

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

    # customer module updates customers.x_DG_E to be the DPV at the beginning of the year,
    # to use end of the year, add customers.x_DG_new_my[reg_year_index, h, m]    
    dpv_pca =
        sum(customers.x_DG_E[h, m] for h in model_data.index_h, m in customers.index_m)
    
    distribution_capex_addition_norm =
        distribution_capex_addition_model.constant +
        distribution_capex_addition_model.saidi_coefficient *
        norm(distribution_utility.SAIDI[reg_year_index], "saidi", norm_inputs) +
        distribution_capex_addition_model.dpv_coefficient * norm(dpv_pca, "dpv", norm_inputs) +
        distribution_capex_addition_model.total_sales_coefficient * norm(total_sale, "total_sales", norm_inputs) +
        distribution_capex_addition_model.residential_customer_coefficient *
        norm(customers.gamma[:Residential], "residential", norm_inputs) +
        distribution_capex_addition_model.commercial_customer_coefficient *
        norm(customers.gamma[:Commercial], "commercial", norm_inputs) +
        distribution_capex_addition_model.industrial_customer_coefficient *
        norm(customers.gamma[:Industrial], "industrial", norm_inputs)

    distribution_om_cost_norm =
        distribution_om_cost_model.constant +
        distribution_om_cost_model.saidi_coefficient *
        norm(distribution_utility.SAIDI[reg_year_index], "saidi", norm_inputs) +
        distribution_om_cost_model.total_sales_coefficient * norm(total_sale, "total_sales", norm_inputs) +
        distribution_om_cost_model.residential_customer_coefficient *
        norm(customers.gamma[:Residential], "residential", norm_inputs) +
        distribution_om_cost_model.commercial_customer_coefficient *
        norm(customers.gamma[:Commercial], "commercial", norm_inputs) +
        distribution_om_cost_model.industrial_customer_coefficient *
        norm(customers.gamma[:Industrial], "industrial", norm_inputs)

    distribution_capex_addition_reverse = reverse_norm(distribution_capex_addition_norm, "distribution_capex_addition", norm_inputs)
    distribution_om_cost_reverse = reverse_norm(distribution_om_cost_norm, "distribution_om_cost", norm_inputs)

    distribution_utility.DistCapExAddition_new_my[reg_year_index] = distribution_capex_addition_reverse
    distribution_utility.DistOMCost_new_my[reg_year_index] = distribution_om_cost_reverse

    DistADITNew = sum(
            distribution_utility.DistCapExAddition_new_my[Symbol(Int(y))] *
            (
                distribution_utility.DistCumuTaxDepre_new_my[Symbol(Int(reg_year - y + 1))] -
                distribution_utility.DistCumuAccoutDepre_new_my[Symbol(Int(reg_year - y + 1))]
            ) *
            distribution_utility.Tax +
            distribution_utility.DistITC_new_my[Symbol(Int(y))] *
            distribution_utility.DistCapExAddition_new_my[Symbol(Int(y))] *
            (1 - distribution_utility.DistCumuITCAmort_new_my[Symbol(Int(reg_year - y + 1))]) for
            y in model_data.year[first(model_data.index_y_fix)]:reg_year
    )

    DistRateBaseNoWC_new = sum(
            distribution_utility.DistCapExAddition_new_my[Symbol(Int(y))] *
            (1 - distribution_utility.DistCumuAccoutDepre_new_my[Symbol(Int(reg_year - y + 1))])
            for y in model_data.year[first(model_data.index_y_fix)]:reg_year
    ) - DistADITNew

    Dist_working_capital = distribution_utility.DaysofWC / 365 * distribution_om_cost_reverse

    DistRateBase_new = DistRateBaseNoWC_new + Dist_working_capital

    distribution_new_annual_accounting_depreciation =
        sum(
            distribution_utility.DistCapExAddition_new_my[Symbol(Int(y))] *
            distribution_utility.DistAnnualAccoutDepre_new_my[Symbol(Int(reg_year - y + 1))] for
            y in model_data.year[first(model_data.index_y_fix)]:reg_year
        )
    
    distribution_new_annual_tax_depreciation =
        sum(
            distribution_utility.DistCapExAddition_new_my[Symbol(Int(y))] *
            distribution_utility.DistAnnualTaxDepre_new_my[Symbol(Int(reg_year - y + 1))] for
            y in model_data.year[first(model_data.index_y_fix)]:reg_year
        )

    return DistRateBase_new, distribution_new_annual_accounting_depreciation, distribution_new_annual_tax_depreciation

end


function solve_agent_problem!(
    distribution_utility::DistributionUtility,
    distribution_utility_opts::AgentOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{<:MarketStructure,<:UseCase},
    agent_store::AgentStore,
    w_iter,
)

    regulator = get_agent(Regulator, agent_store)

    reg_year = model_data.year[first(model_data.index_y)]
    reg_year_index = Symbol(Int(reg_year))

    distribution_cost_before = regulator.distribution_cost

    distribution_existing_balance, distribution_existing_annual_depreciation = 
        existing_distribution_account(distribution_utility, model_data, agent_store, reg_year)                
    
    DistRateBase_new, distribution_new_annual_accounting_depreciation, distribution_new_annual_tax_depreciation = 
        new_distribution_account(distribution_utility, model_data, agent_store, reg_year)

    rate_base = distribution_existing_balance + DistRateBase_new

    debt_interest = rate_base * distribution_utility.DebtRatio * distribution_utility.COD

    return_to_equity = rate_base * (1 - distribution_utility.DebtRatio) * distribution_utility.COE

    income_tax =
        (
            return_to_equity * distribution_utility.Tax +
            (distribution_new_annual_accounting_depreciation - distribution_new_annual_tax_depreciation) * distribution_utility.Tax - 
                distribution_utility.DistITC_new_my[reg_year_index] *
                distribution_utility.DistCapExAddition_new_my[reg_year_index]
        ) / (1 - distribution_utility.Tax)

    operational_cost = distribution_utility.DistOMCost_new_my[reg_year_index]
    
    revenue_requirement =
        debt_interest + return_to_equity + income_tax + operational_cost + distribution_existing_annual_depreciation + distribution_new_annual_accounting_depreciation

    regulator.distribution_cost[reg_year_index] = revenue_requirement

    return compute_difference_percentage_one_norm([
        (distribution_cost_before, regulator.distribution_cost),
    ])

end
