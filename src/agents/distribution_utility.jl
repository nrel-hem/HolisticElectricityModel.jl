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
    current_year::Symbol

    SAIDI::ParamArray
    distribution_capex_balance_model::DistributionCapexBalanceModel
    distribution_capex_addition_model::DistributionCapexAdditionModel
    distribution_om_cost_model::DistributionOMCostModel
    capex_balance_norm_inputs::DataFrame
    capex_addition_norm_inputs::DataFrame
    om_cost_norm_inputs::DataFrame
    beginning_balance_lifetime::Float64
    DistCumuTaxDepre_new_my::ParamArray # cumulative tax depreciation of new capacity (%)
    DistCumuAccoutDepre_new_my::ParamArray # cumulative accounting depreciation of new capacity (%)
    DistITC_new_my::ParamArray # ITC of new capacity (%)
    DistCumuITCAmort_new_my::ParamArray # ITC amortization of new capacity (%)
    DistAnnualAccoutDepre_new_my::ParamArray # annual accounting depreciation of new capacity (%)
    DistAnnualTaxDepre_new_my::ParamArray # annual tax depreciation of new capacity (%) 
    DistCapExAddition_new_my::ParamArray
    DistOMCost_new_my::ParamArray
    Tax::ParamScalar
    DaysofWC::ParamScalar
    DebtRatio::ParamScalar
    COD::ParamScalar
    COE::ParamScalar
end

function DistributionUtility(
    input_filename::AbstractString, # folders where input csv files are stored
    model_data::HEMData; # HEMData
    id = DEFAULT_ID, # currently, id are all default
)

    regression_model_parameters = CSV.read(joinpath(input_filename, "regression_results.csv"), DataFrame)
    regression_inputs_scale = CSV.read(joinpath(input_filename, "regression_inputs_scale.csv"), DataFrame)

    distribution_capex_balance_model = DistributionCapexBalanceModel(
        ParamScalar("constant", regression_model_parameters[(regression_model_parameters.regression .== "capex_balance") .& (regression_model_parameters.variable .== "intercept"), :value][1]),
        ParamScalar("total_sales_coefficient", regression_model_parameters[(regression_model_parameters.regression .== "capex_balance") .& (regression_model_parameters.variable .== "Total Sales MWh"), :value][1]),
        ParamScalar("residential_customer_coefficient", regression_model_parameters[(regression_model_parameters.regression .== "capex_balance") .& (regression_model_parameters.variable .== "Residential Customers"), :value][1]),
        ParamScalar("commercial_customer_coefficient", regression_model_parameters[(regression_model_parameters.regression .== "capex_balance") .& (regression_model_parameters.variable .== "Commercial Customers"), :value][1]),
        ParamScalar("industrial_customer_coefficient", regression_model_parameters[(regression_model_parameters.regression .== "capex_balance") .& (regression_model_parameters.variable .== "Industrial Customers"), :value][1]),
    )

    distribution_capex_addition_model = DistributionCapexAdditionModel(
        ParamScalar("constant", regression_model_parameters[(regression_model_parameters.regression .== "capex_additions") .& (regression_model_parameters.variable .== "intercept"), :value][1]),
        ParamScalar("saidi_coefficient", regression_model_parameters[(regression_model_parameters.regression .== "capex_additions") .& (regression_model_parameters.variable .== "IEEE SAIDI Excluding MED"), :value][1]),
        ParamScalar("dpv_coefficient", regression_model_parameters[(regression_model_parameters.regression .== "capex_additions") .& (regression_model_parameters.variable .== "Solar Capacity MW"), :value][1]),
        ParamScalar("total_sales_coefficient", regression_model_parameters[(regression_model_parameters.regression .== "capex_additions") .& (regression_model_parameters.variable .== "Total Sales MWh"), :value][1]),
        ParamScalar("residential_customer_coefficient", regression_model_parameters[(regression_model_parameters.regression .== "capex_additions") .& (regression_model_parameters.variable .== "Residential Customers"), :value][1]),
        ParamScalar("commercial_customer_coefficient", regression_model_parameters[(regression_model_parameters.regression .== "capex_additions") .& (regression_model_parameters.variable .== "Commercial Customers"), :value][1]),
        ParamScalar("industrial_customer_coefficient", regression_model_parameters[(regression_model_parameters.regression .== "capex_additions") .& (regression_model_parameters.variable .== "Industrial Customers"), :value][1]),
    )

    distribution_om_cost_model = DistributionOMCostModel(
        ParamScalar("constant", regression_model_parameters[(regression_model_parameters.regression .== "om") .& (regression_model_parameters.variable .== "intercept"), :value][1]),
        ParamScalar("saidi_coefficient", regression_model_parameters[(regression_model_parameters.regression .== "om") .& (regression_model_parameters.variable .== "IEEE SAIDI Excluding MED"), :value][1]),
        ParamScalar("total_sales_coefficient", regression_model_parameters[(regression_model_parameters.regression .== "om") .& (regression_model_parameters.variable .== "Total Sales MWh"), :value][1]),
        ParamScalar("residential_customer_coefficient", regression_model_parameters[(regression_model_parameters.regression .== "om") .& (regression_model_parameters.variable .== "Residential Customers"), :value][1]),
        ParamScalar("commercial_customer_coefficient", regression_model_parameters[(regression_model_parameters.regression .== "om") .& (regression_model_parameters.variable .== "Commercial Customers"), :value][1]),
        ParamScalar("industrial_customer_coefficient", regression_model_parameters[(regression_model_parameters.regression .== "om") .& (regression_model_parameters.variable .== "Industrial Customers"), :value][1]),
    )

    # Use the same normalization factors as performing the national-level regression analysis
    capex_balance_norm_inputs = DataFrame(
        id = [
            "total_sales",
            "residential",
            "commercial",
            "industrial",
            "distribution_capex_balance",
        ],
        min = [
            regression_inputs_scale[(regression_inputs_scale.regression .== "capex_balance") .& (regression_inputs_scale.variable .== "Total Sales MWh"), :min][1],
            regression_inputs_scale[(regression_inputs_scale.regression .== "capex_balance") .& (regression_inputs_scale.variable .== "Residential Customers"), :min][1],
            regression_inputs_scale[(regression_inputs_scale.regression .== "capex_balance") .& (regression_inputs_scale.variable .== "Commercial Customers"), :min][1],
            regression_inputs_scale[(regression_inputs_scale.regression .== "capex_balance") .& (regression_inputs_scale.variable .== "Industrial Customers"), :min][1],
            regression_inputs_scale[(regression_inputs_scale.regression .== "capex_balance") .& (regression_inputs_scale.variable .== "Distribution CapEx Balance"), :min][1]
        ],
        max = [
            regression_inputs_scale[(regression_inputs_scale.regression .== "capex_balance") .& (regression_inputs_scale.variable .== "Total Sales MWh"), :max][1],
            regression_inputs_scale[(regression_inputs_scale.regression .== "capex_balance") .& (regression_inputs_scale.variable .== "Residential Customers"), :max][1],
            regression_inputs_scale[(regression_inputs_scale.regression .== "capex_balance") .& (regression_inputs_scale.variable .== "Commercial Customers"), :max][1],
            regression_inputs_scale[(regression_inputs_scale.regression .== "capex_balance") .& (regression_inputs_scale.variable .== "Industrial Customers"), :max][1],
            regression_inputs_scale[(regression_inputs_scale.regression .== "capex_balance") .& (regression_inputs_scale.variable .== "Distribution CapEx Balance"), :max][1]
        ]
    )

    capex_addition_norm_inputs = DataFrame(
        id = [
            "saidi",
            "dpv",
            "total_sales",
            "residential",
            "commercial",
            "industrial",
            "distribution_capex_addition",
        ],
        min = [
            regression_inputs_scale[(regression_inputs_scale.regression .== "capex_additions") .& (regression_inputs_scale.variable .== "IEEE SAIDI Excluding MED"), :min][1],
            regression_inputs_scale[(regression_inputs_scale.regression .== "capex_additions") .& (regression_inputs_scale.variable .== "Solar Capacity MW"), :min][1],
            regression_inputs_scale[(regression_inputs_scale.regression .== "capex_additions") .& (regression_inputs_scale.variable .== "Total Sales MWh"), :min][1],
            regression_inputs_scale[(regression_inputs_scale.regression .== "capex_additions") .& (regression_inputs_scale.variable .== "Residential Customers"), :min][1],
            regression_inputs_scale[(regression_inputs_scale.regression .== "capex_additions") .& (regression_inputs_scale.variable .== "Commercial Customers"), :min][1],
            regression_inputs_scale[(regression_inputs_scale.regression .== "capex_additions") .& (regression_inputs_scale.variable .== "Industrial Customers"), :min][1],
            regression_inputs_scale[(regression_inputs_scale.regression .== "capex_additions") .& (regression_inputs_scale.variable .== "Distribution CapEx Addition"), :min][1]
        ],
        max = [
            regression_inputs_scale[(regression_inputs_scale.regression .== "capex_additions") .& (regression_inputs_scale.variable .== "IEEE SAIDI Excluding MED"), :max][1],
            regression_inputs_scale[(regression_inputs_scale.regression .== "capex_additions") .& (regression_inputs_scale.variable .== "Solar Capacity MW"), :max][1],
            regression_inputs_scale[(regression_inputs_scale.regression .== "capex_additions") .& (regression_inputs_scale.variable .== "Total Sales MWh"), :max][1],
            regression_inputs_scale[(regression_inputs_scale.regression .== "capex_additions") .& (regression_inputs_scale.variable .== "Residential Customers"), :max][1],
            regression_inputs_scale[(regression_inputs_scale.regression .== "capex_additions") .& (regression_inputs_scale.variable .== "Commercial Customers"), :max][1],
            regression_inputs_scale[(regression_inputs_scale.regression .== "capex_additions") .& (regression_inputs_scale.variable .== "Industrial Customers"), :max][1],
            regression_inputs_scale[(regression_inputs_scale.regression .== "capex_additions") .& (regression_inputs_scale.variable .== "Distribution CapEx Addition"), :max][1]
        ],
    )

    om_cost_norm_inputs = DataFrame(
        id = [
            "saidi",
            "total_sales",
            "residential",
            "commercial",
            "industrial",
            "distribution_om_cost",
        ],
        min = [
            regression_inputs_scale[(regression_inputs_scale.regression .== "om") .& (regression_inputs_scale.variable .== "IEEE SAIDI Excluding MED"), :min][1],
            regression_inputs_scale[(regression_inputs_scale.regression .== "om") .& (regression_inputs_scale.variable .== "Total Sales MWh"), :min][1],
            regression_inputs_scale[(regression_inputs_scale.regression .== "om") .& (regression_inputs_scale.variable .== "Residential Customers"), :min][1],
            regression_inputs_scale[(regression_inputs_scale.regression .== "om") .& (regression_inputs_scale.variable .== "Commercial Customers"), :min][1],
            regression_inputs_scale[(regression_inputs_scale.regression .== "om") .& (regression_inputs_scale.variable .== "Industrial Customers"), :min][1],
            regression_inputs_scale[(regression_inputs_scale.regression .== "om") .& (regression_inputs_scale.variable .== "Distribution OM Cost"), :min][1]
        ],
        max = [
            regression_inputs_scale[(regression_inputs_scale.regression .== "om") .& (regression_inputs_scale.variable .== "IEEE SAIDI Excluding MED"), :max][1],
            regression_inputs_scale[(regression_inputs_scale.regression .== "om") .& (regression_inputs_scale.variable .== "Total Sales MWh"), :max][1],
            regression_inputs_scale[(regression_inputs_scale.regression .== "om") .& (regression_inputs_scale.variable .== "Residential Customers"), :max][1],
            regression_inputs_scale[(regression_inputs_scale.regression .== "om") .& (regression_inputs_scale.variable .== "Commercial Customers"), :max][1],
            regression_inputs_scale[(regression_inputs_scale.regression .== "om") .& (regression_inputs_scale.variable .== "Industrial Customers"), :max][1],
            regression_inputs_scale[(regression_inputs_scale.regression .== "om") .& (regression_inputs_scale.variable .== "Distribution OM Cost"), :max][1]
        ]
    )

    return DistributionUtility(
        id,
        first(model_data.index_y),
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
        capex_balance_norm_inputs,
        capex_addition_norm_inputs,
        om_cost_norm_inputs,
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

function normalize(data::Float64, norm_variable::String, norm_inputs::DataFrame)
    data =
        (data - norm_inputs[norm_inputs.id.==norm_variable, "min"][1]) / (
            norm_inputs[norm_inputs.id.==norm_variable, "max"][1] -
            norm_inputs[norm_inputs.id.==norm_variable, "min"][1]
        )
    return data
end

function denormalize(data::Float64, norm_variable::String, norm_inputs::DataFrame)
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
    capex_balance_norm_inputs = distribution_utility.capex_balance_norm_inputs

    total_sale_initial =
        sum(
            customers.gamma(z, h) * model_data.omega(d) * customers.d_my(first(model_data.index_y_fix), h, z, d, t) for
            z in model_data.index_z, h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
        ) +
        # export
        sum(
            model_data.omega(d) * utility.eximport_my(first(model_data.index_y_fix), z, d, t) for
            z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
        ) -
        # DG
        sum(
            model_data.omega(d) * (
                customers.rho_DG(h, m, z, d, t) *
                customers.x_DG_E_my(first(model_data.index_y_fix), h, z, m)
            ) for t in model_data.index_t, h in model_data.index_h, m in customers.index_m, z in model_data.index_z, d in model_data.index_d
        )
    
    distribution_capex_balance_norm =
        distribution_capex_balance_model.constant +
        distribution_capex_balance_model.total_sales_coefficient * normalize(total_sale_initial, "total_sales", capex_balance_norm_inputs) +
        distribution_capex_balance_model.residential_customer_coefficient *
        normalize(sum(customers.gamma(:, :Residential)), "residential", capex_balance_norm_inputs) +
        distribution_capex_balance_model.commercial_customer_coefficient *
        normalize(sum(customers.gamma(:, :Commercial)), "commercial", capex_balance_norm_inputs) +
        distribution_capex_balance_model.industrial_customer_coefficient *
        normalize(sum(customers.gamma(:, :Industrial)), "industrial", capex_balance_norm_inputs)

    distribution_capex_balance_reverse = denormalize(distribution_capex_balance_norm, "distribution_capex_balance", capex_balance_norm_inputs)

    if reg_year - model_data.year_start > distribution_utility.beginning_balance_lifetime
        distribution_existing_balance = 0.0
        distribution_existing_annual_depreciation = 0.0
    else
        distribution_existing_balance = distribution_capex_balance_reverse * (1 - (reg_year - model_data.year_start) / distribution_utility.beginning_balance_lifetime)
        distribution_existing_annual_depreciation = distribution_capex_balance_reverse / distribution_utility.beginning_balance_lifetime
    end

    distribution_existing_balance_per_MWh = distribution_existing_balance / total_sale_initial

    @info "Existing Distribution Balance: $distribution_existing_balance"
    @info "Existing Distribution Annual Depreciation: $distribution_existing_annual_depreciation"
    @info "Existing Distribution Balance per MWh: $distribution_existing_balance_per_MWh"

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
    capex_addition_norm_inputs = distribution_utility.capex_addition_norm_inputs
    om_cost_norm_inputs = distribution_utility.om_cost_norm_inputs

    total_sale =
        sum(
            customers.gamma(z, h) * model_data.omega(d) * customers.d(h, z, d, t) for
            h in model_data.index_h, t in model_data.index_t, z in model_data.index_z, d in model_data.index_d
        ) +
        # export
        sum(
            model_data.omega(d) * utility.eximport_my(reg_year_index, z, d, t) for
            z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
        ) -
        # DG
        sum(
            model_data.omega(d) * (
                customers.rho_DG(h, m, z, d, t) *
                customers.x_DG_E_my(first(model_data.index_y), h, z, m) + sum(
                    customers.rho_DG(h, m, z, d, t) * customers.x_DG_new_my(Symbol(Int(y)), h, z, m)
                    for y = model_data.year(first(model_data.index_y_fix)):reg_year
                )
            ) for z in model_data.index_z, d in model_data.index_d, t in model_data.index_t, h in model_data.index_h, m in customers.index_m
        )

    # customer module updates customers.x_DG_E to be the DPV at the beginning of the year,
    # to use end of the year, add customers.x_DG_new_my(reg_year_index, h, z, m)
    dpv_pca =
        sum(customers.x_DG_E(h, z, m) for h in model_data.index_h, z in model_data.index_z, m in customers.index_m)
    
    distribution_capex_addition_norm =
        distribution_capex_addition_model.constant +
        distribution_capex_addition_model.saidi_coefficient *
        normalize(distribution_utility.SAIDI(reg_year_index), "saidi", capex_addition_norm_inputs) +
        distribution_capex_addition_model.dpv_coefficient * normalize(dpv_pca, "dpv", capex_addition_norm_inputs) +
        distribution_capex_addition_model.total_sales_coefficient * normalize(total_sale, "total_sales", capex_addition_norm_inputs) +
        distribution_capex_addition_model.residential_customer_coefficient *
        normalize(sum(customers.gamma(:, :Residential)), "residential", capex_addition_norm_inputs) +
        distribution_capex_addition_model.commercial_customer_coefficient *
        normalize(sum(customers.gamma(:, :Commercial)), "commercial", capex_addition_norm_inputs) +
        distribution_capex_addition_model.industrial_customer_coefficient *
        normalize(sum(customers.gamma(:, :Industrial)), "industrial", capex_addition_norm_inputs)

    distribution_om_cost_norm =
        distribution_om_cost_model.constant +
        distribution_om_cost_model.saidi_coefficient *
        normalize(distribution_utility.SAIDI(reg_year_index), "saidi", om_cost_norm_inputs) +
        distribution_om_cost_model.total_sales_coefficient * normalize(total_sale, "total_sales", om_cost_norm_inputs) +
        distribution_om_cost_model.residential_customer_coefficient *
        normalize(sum(customers.gamma(:, :Residential)), "residential", om_cost_norm_inputs) +
        distribution_om_cost_model.commercial_customer_coefficient *
        normalize(sum(customers.gamma(:, :Commercial)), "commercial", om_cost_norm_inputs) +
        distribution_om_cost_model.industrial_customer_coefficient *
        normalize(sum(customers.gamma(:, :Industrial)), "industrial", om_cost_norm_inputs)

    distribution_capex_addition_reverse = denormalize(distribution_capex_addition_norm, "distribution_capex_addition", capex_addition_norm_inputs)
    distribution_om_cost_reverse = denormalize(distribution_om_cost_norm, "distribution_om_cost", om_cost_norm_inputs)

    distribution_utility.DistCapExAddition_new_my(reg_year_index, :) .= distribution_capex_addition_reverse
    distribution_utility.DistOMCost_new_my(reg_year_index, :) .= distribution_om_cost_reverse

    DistADITNew = sum(
            distribution_utility.DistCapExAddition_new_my(Symbol(Int(y))) *
            (
                distribution_utility.DistCumuTaxDepre_new_my(Symbol(Int(reg_year - y + 1))) -
                distribution_utility.DistCumuAccoutDepre_new_my(Symbol(Int(reg_year - y + 1)))
            ) *
            distribution_utility.Tax +
            distribution_utility.DistITC_new_my(Symbol(Int(y))) *
            distribution_utility.DistCapExAddition_new_my(Symbol(Int(y))) *
            (1 - distribution_utility.DistCumuITCAmort_new_my(Symbol(Int(reg_year - y + 1)))) for
            y in model_data.year(first(model_data.index_y_fix)):reg_year
    )

    DistRateBaseNoWC_new = sum(
            distribution_utility.DistCapExAddition_new_my(Symbol(Int(y))) *
            (1 - distribution_utility.DistCumuAccoutDepre_new_my(Symbol(Int(reg_year - y + 1))))
            for y in model_data.year(first(model_data.index_y_fix)):reg_year
    ) - DistADITNew

    Dist_working_capital = distribution_utility.DaysofWC / 365 * distribution_om_cost_reverse

    DistRateBase_new = DistRateBaseNoWC_new + Dist_working_capital

    distribution_new_annual_accounting_depreciation =
        sum(
            distribution_utility.DistCapExAddition_new_my(Symbol(Int(y))) *
            distribution_utility.DistAnnualAccoutDepre_new_my(Symbol(Int(reg_year - y + 1))) for
            y in model_data.year(first(model_data.index_y_fix)):reg_year
        )
    
    distribution_new_annual_tax_depreciation =
        sum(
            distribution_utility.DistCapExAddition_new_my(Symbol(Int(y))) *
            distribution_utility.DistAnnualTaxDepre_new_my(Symbol(Int(reg_year - y + 1))) for
            y in model_data.year(first(model_data.index_y_fix)):reg_year
        )

    DistRateBase_new_per_MWh = DistRateBase_new / total_sale

    @info "New Distribution Rate Base: $DistRateBase_new"
    @info "New Distribution Annual Accounting Depreciation: $distribution_new_annual_accounting_depreciation"
    @info "New Distribution Annual Tax Depreciation: $distribution_new_annual_tax_depreciation"
    @info "New Distribution Rate Base per MWh: $DistRateBase_new_per_MWh"

    return DistRateBase_new, distribution_new_annual_accounting_depreciation, distribution_new_annual_tax_depreciation

end


function solve_agent_problem!(
    distribution_utility::DistributionUtility,
    distribution_utility_opts::AgentOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{<:MarketStructure, <:UseCase, <:UseCase, <:UseCase},
    agent_store::AgentStore,
    w_iter,
    window_length,
    jump_model,
    export_file_path,
    update_results::Bool,
    output_intermediate_results::Bool
)

    regulator = get_agent(Regulator, agent_store)

    reg_year, reg_year_index = get_reg_year(model_data)

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
                distribution_utility.DistITC_new_my(reg_year_index) *
                distribution_utility.DistCapExAddition_new_my(reg_year_index)
        ) / (1 - distribution_utility.Tax)

    operational_cost = distribution_utility.DistOMCost_new_my(reg_year_index)
    
    revenue_requirement =
        debt_interest + return_to_equity + income_tax + operational_cost + distribution_existing_annual_depreciation + distribution_new_annual_accounting_depreciation

    regulator.distribution_cost(:, reg_year_index) .= revenue_requirement

    distribution_utility.current_year = reg_year_index

    return compute_difference_percentage_one_norm([
        (distribution_cost_before, regulator.distribution_cost),
    ])

end
