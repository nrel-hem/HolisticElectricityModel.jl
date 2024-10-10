# This module defines the data and functions associated with the regulator

# declare rate designs
abstract type RateDesign end
struct FlatRate <: RateDesign end
struct TOU <: RateDesign end

# declare net metering policies
abstract type NetMeteringPolicy end
struct ExcessRetailRate <: NetMeteringPolicy end
struct ExcessMarginalCost <: NetMeteringPolicy end
struct ExcessZero <: NetMeteringPolicy end

# declare regulator modeling options
abstract type AbstractRegulatorOptions <: AgentOptions end

struct RegulatorOptions{T <: RateDesign, U <: NetMeteringPolicy} <: AbstractRegulatorOptions
    rate_design::T
    net_metering_policy::U
end

function get_file_prefix(options::RegulatorOptions)
    return join(["$(typeof(options.rate_design))", 
                 "$(typeof(options.net_metering_policy))"], "_")
end

abstract type AbstractRegulator <: Agent end

mutable struct Regulator <: AbstractRegulator
    id::String
    index_rate_tou::Dimension
    rep_day_time_tou_mapping::DataFrame
    # Parameters
    "planning reserve (fraction)"
    r::ParamArray
    "allowed return on investment (fraction)"
    z::ParamArray
    distribution_cost::ParamArray
    administration_cost::ParamArray
    transmission_cost::ParamArray
    interconnection_cost::ParamArray
    system_cost::ParamArray
    "other cost not related to the optimization problem"
    othercost::ParamArray
    "Renewable Energy Credits"
    REC::ParamArray

    # Primal Variables
    "retail price"
    p::ParamArray
    "DER excess generation rate"
    p_ex::ParamArray
    "retail price of import/export"
    p_eximport::ParamArray
    "green tariff rate"
    p_green::ParamArray
    "revenue (requirement) of utility company"
    revenue_req::ParamScalar
    "cost of utility company (without return on equity)"
    cost::ParamScalar
    "multi-year retail price"
    p_my::ParamArray
    "multi-year DER excess generation rate"
    p_ex_my::ParamArray
    "multi-year retail price of import/export"
    p_eximport_my::ParamArray
    "revenue (requirement) of utility company by year"
    revenue_req_my::ParamArray
    "cost of utility company (without return on equity) by year"
    cost_my::ParamArray
    "debt interest by year"
    debt_interest_my::ParamArray
    "income tax by year"
    income_tax_my::ParamArray
    "operational cost by year"
    operational_cost_my::ParamArray
    "accounting depreciation by year"
    depreciation_my::ParamArray
    "tax depreciation by year"
    depreciation_tax_my::ParamArray

    p_regression::ParamArray
    p_my_regression::ParamArray
    p_td::ParamArray
    p_my_td::ParamArray
end

function Regulator(input_filename::String, model_data::HEMData; id = DEFAULT_ID)

    index_rate_tou = read_set(
        input_filename,
        "index_rate_tou",
        "index_rate_tou",
        prose_name = "index for time-of-use rates",
    )

    rep_day_time_tou_mapping = CSV.read(joinpath(input_filename, "rep_day_time_tou_mapping.csv"), DataFrame)

    distribution_cost = read_param(
        "distribution_cost",
        input_filename,
        "distribution_cost",
        model_data.index_y,
        [model_data.index_z],
        description = "distribution cost (dollar)",
    )

    administration_cost = read_param(
        "administration_cost",
        input_filename,
        "administration_cost",
        model_data.index_y,
        [model_data.index_z],
        description = "administration cost (dollar)",
    )

    transmission_cost = read_param(
        "transmission_cost",
        input_filename,
        "transmission_cost",
        model_data.index_y,
        [model_data.index_z],
        description = "transmission cost (dollar)",
    )

    interconnection_cost = read_param(
        "interconnection_cost",
        input_filename,
        "interconnection_cost",
        model_data.index_y,
        [model_data.index_z],
        description = "interconnection cost (dollar)",
    )

    system_cost = read_param(
        "system_cost",
        input_filename,
        "system_cost",
        model_data.index_y,
        [model_data.index_z],
        description = "system cost (dollar)",
    )

    return Regulator(
        id,
        index_rate_tou,
        rep_day_time_tou_mapping,
        initialize_param(
            "r",
            model_data.index_z,
            model_data.index_y;
            value = 0.12,
            description = "planning reserve (fraction)",
        ),
        initialize_param(
            "z",
            model_data.index_z,
            model_data.index_y;
            value = 0.112,
            description = "allowed return on investment (fraction)",
        ),
        distribution_cost,
        administration_cost,
        transmission_cost,
        interconnection_cost,
        system_cost,
        initialize_param(
            "othercost",
            model_data.index_z,
            model_data.index_y;
            description = "other cost not related to the optimization problem",
        ),
        initialize_param(
            "REC",
            model_data.index_z,
            model_data.index_y;
            value = 20.0,
            description = "Renewable Energy Credits",
        ),
        initialize_param(
            "p",
            model_data.index_z,
            model_data.index_h,
            model_data.index_d,
            model_data.index_t;
            value = 10.0,
            description = "retail price",
        ),
        initialize_param(
            "p_ex",
            model_data.index_z,
            model_data.index_h,
            model_data.index_d,
            model_data.index_t;
            value = 0.0,
            description = "DER excess generation rate",
        ),
        initialize_param(
            "p_eximport",
            model_data.index_z,
            model_data.index_d,
            model_data.index_t;
            value = 10.0,
            description = "retail price of import/export",
        ),
        initialize_param(
            "p_green",
            model_data.index_z,
            model_data.index_h,
            model_data.index_j;
            value = 20.0,
            description = "green tariff rate",
        ),
        ParamScalar(
            "revenue_req",
            0.00,
            description = "revenue (requirement) of utility company",
        ),
        ParamScalar(
            "cost",
            0.00,
            description = "cost of utility company (without return on equity)",
        ),
        initialize_param(
            "p_my",
            model_data.index_y,
            model_data.index_z,
            model_data.index_h,
            model_data.index_d,
            model_data.index_t;
            value = 10.0,
            description = "multi-year retail price",
        ),
        initialize_param(
            "p_ex_my",
            model_data.index_y,
            model_data.index_z,
            model_data.index_h,
            model_data.index_d,
            model_data.index_t;
            value = 0.0,
            description = "multi-year DER excess generation rate",
        ),
        initialize_param(
            "p_eximport_my",
            model_data.index_y,
            model_data.index_z,
            model_data.index_d,
            model_data.index_t;
            value = 10.0,
            description = "multi-year retail price of import/export",
        ),
        initialize_param(
            "revenue_req_my",
            model_data.index_y,
            model_data.index_z;
            description = "revenue (requirement) of utility company by year",
        ),
        initialize_param(
            "cost_my",
            model_data.index_y,
            model_data.index_z;
            description = "cost of utility company (without return on equity) by year",
        ),
        initialize_param(
            "debt_interest_my",
            model_data.index_y,
            model_data.index_z;
            description = "debt interest by year",
        ),
        initialize_param(
            "income_tax_my",
            model_data.index_y,
            model_data.index_z;
            description = "income tax by year",
        ),
        initialize_param(
            "operational_cost_my",
            model_data.index_y,
            model_data.index_z;
            description = "operational cost by year",
        ),
        initialize_param(
            "depreciation_my",
            model_data.index_y,
            model_data.index_z;
            description = "accounting depreciation by year",
        ),
        initialize_param(
            "depreciation_tax_my",
            model_data.index_y,
            model_data.index_z;
            description = "tax depreciation by year",
        ),
        initialize_param(
            "p_regression",
            model_data.index_h;
            value = 10.0,
            description = "retail price for regression (no T&D cost)",
        ),
        initialize_param(
            "p_my_regression",
            model_data.index_y,
            model_data.index_h;
            value = 10.0,
            description = "multi-year retail price for regression (no T&D cost)",
        ),
        initialize_param(
            "p_td",
            model_data.index_z,
            model_data.index_h;
            value = 0.0,
            description = "T&D component charge",
        ),
        initialize_param(
            "p_my_td",
            model_data.index_y,
            model_data.index_z,
            model_data.index_h;
            value = 0.0,
            description = "multi-year T&D component charge",
        ),
    )
end

get_id(x::Regulator) = x.id

function get_file_prefix(agent::Regulator)
    return "REC$(agent.REC.values[1, 1])"
end

# although Customer is subtype of Agent, 
# Vector{Customer} is not subtype of Vector{Agent}
# But if a vector of customers c1, c2, c3 is defined 
# using the syntax Agent[c1, c2, c3], calling 
# this function will work. Can also:
# Vector{Agent}([c1, c2, c3])

# function solve_agent_problem!(
#     regulator::Regulator,
#     regulator_opts::RegulatorOptions,
#     model_data::HEMData,
#     hem_opts::HEMOptions{VerticallyIntegratedUtility},
#     agent_store::AgentStore,
#     w_iter,
# )

#     for y in model_data.index_y_fix
#         regulator.othercost(y, :) .= regulator.distribution_cost(y) + regulator.administration_cost(y) + regulator.transmission_cost(y) + regulator.interconnection_cost(y) + regulator.system_cost(y)
#     end

#     utility = get_agent(Utility, agent_store)
#     customers = get_agent(CustomerGroup, agent_store)
#     green_developer = get_agent(GreenDeveloper, agent_store)

#     # the year regulator is making a rate case
#     reg_year, reg_year_index = get_reg_year(model_data)
#     # retirement up to the rate case year
#     reg_retirement = KeyedArray(
#         [
#             sum(
#                 utility.x_R_my(Symbol(Int(y_symbol)), k) for
#                 y_symbol in model_data.year(first(model_data.index_y_fix)):reg_year
#             ) for k in utility.index_k_existing
#         ];
#         [get_pair(utility.index_k_existing)]...
#     )

#     # pure volumetric rate
#     energy_cost =
#         sum(
#             model_data.omega(t) *
#             (utility.v_E_my(reg_year_index, k, t) * utility.y_E_my(reg_year_index, k, t))
#             for k in utility.index_k_existing, t in model_data.index_t
#         ) + sum(
#             model_data.omega(t) *
#             (utility.v_C_my(reg_year_index, k, t) * utility.y_C_my(reg_year_index, k, t))
#             for k in utility.index_k_new, t in model_data.index_t
#         )
#     fixed_om =
#         sum(
#             utility.fom_E_my(reg_year_index, k) * (utility.x_E_my(k) - reg_retirement(k))
#             for k in utility.index_k_existing
#         ) + sum(
#             utility.fom_C_my(Symbol(Int(y_symbol)), k) *
#             utility.x_C_my(Symbol(Int(y_symbol)), k) for
#             y_symbol in model_data.year(first(model_data.index_y_fix)):reg_year,
#             k in utility.index_k_new
#         )
#     operational_cost = energy_cost + fixed_om
#     working_capital = utility.DaysofWC / 365 * operational_cost

#     # calculate ADIT and rate base (no working capital) for new builds, "reg_year-y+1" represents the number of years since the new investment is made
#     ADITNew = KeyedArray(
#         [
#             sum(
#                 utility.CapEx_my(Symbol(Int(y)), k) *
#                 utility.x_C_my(Symbol(Int(y)), k) *
#                 (
#                     utility.CumuTaxDepre_new_my(Symbol(Int(reg_year - y + 1)), k) -
#                     utility.CumuAccoutDepre_new_my(Symbol(Int(reg_year - y + 1)), k)
#                 ) *
#                 utility.Tax +
#                 utility.ITC_new_my(Symbol(Int(y)), k) *
#                 utility.CapEx_my(Symbol(Int(y)), k) *
#                 utility.x_C_my(Symbol(Int(y)), k) *
#                 (1 - utility.CumuITCAmort_new_my(Symbol(Int(reg_year - y + 1)), k)) for
#                 y in model_data.year(first(model_data.index_y_fix)):reg_year
#             ) for k in utility.index_k_new
#         ];
#         [get_pair(utility.index_k_new)]...
#     )
#     RateBaseNoWC_new = KeyedArray(
#         [
#             sum(
#                 utility.CapEx_my(Symbol(Int(y)), k) *
#                 utility.x_C_my(Symbol(Int(y)), k) *
#                 (1 - utility.CumuAccoutDepre_new_my(Symbol(Int(reg_year - y + 1)), k))
#                 for y in model_data.year(first(model_data.index_y_fix)):reg_year
#             ) - ADITNew(k) for k in utility.index_k_new
#         ];
#         [get_pair(utility.index_k_new)]...
#     )

#     # calculate total rate base for the year of rate making
#     rate_base =
#         sum(
#             utility.RateBaseNoWC_existing_my(reg_year_index, k) *
#             (utility.x_E_my(k) - reg_retirement(k)) for k in utility.index_k_existing
#         ) +
#         sum(RateBaseNoWC_new(k) for k in utility.index_k_new) +
#         working_capital
#     debt_interest = rate_base * utility.DebtRatio * utility.COD
#     # calculate total depreciation
#     depreciation =
#     # annual depreciation on existing units that have not retired
#         sum(
#             utility.CapEx_existing_my(k) *
#             (utility.x_E_my(k) - reg_retirement(k)) *
#             utility.AnnualAccoutDepre_existing_my(reg_year_index, k) +
#             # existing units that are retired this year will incur their regular annual depreciation, as well as the remaining un-depreciated asset
#             utility.CapEx_existing_my(k) *
#             utility.x_R_my(reg_year_index, k) *
#             (
#                 utility.AnnualAccoutDepre_existing_my(reg_year_index, k) + 1 -
#                 utility.CumuAccoutDepre_existing_my(reg_year_index, k)
#             ) for k in utility.index_k_existing
#         ) +
#         # annual depreciation on new units
#         sum(
#             utility.CapEx_my(Symbol(Int(y)), k) *
#             utility.x_C_my(Symbol(Int(y)), k) *
#             utility.AnnualAccoutDepre_new_my(Symbol(Int(reg_year - y + 1)), k) for
#             y in model_data.year(first(model_data.index_y_fix)):reg_year,
#             k in utility.index_k_new
#         )
#     # calculate total tax depreciation
#     depreciation_tax =
#     # annual depreciation on existing units that have not retired
#         sum(
#             utility.CapEx_existing_my(k) *
#             (utility.x_E_my(k) - reg_retirement(k)) *
#             utility.AnnualTaxDepre_existing_my(reg_year_index, k) +
#             # existing units that are retired this year will incur their regular annual depreciation, as well as the remaining un-depreciated asset
#             utility.CapEx_existing_my(k) *
#             utility.x_R_my(reg_year_index, k) *
#             (
#                 utility.AnnualTaxDepre_existing_my(reg_year_index, k) + 1 -
#                 utility.CumuTaxDepre_existing_my(reg_year_index, k)
#             ) for k in utility.index_k_existing
#         ) +
#         # annual depreciation on new units
#         sum(
#             utility.CapEx_my(Symbol(Int(y)), k) *
#             utility.x_C_my(Symbol(Int(y)), k) *
#             utility.AnnualTaxDepre_new_my(Symbol(Int(reg_year - y + 1)), k) for
#             y in model_data.year(first(model_data.index_y_fix)):reg_year,
#             k in utility.index_k_new
#         )
#     return_to_equity = rate_base * (1 - utility.DebtRatio) * utility.COE
#     #=
#     income_tax = (return_to_equity -
#         sum(utility.CapEx_my(reg_year_index,k)*utility.x_C_my(reg_year_index,k)*utility.ITC_new_my(reg_year_index,k) for k in utility.index_k_new)) /
#         (1-utility.Tax) - return_to_equity
#     =#
#     income_tax =
#         (
#             return_to_equity * utility.Tax +
#             (depreciation - depreciation_tax) * utility.Tax - 
#             sum(
#                 utility.ITC_existing_my(k) *
#                 utility.CapEx_existing_my(k) *
#                 (utility.x_E_my(k) - reg_retirement(k)) *
#                 utility.AnnualITCAmort_existing_my(reg_year_index, k) +
#                 # existing units that are retired this year will incur their regular annual depreciation, as well as the remaining un-depreciated asset
#                 utility.ITC_existing_my(k) *
#                 utility.CapEx_existing_my(k) *
#                 utility.x_R_my(reg_year_index, k) *
#                 (
#                     utility.AnnualITCAmort_existing_my(reg_year_index, k) + 1 -
#                     utility.CumuITCAmort_existing_my(reg_year_index, k)
#                 ) for k in utility.index_k_existing
#             ) -
#             sum(
#                 utility.ITC_new_my(Symbol(Int(y)), k) *
#                 utility.CapEx_my(Symbol(Int(y)), k) *
#                 utility.x_C_my(Symbol(Int(y)), k) *
#                 utility.AnnualITCAmort_new_my(Symbol(Int(reg_year - y + 1)), k) for
#                 y in model_data.year(first(model_data.index_y_fix)):reg_year, k in utility.index_k_new
#             )
#         ) / (1 - utility.Tax)
#     # calculate revenue requirement
#     revenue_requirement =
#         debt_interest + return_to_equity + income_tax + operational_cost + depreciation

#     regulator.revenue_req_my(reg_year_index, :) .= revenue_requirement
#     regulator.cost_my(reg_year_index, :) .=
#         debt_interest + income_tax + operational_cost + depreciation

#     regulator.debt_interest_my(reg_year_index, :) .= debt_interest
#     regulator.income_tax_my(reg_year_index, :) .= income_tax
#     regulator.operational_cost_my(reg_year_index, :) .= operational_cost
#     regulator.depreciation_my(reg_year_index, :) .= depreciation
#     regulator.depreciation_tax_my(reg_year_index, :) .= depreciation_tax

#     # when it comes to net demand, calculate two values: the one with loss is used for cost allocation;
#     # the one without loss is used for rate calculation.
#     net_demand_w_loss = (
#         # demand
#         # when it comes to sharing the revenue requirement (cost), use load including distribution loss
#         # e.g. utility generation is 100 MW, 50 MW to serve load (including distribution loss),
#         #      50 MW for export. It makes sense to allocate the same cost for internal load and export.
#         sum(
#             customers.gamma(h) * model_data.omega(t) * customers.d(h, t) for
#             h in model_data.index_h, t in model_data.index_t
#         ) +
#         # export
#         sum(
#             model_data.omega(t) * utility.eximport_my(reg_year_index, t) for
#             t in model_data.index_t
#         ) -
#         # DG
#         # when it comes to cost allocation, simply use total DER generation to offset total load;
#         # this does not consider enery offset at the household level, which is only considered when 
#         # calculating retail rates.
#         sum(
#             model_data.omega(t) * (
#                 customers.rho_DG(h, m, t) *
#                 customers.x_DG_E_my(first(model_data.index_y), h, m) + sum(
#                     customers.rho_DG(h, m, t) * customers.x_DG_new_my(Symbol(Int(y)), h, m) for
#                     y in model_data.year(first(model_data.index_y_fix)):reg_year
#                 )
#             ) for t in model_data.index_t, h in model_data.index_h,
#             m in customers.index_m
#         ) -
#         # green technology subscription
#         sum(
#             model_data.omega(t) * utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
#             model_data.year(first(model_data.index_y_fix)):reg_year)
#             for t in model_data.index_t, j in model_data.index_j, h in model_data.index_h
#         )
#     )

#     # In the presence of distribution loss, multiply load (including distribution loss) by (1 - loss factor),
#     # DER generation above this value shall be compansated for DER excess credits
#     # e.g. DER generation is 100 MW, load (without loss) is 95 MW, receive 5 MW excess credits
#     der_excess_cost_h = KeyedArray(
#         [
#             sum(
#                 model_data.omega(t) *
#                 regulator.p_ex(h, t) *
#                 (
#                     max(
#                         0,
#                         customers.rho_DG(h, m, t) * customers.Opti_DG_E(h, m) -
#                         customers.d(h, t) * (1 - utility.loss_dist),
#                     ) * customers.x_DG_E_my(first(model_data.index_y), h, m) /
#                     customers.Opti_DG_E(h, m) + sum(
#                         max(
#                             0,
#                             customers.rho_DG(h, m, t) *
#                             customers.Opti_DG_my(Symbol(Int(y)), h, m) -
#                             customers.d(h, t) * (1 - utility.loss_dist),
#                         ) * customers.x_DG_new_my(Symbol(Int(y)), h, m) /
#                         customers.Opti_DG_my(Symbol(Int(y)), h, m) for
#                         y in model_data.year(first(model_data.index_y_fix)):reg_year
#                     )
#                 ) for t in model_data.index_t, m in customers.index_m
#             ) for h in model_data.index_h
#         ];
#         [get_pair(model_data.index_h)]...
#     )

#     net_demand_h_w_loss = KeyedArray(
#         # demand
#         [
#             sum(
#                 customers.gamma(h) * model_data.omega(t) * customers.d(h, t) for
#                 t in model_data.index_t
#             ) -
#             # DG
#             sum(
#                 model_data.omega(t) * (
#                     customers.rho_DG(h, m, t) *
#                     customers.x_DG_E_my(first(model_data.index_y), h, m) + sum(
#                         customers.rho_DG(h, m, t) *
#                         customers.x_DG_new_my(Symbol(Int(y)), h, m) for
#                         y in model_data.year(first(model_data.index_y_fix)):reg_year
#                     )
#                 ) for t in model_data.index_t, m in customers.index_m
#             ) -
#             # green technology subscription
#             sum(
#                 model_data.omega(t) * utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
#                 model_data.year(first(model_data.index_y_fix)):reg_year)
#                 for t in model_data.index_t, j in model_data.index_j
#             )
#             for h in model_data.index_h
#         ];
#         [get_pair(model_data.index_h)]...
#     )

#     net_demand_h_wo_loss = KeyedArray(
#         # demand
#         [
#             sum(
#                 customers.gamma(h) *
#                 model_data.omega(t) *
#                 customers.d(h, t) *
#                 (1 - utility.loss_dist) for t in model_data.index_t
#             ) -
#             # DG
#             # since this net demand is for rate calculation, we need to consider energy offset at the household level.
#             # e.g. two household, with 100 MW load each (without loss), if one of them has DER and generated 
#             # 120 MW of energy, he/she does not need to pay for energy, but the other one still have to pay 
#             # 100 MW instead of 80 MW.
#             sum(
#                 model_data.omega(t) * (
#                     min(
#                         customers.rho_DG(h, m, t) * customers.Opti_DG_E(h, m),
#                         customers.d(h, t) * (1 - utility.loss_dist),
#                     ) * customers.x_DG_E_my(first(model_data.index_y), h, m) /
#                     customers.Opti_DG_E(h, m) + sum(
#                         min(
#                             customers.rho_DG(h, m, t) *
#                             customers.Opti_DG_my(Symbol(Int(y)), h, m),
#                             customers.d(h, t) * (1 - utility.loss_dist),
#                         ) * customers.x_DG_new_my(Symbol(Int(y)), h, m) /
#                         customers.Opti_DG_my(Symbol(Int(y)), h, m) for
#                         y in model_data.year(first(model_data.index_y_fix)):reg_year
#                     )
#                 ) for t in model_data.index_t, m in customers.index_m
#             ) -
#             # green technology subscription
#             sum(
#                 model_data.omega(t) * utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
#                 model_data.year(first(model_data.index_y_fix)):reg_year)
#                 for t in model_data.index_t, j in model_data.index_j
#             )
#             for h in model_data.index_h
#         ];
#         [get_pair(model_data.index_h)]...
#     )

#     net_demand_wo_green_tech_h_wo_loss = KeyedArray(
#         [
#             # demand
#             sum(
#                 customers.gamma(h) *
#                 model_data.omega(t) *
#                 customers.d(h, t) *
#                 (1 - utility.loss_dist) for t in model_data.index_t
#             ) -
#             # DG
#             # since this net demand is for rate calculation, we need to consider energy offset at the household level.
#             # e.g. two household, with 100 MW load each (without loss), if one of them has DER and generated 
#             # 120 MW of energy, he/she does not need to pay for energy, but the other one still have to pay 
#             # 100 MW instead of 80 MW.
#             sum(
#                 model_data.omega(t) * (
#                     min(
#                         customers.rho_DG(h, m, t) * customers.Opti_DG_E(h, m),
#                         customers.d(h, t) * (1 - utility.loss_dist),
#                     ) * customers.x_DG_E_my(first(model_data.index_y), h, m) /
#                     customers.Opti_DG_E(h, m) + sum(
#                         min(
#                             customers.rho_DG(h, m, t) *
#                             customers.Opti_DG_my(Symbol(Int(y)), h, m),
#                             customers.d(h, t) * (1 - utility.loss_dist),
#                         ) * customers.x_DG_new_my(Symbol(Int(y)), h, m) /
#                         customers.Opti_DG_my(Symbol(Int(y)), h, m) for
#                         y in model_data.year(first(model_data.index_y_fix)):reg_year
#                     )
#                 ) for t in model_data.index_t, m in customers.index_m
#             )
#             for h in model_data.index_h
#         ];
#         [get_pair(model_data.index_h)]...
#     )

#     #=
#     energy_cost_t = Dict(t => sum(utility.v_E(k, t)*utility.y_E(k, t) for k in utility.index_k_existing) +
#         sum(utility.v_C(k, t)*utility.y_C(k, t) for k in utility.index_k_new) for t in model_data.index_t)
#     =#
#     energy_cost_t = KeyedArray(
#         [
#             sum(
#                 utility.v_E_my(reg_year_index, k, t) * utility.y_E_my(reg_year_index, k, t) for k in utility.index_k_existing
#             ) + sum(
#                 utility.v_C_my(reg_year_index, k, t) * utility.y_C_my(reg_year_index, k, t) for k in utility.index_k_new
#             ) for t in model_data.index_t
#         ];
#         [get_pair(model_data.index_t)]...
#     )

#     net_demand_t_w_loss = KeyedArray(
#         # demand
#         [
#             sum(customers.gamma(h) * customers.d(h, t) for h in model_data.index_h) +
#             # export
#             utility.eximport_my(reg_year_index, t) -
#             # DG
#             sum(
#                 customers.rho_DG(h, m, t) *
#                 customers.x_DG_E_my(first(model_data.index_y), h, m) + sum(
#                     customers.rho_DG(h, m, t) * customers.x_DG_new_my(Symbol(Int(y)), h, m) for
#                     y in model_data.year(first(model_data.index_y_fix)):reg_year
#                 ) for h in model_data.index_h, m in customers.index_m
#             ) -
#             # green technology subscription
#             sum(
#                 utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
#                 model_data.year(first(model_data.index_y_fix)):reg_year)
#                 for h in model_data.index_h, j in model_data.index_j
#             )
#             for t in model_data.index_t
#         ];
#         [get_pair(model_data.index_t)]...
#     )

#     der_excess_cost_h_t = make_keyed_array(model_data.index_h, model_data.index_t)
#     for h in model_data.index_h, t in model_data.index_t
#         der_excess_cost_h_t(h, t, :) .= sum(
#             regulator.p_ex(h, t) * (
#                 max(
#                     0,
#                     customers.rho_DG(h, m, t) * customers.Opti_DG_E(h, m) -
#                     customers.d(h, t) * (1 - utility.loss_dist),
#                 ) * customers.x_DG_E_my(first(model_data.index_y), h, m) /
#                 customers.Opti_DG_E(h, m) + sum(
#                     max(
#                         0,
#                         customers.rho_DG(h, m, t) *
#                         customers.Opti_DG_my(Symbol(Int(y)), h, m) -
#                         customers.d(h, t) * (1 - utility.loss_dist),
#                     ) * customers.x_DG_new_my(Symbol(Int(y)), h, m) /
#                     customers.Opti_DG_my(Symbol(Int(y)), h, m) for
#                     y in model_data.year(first(model_data.index_y_fix)):reg_year
#                 )
#             ) for m in customers.index_m
#         )
#     end

#     net_demand_h_t_w_loss = make_keyed_array(model_data.index_h, model_data.index_t)
#     for h in model_data.index_h, t in model_data.index_t
#         net_demand_h_t_w_loss(h, t, :) .=
#         # demand
#             customers.gamma(h) * customers.d(h, t) -
#             # DG
#             sum(
#                 customers.rho_DG(h, m, t) *
#                 customers.x_DG_E_my(first(model_data.index_y), h, m) + sum(
#                     customers.rho_DG(h, m, t) * customers.x_DG_new_my(Symbol(Int(y)), h, m)
#                     for y in model_data.year(first(model_data.index_y_fix)):reg_year
#                 ) for m in customers.index_m
#             ) -
#             sum(
#                 utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
#                 model_data.year(first(model_data.index_y_fix)):reg_year)
#                 for j in model_data.index_j
#             )
#     end

#     net_demand_h_t_wo_loss = make_keyed_array(model_data.index_h, model_data.index_t)
#     for h in model_data.index_h, t in model_data.index_t
#         net_demand_h_t_wo_loss(h, t, :) .=
#         # demand
#             customers.gamma(h) * customers.d(h, t) * (1 - utility.loss_dist) -
#             # DG
#             sum(
#                 min(
#                     customers.rho_DG(h, m, t) * customers.Opti_DG_E(h, m),
#                     customers.d(h, t) * (1 - utility.loss_dist),
#                 ) * customers.x_DG_E_my(first(model_data.index_y), h, m) /
#                 customers.Opti_DG_E(h, m) + sum(
#                     min(
#                         customers.rho_DG(h, m, t) *
#                         customers.Opti_DG_my(Symbol(Int(y)), h, m),
#                         customers.d(h, t) * (1 - utility.loss_dist),
#                     ) * customers.x_DG_new_my(Symbol(Int(y)), h, m) /
#                     customers.Opti_DG_my(Symbol(Int(y)), h, m) for
#                     y in model_data.year(first(model_data.index_y_fix)):reg_year
#                 ) for m in customers.index_m
#             ) -
#             sum(
#                 utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
#                 model_data.year(first(model_data.index_y_fix)):reg_year)
#                 for j in model_data.index_j
#             )
#     end

#     # for the purpose of calculating net peak load, use load including distribution loss
#     net_demand_for_peak_h_t = make_keyed_array(model_data.index_h, model_data.index_t)
#     for h in model_data.index_h, t in model_data.index_t
#         net_demand_for_peak_h_t(h, t, :) .=
#             customers.gamma(h) * customers.d(h, t) - sum(
#                 customers.rho_DG(h, m, t) *
#                 customers.x_DG_E_my(first(model_data.index_y), h, m) + sum(
#                     customers.rho_DG(h, m, t) * customers.x_DG_new_my(Symbol(Int(y)), h, m)
#                     for y in model_data.year(first(model_data.index_y_fix)):reg_year
#                 ) for m in customers.index_m
#             ) -
#             sum(
#                 utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
#                 model_data.year(first(model_data.index_y_fix)):reg_year)
#                 for j in model_data.index_j
#             )
#     end

#     # excluding green tech generation offset when it comes to sharing T&D costs
#     net_demand_for_peak_wo_green_tech_h_t = make_keyed_array(model_data.index_h, model_data.index_t)
#     for h in model_data.index_h, t in model_data.index_t
#         net_demand_for_peak_wo_green_tech_h_t(h, t, :) .=
#             customers.gamma(h) * customers.d(h, t) - sum(
#                 customers.rho_DG(h, m, t) *
#                 customers.x_DG_E_my(first(model_data.index_y), h, m) + sum(
#                     customers.rho_DG(h, m, t) * customers.x_DG_new_my(Symbol(Int(y)), h, m) for
#                     y in model_data.year(first(model_data.index_y_fix)):reg_year
#                 ) for m in customers.index_m
#             )
#     end

#     net_peak_load_h = KeyedArray(
#         [
#             findmax(Dict(t => net_demand_for_peak_h_t(h, t) for t in model_data.index_t))[1] for h in model_data.index_h
#         ];
#         [get_pair(model_data.index_h)]...
#     )

#     net_peak_load_wo_green_tech_h = KeyedArray(
#         [ 
#             findmax(Dict(t => net_demand_for_peak_wo_green_tech_h_t(h, t) for t in model_data.index_t))[1] for h in model_data.index_h
#         ];
#         [get_pair(model_data.index_h)]...
#     )

#     # Cost Classification/Allocation
#     energy_cost_allocation_h = KeyedArray(
#         [
#             energy_cost * net_demand_h_w_loss(h) / net_demand_w_loss + der_excess_cost_h(h) for h in model_data.index_h
#         ];
#         [get_pair(model_data.index_h)]...
#     )
#     energy_cost_allocation_eximport =
#         energy_cost * sum(
#             model_data.omega(t) * utility.eximport_my(reg_year_index, t) for
#             t in model_data.index_t
#         ) / net_demand_w_loss

#     # allocate (revenue_requirement - energy_cost) by net peak load with green tech generation offset;
#     # allocate regulator.othercost (T&D costs) by net peak load without green tech generation offset;
#     demand_cost_allocation_h = KeyedArray(
#         [
#             (revenue_requirement - energy_cost) * net_peak_load_h(h) / (
#                 sum(net_peak_load_h(h) for h in model_data.index_h) +
#                 utility.Peak_eximport_my(reg_year_index)
#             ) + 
#             regulator.othercost(reg_year_index) * net_peak_load_wo_green_tech_h(h) / (
#                 sum(net_peak_load_wo_green_tech_h(h) for h in model_data.index_h) +
#                 utility.Peak_eximport_my(reg_year_index)
#             )
#             for h in model_data.index_h
#         ];
#         [get_pair(model_data.index_h)]...
#     )
#     demand_cost_allocation_capacity_h = KeyedArray(
#         [
#             (revenue_requirement - energy_cost) * net_peak_load_h(h) / (
#                 sum(net_peak_load_h(h) for h in model_data.index_h) +
#                 utility.Peak_eximport_my(reg_year_index)
#             )
#             for h in model_data.index_h
#         ];
#         [get_pair(model_data.index_h)]...
#     )
#     demand_cost_allocation_othercost_h = KeyedArray(
#         [
#             regulator.othercost(reg_year_index) * net_peak_load_wo_green_tech_h(h) / (
#                 sum(net_peak_load_wo_green_tech_h(h) for h in model_data.index_h) +
#                 utility.Peak_eximport_my(reg_year_index)
#             )
#             for h in model_data.index_h
#         ];
#         [get_pair(model_data.index_h)]...
#     )

#     demand_cost_allocation_eximport =
#         (revenue_requirement - energy_cost) *
#         utility.Peak_eximport_my(reg_year_index) / (
#             sum(net_peak_load_h(h) for h in model_data.index_h) +
#             utility.Peak_eximport_my(reg_year_index)
#         ) +
#         regulator.othercost(reg_year_index) *
#         utility.Peak_eximport_my(reg_year_index) / (
#             sum(net_peak_load_wo_green_tech_h(h) for h in model_data.index_h) +
#             utility.Peak_eximport_my(reg_year_index)
#         )

#     energy_cost_allocation_h_t = make_keyed_array(model_data.index_h, model_data.index_t)
#     for h in model_data.index_h, t in model_data.index_t
#         energy_cost_allocation_h_t(h, t, :) .=
#             energy_cost_t(t) * net_demand_h_t_w_loss(h, t) / net_demand_t_w_loss(t) +
#             der_excess_cost_h_t(h, t)
#     end
#     energy_cost_allocation_eximport_t = KeyedArray(
#         [
#             energy_cost_t(t) * utility.eximport_my(reg_year_index, t) /
#             net_demand_t_w_loss(t) for t in model_data.index_t
#         ];
#         [get_pair(model_data.index_t)]...
#     )

#     # compute the retail price
#     p_before = ParamArray(regulator.p, "p_before")
#     fill!(p_before, NaN)  # TODO DT: debug only
#     for h in model_data.index_h, t in model_data.index_t
#         p_before(h, t, :) .= regulator.p_my(reg_year_index, h, t)
#     end

#     p_before_wavg = ParamArray(regulator.p_td, "p_before_wavg")
#     fill!(p_before_wavg, NaN)  # TODO DT: debug only
#     for h in model_data.index_h
#         p_before_wavg(h, :) .= 
#             sum(regulator.p_my(reg_year_index, h, t) * model_data.omega(t) * customers.d(h, t) for t in model_data.index_t) / 
#             sum(model_data.omega(t) * customers.d(h, t) for t in model_data.index_t)
#     end

#     p_ex_before = ParamArray(regulator.p_ex, "p_ex_before")
#     fill!(p_ex_before, NaN)
#     for h in model_data.index_h, t in model_data.index_t
#         p_ex_before(h, t, :) .= regulator.p_ex_my(reg_year_index, h, t)
#     end

#     # TODO: Call a function instead of using if-then
#     if regulator_opts.rate_design isa FlatRate
#         fill!(regulator.p, NaN)
#         for h in model_data.index_h, t in model_data.index_t
#             regulator.p(h, t, :) .=
#                 (energy_cost_allocation_h(h) + demand_cost_allocation_capacity_h(h)) /
#                 net_demand_h_wo_loss(h) +
#                 demand_cost_allocation_othercost_h(h) / net_demand_wo_green_tech_h_wo_loss(h)
#         end
#         fill!(regulator.p_eximport, NaN)
#         for t in model_data.index_t
#             regulator.p_eximport(t, :) .=
#                 (energy_cost_allocation_eximport + demand_cost_allocation_eximport) /
#                 sum(
#                     model_data.omega(t) * utility.eximport_my(reg_year_index, t) for
#                     t in model_data.index_t
#                 )
#         end
#     elseif regulator_opts.rate_design isa TOU
#         fill!(regulator.p, NaN)
#         for h in model_data.index_h, t in model_data.index_t
#             regulator.p(h, t, :) .=
#                 energy_cost_allocation_h_t(h, t) / net_demand_h_t_wo_loss(h, t) +
#                 demand_cost_allocation_capacity_h(h) / net_demand_h_wo_loss(h) +
#                 demand_cost_allocation_othercost_h(h) / net_demand_wo_green_tech_h_wo_loss(h)
#         end
#         fill!(regulator.p_eximport, NaN)
#         for t in model_data.index_t
#             regulator.p_eximport(t, :) .=
#                 energy_cost_allocation_eximport_t(t) /
#                 utility.eximport_my(reg_year_index, t) +
#                 demand_cost_allocation_eximport / sum(
#                     model_data.omega(t) * utility.eximport_my(reg_year_index, t) for
#                     t in model_data.index_t
#                 )
#         end
#     end

#     fill!(regulator.p_regression, NaN)
#     for h in model_data.index_h
#         regulator.p_regression(h, :) .=
#             (energy_cost_allocation_h(h) + demand_cost_allocation_capacity_h(h)) /
#             net_demand_h_wo_loss(h)
#     end

#     fill!(regulator.p_td, NaN)
#     for h in model_data.index_h
#         regulator.p_td(h, :) .=
#             demand_cost_allocation_othercost_h(h) / net_demand_wo_green_tech_h_wo_loss(h)
#     end

#     # TODO: Call a function instead of using if-then
#     if regulator_opts.net_metering_policy isa ExcessRetailRate
#         regulator.p_ex = ParamArray(regulator.p)
#     elseif regulator_opts.net_metering_policy isa ExcessMarginalCost
#         fill!(regulator.p_ex, NaN)
#         for h in model_data.index_h, t in model_data.index_t
#             regulator.p_ex(h, t, :) .=
#                 utility.miu_my(reg_year_index, t) /
#                 (model_data.omega(t) * utility.pvf_onm(reg_year_index))
#         end
#     elseif regulator_opts.net_metering_policy isa ExcessZero
#         fill!(regulator.p_ex, 0.0)
#     end

#     for h in model_data.index_h, t in model_data.index_t
#         regulator.p_my(reg_year_index, h, t, :) .= regulator.p(h, t)
#         regulator.p_ex_my(reg_year_index, h, t, :) .= regulator.p_ex(h, t)
#         regulator.p_eximport_my(reg_year_index, t, :) .= regulator.p_eximport(t)
#     end

#     for h in model_data.index_h
#         regulator.p_my_regression(reg_year_index, h, :) .= regulator.p_regression(h)
#         regulator.p_my_td(reg_year_index, h, :) .= regulator.p_td(h)
#     end

#     p_after_wavg = ParamArray(regulator.p_td, "p_after_wavg")
#     fill!(p_after_wavg, NaN)  # TODO DT: debug only
#     for h in model_data.index_h
#         p_after_wavg(h, :) .= 
#             sum(regulator.p_my(reg_year_index, h, t) * model_data.omega(t) * customers.d(h, t) for t in model_data.index_t) / 
#             sum(model_data.omega(t) * customers.d(h, t) for t in model_data.index_t)
#     end

#     @info "Original retail price" p_before
#     @info "Original DER excess rate" p_ex_before
#     @info "New retail price" regulator.p
#     @info "New DER excess rate" regulator.p_ex

#     return compute_difference_percentage_one_norm([
#         (p_before_wavg.values, p_after_wavg.values),
#     ])
# end

# regulator's module with updated dimensions

function solve_agent_problem!(
    regulator::Regulator,
    regulator_opts::RegulatorOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{VerticallyIntegratedUtility},
    agent_store::AgentStore,
    w_iter,
    jump_model,
    export_file_path,
    update_results::Bool
)

    delta_t = get_delta_t(model_data)

    utility = get_agent(Utility, agent_store)
    customers = get_agent(CustomerGroup, agent_store)
    green_developer = get_agent(GreenDeveloper, agent_store)
    der_aggregator = get_agent(DERAggregator, agent_store)

    # the year regulator is making a rate case
    reg_year, reg_year_index = get_reg_year(model_data)
    # retirement up to the rate case year
    reg_retirement = sum(
        utility.x_R_my(Symbol(Int(y_symbol)), :, :) for
        y_symbol in model_data.year(first(model_data.index_y_fix)):reg_year
    )
    reg_retirement_stor = sum(
        utility.x_stor_R_my(Symbol(Int(y_symbol)), :, :) for
        y_symbol in model_data.year(first(model_data.index_y_fix)):reg_year
    )

    for h in model_data.index_h, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
        customers.d(h, z, d, t, :) .= customers.d_my(reg_year_index, h, z, d, t)
        # customers.DERGen(h, t, :) .= customers.DERGen_my(reg_year_index, h, t)
    end

    # since regulator problem is ahead of DERAggregator probelm, use previous year's aggregation results.
    reg_year_dera, reg_year_index_dera = get_prev_reg_year(model_data, w_iter)

    for z in model_data.index_z
        regulator.othercost(z, reg_year_index, :) .= regulator.distribution_cost(z, reg_year_index) + regulator.administration_cost(z, reg_year_index) + regulator.transmission_cost(z, reg_year_index) + regulator.interconnection_cost(z, reg_year_index) + regulator.system_cost(z, reg_year_index) + der_aggregator.revenue(reg_year_index_dera, z)
    end

    total_der_stor_capacity = make_keyed_array(model_data.index_z, model_data.index_h)
    # this total_der_pv_capacity is the approximate capacity of pv portion of pv+storage tech, not all dpv capacity
    total_der_pv_capacity = make_keyed_array(model_data.index_z, model_data.index_h)
    for z in model_data.index_z, h in model_data.index_h
        if w_iter >= 2
            total_der_stor_capacity(z, h, :) .=
                customers.x_DG_E_my(reg_year_index_dera, h, z, :BTMStorage) + sum(
                    customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMStorage) for
                    y in model_data.year(first(model_data.index_y_fix)):reg_year_dera
                )
        else
            total_der_stor_capacity(z, h, :) .= customers.x_DG_E_my(reg_year_index_dera, h, z, :BTMStorage)
        end
        total_der_pv_capacity(z, h, :) .= total_der_stor_capacity(z, h) / customers.Opti_DG_E(z, h, :BTMStorage) * customers.Opti_DG_E(z, h, :BTMPV)
    end
    
    # pure volumetric rate
    energy_cost = make_keyed_array(model_data.index_z)
    for z in model_data.index_z
        energy_cost(z, :) .=
            sum(
                model_data.omega(d) * delta_t *
                (utility.v_E_my(reg_year_index, k, z, d, t) * utility.y_E_my(reg_year_index, k, z, d, t))
                for k in utility.index_k_existing, d in model_data.index_d, t in model_data.index_t
            ) + sum(
                model_data.omega(d) * delta_t *
                (utility.v_C_my(reg_year_index, k, z, d, t) * utility.y_C_my(reg_year_index, k, z, d, t))
                for k in utility.index_k_new, d in model_data.index_d, t in model_data.index_t
            )
    end
    fixed_om = make_keyed_array(model_data.index_z)
    for z in model_data.index_z
        fixed_om(z, :) .=
            sum(
                utility.fom_E_my(reg_year_index, z, k) * (utility.x_E_my(z, k) - reg_retirement(k, z))
                for k in utility.index_k_existing
            ) + sum(
                utility.fom_C_my(Symbol(Int(y_symbol)), z, k) *
                utility.x_C_my(Symbol(Int(y_symbol)), k, z) for
                y_symbol in model_data.year(first(model_data.index_y_fix)):reg_year,
                k in utility.index_k_new
            ) + sum(
                utility.fom_stor_E_my(reg_year_index, z, s) * (utility.x_stor_E_my(z, s) - reg_retirement_stor(s, z))
                for s in utility.index_stor_existing
            ) + sum(
                utility.fom_stor_C_my(Symbol(Int(y_symbol)), z, s) *
                utility.x_stor_C_my(Symbol(Int(y_symbol)), s, z) for
                y_symbol in model_data.year(first(model_data.index_y_fix)):reg_year,
                s in utility.index_stor_new
            )
    end

    eximport_cost = make_keyed_array(model_data.index_z, utility.index_l, model_data.index_d, model_data.index_t)
    for z in model_data.index_z, l in utility.index_l, d in model_data.index_d, t in model_data.index_t
        # import is charged at the sink (receiving node); export is credited at the sink node (receiving node)
        # cost here shows export cost being positive, import cost being negative
        if utility.trans_topology(l, z) * utility.flow_my(reg_year_index, l, d, t) >= 0.0 # export from zone z through line l, l is defined from zone z to somewhere else
            l_sink = utility.trans_topology.dims[2].elements[findall(x -> x == -1, utility.trans_topology(l, :))[1]]    # get the sink node of line l
            eximport_cost(z, l, d, t, :) .= utility.p_energy_cem_my(reg_year_index, l_sink, d, t) * utility.trans_topology(l, z) * utility.flow_my(reg_year_index, l, d, t) # if price is positive, this should be a positive value
        else    # import into zone z through l
            eximport_cost(z, l, d, t, :) .= utility.p_energy_cem_my(reg_year_index, z, d, t) * utility.trans_topology(l, z) * utility.flow_my(reg_year_index, l, d, t) # if price is positive, this should be a negative value
        end
    end

    net_eximport_cost = make_keyed_array(model_data.index_z)
    for z in model_data.index_z
        net_eximport_cost(z, :) .=
            # export and import cost of modeled areas
            sum(
                model_data.omega(d) * delta_t * eximport_cost(z, l, d, t)
                for l in utility.index_l, d in model_data.index_d, t in model_data.index_t
            ) + 
            # export and import cost outside of modeled areas
            sum(
                model_data.omega(d) * delta_t *
                utility.p_energy_cem_my(reg_year_index, z, d, t) * utility.eximport_my(reg_year_index, z, d, t)
                for d in model_data.index_d, t in model_data.index_t
            )
    end

    eximport_cap_cost = make_keyed_array(model_data.index_z, utility.index_l)
    for z in model_data.index_z, l in utility.index_l
        # import is charged at the sink (receiving node); export is credited at the sink node (receiving node)
        # cost here shows export cost being positive, import cost being negative
        if utility.trans_topology(l, z) * utility.flow_cap_my(reg_year_index, l) >= 0.0 # export from zone z through line l, l is defined from zone z to somewhere else
            l_sink = utility.trans_topology.dims[2].elements[findall(x -> x == -1, utility.trans_topology(l, :))[1]]    # get the sink node of line l
            eximport_cap_cost(z, l, :) .= utility.p_cap_cem_my(reg_year_index, l_sink) * utility.trans_topology(l, z) * utility.flow_cap_my(reg_year_index, l) # if price is positive, this should be a positive value
        else    # import into zone z through l
            eximport_cap_cost(z, l, :) .= utility.p_cap_cem_my(reg_year_index, z) * utility.trans_topology(l, z) * utility.flow_cap_my(reg_year_index, l) # if price is positive, this should be a negative value
        end
    end

    net_eximport_cap_cost = make_keyed_array(model_data.index_z)
    for z in model_data.index_z
        net_eximport_cap_cost(z, :) .=
            # export and import cost of modeled areas
            sum(
                eximport_cap_cost(z, l) for l in utility.index_l
            ) + 
            # export and import cost outside of modeled areas
            utility.p_cap_cem_my(reg_year_index, z) * utility.eximport_my(reg_year_index, z, utility.Max_Net_Load_my_dict[reg_year_index, z][1], utility.Max_Net_Load_my_dict[reg_year_index, z][2])
    end

    operational_cost = energy_cost .- net_eximport_cost .+ fixed_om .- net_eximport_cap_cost
    working_capital = utility.DaysofWC * operational_cost / 365

    # calculate ADIT and rate base (no working capital) for new builds, "reg_year-y+1" represents the number of years since the new investment is made
    ADITNew = make_keyed_array(model_data.index_z, utility.index_k_new)
    for z in model_data.index_z
        for k in utility.index_k_new
            ADITNew(z, k, :) .= sum(
                utility.CapEx_my(Symbol(Int(y)), z, k) *
                utility.x_C_my(Symbol(Int(y)), k, z) *
                (
                    utility.CumuTaxDepre_new_my(Symbol(Int(reg_year - y + 1)), k) -
                    utility.CumuAccoutDepre_new_my(Symbol(Int(reg_year - y + 1)), k)
                ) *
                utility.Tax +
                utility.ITC_new_my(Symbol(Int(y)), k) *
                utility.CapEx_my(Symbol(Int(y)), z, k) *
                utility.x_C_my(Symbol(Int(y)), k, z) *
                (1 - utility.CumuITCAmort_new_my(Symbol(Int(reg_year - y + 1)), k)) for
                y in model_data.year(first(model_data.index_y_fix)):reg_year
            )
        end
    end
    RateBaseNoWC_new = make_keyed_array(model_data.index_z, utility.index_k_new)
    for z in model_data.index_z
        for k in utility.index_k_new
            RateBaseNoWC_new(z, k, :) .= sum(
                utility.CapEx_my(Symbol(Int(y)), z, k) *
                utility.x_C_my(Symbol(Int(y)), k, z) *
                (1 - utility.CumuAccoutDepre_new_my(Symbol(Int(reg_year - y + 1)), k))
                for y in model_data.year(first(model_data.index_y_fix)):reg_year
            ) - ADITNew(z, k)
        end
    end

    ADITStorNew = make_keyed_array(model_data.index_z, utility.index_stor_new)
    for z in model_data.index_z
        for s in utility.index_stor_new
            ADITStorNew(z, s, :) .= sum(
                utility.CapEx_stor_my(Symbol(Int(y)), z, s) *
                utility.x_stor_C_my(Symbol(Int(y)), s, z) *
                (
                    utility.CumuTaxDepreStor_new_my(Symbol(Int(reg_year - y + 1)), s) -
                    utility.CumuAccoutDepreStor_new_my(Symbol(Int(reg_year - y + 1)), s)
                ) *
                utility.Tax +
                utility.ITCStor_new_my(Symbol(Int(y)), s) *
                utility.CapEx_stor_my(Symbol(Int(y)), z, s) *
                utility.x_stor_C_my(Symbol(Int(y)), s, z) *
                (1 - utility.CumuITCAmortStor_new_my(Symbol(Int(reg_year - y + 1)), s)) for
                y in model_data.year(first(model_data.index_y_fix)):reg_year
            )
        end
    end
    RateBaseNoWCStor_new = make_keyed_array(model_data.index_z, utility.index_stor_new)
    for z in model_data.index_z
        for s in utility.index_stor_new
            RateBaseNoWCStor_new(z, s, :) .= sum(
                utility.CapEx_stor_my(Symbol(Int(y)), z, s) *
                utility.x_stor_C_my(Symbol(Int(y)), s, z) *
                (1 - utility.CumuAccoutDepreStor_new_my(Symbol(Int(reg_year - y + 1)), s))
                for y in model_data.year(first(model_data.index_y_fix)):reg_year
            ) - ADITStorNew(z, s)
        end
    end

    # calculate total rate base for the year of rate making
    rate_base = make_keyed_array(model_data.index_z)
    for z in model_data.index_z
        rate_base(z, :) .= 
        sum(
            utility.RateBaseNoWC_existing_my(reg_year_index, z, k) *
            (utility.x_E_my(z, k) - reg_retirement(k, z)) for k in utility.index_k_existing
        ) +
        sum(RateBaseNoWC_new(z, k) for k in utility.index_k_new) +
        sum(
            utility.RateBaseNoWCStor_existing_my(reg_year_index, z, s) *
            (utility.x_stor_E_my(z, s) - reg_retirement_stor(s, z)) for s in utility.index_stor_existing
        ) +
        sum(RateBaseNoWCStor_new(z, s) for s in utility.index_stor_new) +
        working_capital(z)
    end
    debt_interest = rate_base * utility.DebtRatio.value * utility.COD.value
    # calculate total depreciation
    depreciation = make_keyed_array(model_data.index_z)
    for z in model_data.index_z
        depreciation(z, :) .= 
        # annual depreciation on existing units that have not retired
            sum(
                utility.CapEx_existing_my(z, k) *
                (utility.x_E_my(z, k) - reg_retirement(k, z)) *
                utility.AnnualAccoutDepre_existing_my(reg_year_index, k) +
                # existing units that are retired this year will incur their regular annual depreciation, as well as the remaining un-depreciated asset
                utility.CapEx_existing_my(z, k) *
                utility.x_R_my(reg_year_index, k, z) *
                (
                    utility.AnnualAccoutDepre_existing_my(reg_year_index, k) + 1 -
                    utility.CumuAccoutDepre_existing_my(reg_year_index, k)
                ) for k in utility.index_k_existing
            ) +
            # annual depreciation on new units
            sum(
                utility.CapEx_my(Symbol(Int(y)), z, k) *
                utility.x_C_my(Symbol(Int(y)), k, z) *
                utility.AnnualAccoutDepre_new_my(Symbol(Int(reg_year - y + 1)), k) for
                y in model_data.year(first(model_data.index_y_fix)):reg_year,
                k in utility.index_k_new
            ) +
            sum(
                utility.CapExStor_existing_my(z, s) *
                (utility.x_stor_E_my(z, s) - reg_retirement_stor(s, z)) *
                utility.AnnualAccoutDepreStor_existing_my(reg_year_index, s) +
                # existing units that are retired this year will incur their regular annual depreciation, as well as the remaining un-depreciated asset
                utility.CapExStor_existing_my(z, s) *
                utility.x_stor_R_my(reg_year_index, s, z) *
                (
                    utility.AnnualAccoutDepreStor_existing_my(reg_year_index, s) + 1 -
                    utility.CumuAccoutDepreStor_existing_my(reg_year_index, s)
                ) for s in utility.index_stor_existing
            ) +
            # annual depreciation on new units
            sum(
                utility.CapEx_stor_my(Symbol(Int(y)), z, s) *
                utility.x_stor_C_my(Symbol(Int(y)), s, z) *
                utility.AnnualAccoutDepreStor_new_my(Symbol(Int(reg_year - y + 1)), s) for
                y in model_data.year(first(model_data.index_y_fix)):reg_year,
                s in utility.index_stor_new
            )
    end
    # calculate total tax depreciation
    depreciation_tax = make_keyed_array(model_data.index_z)
    for z in model_data.index_z
        depreciation_tax(z, :) .= 
        # annual depreciation on existing units that have not retired
            sum(
                utility.CapEx_existing_my(z, k) *
                (utility.x_E_my(z, k) - reg_retirement(k, z)) *
                utility.AnnualTaxDepre_existing_my(reg_year_index, k) +
                # existing units that are retired this year will incur their regular annual depreciation, as well as the remaining un-depreciated asset
                utility.CapEx_existing_my(z, k) *
                utility.x_R_my(reg_year_index, k, z) *
                (
                    utility.AnnualTaxDepre_existing_my(reg_year_index, k) + 1 -
                    utility.CumuTaxDepre_existing_my(reg_year_index, k)
                ) for k in utility.index_k_existing
            ) +
            # annual depreciation on new units
            sum(
                utility.CapEx_my(Symbol(Int(y)), z, k) *
                utility.x_C_my(Symbol(Int(y)), k, z) *
                utility.AnnualTaxDepre_new_my(Symbol(Int(reg_year - y + 1)), k) for
                y in model_data.year(first(model_data.index_y_fix)):reg_year,
                k in utility.index_k_new
            ) +
            sum(
                utility.CapExStor_existing_my(z, s) *
                (utility.x_stor_E_my(z, s) - reg_retirement_stor(s, z)) *
                utility.AnnualTaxDepreStor_existing_my(reg_year_index, s) +
                # existing units that are retired this year will incur their regular annual depreciation, as well as the remaining un-depreciated asset
                utility.CapExStor_existing_my(z, s) *
                utility.x_stor_R_my(reg_year_index, s, z) *
                (
                    utility.AnnualTaxDepreStor_existing_my(reg_year_index, s) + 1 -
                    utility.CumuTaxDepreStor_existing_my(reg_year_index, s)
                ) for s in utility.index_stor_existing
            ) +
            # annual depreciation on new units
            sum(
                utility.CapEx_stor_my(Symbol(Int(y)), z, s) *
                utility.x_stor_C_my(Symbol(Int(y)), s, z) *
                utility.AnnualTaxDepreStor_new_my(Symbol(Int(reg_year - y + 1)), s) for
                y in model_data.year(first(model_data.index_y_fix)):reg_year,
                s in utility.index_stor_new
            )
    end
    return_to_equity = rate_base * (1 - utility.DebtRatio.value) * utility.COE.value
    #=
    income_tax = (return_to_equity -
        sum(utility.CapEx_my(reg_year_index,k)*utility.x_C_my(reg_year_index,k)*utility.ITC_new_my(reg_year_index,k) for k in utility.index_k_new)) /
        (1-utility.Tax) - return_to_equity
    =#
    income_tax = make_keyed_array(model_data.index_z)
    for z in model_data.index_z
        income_tax(z, :) .= 
            (
                return_to_equity(z) * utility.Tax +
                (depreciation(z) - depreciation_tax(z)) * utility.Tax - 
                sum(
                    utility.ITC_existing_my(k) *
                    utility.CapEx_existing_my(z, k) *
                    (utility.x_E_my(z, k) - reg_retirement(k, z)) *
                    utility.AnnualITCAmort_existing_my(reg_year_index, k) +
                    # existing units that are retired this year will incur their regular annual depreciation, as well as the remaining un-depreciated asset
                    utility.ITC_existing_my(k) *
                    utility.CapEx_existing_my(z, k) *
                    utility.x_R_my(reg_year_index, k, z) *
                    (
                        utility.AnnualITCAmort_existing_my(reg_year_index, k) + 1 -
                        utility.CumuITCAmort_existing_my(reg_year_index, k)
                    ) for k in utility.index_k_existing
                ) -
                sum(
                    utility.ITC_new_my(Symbol(Int(y)), k) *
                    utility.CapEx_my(Symbol(Int(y)), z, k) *
                    utility.x_C_my(Symbol(Int(y)), k, z) *
                    utility.AnnualITCAmort_new_my(Symbol(Int(reg_year - y + 1)), k) for
                    y in model_data.year(first(model_data.index_y_fix)):reg_year, k in utility.index_k_new
                ) -
                sum(
                    utility.ITCStor_existing_my(s) *
                    utility.CapExStor_existing_my(z, s) *
                    (utility.x_stor_E_my(z, s) - reg_retirement_stor(s, z)) *
                    utility.AnnualITCAmortStor_existing_my(reg_year_index, s) +
                    # existing units that are retired this year will incur their regular annual depreciation, as well as the remaining un-depreciated asset
                    utility.ITCStor_existing_my(s) *
                    utility.CapExStor_existing_my(z, s) *
                    utility.x_stor_R_my(reg_year_index, s, z) *
                    (
                        utility.AnnualITCAmortStor_existing_my(reg_year_index, s) + 1 -
                        utility.CumuITCAmortStor_existing_my(reg_year_index, s)
                    ) for s in utility.index_stor_existing
                ) -
                sum(
                    utility.ITCStor_new_my(Symbol(Int(y)), s) *
                    utility.CapEx_stor_my(Symbol(Int(y)), z, s) *
                    utility.x_stor_C_my(Symbol(Int(y)), s, z) *
                    utility.AnnualITCAmortStor_new_my(Symbol(Int(reg_year - y + 1)), s) for
                    y in model_data.year(first(model_data.index_y_fix)):reg_year, s in utility.index_stor_new
                )
            ) / (1 - utility.Tax)
    end
    # calculate revenue requirement
    revenue_requirement =
        debt_interest .+ return_to_equity .+ income_tax .+ operational_cost .+ depreciation

    regulator.revenue_req_my(reg_year_index, :) .= revenue_requirement
    regulator.cost_my(reg_year_index, :) .=
        debt_interest .+ income_tax .+ operational_cost .+ depreciation

    regulator.debt_interest_my(reg_year_index, :) .= debt_interest
    regulator.income_tax_my(reg_year_index, :) .= income_tax
    regulator.operational_cost_my(reg_year_index, :) .= operational_cost
    regulator.depreciation_my(reg_year_index, :) .= depreciation
    regulator.depreciation_tax_my(reg_year_index, :) .= depreciation_tax

    # when it comes to net demand, calculate two values: the one with loss is used for cost allocation;
    # the one without loss is used for rate calculation.
    net_demand_w_loss = make_keyed_array(model_data.index_z)
    for z in model_data.index_z
        net_demand_w_loss(z, :) .= 
            # demand
            # when it comes to sharing the revenue requirement (cost), use load including distribution loss
            # e.g. utility generation is 100 MW, 50 MW to serve load (including distribution loss),
            #      50 MW for export. It makes sense to allocate the same cost for internal load and export.
            sum(
                customers.gamma(z, h) * model_data.omega(d) * delta_t * customers.d(h, z, d, t) for
                h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
            ) +
            # exogenous export/import
            sum(
                model_data.omega(d) * delta_t * utility.eximport_my(reg_year_index, z, d, t) for
                d in model_data.index_d, t in model_data.index_t
            ) +
            # endogenous export/import (flow out of zone z)
            sum(
                model_data.omega(d) * delta_t * utility.trans_topology(l, z) * utility.flow_my(reg_year_index, l, d, t) for
                l in utility.index_l, d in model_data.index_d, t in model_data.index_t
            ) -
            # DG
            # when it comes to cost allocation, simply use total DER generation to offset total load;
            # this does not consider enery offset at the household level, which is only considered when 
            # calculating retail rates (e.g., DER excess compensation).
            sum(
                model_data.omega(d) * delta_t * (
                    customers.rho_DG(h, m, z, d, t) *
                    customers.x_DG_E_my(first(model_data.index_y), h, z, m) + sum(
                        customers.rho_DG(h, m, z, d, t) * customers.x_DG_new_my(Symbol(Int(y)), h, z, m) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    )
                ) for h in model_data.index_h, m in customers.index_m, d in model_data.index_d, t in model_data.index_t
            ) + 
            # remove aggregated behind-the-meter storage/pv generation/consumption since they're front-of-the-meter now
            sum(
                model_data.omega(d) * delta_t * customers.rho_DG(h, :BTMStorage, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_stor_capacity(z, h) 
                for h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
            ) + 
            sum(
                model_data.omega(d) * delta_t * customers.rho_DG(h, :BTMPV, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_pv_capacity(z, h) 
                for h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
            ) -
            # green technology subscription
            sum(
                model_data.omega(d) * delta_t * utility.rho_C_my(j, z, d, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, z, h) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):reg_year)
                for j in model_data.index_j, h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
            )
    end

    net_demand_w_loss_no_eximport = make_keyed_array(model_data.index_z)
    for z in model_data.index_z
        net_demand_w_loss_no_eximport(z, :) .= 
            # demand
            # when it comes to sharing the revenue requirement (cost), use load including distribution loss
            # e.g. utility generation is 100 MW, 50 MW to serve load (including distribution loss),
            #      50 MW for export. It makes sense to allocate the same cost for internal load and export.
            sum(
                customers.gamma(z, h) * model_data.omega(d) * delta_t * customers.d(h, z, d, t) for
                h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
            ) -
            # DG
            # when it comes to cost allocation, simply use total DER generation to offset total load;
            # this does not consider enery offset at the household level, which is only considered when 
            # calculating retail rates (e.g., DER excess compensation).
            sum(
                model_data.omega(d) * delta_t * (
                    customers.rho_DG(h, m, z, d, t) *
                    customers.x_DG_E_my(first(model_data.index_y), h, z, m) + sum(
                        customers.rho_DG(h, m, z, d, t) * customers.x_DG_new_my(Symbol(Int(y)), h, z, m) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    )
                ) for h in model_data.index_h, m in customers.index_m, d in model_data.index_d, t in model_data.index_t
            ) + 
            # remove aggregated behind-the-meter storage/pv generation/consumption since they're front-of-the-meter now
            sum(
                model_data.omega(d) * delta_t * customers.rho_DG(h, :BTMStorage, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_stor_capacity(z, h) 
                for h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
            ) + 
            sum(
                model_data.omega(d) * delta_t * customers.rho_DG(h, :BTMPV, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_pv_capacity(z, h) 
                for h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
            ) -
            # green technology subscription
            sum(
                model_data.omega(d) * delta_t * utility.rho_C_my(j, z, d, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, z, h) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):reg_year)
                for j in model_data.index_j, h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
            )
    end

    # In the presence of distribution loss, multiply load (including distribution loss) by (1 - loss factor),
    # DER generation above this value shall be compansated for DER excess credits
    # e.g. DER generation is 100 MW, load (without loss) is 95 MW, receive 5 MW excess credits
    # Another important issue not considered here is that, if there are multiple DER technologies m, whether they are adopted by the same customer?
    # this will affect DER excess compensation
    # the script here assumes these DER tech m are adopted by different customers
    # need information on how many households adopt stand-alone tech, and how many adopt multiple techs?
    # TODO: we might want to iterate on regulator's module a few times since regulator.p_ex is used as inputs to calculate the rates.
    der_excess_cost_h = make_keyed_array(model_data.index_z, model_data.index_h)
    for z in model_data.index_z, h in model_data.index_h
        der_excess_cost_h(z, h, :) .= 
        # der excess for pv-only customers
        sum(
            model_data.omega(d) * delta_t *
            regulator.p_ex(z, h, d, t) *
            (
                max(
                    0,
                    customers.rho_DG(h, :BTMPV, z, d, t) * customers.Opti_DG_E(z, h, :BTMPV) -
                    customers.d(h, z, d, t) * (1 - utility.loss_dist),
                ) * customers.x_DG_E_my(first(model_data.index_y), h, z, :BTMPV) /
                customers.Opti_DG_E(z, h, :BTMPV) + sum(
                    max(
                        0,
                        customers.rho_DG(h, :BTMPV, z, d, t) *
                        customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMPV) -
                        customers.d(h, z, d, t) * (1 - utility.loss_dist),
                    ) * customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMPV) /
                    customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMPV) for
                    y in model_data.year(first(model_data.index_y_fix)):reg_year
                ) - 
                # minus pv portion of pv_storage tech
                max(
                    0,
                    customers.rho_DG(h, :BTMPV, z, d, t) * customers.Opti_DG_E(z, h, :BTMPV) -
                    customers.d(h, z, d, t) * (1 - utility.loss_dist),
                ) * total_der_pv_capacity(z, h) /
                customers.Opti_DG_E(z, h, :BTMPV)
            )
            for d in model_data.index_d, t in model_data.index_t
        ) + 
        # der excess for pv+storage customers who did not participate in aggregation
        sum(
            model_data.omega(d) * delta_t *
            regulator.p_ex(z, h, d, t) *
            (
                max(
                    0,
                    sum(customers.rho_DG(h, m, z, d, t) * customers.Opti_DG_E(z, h, m) for m in customers.index_m) -
                    customers.d(h, z, d, t) * (1 - utility.loss_dist),
                ) * customers.x_DG_E_my(first(model_data.index_y), h, z, :BTMStorage) /
                customers.Opti_DG_E(z, h, :BTMStorage) + sum(
                    max(
                        0,
                        sum(customers.rho_DG(h, m, z, d, t) *
                        customers.Opti_DG_my(Symbol(Int(y)), z, h, m) for m in customers.index_m) -
                        customers.d(h, z, d, t) * (1 - utility.loss_dist),
                    ) * customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMStorage) /
                    customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMStorage) for
                    y in model_data.year(first(model_data.index_y_fix)):reg_year
                )
            ) * (1 - der_aggregator.aggregation_level(reg_year_index_dera, z))
            for d in model_data.index_d, t in model_data.index_t
        )
    end

    net_demand_h_w_loss = make_keyed_array(model_data.index_z, model_data.index_h)
    for z in model_data.index_z, h in model_data.index_h
        # Demand
        net_demand_h_w_loss(z, h, :) .= sum(
                customers.gamma(z, h) * model_data.omega(d) * delta_t * customers.d(h, z, d, t) for
                d in model_data.index_d, t in model_data.index_t
            ) -
            # DG
            sum(
                model_data.omega(d) * delta_t * (
                    customers.rho_DG(h, m, z, d, t) *
                    customers.x_DG_E_my(first(model_data.index_y), h, z, m) + sum(
                        customers.rho_DG(h, m, z, d, t) * customers.x_DG_new_my(Symbol(Int(y)), h, z, m) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    )
                ) for m in customers.index_m, d in model_data.index_d, t in model_data.index_t
            ) + 
            # remove aggregated behind-the-meter storage/pv generation/consumption since they're front-of-the-meter now
            sum(model_data.omega(d) * delta_t * customers.rho_DG(h, :BTMStorage, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_stor_capacity(z, h) for d in model_data.index_d, t in model_data.index_t) + 
            sum(model_data.omega(d) * delta_t * customers.rho_DG(h, :BTMPV, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_pv_capacity(z, h) for d in model_data.index_d, t in model_data.index_t) -
            # green technology subscription
            sum(
                model_data.omega(d) * delta_t * utility.rho_C_my(j, z, d, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, z, h) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):reg_year)
                for j in model_data.index_j, d in model_data.index_d, t in model_data.index_t
            )
    end

    net_demand_h_wo_loss = make_keyed_array(model_data.index_z, model_data.index_h)
    for z in model_data.index_z, h in model_data.index_h
        net_demand_h_wo_loss(z, h, :) .= 
            # Demand without loss
            sum(
                customers.gamma(z, h) * model_data.omega(d) * delta_t * customers.d(h, z, d, t) * (1 - utility.loss_dist) for
                d in model_data.index_d, t in model_data.index_t
            ) -
            # DG
            # since this net demand is for rate calculation, we need to consider energy offset at the household level.
            # e.g. two household, with 100 MW load each (without loss), if one of them has DER and generated 
            # 120 MW of energy, he/she does not need to pay for energy, but the other one still have to pay 
            # 100 MW instead of 80 MW.

            # behind the meter generation for pv-only customers
            sum(
                model_data.omega(d) * delta_t *
                (
                    min(
                        customers.rho_DG(h, :BTMPV, z, d, t) * customers.Opti_DG_E(z, h, :BTMPV),
                        customers.d(h, z, d, t) * (1 - utility.loss_dist),
                    ) * customers.x_DG_E_my(first(model_data.index_y), h, z, :BTMPV) /
                    customers.Opti_DG_E(z, h, :BTMPV) + sum(
                        min(
                            customers.rho_DG(h, :BTMPV, z, d, t) *
                            customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMPV),
                            customers.d(h, z, d, t) * (1 - utility.loss_dist)
                        ) * customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMPV) /
                        customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMPV) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    ) - 
                    # minus pv portion of pv_storage tech
                    min(
                        customers.rho_DG(h, :BTMPV, z, d, t) * customers.Opti_DG_E(z, h, :BTMPV),
                        customers.d(h, z, d, t) * (1 - utility.loss_dist)
                    ) * total_der_pv_capacity(z, h) /
                    customers.Opti_DG_E(z, h, :BTMPV)
                )
                for d in model_data.index_d, t in model_data.index_t
            ) - 
            # behind the meter generation for pv+storage customers who did not participate in aggregation
            sum(
                model_data.omega(d) * delta_t *
                (
                    min(
                        sum(customers.rho_DG(h, m, z, d, t) * customers.Opti_DG_E(z, h, m) for m in customers.index_m),
                        customers.d(h, z, d, t) * (1 - utility.loss_dist)
                    ) * customers.x_DG_E_my(first(model_data.index_y), h, z, :BTMStorage) /
                    customers.Opti_DG_E(z, h, :BTMStorage) + sum(
                        min(
                            sum(customers.rho_DG(h, m, z, d, t) *
                            customers.Opti_DG_my(Symbol(Int(y)), z, h, m) for m in customers.index_m),
                            customers.d(h, z, d, t) * (1 - utility.loss_dist)
                        ) * customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMStorage) /
                        customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMStorage) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    )
                ) * (1 - der_aggregator.aggregation_level(reg_year_index_dera, z))
                for d in model_data.index_d, t in model_data.index_t
            ) -
            # green technology subscription
            sum(
                model_data.omega(d) * delta_t * utility.rho_C_my(j, z, d, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, z, h) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):reg_year)
                for j in model_data.index_j, d in model_data.index_d, t in model_data.index_t
            )
    end

    net_demand_wo_green_tech_h_wo_loss = make_keyed_array(model_data.index_z, model_data.index_h)
    for z in model_data.index_z, h in model_data.index_h
        net_demand_wo_green_tech_h_wo_loss(z, h, :) .= 
            # Demand with loss
            sum(
                customers.gamma(z, h) * model_data.omega(d) * delta_t * customers.d(h, z, d, t) * (1 - utility.loss_dist) for
                d in model_data.index_d, t in model_data.index_t
            ) -
            # DG
            # since this net demand is for rate calculation, we need to consider energy offset at the household level.
            # e.g. two household, with 100 MW load each (without loss), if one of them has DER and generated 
            # 120 MW of energy, he/she does not need to pay for energy, but the other one still have to pay 
            # 100 MW instead of 80 MW.
            # behind the meter generation for pv-only customers
            sum(
                model_data.omega(d) * delta_t *
                (
                    min(
                        customers.rho_DG(h, :BTMPV, z, d, t) * customers.Opti_DG_E(z, h, :BTMPV),
                        customers.d(h, z, d, t) * (1 - utility.loss_dist),
                    ) * customers.x_DG_E_my(first(model_data.index_y), h, z, :BTMPV) /
                    customers.Opti_DG_E(z, h, :BTMPV) + sum(
                        min(
                            customers.rho_DG(h, :BTMPV, z, d, t) *
                            customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMPV),
                            customers.d(h, z, d, t) * (1 - utility.loss_dist)
                        ) * customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMPV) /
                        customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMPV) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    ) - 
                    # minus pv portion of pv_storage tech
                    min(
                        customers.rho_DG(h, :BTMPV, z, d, t) * customers.Opti_DG_E(z, h, :BTMPV),
                        customers.d(h, z, d, t) * (1 - utility.loss_dist)
                    ) * total_der_pv_capacity(z, h) /
                    customers.Opti_DG_E(z, h, :BTMPV)
                )
                for d in model_data.index_d, t in model_data.index_t
            ) - 
            # behind the meter generation for pv+storage customers who did not participate in aggregation
            sum(
                model_data.omega(d) * delta_t *
                (
                    min(
                        sum(customers.rho_DG(h, m, z, d, t) * customers.Opti_DG_E(z, h, m) for m in customers.index_m),
                        customers.d(h, z, d, t) * (1 - utility.loss_dist)
                    ) * customers.x_DG_E_my(first(model_data.index_y), h, z, :BTMStorage) /
                    customers.Opti_DG_E(z, h, :BTMStorage) + sum(
                        min(
                            sum(customers.rho_DG(h, m, z, d, t) *
                            customers.Opti_DG_my(Symbol(Int(y)), z, h, m) for m in customers.index_m),
                            customers.d(h, z, d, t) * (1 - utility.loss_dist)
                        ) * customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMStorage) /
                        customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMStorage) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    )
                ) * (1 - der_aggregator.aggregation_level(reg_year_index_dera, z))
                for d in model_data.index_d, t in model_data.index_t
            )
    end

    #=
    energy_cost_t = Dict(t => sum(utility.v_E(k, t)*utility.y_E(k, t) for k in utility.index_k_existing) +
        sum(utility.v_C(k, t)*utility.y_C(k, t) for k in utility.index_k_new) for t in model_data.index_t)
    =#
    energy_cost_t = make_keyed_array(model_data.index_z, regulator.index_rate_tou)
    for z in model_data.index_z, tou in regulator.index_rate_tou
        energy_cost_t_temp = 0.0
        for i in 1:size(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :])[1]
            d = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_d][i])
            t = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_t][i])
            energy_cost_t_temp = energy_cost_t_temp + 
                sum(
                    model_data.omega(d) * delta_t *
                    (utility.v_E_my(reg_year_index, k, z, d, t) * utility.y_E_my(reg_year_index, k, z, d, t))
                    for k in utility.index_k_existing
                ) + sum(
                    model_data.omega(d) * delta_t *
                    (utility.v_C_my(reg_year_index, k, z, d, t) * utility.y_C_my(reg_year_index, k, z, d, t))
                    for k in utility.index_k_new
                )
        end
        energy_cost_t(z, tou, :) .= energy_cost_t_temp
    end

    net_eximport_cost_t = make_keyed_array(model_data.index_z, regulator.index_rate_tou)
    for z in model_data.index_z, tou in regulator.index_rate_tou
        net_eximport_cost_t_temp = 0.0
        for i in 1:size(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :])[1]
            d = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_d][i])
            t = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_t][i])
            net_eximport_cost_t_temp = net_eximport_cost_t_temp +
                # export and import cost of modeled areas
                sum(
                    model_data.omega(d) * delta_t * eximport_cost(z, l, d, t)
                    for l in utility.index_l
                ) + 
                # export and import cost outside of modeled areas
                sum(
                    model_data.omega(d) * delta_t *
                    utility.p_energy_cem_my(reg_year_index, z, d, t) * utility.eximport_my(reg_year_index, z, d, t)
                )
        end
        net_eximport_cost_t(z, tou, :) .= net_eximport_cost_t_temp
    end

    net_demand_t_w_loss = make_keyed_array(model_data.index_z, regulator.index_rate_tou)
    for z in model_data.index_z, tou in regulator.index_rate_tou
        net_demand_t_w_loss_temp = 0.0
        for i in 1:size(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :])[1]
            d = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_d][i])
            t = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_t][i])
            net_demand_t_w_loss_temp = net_demand_t_w_loss_temp + 
                # demand
                # when it comes to sharing the revenue requirement (cost), use load including distribution loss
                # e.g. utility generation is 100 MW, 50 MW to serve load (including distribution loss),
                #      50 MW for export. It makes sense to allocate the same cost for internal load and export.
                sum(
                    customers.gamma(z, h) * model_data.omega(d) * delta_t * customers.d(h, z, d, t) for
                    h in model_data.index_h
                ) +
                # exogenous export/import
                sum(
                    model_data.omega(d) * delta_t * utility.eximport_my(reg_year_index, z, d, t)
                ) +
                # endogenous export/import (flow out of zone z)
                sum(
                    model_data.omega(d) * delta_t * utility.trans_topology(l, z) * utility.flow_my(reg_year_index, l, d, t) for
                    l in utility.index_l
                ) -
                # DG
                # when it comes to cost allocation, simply use total DER generation to offset total load;
                # this does not consider enery offset at the household level, which is only considered when 
                # calculating retail rates (e.g., DER excess compensation).
                sum(
                    model_data.omega(d) * delta_t * (
                        customers.rho_DG(h, m, z, d, t) *
                        customers.x_DG_E_my(first(model_data.index_y), h, z, m) + sum(
                            customers.rho_DG(h, m, z, d, t) * customers.x_DG_new_my(Symbol(Int(y)), h, z, m) for
                            y in model_data.year(first(model_data.index_y_fix)):reg_year
                        )
                    ) for h in model_data.index_h, m in customers.index_m
                ) +
                # remove aggregated behind-the-meter storage/pv generation/consumption since they're front-of-the-meter now
                sum(
                    model_data.omega(d) * delta_t * customers.rho_DG(h, :BTMStorage, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_stor_capacity(z, h) 
                    for h in model_data.index_h
                ) + 
                sum(
                    model_data.omega(d) * delta_t * customers.rho_DG(h, :BTMPV, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_pv_capacity(z, h) 
                    for h in model_data.index_h
                ) -
                # green technology subscription
                sum(
                    model_data.omega(d) * delta_t * utility.rho_C_my(j, z, d, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, z, h) for y_symbol in
                    model_data.year(first(model_data.index_y_fix)):reg_year)
                    for j in model_data.index_j, h in model_data.index_h
                )
        end
        net_demand_t_w_loss(z, tou, :) .= net_demand_t_w_loss_temp
    end

    net_demand_t_w_loss_no_eximport = make_keyed_array(model_data.index_z, regulator.index_rate_tou)
    for z in model_data.index_z, tou in regulator.index_rate_tou
        net_demand_t_w_loss_no_eximport_temp = 0.0
        for i in 1:size(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :])[1]
            d = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_d][i])
            t = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_t][i])
            net_demand_t_w_loss_no_eximport_temp = net_demand_t_w_loss_no_eximport_temp + 
                # demand
                # when it comes to sharing the revenue requirement (cost), use load including distribution loss
                # e.g. utility generation is 100 MW, 50 MW to serve load (including distribution loss),
                #      50 MW for export. It makes sense to allocate the same cost for internal load and export.
                sum(
                    customers.gamma(z, h) * model_data.omega(d) * delta_t * customers.d(h, z, d, t) for
                    h in model_data.index_h
                ) -
                # DG
                # when it comes to cost allocation, simply use total DER generation to offset total load;
                # this does not consider enery offset at the household level, which is only considered when 
                # calculating retail rates (e.g., DER excess compensation).
                sum(
                    model_data.omega(d) * delta_t * (
                        customers.rho_DG(h, m, z, d, t) *
                        customers.x_DG_E_my(first(model_data.index_y), h, z, m) + sum(
                            customers.rho_DG(h, m, z, d, t) * customers.x_DG_new_my(Symbol(Int(y)), h, z, m) for
                            y in model_data.year(first(model_data.index_y_fix)):reg_year
                        )
                    ) for h in model_data.index_h, m in customers.index_m
                ) + 
                # remove aggregated behind-the-meter storage/pv generation/consumption since they're front-of-the-meter now
                sum(
                    model_data.omega(d) * delta_t * customers.rho_DG(h, :BTMStorage, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_stor_capacity(z, h) 
                    for h in model_data.index_h
                ) + 
                sum(
                    model_data.omega(d) * delta_t * customers.rho_DG(h, :BTMPV, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_pv_capacity(z, h) 
                    for h in model_data.index_h
                ) -
                # green technology subscription
                sum(
                    model_data.omega(d) * delta_t * utility.rho_C_my(j, z, d, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, z, h) for y_symbol in
                    model_data.year(first(model_data.index_y_fix)):reg_year)
                    for j in model_data.index_j, h in model_data.index_h
                )
        end
        net_demand_t_w_loss_no_eximport(z, tou, :) .= net_demand_t_w_loss_no_eximport_temp
    end

    exo_eximport_demand_t = make_keyed_array(model_data.index_z, regulator.index_rate_tou)
    for z in model_data.index_z, tou in regulator.index_rate_tou
        exo_eximport_demand_t_temp = 0.0
        for i in 1:size(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :])[1]
            d = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_d][i])
            t = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_t][i])
            exo_eximport_demand_t_temp = exo_eximport_demand_t_temp + 
                # exogenous export/import
                sum(
                    model_data.omega(d) * delta_t * utility.eximport_my(reg_year_index, z, d, t)
                )
        end
        exo_eximport_demand_t(z, tou, :) .= exo_eximport_demand_t_temp
    end

    edo_eximport_demand_t = make_keyed_array(model_data.index_z, utility.index_l, regulator.index_rate_tou)
    for z in model_data.index_z, l in utility.index_l, tou in regulator.index_rate_tou
        edo_eximport_demand_t_temp = 0.0
        for i in 1:size(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :])[1]
            d = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_d][i])
            t = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_t][i])
            edo_eximport_demand_t_temp = edo_eximport_demand_t_temp + 
                # endogenous export/import (flow out of zone z)
                sum(
                    model_data.omega(d) * delta_t * utility.trans_topology(l, z) * utility.flow_my(reg_year_index, l, d, t)
                )
        end
        edo_eximport_demand_t(z, l, tou, :) .= edo_eximport_demand_t_temp
    end

    der_excess_cost_h_t = make_keyed_array(model_data.index_z, model_data.index_h, regulator.index_rate_tou)
    for z in model_data.index_z, h in model_data.index_h, tou in regulator.index_rate_tou
        der_excess_cost_h_t_temp = 0.0
        for i in 1:size(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :])[1]
            d = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_d][i])
            t = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_t][i])
            der_excess_cost_h_t_temp = der_excess_cost_h_t_temp + 
            # der excess for pv-only customers
            model_data.omega(d) * delta_t *
            regulator.p_ex(z, h, d, t) *
            (
                max(
                    0,
                    customers.rho_DG(h, :BTMPV, z, d, t) * customers.Opti_DG_E(z, h, :BTMPV) -
                    customers.d(h, z, d, t) * (1 - utility.loss_dist),
                ) * customers.x_DG_E_my(first(model_data.index_y), h, z, :BTMPV) /
                customers.Opti_DG_E(z, h, :BTMPV) + sum(
                    max(
                        0,
                        customers.rho_DG(h, :BTMPV, z, d, t) *
                        customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMPV) -
                        customers.d(h, z, d, t) * (1 - utility.loss_dist),
                    ) * customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMPV) /
                    customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMPV) for
                    y in model_data.year(first(model_data.index_y_fix)):reg_year
                ) - 
                # minus pv portion of pv_storage tech
                max(
                    0,
                    customers.rho_DG(h, :BTMPV, z, d, t) * customers.Opti_DG_E(z, h, :BTMPV) -
                    customers.d(h, z, d, t) * (1 - utility.loss_dist),
                ) * total_der_pv_capacity(z, h) /
                customers.Opti_DG_E(z, h, :BTMPV)
            ) + 
            # der excess for pv+storage customers who did not participate in aggregation
            model_data.omega(d) * delta_t *
            regulator.p_ex(z, h, d, t) *
            (
                max(
                    0,
                    sum(customers.rho_DG(h, m, z, d, t) * customers.Opti_DG_E(z, h, m) for m in customers.index_m) -
                    customers.d(h, z, d, t) * (1 - utility.loss_dist),
                ) * customers.x_DG_E_my(first(model_data.index_y), h, z, :BTMStorage) /
                customers.Opti_DG_E(z, h, :BTMStorage) + sum(
                    max(
                        0,
                        sum(customers.rho_DG(h, m, z, d, t) *
                        customers.Opti_DG_my(Symbol(Int(y)), z, h, m) for m in customers.index_m) -
                        customers.d(h, z, d, t) * (1 - utility.loss_dist),
                    ) * customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMStorage) /
                    customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMStorage) for
                    y in model_data.year(first(model_data.index_y_fix)):reg_year
                )
            ) * (1 - der_aggregator.aggregation_level(reg_year_index_dera, z))
        end
        der_excess_cost_h_t(z, h, tou, :) .= der_excess_cost_h_t_temp
    end

    net_demand_h_t_w_loss = make_keyed_array(model_data.index_z, model_data.index_h, regulator.index_rate_tou)
    for z in model_data.index_z, h in model_data.index_h, tou in regulator.index_rate_tou
        net_demand_h_t_w_loss_temp = 0.0
        for i in 1:size(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :])[1]
            d = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_d][i])
            t = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_t][i])
            net_demand_h_t_w_loss_temp = net_demand_h_t_w_loss_temp +
                sum(
                    customers.gamma(z, h) * model_data.omega(d) * delta_t * customers.d(h, z, d, t)
                ) -
                # DG
                sum(
                    model_data.omega(d) * delta_t * (
                        customers.rho_DG(h, m, z, d, t) *
                        customers.x_DG_E_my(first(model_data.index_y), h, z, m) + sum(
                            customers.rho_DG(h, m, z, d, t) * customers.x_DG_new_my(Symbol(Int(y)), h, z, m) for
                            y in model_data.year(first(model_data.index_y_fix)):reg_year
                        )
                    ) for m in customers.index_m
                ) +
                # remove aggregated behind-the-meter storage/pv generation/consumption since they're front-of-the-meter now
                model_data.omega(d) * delta_t * customers.rho_DG(h, :BTMStorage, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_stor_capacity(z, h) + 
                model_data.omega(d) * delta_t * customers.rho_DG(h, :BTMPV, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_pv_capacity(z, h) -
                # green technology subscription
                sum(
                    model_data.omega(d) * delta_t * utility.rho_C_my(j, z, d, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, z, h) for y_symbol in
                    model_data.year(first(model_data.index_y_fix)):reg_year)
                    for j in model_data.index_j
                )
        end
        net_demand_h_t_w_loss(z, h, tou, :) .= net_demand_h_t_w_loss_temp
    end

    net_demand_h_t_wo_loss = make_keyed_array(model_data.index_z, model_data.index_h, regulator.index_rate_tou)
    for z in model_data.index_z, h in model_data.index_h, tou in regulator.index_rate_tou
        net_demand_h_t_wo_loss_temp = 0.0
        for i in 1:size(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :])[1]
            d = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_d][i])
            t = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_t][i])
            net_demand_h_t_wo_loss_temp = net_demand_h_t_wo_loss_temp +
                # Demand without loss
                sum(
                    customers.gamma(z, h) * model_data.omega(d) * delta_t * customers.d(h, z, d, t) * (1 - utility.loss_dist)
                ) -
                # DG
                # since this net demand is for rate calculation, we need to consider energy offset at the household level.
                # e.g. two household, with 100 MW load each (without loss), if one of them has DER and generated 
                # 120 MW of energy, he/she does not need to pay for energy, but the other one still have to pay 
                # 100 MW instead of 80 MW.
                # behind the meter generation for pv-only customers
                model_data.omega(d) * delta_t *
                (
                    min(
                        customers.rho_DG(h, :BTMPV, z, d, t) * customers.Opti_DG_E(z, h, :BTMPV),
                        customers.d(h, z, d, t) * (1 - utility.loss_dist),
                    ) * customers.x_DG_E_my(first(model_data.index_y), h, z, :BTMPV) /
                    customers.Opti_DG_E(z, h, :BTMPV) + sum(
                        min(
                            customers.rho_DG(h, :BTMPV, z, d, t) *
                            customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMPV),
                            customers.d(h, z, d, t) * (1 - utility.loss_dist)
                        ) * customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMPV) /
                        customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMPV) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    ) - 
                    # minus pv portion of pv_storage tech
                    min(
                        customers.rho_DG(h, :BTMPV, z, d, t) * customers.Opti_DG_E(z, h, :BTMPV),
                        customers.d(h, z, d, t) * (1 - utility.loss_dist)
                    ) * total_der_pv_capacity(z, h) /
                    customers.Opti_DG_E(z, h, :BTMPV)
                ) - 
                # behind the meter generation for pv+storage customers who did not participate in aggregation
                model_data.omega(d) * delta_t *
                (
                    min(
                        sum(customers.rho_DG(h, m, z, d, t) * customers.Opti_DG_E(z, h, m) for m in customers.index_m),
                        customers.d(h, z, d, t) * (1 - utility.loss_dist)
                    ) * customers.x_DG_E_my(first(model_data.index_y), h, z, :BTMStorage) /
                    customers.Opti_DG_E(z, h, :BTMStorage) + sum(
                        min(
                            sum(customers.rho_DG(h, m, z, d, t) *
                            customers.Opti_DG_my(Symbol(Int(y)), z, h, m) for m in customers.index_m),
                            customers.d(h, z, d, t) * (1 - utility.loss_dist)
                        ) * customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMStorage) /
                        customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMStorage) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    )
                ) * (1 - der_aggregator.aggregation_level(reg_year_index_dera, z)) -
                # green technology subscription
                sum(
                    model_data.omega(d) * delta_t * utility.rho_C_my(j, z, d, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, z, h) for y_symbol in
                    model_data.year(first(model_data.index_y_fix)):reg_year)
                    for j in model_data.index_j
                )
        end
        net_demand_h_t_wo_loss(z, h, tou, :) .= net_demand_h_t_wo_loss_temp
    end

    # for the purpose of calculating net peak load, use load including distribution loss
    # green tech offset is included here because these investments are paid for by corresponding load (treated in the same way as DPV)
    net_demand_for_peak_h_t = make_keyed_array(model_data.index_z, model_data.index_h, model_data.index_d, model_data.index_t)
    for z in model_data.index_z, h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
        net_demand_for_peak_h_t(z, h, d, t, :) .=
            sum(
                customers.gamma(z, h) * customers.d(h, z, d, t)
            ) -
            # DG
            sum(
                (
                    customers.rho_DG(h, m, z, d, t) *
                    customers.x_DG_E_my(first(model_data.index_y), h, z, m) + sum(
                        customers.rho_DG(h, m, z, d, t) * customers.x_DG_new_my(Symbol(Int(y)), h, z, m) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    )
                ) for m in customers.index_m
            ) + 
            # remove aggregated behind-the-meter storage/pv generation/consumption since they're front-of-the-meter now
            customers.rho_DG(h, :BTMStorage, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_stor_capacity(z, h) + 
            customers.rho_DG(h, :BTMPV, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_pv_capacity(z, h) -
            # green technology subscription
            sum(
                utility.rho_C_my(j, z, d, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, z, h) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):reg_year)
                for j in model_data.index_j
            )
    end

    # excluding green tech generation offset when it comes to sharing T&D costs
    net_demand_for_peak_wo_green_tech_h_t = make_keyed_array(model_data.index_z, model_data.index_h, model_data.index_d, model_data.index_t)
    for z in model_data.index_z, h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
        net_demand_for_peak_wo_green_tech_h_t(z, h, d, t, :) .=
            sum(
                customers.gamma(z, h) * customers.d(h, z, d, t)
            ) -
            # DG
            sum(
                (
                    customers.rho_DG(h, m, z, d, t) *
                    customers.x_DG_E_my(first(model_data.index_y), h, z, m) + sum(
                        customers.rho_DG(h, m, z, d, t) * customers.x_DG_new_my(Symbol(Int(y)), h, z, m) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    )
                ) for m in customers.index_m
            ) + 
            # remove aggregated behind-the-meter storage/pv generation/consumption since they're front-of-the-meter now
            customers.rho_DG(h, :BTMStorage, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_stor_capacity(z, h) + 
            customers.rho_DG(h, :BTMPV, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_pv_capacity(z, h)
    end

    # non-coincident peak load
    net_peak_load_h = make_keyed_array(model_data.index_z, model_data.index_h)
    for z in model_data.index_z, h in model_data.index_h
        net_peak_load_h(z, h, :) .=
            findmax(Dict((d, t) => net_demand_for_peak_h_t(z, h, d, t) for d in model_data.index_d, t in model_data.index_t))[1]
    end

    net_peak_load_wo_green_tech_h = make_keyed_array(model_data.index_z, model_data.index_h)
    for z in model_data.index_z, h in model_data.index_h
        net_peak_load_wo_green_tech_h(z, h, :) .=
            findmax(Dict((d, t) => net_demand_for_peak_wo_green_tech_h_t(z, h, d, t) for d in model_data.index_d, t in model_data.index_t))[1]
    end

    peak_load_edo_eximport = make_keyed_array(model_data.index_z, utility.index_l)
    for z in model_data.index_z, l in utility.index_l
        peak_load_edo_eximport(z, l, :) .=
            findmax(Dict((d, t) => utility.trans_topology(l, z) * utility.flow_my(reg_year_index, l, d, t) for d in model_data.index_d, t in model_data.index_t))[1]
    end

    # Cost Classification/Allocation
    energy_cost_allocation_h = make_keyed_array(model_data.index_z, model_data.index_h)
    for z in model_data.index_z, h in model_data.index_h
        energy_cost_allocation_h(z, h, :) .=
            (energy_cost(z) - net_eximport_cost(z)) * net_demand_h_w_loss(z, h) / net_demand_w_loss_no_eximport(z) + der_excess_cost_h(z, h)
    end
    replace!(energy_cost_allocation_h, NaN => 0.0)

    # allocate (revenue_requirement - energy_cost + net_eximport_cost) by net peak load with green tech generation offset;
    # allocate regulator.othercost (T&D costs) by net peak load without green tech generation offset;
    demand_cost_allocation_h = make_keyed_array(model_data.index_z, model_data.index_h)
    for z in model_data.index_z, h in model_data.index_h
        demand_cost_allocation_h(z, h, :) .=
            (revenue_requirement(z) - energy_cost(z) + net_eximport_cost(z)) * net_peak_load_h(z, h) / (
                sum(net_peak_load_h(z, h) for h in model_data.index_h)
            ) + 
            regulator.othercost(z, reg_year_index) * net_peak_load_wo_green_tech_h(z, h) / (
                sum(net_peak_load_wo_green_tech_h(z, h) for h in model_data.index_h)
            )
    end
    replace!(demand_cost_allocation_h, NaN => 0.0)

    demand_cost_allocation_capacity_h = make_keyed_array(model_data.index_z, model_data.index_h)
    for z in model_data.index_z, h in model_data.index_h
        demand_cost_allocation_capacity_h(z, h, :) .=
            (revenue_requirement(z) - energy_cost(z) + net_eximport_cost(z)) * net_peak_load_h(z, h) / (
                sum(net_peak_load_h(z, h) for h in model_data.index_h)
            )
    end
    replace!(demand_cost_allocation_capacity_h, NaN => 0.0)

    demand_cost_allocation_othercost_h = make_keyed_array(model_data.index_z, model_data.index_h)
    for z in model_data.index_z, h in model_data.index_h
        demand_cost_allocation_othercost_h(z, h, :) .=
            regulator.othercost(z, reg_year_index) * net_peak_load_wo_green_tech_h(z, h) / (
                sum(net_peak_load_wo_green_tech_h(z, h) for h in model_data.index_h)
            )
    end
    replace!(demand_cost_allocation_othercost_h, NaN => 0.0)

    energy_cost_allocation_h_t = make_keyed_array(model_data.index_z, model_data.index_h, regulator.index_rate_tou)
    for z in model_data.index_z, h in model_data.index_h, tou in regulator.index_rate_tou
        energy_cost_allocation_h_t(z, h, tou, :) .=
            (energy_cost_t(z, tou) - net_eximport_cost_t(z, tou)) * net_demand_h_t_w_loss(z, h, tou) / net_demand_t_w_loss_no_eximport(z, tou) +
            der_excess_cost_h_t(z, h, tou)
    end
    replace!(energy_cost_allocation_h_t, NaN => 0.0)

    # compute the retail price
    p_before = ParamArray(regulator.p, "p_before")
    fill!(p_before, NaN)  # TODO DT: debug only
    for z in model_data.index_z, h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
        p_before(z, h, d, t, :) .= regulator.p_my(reg_year_index, z, h, d, t)
    end

    # p_before_wavg = initialize_param("p_before_wavg", model_data.index_h)
    # for h in model_data.index_h
    #     p_before_wavg(h, :) .= 
    #         sum(regulator.p_my(reg_year_index, z, h, d, t) * model_data.omega(d) * delta_t * customers.d(h, z, d, t) for z in model_data.index_z, d in model_data.index_d, t in model_data.index_t) / 
    #         sum(model_data.omega(d) * delta_t * customers.d(h, z, d, t) for z in model_data.index_z, d in model_data.index_d, t in model_data.index_t)
    # end
    # replace!(p_before_wavg.values, NaN => 0.0)

    # TODO: Call a function instead of using if-then
    if regulator_opts.rate_design isa FlatRate
        fill!(regulator.p, NaN)
        for z in model_data.index_z, h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
            regulator.p(z, h, d, t, :) .=
                (energy_cost_allocation_h(z, h) + demand_cost_allocation_capacity_h(z, h)) /
                net_demand_h_wo_loss(z, h) +
                demand_cost_allocation_othercost_h(z, h) / net_demand_wo_green_tech_h_wo_loss(z, h)
        end
        replace!(regulator.p.values, NaN => 0.0)
    elseif regulator_opts.rate_design isa TOU
        fill!(regulator.p, NaN)
        for z in model_data.index_z, h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
            tou = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_d.==String(d)) .& (regulator.rep_day_time_tou_mapping.index_t.==String(t)), :index_rate_tou][1])
            regulator.p(z, h, d, t, :) .=
                energy_cost_allocation_h_t(z, h, tou) / net_demand_h_t_wo_loss(z, h, tou) +
                demand_cost_allocation_capacity_h(z, h) / net_demand_h_wo_loss(z, h) +
                demand_cost_allocation_othercost_h(z, h) / net_demand_wo_green_tech_h_wo_loss(z, h)
        end
        replace!(regulator.p.values, NaN => 0.0)
    end

    # fill!(regulator.p_regression, NaN)
    # for h in model_data.index_h
    #     regulator.p_regression(h, :) .=
    #         (energy_cost_allocation_h(h) + demand_cost_allocation_capacity_h(h)) /
    #         net_demand_h_wo_loss(h)
    # end

    # fill!(regulator.p_td, NaN)
    # for h in model_data.index_h
    #     regulator.p_td(h, :) .=
    #         demand_cost_allocation_othercost_h(h) / net_demand_wo_green_tech_h_wo_loss(h)
    # end

    # TODO: Call a function instead of using if-then
    if regulator_opts.net_metering_policy isa ExcessRetailRate
        regulator.p_ex = ParamArray(regulator.p)
    elseif regulator_opts.net_metering_policy isa ExcessMarginalCost
        fill!(regulator.p_ex, NaN)
        for z in model_data.index_z, h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
            regulator.p_ex(z, h, d, t, :) .=
                utility.p_energy_cem_my(reg_year_index, z, d, t)
        end
    elseif regulator_opts.net_metering_policy isa ExcessZero
        fill!(regulator.p_ex, 0.0)
    end

    for z in model_data.index_z, h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
        regulator.p_my(reg_year_index, z, h, d, t, :) .= regulator.p(z, h, d, t)
        regulator.p_ex_my(reg_year_index, z, h, d, t, :) .= regulator.p_ex(z, h, d, t)
    end

    # for h in model_data.index_h
    #     regulator.p_my_regression(reg_year_index, h, :) .= regulator.p_regression(h)
    #     regulator.p_my_td(reg_year_index, h, :) .= regulator.p_td(h)
    # end

    # p_after_wavg = initialize_param("p_after_wavg", model_data.index_h)
    # for h in model_data.index_h
    #     p_after_wavg(h, :) .= 
    #         sum(regulator.p_my(reg_year_index, z, h, d, t) * model_data.omega(d) * delta_t * customers.d(h, z, d, t) * customers.gamma(z, h) for z in model_data.index_z, d in model_data.index_d, t in model_data.index_t) / 
    #         sum(model_data.omega(d) * delta_t * customers.d(h, z, d, t) * customers.gamma(z, h) for z in model_data.index_z, d in model_data.index_d, t in model_data.index_t)
    # end
    # replace!(p_after_wavg.values, NaN => 0.0)

    # @info "Original retail price" p_before
    # @info "Original DER excess rate" p_ex_before
    # @info "New retail price" regulator.p
    # @info "New DER excess rate" regulator.p_ex

    return compute_difference_percentage_maximum_one_norm([
        (p_before.values, regulator.p.values),
    ])
end

# function solve_agent_problem!(
#     regulator::Regulator,
#     regulator_opts::RegulatorOptions,
#     model_data::HEMData,
#     hem_opts::HEMOptions{WholesaleMarket},
#     agent_store::AgentStore,
#     w_iter,
# )

#     for y in model_data.index_y_fix
#         regulator.othercost(y, :) .= (
#             regulator.distribution_cost(y) + 
#             regulator.administration_cost(y) + 
#             regulator.transmission_cost(y) + 
#             regulator.interconnection_cost(y) + 
#             regulator.system_cost(y)
#         )
#     end

#     customers = get_agent(CustomerGroup, agent_store)
#     ipp = get_agent(IPPGroup, agent_store)
#     utility = get_agent(Utility, agent_store)
#     green_developer = get_agent(GreenDeveloper, agent_store)

#     # the year regulator is making a rate case
#     reg_year, reg_year_index = get_reg_year(model_data)

#     net_demand_w_loss = (
#         # demand
#         # when it comes to sharing the revenue requirement (cost), use load including distribution loss
#         # e.g. utility generation is 100 MW, 50 MW to serve load (including distribution loss),
#         #      50 MW for export. It makes sense to allocate the same cost for internal load and export.
#         sum(
#             customers.gamma(h) * model_data.omega(t) * customers.d(h, t) for
#             h in model_data.index_h, t in model_data.index_t
#         ) +
#         # export
#         sum(
#             model_data.omega(t) * ipp.eximport_my(reg_year_index, t) for
#             t in model_data.index_t
#         ) -
#         # DG
#         # when it comes to cost allocation, simply use total DER generation to offset total load;
#         # this does not consider enery offset at the household level, which is only considered when 
#         # calculating retail rates.
#         sum(
#             model_data.omega(t) * (
#                 customers.rho_DG(h, m, t) *
#                 customers.x_DG_E_my(first(model_data.index_y), h, m) + sum(
#                     customers.rho_DG(h, m, t) * customers.x_DG_new_my(Symbol(Int(y)), h, m) for
#                     y in model_data.year(first(model_data.index_y_fix)):reg_year
#                 )
#             ) for t in model_data.index_t, h in model_data.index_h,
#             m in customers.index_m
#         ) -
#         # green technology subscription
#         sum(
#             model_data.omega(t) * utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
#             model_data.year(first(model_data.index_y_fix)):reg_year)
#             for t in model_data.index_t, j in model_data.index_j, h in model_data.index_h
#         )
#     )

#     der_excess_cost_h = KeyedArray(
#         [
#             sum(
#                 model_data.omega(t) *
#                 regulator.p_ex(h, t) *
#                 (
#                     max(
#                         0,
#                         customers.rho_DG(h, m, t) * customers.Opti_DG_E(h, m) -
#                         customers.d(h, t) * (1 - utility.loss_dist),
#                     ) * customers.x_DG_E_my(first(model_data.index_y), h, m) /
#                     customers.Opti_DG_E(h, m) + sum(
#                         max(
#                             0,
#                             customers.rho_DG(h, m, t) *
#                             customers.Opti_DG_my(Symbol(Int(y)), h, m) -
#                             customers.d(h, t) * (1 - utility.loss_dist),
#                         ) * customers.x_DG_new_my(Symbol(Int(y)), h, m) /
#                         customers.Opti_DG_my(Symbol(Int(y)), h, m) for
#                         y in model_data.year(first(model_data.index_y_fix)):reg_year
#                     )
#                 ) for t in model_data.index_t, m in customers.index_m
#             ) for h in model_data.index_h
#         ];
#         [get_pair(model_data.index_h)]...
#     )

#     net_demand_h_w_loss = KeyedArray(
#         # demand
#         [
#             sum(
#                 customers.gamma(h) * model_data.omega(t) * customers.d(h, t) for
#                 t in model_data.index_t
#             ) -
#             # DG
#             sum(
#                 model_data.omega(t) * (
#                     customers.rho_DG(h, m, t) *
#                     customers.x_DG_E_my(first(model_data.index_y), h, m) + sum(
#                         customers.rho_DG(h, m, t) *
#                         customers.x_DG_new_my(Symbol(Int(y)), h, m) for
#                         y in model_data.year(first(model_data.index_y_fix)):reg_year
#                     )
#                 ) for t in model_data.index_t, m in customers.index_m
#             ) -
#             # green technology subscription
#             sum(
#                 model_data.omega(t) * utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
#                 model_data.year(first(model_data.index_y_fix)):reg_year)
#                 for t in model_data.index_t, j in model_data.index_j
#             )
#             for h in model_data.index_h
#         ];
#         [get_pair(model_data.index_h)]...
#     )

#     net_demand_h_wo_loss = KeyedArray(
#         # demand
#         [
#             sum(
#                 customers.gamma(h) *
#                 model_data.omega(t) *
#                 customers.d(h, t) *
#                 (1 - utility.loss_dist) for t in model_data.index_t
#             ) -
#             # DG
#             # since this net demand is for rate calculation, we need to consider energy offset at the household level.
#             # e.g. two household, with 100 MW load each (without loss), if one of them has DER and generated 
#             # 120 MW of energy, he/she does not need to pay for energy, but the other one still have to pay 
#             # 100 MW instead of 80 MW.
#             sum(
#                 model_data.omega(t) * (
#                     min(
#                         customers.rho_DG(h, m, t) * customers.Opti_DG_E(h, m),
#                         customers.d(h, t) * (1 - utility.loss_dist),
#                     ) * customers.x_DG_E_my(first(model_data.index_y), h, m) /
#                     customers.Opti_DG_E(h, m) + sum(
#                         min(
#                             customers.rho_DG(h, m, t) *
#                             customers.Opti_DG_my(Symbol(Int(y)), h, m),
#                             customers.d(h, t) * (1 - utility.loss_dist),
#                         ) * customers.x_DG_new_my(Symbol(Int(y)), h, m) /
#                         customers.Opti_DG_my(Symbol(Int(y)), h, m) for
#                         y in model_data.year(first(model_data.index_y_fix)):reg_year
#                     )
#                 ) for t in model_data.index_t, m in customers.index_m
#             ) -
#             # green technology subscription
#             sum(
#                 model_data.omega(t) * utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
#                 model_data.year(first(model_data.index_y_fix)):reg_year)
#                 for t in model_data.index_t, j in model_data.index_j
#             )
#             for h in model_data.index_h
#         ];
#         [get_pair(model_data.index_h)]...
#     )

#     net_demand_wo_green_tech_h_wo_loss = KeyedArray(
#         [
#         # demand
#             sum(
#                 customers.gamma(h) *
#                 model_data.omega(t) *
#                 customers.d(h, t) *
#                 (1 - utility.loss_dist) for t in model_data.index_t
#             ) -
#             # DG
#             # since this net demand is for rate calculation, we need to consider energy offset at the household level.
#             # e.g. two household, with 100 MW load each (without loss), if one of them has DER and generated 
#             # 120 MW of energy, he/she does not need to pay for energy, but the other one still have to pay 
#             # 100 MW instead of 80 MW.
#             sum(
#                 model_data.omega(t) * (
#                     min(
#                         customers.rho_DG(h, m, t) * customers.Opti_DG_E(h, m),
#                         customers.d(h, t) * (1 - utility.loss_dist),
#                     ) * customers.x_DG_E_my(first(model_data.index_y), h, m) /
#                     customers.Opti_DG_E(h, m) + sum(
#                         min(
#                             customers.rho_DG(h, m, t) *
#                             customers.Opti_DG_my(Symbol(Int(y)), h, m),
#                             customers.d(h, t) * (1 - utility.loss_dist),
#                         ) * customers.x_DG_new_my(Symbol(Int(y)), h, m) /
#                         customers.Opti_DG_my(Symbol(Int(y)), h, m) for
#                         y in model_data.year(first(model_data.index_y_fix)):reg_year
#                     )
#                 ) for t in model_data.index_t, m in customers.index_m
#             )
#             for h in model_data.index_h
#         ];
#         [get_pair(model_data.index_h)]...
#     )

#     net_demand_t_w_loss = KeyedArray(
#         # demand
#         [
#             sum(customers.gamma(h) * customers.d(h, t) for h in model_data.index_h) +
#             # export
#             utility.eximport_my(reg_year_index, t) -
#             # DG
#             sum(
#                 customers.rho_DG(h, m, t) *
#                 customers.x_DG_E_my(first(model_data.index_y), h, m) + sum(
#                     customers.rho_DG(h, m, t) * customers.x_DG_new_my(Symbol(Int(y)), h, m) for
#                     y in model_data.year(first(model_data.index_y_fix)):reg_year
#                 ) for h in model_data.index_h, m in customers.index_m
#             ) -
#             # green technology subscription
#             sum(
#                 utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
#                 model_data.year(first(model_data.index_y_fix)):reg_year)
#                 for h in model_data.index_h, j in model_data.index_j
#             )
#             for t in model_data.index_t
#         ];
#         [get_pair(model_data.index_t)]...
#     )

#     der_excess_cost_h_t = make_keyed_array(model_data.index_h, model_data.index_t)
#     for h in model_data.index_h, t in model_data.index_t
#         der_excess_cost_h_t(h, t, :) .= sum(
#             regulator.p_ex(h, t) * (
#                 max(
#                     0,
#                     customers.rho_DG(h, m, t) * customers.Opti_DG_E(h, m) -
#                     customers.d(h, t) * (1 - utility.loss_dist),
#                 ) * customers.x_DG_E_my(first(model_data.index_y), h, m) /
#                 customers.Opti_DG_E(h, m) + sum(
#                     max(
#                         0,
#                         customers.rho_DG(h, m, t) *
#                         customers.Opti_DG_my(Symbol(Int(y)), h, m) -
#                         customers.d(h, t) * (1 - utility.loss_dist),
#                     ) * customers.x_DG_new_my(Symbol(Int(y)), h, m) /
#                     customers.Opti_DG_my(Symbol(Int(y)), h, m) for
#                     y in model_data.year(first(model_data.index_y_fix)):reg_year
#                 )
#             ) for m in customers.index_m
#         )
#     end

#     net_demand_h_t_w_loss = make_keyed_array(model_data.index_h, model_data.index_t)
#     for h in model_data.index_h, t in model_data.index_t
#         net_demand_h_t_w_loss(h, t, :) .=
#         # demand
#             customers.gamma(h) * customers.d(h, t) -
#             # DG
#             sum(
#                 customers.rho_DG(h, m, t) *
#                 customers.x_DG_E_my(first(model_data.index_y), h, m) + sum(
#                     customers.rho_DG(h, m, t) * customers.x_DG_new_my(Symbol(Int(y)), h, m)
#                     for y in model_data.year(first(model_data.index_y_fix)):reg_year
#                 ) for m in customers.index_m
#             ) -
#             sum(
#                 utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
#                 model_data.year(first(model_data.index_y_fix)):reg_year)
#                 for j in model_data.index_j
#             )
#     end

#     net_demand_h_t_wo_loss = make_keyed_array(model_data.index_h, model_data.index_t)
#     for h in model_data.index_h, t in model_data.index_t
#         net_demand_h_t_wo_loss(h, t, :) .=
#         # demand
#             customers.gamma(h) * customers.d(h, t) * (1 - utility.loss_dist) -
#             # DG
#             sum(
#                 min(
#                     customers.rho_DG(h, m, t) * customers.Opti_DG_E(h, m),
#                     customers.d(h, t) * (1 - utility.loss_dist),
#                 ) * customers.x_DG_E_my(first(model_data.index_y), h, m) /
#                 customers.Opti_DG_E(h, m) + sum(
#                     min(
#                         customers.rho_DG(h, m, t) *
#                         customers.Opti_DG_my(Symbol(Int(y)), h, m),
#                         customers.d(h, t) * (1 - utility.loss_dist),
#                     ) * customers.x_DG_new_my(Symbol(Int(y)), h, m) /
#                     customers.Opti_DG_my(Symbol(Int(y)), h, m) for
#                     y in model_data.year(first(model_data.index_y_fix)):reg_year
#                 ) for m in customers.index_m
#             ) -
#             sum(
#                 utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
#                 model_data.year(first(model_data.index_y_fix)):reg_year)
#                 for j in model_data.index_j
#             )
#     end

#     # for the purpose of calculating net peak load, use load including distribution loss
#     net_demand_for_peak_h_t = make_keyed_array(model_data.index_h, model_data.index_t)
#     for h in model_data.index_h, t in model_data.index_t
#         net_demand_for_peak_h_t(h, t, :) .=
#             customers.gamma(h) * customers.d(h, t) - sum(
#                 customers.rho_DG(h, m, t) *
#                 customers.x_DG_E_my(first(model_data.index_y), h, m) + sum(
#                     customers.rho_DG(h, m, t) * customers.x_DG_new_my(Symbol(Int(y)), h, m)
#                     for y in model_data.year(first(model_data.index_y_fix)):reg_year
#                 ) for m in customers.index_m
#             ) -
#             sum(
#                 utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
#                 model_data.year(first(model_data.index_y_fix)):reg_year)
#                 for j in model_data.index_j
#             )
#     end

#     # excluding green tech generation offset when it comes to sharing T&D costs
#     net_demand_for_peak_wo_green_tech_h_t = make_keyed_array(model_data.index_h, model_data.index_t)
#     for h in model_data.index_h, t in model_data.index_t
#         net_demand_for_peak_wo_green_tech_h_t(h, t, :) .=
#                 customers.gamma(h) * customers.d(h, t) - sum(
#                     customers.rho_DG(h, m, t) *
#                     customers.x_DG_E_my(first(model_data.index_y), h, m) + sum(
#                         customers.rho_DG(h, m, t) * customers.x_DG_new_my(Symbol(Int(y)), h, m) for
#                         y in model_data.year(first(model_data.index_y_fix)):reg_year
#                     ) for m in customers.index_m
#                 )
#     end

#     net_peak_load_h = KeyedArray(
#         [
#             findmax(Dict(t => net_demand_for_peak_h_t(h, t) for t in model_data.index_t))[1] for h in model_data.index_h
#         ];
#         [get_pair(model_data.index_h)]...
#     )

#     net_peak_load_wo_green_tech_h = KeyedArray(
#         [ 
#             findmax(Dict(t => net_demand_for_peak_wo_green_tech_h_t(h, t) for t in model_data.index_t))[1] for h in model_data.index_h
#         ];
#         [get_pair(model_data.index_h)]...
#     )

#     # Cost Classification/Allocation
#     energy_purchase_cost =
#         sum(
#             ipp.miu_my(reg_year_index, t) * (
#                 sum(
#                     ipp.y_E_my(reg_year_index, p, k, t) for k in ipp.index_k_existing,
#                     p in ipp.index_p
#                 ) + sum(
#                     ipp.y_C_my(reg_year_index, p, k, t) for k in ipp.index_k_new,
#                     p in ipp.index_p
#                 )
#             ) for t in model_data.index_t
#         ) +
#         # incorporate REC cost into energy purchase cost
#         sum(
#             regulator.REC *
#             model_data.omega(t) *
#             (
#                 sum(
#                     ipp.y_E_my(reg_year_index, p, rps, t) for rps in ipp.index_rps,
#                     p in ipp.index_p
#                 ) + sum(
#                     ipp.y_C_my(reg_year_index, p, rps, t) for rps in ipp.index_rps,
#                     p in ipp.index_p
#                 )
#             ) for t in model_data.index_t
#         )
#     energy_purchase_cost_t = KeyedArray(
#         [
#             ipp.miu_my(reg_year_index, t) / model_data.omega(t) * (
#                 sum(
#                     ipp.y_E_my(reg_year_index, p, k, t) for k in ipp.index_k_existing,
#                     p in ipp.index_p
#                 ) + sum(
#                     ipp.y_C_my(reg_year_index, p, k, t) for k in ipp.index_k_new,
#                     p in ipp.index_p
#                 )
#             ) +
#             # incorporate REC cost into energy purchase cost
#             regulator.REC * (
#                 sum(
#                     ipp.y_E_my(reg_year_index, p, rps, t) for rps in ipp.index_rps,
#                     p in ipp.index_p
#                 ) + sum(
#                     ipp.y_C_my(reg_year_index, p, rps, t) for rps in ipp.index_rps,
#                     p in ipp.index_p
#                 )
#             ) for t in model_data.index_t
#         ];
#         [get_pair(model_data.index_t)]...
#     )

#     capacity_purchase_cost = sum(
#         ipp.ucap(reg_year_index, p) * (
#             ipp.Capacity_intercept_my(reg_year_index) +
#             ipp.Capacity_slope_my(reg_year_index) *
#             sum(ipp.ucap(reg_year_index, p) for p in ipp.index_p)
#         ) for p in ipp.index_p
#     )

#     energy_cost_allocation_h = KeyedArray(
#         [
#             energy_purchase_cost * net_demand_h_w_loss(h) / net_demand_w_loss +
#             der_excess_cost_h(h) for h in model_data.index_h
#         ];
#         [get_pair(model_data.index_h)]...
#     )
#     energy_cost_allocation_eximport =
#         energy_purchase_cost * sum(
#             model_data.omega(t) * utility.eximport_my(reg_year_index, t) for
#             t in model_data.index_t
#         ) / net_demand_w_loss

#     # allocate capacity_purchase_cost by net peak load with green tech generation offset;
#     # allocate regulator.othercost (T&D costs) by net peak load without green tech generation offset;
#     demand_cost_allocation_h = KeyedArray(
#         [
#             capacity_purchase_cost * net_peak_load_h(h) / (
#                 sum(net_peak_load_h(h) for h in model_data.index_h) +
#                 utility.Peak_eximport_my(reg_year_index)
#             ) +
#             regulator.othercost(reg_year_index) * net_peak_load_wo_green_tech_h(h) / (
#                 sum(net_peak_load_wo_green_tech_h(h) for h in model_data.index_h) +
#                 utility.Peak_eximport_my(reg_year_index)
#             )
#             for h in model_data.index_h
#         ];
#         [get_pair(model_data.index_h)]...
#     )
#     demand_cost_allocation_capacity_h = KeyedArray(
#         [
#             capacity_purchase_cost * net_peak_load_h(h) / (
#                 sum(net_peak_load_h(h) for h in model_data.index_h) +
#                 utility.Peak_eximport_my(reg_year_index)
#             )
#             for h in model_data.index_h
#         ];
#         [get_pair(model_data.index_h)]...
#     )
#     demand_cost_allocation_othercost_h = KeyedArray(
#         [
#             regulator.othercost(reg_year_index) * net_peak_load_wo_green_tech_h(h) / (
#                 sum(net_peak_load_wo_green_tech_h(h) for h in model_data.index_h) +
#                 utility.Peak_eximport_my(reg_year_index)
#             )
#             for h in model_data.index_h
#         ];
#         [get_pair(model_data.index_h)]...
#     )

#     demand_cost_allocation_eximport =
#         capacity_purchase_cost *
#         utility.Peak_eximport_my(reg_year_index) / (
#             sum(net_peak_load_h(h) for h in model_data.index_h) +
#             utility.Peak_eximport_my(reg_year_index)
#         ) +
#         regulator.othercost(reg_year_index) *
#         utility.Peak_eximport_my(reg_year_index) / (
#             sum(net_peak_load_wo_green_tech_h(h) for h in model_data.index_h) +
#             utility.Peak_eximport_my(reg_year_index)
#         )
#     energy_cost_allocation_h_t = make_keyed_array(model_data.index_h, model_data.index_t)
#     for h in model_data.index_h, t in model_data.index_t
#         energy_cost_allocation_h_t(h, t, :) .=
#             energy_purchase_cost_t(t) * net_demand_h_t_w_loss(h, t) /
#             net_demand_t_w_loss(t) + der_excess_cost_h_t(h, t)
#     end

#     energy_cost_allocation_eximport_t = KeyedArray(
#         [
#             energy_purchase_cost_t(t) * utility.eximport_my(reg_year_index, t) /
#             net_demand_t_w_loss(t) for t in model_data.index_t
#         ];
#         [get_pair(model_data.index_t)]...
#     )

#     p_before = ParamArray(regulator.p, "p_before")
#     fill!(p_before, NaN)
#     for h in model_data.index_h, t in model_data.index_t
#         p_before(h, t, :) .= regulator.p_my(reg_year_index, h, t)
#     end

#     p_before_wavg = ParamArray(regulator.p_td, "p_before_wavg")
#     fill!(p_before_wavg, NaN)  # TODO DT: debug only
#     for h in model_data.index_h
#         p_before_wavg(h, :) .= 
#             sum(regulator.p_my(reg_year_index, h, t) * model_data.omega(t) * customers.d(h, t) for t in model_data.index_t) / 
#             sum(model_data.omega(t) * customers.d(h, t) for t in model_data.index_t)
#     end

#     p_ex_before = ParamArray(regulator.p_ex, "p_ex_before")
#     fill!(p_ex_before, NaN)
#     for h in model_data.index_h, t in model_data.index_t
#         p_ex_before(h, t, :) .= regulator.p_ex_my(reg_year_index, h, t)
#     end

#     # TODO: Call a function instead of using if-then
#     if regulator_opts.rate_design isa FlatRate
#         fill!(regulator.p, NaN)
#         for h in model_data.index_h, t in model_data.index_t
#             regulator.p(h, t, :) .=
#                 (energy_cost_allocation_h(h) + demand_cost_allocation_capacity_h(h)) /
#                 net_demand_h_wo_loss(h) +
#                 demand_cost_allocation_othercost_h(h) / net_demand_wo_green_tech_h_wo_loss(h)
#         end
#         fill!(regulator.p_eximport, NaN)
#         for t in model_data.index_t
#             regulator.p_eximport(t, :) .=
#                 (energy_cost_allocation_eximport + demand_cost_allocation_eximport) /
#                 sum(
#                     model_data.omega(t) * utility.eximport_my(reg_year_index, t) for
#                     t in model_data.index_t
#                 )
#         end
#     elseif regulator_opts.rate_design isa TOU
#         fill!(regulator.p, NaN)
#         for h in model_data.index_h, t in model_data.index_t
#             regulator.p(h, t, :) .=
#                 energy_cost_allocation_h_t(h, t) / net_demand_h_t_wo_loss(h, t) +
#                 demand_cost_allocation_capacity_h(h) / net_demand_h_wo_loss(h) +
#                 demand_cost_allocation_othercost_h(h) / net_demand_wo_green_tech_h_wo_loss(h)
#         end
#         fill!(regulator.p_eximport, NaN)
#         for t in model_data.index_t
#             regulator.p_eximport(t, :) .=
#                 energy_cost_allocation_eximport_t(t) /
#                 utility.eximport_my(reg_year_index, t) +
#                 demand_cost_allocation_eximport / sum(
#                     model_data.omega(t) * utility.eximport_my(reg_year_index, t) for
#                     t in model_data.index_t
#                 )
#         end
#     end

#     fill!(regulator.p_regression, NaN)
#     for h in model_data.index_h
#         regulator.p_regression(h, :) .=
#             (energy_cost_allocation_h(h) + demand_cost_allocation_capacity_h(h)) /
#             net_demand_h_wo_loss(h)
#     end

#     fill!(regulator.p_td, NaN)
#     for h in model_data.index_h
#         regulator.p_td(h, :) .=
#             demand_cost_allocation_othercost_h(h) / net_demand_wo_green_tech_h_wo_loss(h)
#     end

#     # TODO: Call a function instead of using if-then
#     if regulator_opts.net_metering_policy isa ExcessRetailRate
#         regulator.p_ex = ParamArray(regulator.p)
#     elseif regulator_opts.net_metering_policy isa ExcessMarginalCost
#         fill!(regulator.p_ex, NaN)
#         for h in model_data.index_h, t in model_data.index_t
#             regulator.p_ex(h, t, :) .= ipp.miu_my(reg_year_index, t) / model_data.omega(t)
#         end
#     elseif regulator_opts.net_metering_policy isa ExcessZero
#         fill!(regulator.p_ex, 0.0)
#     end

#     for h in model_data.index_h, t in model_data.index_t
#         regulator.p_my(reg_year_index, h, t, :) .= regulator.p(h, t)
#         regulator.p_ex_my(reg_year_index, h, t, :) .= regulator.p_ex(h, t)
#         regulator.p_eximport_my(reg_year_index, t, :) .= regulator.p_eximport(t)
#     end

#     for h in model_data.index_h
#         regulator.p_my_regression(reg_year_index, h, :) .= regulator.p_regression(h)
#         regulator.p_my_td(reg_year_index, h, :) .= regulator.p_td(h)
#     end

#     p_after_wavg = ParamArray(regulator.p_td, "p_after_wavg")
#     fill!(p_after_wavg, NaN)  # TODO DT: debug only
#     for h in model_data.index_h
#         p_after_wavg(h, :) .= 
#             sum(regulator.p_my(reg_year_index, h, t) * model_data.omega(t) * customers.d(h, t) for t in model_data.index_t) / 
#             sum(model_data.omega(t) * customers.d(h, t) for t in model_data.index_t)
#     end

#     @info "Original retail price" p_before
#     @info "Original DER excess rate" p_ex_before
#     @info "New retail price" regulator.p
#     @info "New DER excess rate" regulator.p_ex

#     return compute_difference_percentage_one_norm([
#         (p_before_wavg.values, p_after_wavg.values),
#     ])
# end

# regulator's module with updated dimensions

function solve_agent_problem!(
    regulator::Regulator,
    regulator_opts::RegulatorOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{WholesaleMarket},
    agent_store::AgentStore,
    w_iter,
    jump_model,
    export_file_path,
    update_results::Bool
)

    for y in model_data.index_y_fix
        for z in model_data.index_z
            regulator.othercost(z, y, :) .= regulator.distribution_cost(z, y) + regulator.administration_cost(z, y) + regulator.transmission_cost(z, y) + regulator.interconnection_cost(z, y) + regulator.system_cost(z, y)
        end
    end

    delta_t = get_delta_t(model_data)

    customers = get_agent(CustomerGroup, agent_store)
    ipp = get_agent(IPPGroup, agent_store)
    utility = get_agent(Utility, agent_store)
    green_developer = get_agent(GreenDeveloper, agent_store)
    der_aggregator = get_agent(DERAggregator, agent_store)

    # the year regulator is making a rate case
    reg_year, reg_year_index = get_reg_year(model_data)

    for h in model_data.index_h, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
        customers.d(h, z, d, t, :) .= customers.d_my(reg_year_index, h, z, d, t)
        # customers.DERGen(h, t, :) .= customers.DERGen_my(reg_year_index, h, t)
    end

    # since regulator problem is ahead of DERAggregator probelm, use previous year's aggregation results.
    reg_year_dera, reg_year_index_dera = get_prev_reg_year(model_data, w_iter)

    total_der_stor_capacity = make_keyed_array(model_data.index_z, model_data.index_h)
    # this total_der_pv_capacity is the approximate capacity of pv portion of pv+storage tech, not all dpv capacity
    total_der_pv_capacity = make_keyed_array(model_data.index_z, model_data.index_h)
    for z in model_data.index_z, h in model_data.index_h
        if w_iter >= 2
            total_der_stor_capacity(z, h, :) .=
                customers.x_DG_E_my(reg_year_index_dera, h, z, :BTMStorage) + sum(
                    customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMStorage) for
                    y in model_data.year(first(model_data.index_y_fix)):reg_year_dera
                )
        else
            total_der_stor_capacity(z, h, :) .= customers.x_DG_E_my(reg_year_index_dera, h, z, :BTMStorage)
        end
        total_der_pv_capacity(z, h, :) .= total_der_stor_capacity(z, h) / customers.Opti_DG_E(z, h, :BTMStorage) * customers.Opti_DG_E(z, h, :BTMPV)
    end

    net_demand_w_loss = make_keyed_array(model_data.index_z)
    for z in model_data.index_z
        net_demand_w_loss(z, :) .= 
            # demand
            # when it comes to sharing the revenue requirement (cost), use load including distribution loss
            # e.g. utility generation is 100 MW, 50 MW to serve load (including distribution loss),
            #      50 MW for export. It makes sense to allocate the same cost for internal load and export.
            sum(
                customers.gamma(z, h) * model_data.omega(d) * delta_t * customers.d(h, z, d, t) for
                h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
            ) +
            # exogenous export/import
            sum(
                model_data.omega(d) * delta_t * ipp.eximport_my(reg_year_index, z, d, t) for
                d in model_data.index_d, t in model_data.index_t
            ) +
            # endogenous export/import (flow out of zone z)
            sum(
                model_data.omega(d) * delta_t * ipp.trans_topology(l, z) * ipp.flow_my(reg_year_index, l, d, t) for
                l in utility.index_l, d in model_data.index_d, t in model_data.index_t
            ) -
            # DG
            # when it comes to cost allocation, simply use total DER generation to offset total load;
            # this does not consider enery offset at the household level, which is only considered when 
            # calculating retail rates (e.g., DER excess compensation).
            sum(
                model_data.omega(d) * delta_t * (
                    customers.rho_DG(h, m, z, d, t) *
                    customers.x_DG_E_my(first(model_data.index_y), h, z, m) + sum(
                        customers.rho_DG(h, m, z, d, t) * customers.x_DG_new_my(Symbol(Int(y)), h, z, m) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    )
                ) for h in model_data.index_h, m in customers.index_m, d in model_data.index_d, t in model_data.index_t
            ) + 
            # remove aggregated behind-the-meter storage/pv generation/consumption since they're front-of-the-meter now
            sum(
                model_data.omega(d) * delta_t * customers.rho_DG(h, :BTMStorage, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_stor_capacity(z, h) 
                for h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
            ) + 
            sum(
                model_data.omega(d) * delta_t * customers.rho_DG(h, :BTMPV, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_pv_capacity(z, h) 
                for h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
            ) -
            # green technology subscription
            sum(
                model_data.omega(d) * delta_t * utility.rho_C_my(j, z, d, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, z, h) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):reg_year)
                for j in model_data.index_j, h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
            )
    end

    net_demand_w_loss_no_eximport = make_keyed_array(model_data.index_z)
    for z in model_data.index_z
        net_demand_w_loss_no_eximport(z, :) .= 
            # demand
            # when it comes to sharing the revenue requirement (cost), use load including distribution loss
            # e.g. utility generation is 100 MW, 50 MW to serve load (including distribution loss),
            #      50 MW for export. It makes sense to allocate the same cost for internal load and export.
            sum(
                customers.gamma(z, h) * model_data.omega(d) * delta_t * customers.d(h, z, d, t) for
                h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
            ) -
            # DG
            # when it comes to cost allocation, simply use total DER generation to offset total load;
            # this does not consider enery offset at the household level, which is only considered when 
            # calculating retail rates (e.g., DER excess compensation).
            sum(
                model_data.omega(d) * delta_t * (
                    customers.rho_DG(h, m, z, d, t) *
                    customers.x_DG_E_my(first(model_data.index_y), h, z, m) + sum(
                        customers.rho_DG(h, m, z, d, t) * customers.x_DG_new_my(Symbol(Int(y)), h, z, m) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    )
                ) for h in model_data.index_h, m in customers.index_m, d in model_data.index_d, t in model_data.index_t
            ) + 
            # remove aggregated behind-the-meter storage generation/consumption since they're front-of-the-meter now
            sum(
                model_data.omega(d) * delta_t * customers.rho_DG(h, :BTMStorage, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_stor_capacity(z, h) 
                for h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
            ) + 
            sum(
                model_data.omega(d) * delta_t * customers.rho_DG(h, :BTMPV, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_pv_capacity(z, h) 
                for h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
            ) -
            # green technology subscription
            sum(
                model_data.omega(d) * delta_t * utility.rho_C_my(j, z, d, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, z, h) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):reg_year)
                for j in model_data.index_j, h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
            )
    end

    # In the presence of distribution loss, multiply load (including distribution loss) by (1 - loss factor),
    # DER generation above this value shall be compansated for DER excess credits
    # e.g. DER generation is 100 MW, load (without loss) is 95 MW, receive 5 MW excess credits
    # this formula is adapted for BTM PV+Storage only, with sum(customers.rho_DG(h, m, z, d, t) * customers.Opti_DG_E(z, h, m) for m in customers.index_m)
    # and sum(customers.rho_DG(h, m, z, d, t) * customers.Opti_DG_my(Symbol(Int(y)), z, h, m) for m in customers.index_m) indicating generation
    # from both BTM PV and BTM Storage, and customers.x_DG_E_my(first(model_data.index_y), h, z, :BTMPV) / customers.Opti_DG_E(z, h, :BTMPV) and 
    # customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMPV) / customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMPV) denoting how many households are there.
    # need to think more broadly about number of households when adopting multiple DERs.
    der_excess_cost_h = make_keyed_array(model_data.index_z, model_data.index_h)
    for z in model_data.index_z, h in model_data.index_h
        der_excess_cost_h(z, h, :) .= 
        # der excess for pv-only customers
        sum(
            model_data.omega(d) * delta_t *
            regulator.p_ex(z, h, d, t) *
            (
                max(
                    0,
                    customers.rho_DG(h, :BTMPV, z, d, t) * customers.Opti_DG_E(z, h, :BTMPV) -
                    customers.d(h, z, d, t) * (1 - utility.loss_dist),
                ) * customers.x_DG_E_my(first(model_data.index_y), h, z, :BTMPV) /
                customers.Opti_DG_E(z, h, :BTMPV) + sum(
                    max(
                        0,
                        customers.rho_DG(h, :BTMPV, z, d, t) *
                        customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMPV) -
                        customers.d(h, z, d, t) * (1 - utility.loss_dist),
                    ) * customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMPV) /
                    customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMPV) for
                    y in model_data.year(first(model_data.index_y_fix)):reg_year
                ) - 
                # minus pv portion of pv_storage tech
                max(
                    0,
                    customers.rho_DG(h, :BTMPV, z, d, t) * customers.Opti_DG_E(z, h, :BTMPV) -
                    customers.d(h, z, d, t) * (1 - utility.loss_dist),
                ) * total_der_pv_capacity(z, h) /
                customers.Opti_DG_E(z, h, :BTMPV)
            )
            for d in model_data.index_d, t in model_data.index_t
        ) + 
        # der excess for pv+storage customers who did not participate in aggregation
        sum(
            model_data.omega(d) * delta_t *
            regulator.p_ex(z, h, d, t) *
            (
                max(
                    0,
                    sum(customers.rho_DG(h, m, z, d, t) * customers.Opti_DG_E(z, h, m) for m in customers.index_m) -
                    customers.d(h, z, d, t) * (1 - utility.loss_dist),
                ) * customers.x_DG_E_my(first(model_data.index_y), h, z, :BTMStorage) /
                customers.Opti_DG_E(z, h, :BTMStorage) + sum(
                    max(
                        0,
                        sum(customers.rho_DG(h, m, z, d, t) *
                        customers.Opti_DG_my(Symbol(Int(y)), z, h, m) for m in customers.index_m) -
                        customers.d(h, z, d, t) * (1 - utility.loss_dist),
                    ) * customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMStorage) /
                    customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMStorage) for
                    y in model_data.year(first(model_data.index_y_fix)):reg_year
                )
            ) * (1 - der_aggregator.aggregation_level(reg_year_index_dera, z))
            for d in model_data.index_d, t in model_data.index_t
        )
    end

    net_demand_h_w_loss = make_keyed_array(model_data.index_z, model_data.index_h)
    for z in model_data.index_z, h in model_data.index_h
        # Demand
        net_demand_h_w_loss(z, h, :) .= sum(
                customers.gamma(z, h) * model_data.omega(d) * delta_t * customers.d(h, z, d, t) for
                d in model_data.index_d, t in model_data.index_t
            ) -
            # DG
            sum(
                model_data.omega(d) * delta_t * (
                    customers.rho_DG(h, m, z, d, t) *
                    customers.x_DG_E_my(first(model_data.index_y), h, z, m) + sum(
                        customers.rho_DG(h, m, z, d, t) * customers.x_DG_new_my(Symbol(Int(y)), h, z, m) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    )
                ) for m in customers.index_m, d in model_data.index_d, t in model_data.index_t
            ) + 
            # remove aggregated behind-the-meter storage/pv generation/consumption since they're front-of-the-meter now
            sum(model_data.omega(d) * delta_t * customers.rho_DG(h, :BTMStorage, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_stor_capacity(z, h) for d in model_data.index_d, t in model_data.index_t) + 
            sum(model_data.omega(d) * delta_t * customers.rho_DG(h, :BTMPV, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_pv_capacity(z, h) for d in model_data.index_d, t in model_data.index_t) -
            # green technology subscription
            sum(
                model_data.omega(d) * delta_t * utility.rho_C_my(j, z, d, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, z, h) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):reg_year)
                for j in model_data.index_j, d in model_data.index_d, t in model_data.index_t
            )
    end

    net_demand_h_wo_loss = make_keyed_array(model_data.index_z, model_data.index_h)
    for z in model_data.index_z, h in model_data.index_h
        net_demand_h_wo_loss(z, h, :) .= 
            # Demand without loss
            sum(
                customers.gamma(z, h) * model_data.omega(d) * delta_t * customers.d(h, z, d, t) * (1 - utility.loss_dist) for
                d in model_data.index_d, t in model_data.index_t
            ) -
            # DG
            # since this net demand is for rate calculation, we need to consider energy offset at the household level.
            # e.g. two household, with 100 MWh load each (without loss), if one of them has DER and generated 
            # 120 MWh of energy, he/she does not need to pay for energy, but the other one still have to pay 
            # 100 MWh instead of 80 MWh.
            
            # behind the meter generation for pv-only customers
            sum(
                model_data.omega(d) * delta_t *
                (
                    min(
                        customers.rho_DG(h, :BTMPV, z, d, t) * customers.Opti_DG_E(z, h, :BTMPV),
                        customers.d(h, z, d, t) * (1 - utility.loss_dist),
                    ) * customers.x_DG_E_my(first(model_data.index_y), h, z, :BTMPV) /
                    customers.Opti_DG_E(z, h, :BTMPV) + sum(
                        min(
                            customers.rho_DG(h, :BTMPV, z, d, t) *
                            customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMPV),
                            customers.d(h, z, d, t) * (1 - utility.loss_dist)
                        ) * customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMPV) /
                        customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMPV) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    ) - 
                    # minus pv portion of pv_storage tech
                    min(
                        customers.rho_DG(h, :BTMPV, z, d, t) * customers.Opti_DG_E(z, h, :BTMPV),
                        customers.d(h, z, d, t) * (1 - utility.loss_dist)
                    ) * total_der_pv_capacity(z, h) /
                    customers.Opti_DG_E(z, h, :BTMPV)
                )
                for d in model_data.index_d, t in model_data.index_t
            ) - 
            # behind the meter generation for pv+storage customers who did not participate in aggregation
            sum(
                model_data.omega(d) * delta_t *
                (
                    min(
                        sum(customers.rho_DG(h, m, z, d, t) * customers.Opti_DG_E(z, h, m) for m in customers.index_m),
                        customers.d(h, z, d, t) * (1 - utility.loss_dist)
                    ) * customers.x_DG_E_my(first(model_data.index_y), h, z, :BTMStorage) /
                    customers.Opti_DG_E(z, h, :BTMStorage) + sum(
                        min(
                            sum(customers.rho_DG(h, m, z, d, t) *
                            customers.Opti_DG_my(Symbol(Int(y)), z, h, m) for m in customers.index_m),
                            customers.d(h, z, d, t) * (1 - utility.loss_dist)
                        ) * customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMStorage) /
                        customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMStorage) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    )
                ) * (1 - der_aggregator.aggregation_level(reg_year_index_dera, z))
                for d in model_data.index_d, t in model_data.index_t
            ) - 
            # sum(
            #     model_data.omega(d) * delta_t * (
            #         min(
            #             sum(customers.rho_DG(h, m, z, d, t) * customers.Opti_DG_E(z, h, m) for m in customers.index_m),
            #             customers.d(h, z, d, t) * (1 - utility.loss_dist),
            #         ) * customers.x_DG_E_my(first(model_data.index_y), h, z, :BTMPV) /
            #         customers.Opti_DG_E(z, h, :BTMPV) + sum(
            #             min(
            #                 sum(customers.rho_DG(h, m, z, d, t) *
            #                 customers.Opti_DG_my(Symbol(Int(y)), z, h, m) for m in customers.index_m),
            #                 customers.d(h, z, d, t) * (1 - utility.loss_dist),
            #             ) * customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMPV) /
            #             customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMPV) for
            #             y in model_data.year(first(model_data.index_y_fix)):reg_year
            #         )
            #     ) for d in model_data.index_d, t in model_data.index_t
            # ) -
            # green technology subscription
            sum(
                model_data.omega(d) * delta_t * utility.rho_C_my(j, z, d, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, z, h) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):reg_year)
                for j in model_data.index_j, d in model_data.index_d, t in model_data.index_t
            )
    end

    net_demand_wo_green_tech_h_wo_loss = make_keyed_array(model_data.index_z, model_data.index_h)
    for z in model_data.index_z, h in model_data.index_h
        net_demand_wo_green_tech_h_wo_loss(z, h, :) .= 
            # Demand with loss
            sum(
                customers.gamma(z, h) * model_data.omega(d) * delta_t * customers.d(h, z, d, t) * (1 - utility.loss_dist) for
                d in model_data.index_d, t in model_data.index_t
            ) -
            # DG
            # since this net demand is for rate calculation, we need to consider energy offset at the household level.
            # e.g. two household, with 100 MW load each (without loss), if one of them has DER and generated 
            # 120 MW of energy, he/she does not need to pay for energy, but the other one still have to pay 
            # 100 MW instead of 80 MW.
            # behind the meter generation for pv-only customers
            sum(
                model_data.omega(d) * delta_t *
                (
                    min(
                        customers.rho_DG(h, :BTMPV, z, d, t) * customers.Opti_DG_E(z, h, :BTMPV),
                        customers.d(h, z, d, t) * (1 - utility.loss_dist),
                    ) * customers.x_DG_E_my(first(model_data.index_y), h, z, :BTMPV) /
                    customers.Opti_DG_E(z, h, :BTMPV) + sum(
                        min(
                            customers.rho_DG(h, :BTMPV, z, d, t) *
                            customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMPV),
                            customers.d(h, z, d, t) * (1 - utility.loss_dist)
                        ) * customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMPV) /
                        customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMPV) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    ) - 
                    # minus pv portion of pv_storage tech
                    min(
                        customers.rho_DG(h, :BTMPV, z, d, t) * customers.Opti_DG_E(z, h, :BTMPV),
                        customers.d(h, z, d, t) * (1 - utility.loss_dist)
                    ) * total_der_pv_capacity(z, h) /
                    customers.Opti_DG_E(z, h, :BTMPV)
                )
                for d in model_data.index_d, t in model_data.index_t
            ) - 
            # behind the meter generation for pv+storage customers who did not participate in aggregation
            sum(
                model_data.omega(d) * delta_t *
                (
                    min(
                        sum(customers.rho_DG(h, m, z, d, t) * customers.Opti_DG_E(z, h, m) for m in customers.index_m),
                        customers.d(h, z, d, t) * (1 - utility.loss_dist)
                    ) * customers.x_DG_E_my(first(model_data.index_y), h, z, :BTMStorage) /
                    customers.Opti_DG_E(z, h, :BTMStorage) + sum(
                        min(
                            sum(customers.rho_DG(h, m, z, d, t) *
                            customers.Opti_DG_my(Symbol(Int(y)), z, h, m) for m in customers.index_m),
                            customers.d(h, z, d, t) * (1 - utility.loss_dist)
                        ) * customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMStorage) /
                        customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMStorage) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    )
                ) * (1 - der_aggregator.aggregation_level(reg_year_index_dera, z))
                for d in model_data.index_d, t in model_data.index_t
            )
            # sum(
            #     model_data.omega(d) * delta_t * (
            #         min(
            #             sum(customers.rho_DG(h, m, z, d, t) * customers.Opti_DG_E(z, h, m) for m in customers.index_m),
            #             customers.d(h, z, d, t) * (1 - utility.loss_dist),
            #         ) * customers.x_DG_E_my(first(model_data.index_y), h, z, :BTMPV) /
            #         customers.Opti_DG_E(z, h, :BTMPV) + sum(
            #             min(
            #                 sum(customers.rho_DG(h, m, z, d, t) *
            #                 customers.Opti_DG_my(Symbol(Int(y)), z, h, m) for m in customers.index_m),
            #                 customers.d(h, z, d, t) * (1 - utility.loss_dist),
            #             ) * customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMPV) /
            #             customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMPV) for
            #             y in model_data.year(first(model_data.index_y_fix)):reg_year
            #         )
            #     ) for d in model_data.index_d, t in model_data.index_t
            # )
    end

    net_demand_t_w_loss = make_keyed_array(model_data.index_z, regulator.index_rate_tou)
    for z in model_data.index_z, tou in regulator.index_rate_tou
        net_demand_t_w_loss_temp = 0.0
        for i in 1:size(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :])[1]
            d = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_d][i])
            t = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_t][i])
            net_demand_t_w_loss_temp = net_demand_t_w_loss_temp + 
                # demand
                # when it comes to sharing the revenue requirement (cost), use load including distribution loss
                # e.g. utility generation is 100 MW, 50 MW to serve load (including distribution loss),
                #      50 MW for export. It makes sense to allocate the same cost for internal load and export.
                sum(
                    customers.gamma(z, h) * model_data.omega(d) * delta_t * customers.d(h, z, d, t) for
                    h in model_data.index_h
                ) +
                # exogenous export/import
                sum(
                    model_data.omega(d) * delta_t * ipp.eximport_my(reg_year_index, z, d, t)
                ) +
                # endogenous export/import (flow out of zone z)
                sum(
                    model_data.omega(d) * delta_t * ipp.trans_topology(l, z) * ipp.flow_my(reg_year_index, l, d, t) for
                    l in utility.index_l
                ) -
                # DG
                # when it comes to cost allocation, simply use total DER generation to offset total load;
                # this does not consider enery offset at the household level, which is only considered when 
                # calculating retail rates (e.g., DER excess compensation).
                sum(
                    model_data.omega(d) * delta_t * (
                        customers.rho_DG(h, m, z, d, t) *
                        customers.x_DG_E_my(first(model_data.index_y), h, z, m) + sum(
                            customers.rho_DG(h, m, z, d, t) * customers.x_DG_new_my(Symbol(Int(y)), h, z, m) for
                            y in model_data.year(first(model_data.index_y_fix)):reg_year
                        )
                    ) for h in model_data.index_h, m in customers.index_m
                ) +
                # remove aggregated behind-the-meter storage/pv generation/consumption since they're front-of-the-meter now
                sum(
                    model_data.omega(d) * delta_t * customers.rho_DG(h, :BTMStorage, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_stor_capacity(z, h) 
                    for h in model_data.index_h
                ) + 
                sum(
                    model_data.omega(d) * delta_t * customers.rho_DG(h, :BTMPV, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_pv_capacity(z, h) 
                    for h in model_data.index_h
                ) -
                # green technology subscription
                sum(
                    model_data.omega(d) * delta_t * utility.rho_C_my(j, z, d, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, z, h) for y_symbol in
                    model_data.year(first(model_data.index_y_fix)):reg_year)
                    for j in model_data.index_j, h in model_data.index_h
                )
        end
        net_demand_t_w_loss(z, tou, :) .= net_demand_t_w_loss_temp
    end

    net_demand_t_w_loss_no_eximport = make_keyed_array(model_data.index_z, regulator.index_rate_tou)
    for z in model_data.index_z, tou in regulator.index_rate_tou
        net_demand_t_w_loss_no_eximport_temp = 0.0
        for i in 1:size(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :])[1]
            d = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_d][i])
            t = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_t][i])
            net_demand_t_w_loss_no_eximport_temp = net_demand_t_w_loss_no_eximport_temp + 
                # demand
                # when it comes to sharing the revenue requirement (cost), use load including distribution loss
                # e.g. utility generation is 100 MW, 50 MW to serve load (including distribution loss),
                #      50 MW for export. It makes sense to allocate the same cost for internal load and export.
                sum(
                    customers.gamma(z, h) * model_data.omega(d) * delta_t * customers.d(h, z, d, t) for
                    h in model_data.index_h
                ) -
                # DG
                # when it comes to cost allocation, simply use total DER generation to offset total load;
                # this does not consider enery offset at the household level, which is only considered when 
                # calculating retail rates (e.g., DER excess compensation).
                sum(
                    model_data.omega(d) * delta_t * (
                        customers.rho_DG(h, m, z, d, t) *
                        customers.x_DG_E_my(first(model_data.index_y), h, z, m) + sum(
                            customers.rho_DG(h, m, z, d, t) * customers.x_DG_new_my(Symbol(Int(y)), h, z, m) for
                            y in model_data.year(first(model_data.index_y_fix)):reg_year
                        )
                    ) for h in model_data.index_h, m in customers.index_m
                ) + 
                # remove aggregated behind-the-meter storage/pv generation/consumption since they're front-of-the-meter now
                sum(
                    model_data.omega(d) * delta_t * customers.rho_DG(h, :BTMStorage, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_stor_capacity(z, h) 
                    for h in model_data.index_h
                ) + 
                sum(
                    model_data.omega(d) * delta_t * customers.rho_DG(h, :BTMPV, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_pv_capacity(z, h) 
                    for h in model_data.index_h
                ) -
                # green technology subscription
                sum(
                    model_data.omega(d) * delta_t * utility.rho_C_my(j, z, d, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, z, h) for y_symbol in
                    model_data.year(first(model_data.index_y_fix)):reg_year)
                    for j in model_data.index_j, h in model_data.index_h
                )
        end
        net_demand_t_w_loss_no_eximport(z, tou, :) .= net_demand_t_w_loss_no_eximport_temp
    end

    exo_eximport_demand_t = make_keyed_array(model_data.index_z, regulator.index_rate_tou)
    for z in model_data.index_z, tou in regulator.index_rate_tou
        exo_eximport_demand_t_temp = 0.0
        for i in 1:size(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :])[1]
            d = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_d][i])
            t = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_t][i])
            exo_eximport_demand_t_temp = exo_eximport_demand_t_temp + 
                # exogenous export/import
                sum(
                    model_data.omega(d) * delta_t * ipp.eximport_my(reg_year_index, z, d, t)
                )
        end
        exo_eximport_demand_t(z, tou, :) .= exo_eximport_demand_t_temp
    end

    edo_eximport_demand_t = make_keyed_array(model_data.index_z, utility.index_l, regulator.index_rate_tou)
    for z in model_data.index_z, l in utility.index_l, tou in regulator.index_rate_tou
        edo_eximport_demand_t_temp = 0.0
        for i in 1:size(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :])[1]
            d = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_d][i])
            t = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_t][i])
            edo_eximport_demand_t_temp = edo_eximport_demand_t_temp + 
                # endogenous export/import (flow out of zone z)
                sum(
                    model_data.omega(d) * delta_t * ipp.trans_topology(l, z) * ipp.flow_my(reg_year_index, l, d, t)
                )
        end
        edo_eximport_demand_t(z, l, tou, :) .= edo_eximport_demand_t_temp
    end

    der_excess_cost_h_t = make_keyed_array(model_data.index_z, model_data.index_h, regulator.index_rate_tou)
    for z in model_data.index_z, h in model_data.index_h, tou in regulator.index_rate_tou
        der_excess_cost_h_t_temp = 0.0
        for i in 1:size(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :])[1]
            d = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_d][i])
            t = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_t][i])
            der_excess_cost_h_t_temp = der_excess_cost_h_t_temp + 
            # der excess for pv-only customers
            model_data.omega(d) * delta_t *
            regulator.p_ex(z, h, d, t) *
            (
                max(
                    0,
                    customers.rho_DG(h, :BTMPV, z, d, t) * customers.Opti_DG_E(z, h, :BTMPV) -
                    customers.d(h, z, d, t) * (1 - utility.loss_dist),
                ) * customers.x_DG_E_my(first(model_data.index_y), h, z, :BTMPV) /
                customers.Opti_DG_E(z, h, :BTMPV) + sum(
                    max(
                        0,
                        customers.rho_DG(h, :BTMPV, z, d, t) *
                        customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMPV) -
                        customers.d(h, z, d, t) * (1 - utility.loss_dist),
                    ) * customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMPV) /
                    customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMPV) for
                    y in model_data.year(first(model_data.index_y_fix)):reg_year
                ) - 
                # minus pv portion of pv_storage tech
                max(
                    0,
                    customers.rho_DG(h, :BTMPV, z, d, t) * customers.Opti_DG_E(z, h, :BTMPV) -
                    customers.d(h, z, d, t) * (1 - utility.loss_dist),
                ) * total_der_pv_capacity(z, h) /
                customers.Opti_DG_E(z, h, :BTMPV)
            ) + 
            # der excess for pv+storage customers who did not participate in aggregation
            model_data.omega(d) * delta_t *
            regulator.p_ex(z, h, d, t) *
            (
                max(
                    0,
                    sum(customers.rho_DG(h, m, z, d, t) * customers.Opti_DG_E(z, h, m) for m in customers.index_m) -
                    customers.d(h, z, d, t) * (1 - utility.loss_dist),
                ) * customers.x_DG_E_my(first(model_data.index_y), h, z, :BTMStorage) /
                customers.Opti_DG_E(z, h, :BTMStorage) + sum(
                    max(
                        0,
                        sum(customers.rho_DG(h, m, z, d, t) *
                        customers.Opti_DG_my(Symbol(Int(y)), z, h, m) for m in customers.index_m) -
                        customers.d(h, z, d, t) * (1 - utility.loss_dist),
                    ) * customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMStorage) /
                    customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMStorage) for
                    y in model_data.year(first(model_data.index_y_fix)):reg_year
                )
            ) * (1 - der_aggregator.aggregation_level(reg_year_index_dera, z))
            # der_excess_cost_h_t_temp = der_excess_cost_h_t_temp + 
            #     model_data.omega(d) * delta_t *
            #     regulator.p_ex(z, h, d, t) *
            #     (
            #         max(
            #             0,
            #             sum(customers.rho_DG(h, m, z, d, t) * customers.Opti_DG_E(z, h, m) for m in customers.index_m) -
            #             customers.d(h, z, d, t) * (1 - utility.loss_dist),
            #         ) * customers.x_DG_E_my(first(model_data.index_y), h, z, :BTMPV) /
            #         customers.Opti_DG_E(z, h, :BTMPV) + sum(
            #             max(
            #                 0,
            #                 sum(customers.rho_DG(h, m, z, d, t) *
            #                 customers.Opti_DG_my(Symbol(Int(y)), z, h, m) for m in customers.index_m) -
            #                 customers.d(h, z, d, t) * (1 - utility.loss_dist),
            #             ) * customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMPV) /
            #             customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMPV) for
            #             y in model_data.year(first(model_data.index_y_fix)):reg_year
            #         )
            #     )
        end
        der_excess_cost_h_t(z, h, tou, :) .= der_excess_cost_h_t_temp
    end

    net_demand_h_t_w_loss = make_keyed_array(model_data.index_z, model_data.index_h, regulator.index_rate_tou)
    for z in model_data.index_z, h in model_data.index_h, tou in regulator.index_rate_tou
        net_demand_h_t_w_loss_temp = 0.0
        for i in 1:size(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :])[1]
            d = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_d][i])
            t = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_t][i])
            net_demand_h_t_w_loss_temp = net_demand_h_t_w_loss_temp +
                sum(
                    customers.gamma(z, h) * model_data.omega(d) * delta_t * customers.d(h, z, d, t)
                ) -
                # DG
                sum(
                    model_data.omega(d) * delta_t * (
                        customers.rho_DG(h, m, z, d, t) *
                        customers.x_DG_E_my(first(model_data.index_y), h, z, m) + sum(
                            customers.rho_DG(h, m, z, d, t) * customers.x_DG_new_my(Symbol(Int(y)), h, z, m) for
                            y in model_data.year(first(model_data.index_y_fix)):reg_year
                        )
                    ) for m in customers.index_m
                ) +
                # remove aggregated behind-the-meter storage/pv generation/consumption since they're front-of-the-meter now
                model_data.omega(d) * delta_t * customers.rho_DG(h, :BTMStorage, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_stor_capacity(z, h) + 
                model_data.omega(d) * delta_t * customers.rho_DG(h, :BTMPV, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_pv_capacity(z, h) -
                # green technology subscription
                sum(
                    model_data.omega(d) * delta_t * utility.rho_C_my(j, z, d, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, z, h) for y_symbol in
                    model_data.year(first(model_data.index_y_fix)):reg_year)
                    for j in model_data.index_j
                )
        end
        net_demand_h_t_w_loss(z, h, tou, :) .= net_demand_h_t_w_loss_temp
    end

    net_demand_h_t_wo_loss = make_keyed_array(model_data.index_z, model_data.index_h, regulator.index_rate_tou)
    for z in model_data.index_z, h in model_data.index_h, tou in regulator.index_rate_tou
        net_demand_h_t_wo_loss_temp = 0.0
        for i in 1:size(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :])[1]
            d = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_d][i])
            t = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_t][i])
            net_demand_h_t_wo_loss_temp = net_demand_h_t_wo_loss_temp +
                # Demand without loss
                sum(
                    customers.gamma(z, h) * model_data.omega(d) * delta_t * customers.d(h, z, d, t) * (1 - utility.loss_dist)
                ) -
                # DG
                # since this net demand is for rate calculation, we need to consider energy offset at the household level.
                # e.g. two household, with 100 MW load each (without loss), if one of them has DER and generated 
                # 120 MW of energy, he/she does not need to pay for energy, but the other one still have to pay 
                # 100 MW instead of 80 MW.

                # behind the meter generation for pv-only customers
                model_data.omega(d) * delta_t *
                (
                    min(
                        customers.rho_DG(h, :BTMPV, z, d, t) * customers.Opti_DG_E(z, h, :BTMPV),
                        customers.d(h, z, d, t) * (1 - utility.loss_dist),
                    ) * customers.x_DG_E_my(first(model_data.index_y), h, z, :BTMPV) /
                    customers.Opti_DG_E(z, h, :BTMPV) + sum(
                        min(
                            customers.rho_DG(h, :BTMPV, z, d, t) *
                            customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMPV),
                            customers.d(h, z, d, t) * (1 - utility.loss_dist)
                        ) * customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMPV) /
                        customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMPV) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    ) - 
                    # minus pv portion of pv_storage tech
                    min(
                        customers.rho_DG(h, :BTMPV, z, d, t) * customers.Opti_DG_E(z, h, :BTMPV),
                        customers.d(h, z, d, t) * (1 - utility.loss_dist)
                    ) * total_der_pv_capacity(z, h) /
                    customers.Opti_DG_E(z, h, :BTMPV)
                ) - 
                # behind the meter generation for pv+storage customers who did not participate in aggregation
                model_data.omega(d) * delta_t *
                (
                    min(
                        sum(customers.rho_DG(h, m, z, d, t) * customers.Opti_DG_E(z, h, m) for m in customers.index_m),
                        customers.d(h, z, d, t) * (1 - utility.loss_dist)
                    ) * customers.x_DG_E_my(first(model_data.index_y), h, z, :BTMStorage) /
                    customers.Opti_DG_E(z, h, :BTMStorage) + sum(
                        min(
                            sum(customers.rho_DG(h, m, z, d, t) *
                            customers.Opti_DG_my(Symbol(Int(y)), z, h, m) for m in customers.index_m),
                            customers.d(h, z, d, t) * (1 - utility.loss_dist)
                        ) * customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMStorage) /
                        customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMStorage) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    )
                ) * (1 - der_aggregator.aggregation_level(reg_year_index_dera, z)) -
                # model_data.omega(d) * delta_t * (
                #     min(
                #         sum(customers.rho_DG(h, m, z, d, t) * customers.Opti_DG_E(z, h, m) for m in customers.index_m),
                #         customers.d(h, z, d, t) * (1 - utility.loss_dist),
                #     ) * customers.x_DG_E_my(first(model_data.index_y), h, z, :BTMPV) /
                #     customers.Opti_DG_E(z, h, :BTMPV) + sum(
                #         min(
                #             sum(customers.rho_DG(h, m, z, d, t) *
                #             customers.Opti_DG_my(Symbol(Int(y)), z, h, m) for m in customers.index_m),
                #             customers.d(h, z, d, t) * (1 - utility.loss_dist),
                #         ) * customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMPV) /
                #         customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMPV) for
                #         y in model_data.year(first(model_data.index_y_fix)):reg_year
                #     )
                # ) -
                # green technology subscription
                sum(
                    model_data.omega(d) * delta_t * utility.rho_C_my(j, z, d, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, z, h) for y_symbol in
                    model_data.year(first(model_data.index_y_fix)):reg_year)
                    for j in model_data.index_j
                )
        end
        net_demand_h_t_wo_loss(z, h, tou, :) .= net_demand_h_t_wo_loss_temp
    end

    net_demand_wo_green_tech_h_t_wo_loss = make_keyed_array(model_data.index_z, model_data.index_h, regulator.index_rate_tou)
    for z in model_data.index_z, h in model_data.index_h, tou in regulator.index_rate_tou
        net_demand_wo_green_tech_h_t_wo_loss_temp = 0.0
        for i in 1:size(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :])[1]
            d = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_d][i])
            t = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_t][i])
            net_demand_wo_green_tech_h_t_wo_loss_temp = net_demand_wo_green_tech_h_t_wo_loss_temp +
                # Demand without loss
                sum(
                    customers.gamma(z, h) * model_data.omega(d) * delta_t * customers.d(h, z, d, t) * (1 - utility.loss_dist)
                ) -
                # DG
                # since this net demand is for rate calculation, we need to consider energy offset at the household level.
                # e.g. two household, with 100 MW load each (without loss), if one of them has DER and generated 
                # 120 MW of energy, he/she does not need to pay for energy, but the other one still have to pay 
                # 100 MW instead of 80 MW.

                # behind the meter generation for pv-only customers
                model_data.omega(d) * delta_t *
                (
                    min(
                        customers.rho_DG(h, :BTMPV, z, d, t) * customers.Opti_DG_E(z, h, :BTMPV),
                        customers.d(h, z, d, t) * (1 - utility.loss_dist),
                    ) * customers.x_DG_E_my(first(model_data.index_y), h, z, :BTMPV) /
                    customers.Opti_DG_E(z, h, :BTMPV) + sum(
                        min(
                            customers.rho_DG(h, :BTMPV, z, d, t) *
                            customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMPV),
                            customers.d(h, z, d, t) * (1 - utility.loss_dist)
                        ) * customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMPV) /
                        customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMPV) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    ) - 
                    # minus pv portion of pv_storage tech
                    min(
                        customers.rho_DG(h, :BTMPV, z, d, t) * customers.Opti_DG_E(z, h, :BTMPV),
                        customers.d(h, z, d, t) * (1 - utility.loss_dist)
                    ) * total_der_pv_capacity(z, h) /
                    customers.Opti_DG_E(z, h, :BTMPV)
                ) - 
                # behind the meter generation for pv+storage customers who did not participate in aggregation
                model_data.omega(d) * delta_t *
                (
                    min(
                        sum(customers.rho_DG(h, m, z, d, t) * customers.Opti_DG_E(z, h, m) for m in customers.index_m),
                        customers.d(h, z, d, t) * (1 - utility.loss_dist)
                    ) * customers.x_DG_E_my(first(model_data.index_y), h, z, :BTMStorage) /
                    customers.Opti_DG_E(z, h, :BTMStorage) + sum(
                        min(
                            sum(customers.rho_DG(h, m, z, d, t) *
                            customers.Opti_DG_my(Symbol(Int(y)), z, h, m) for m in customers.index_m),
                            customers.d(h, z, d, t) * (1 - utility.loss_dist)
                        ) * customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMStorage) /
                        customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMStorage) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    )
                ) * (1 - der_aggregator.aggregation_level(reg_year_index_dera, z))

                # model_data.omega(d) * delta_t * (
                #     min(
                #         sum(customers.rho_DG(h, m, z, d, t) * customers.Opti_DG_E(z, h, m) for m in customers.index_m),
                #         customers.d(h, z, d, t) * (1 - utility.loss_dist),
                #     ) * customers.x_DG_E_my(first(model_data.index_y), h, z, :BTMPV) /
                #     customers.Opti_DG_E(z, h, :BTMPV) + sum(
                #         min(
                #             sum(customers.rho_DG(h, m, z, d, t) *
                #             customers.Opti_DG_my(Symbol(Int(y)), z, h, m) for m in customers.index_m),
                #             customers.d(h, z, d, t) * (1 - utility.loss_dist),
                #         ) * customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMPV) /
                #         customers.Opti_DG_my(Symbol(Int(y)), z, h, :BTMPV) for
                #         y in model_data.year(first(model_data.index_y_fix)):reg_year
                #     )
                # )
        end
        net_demand_wo_green_tech_h_t_wo_loss(z, h, tou, :) .= net_demand_wo_green_tech_h_t_wo_loss_temp
    end

    # for the purpose of calculating net peak load, use load including distribution loss
    net_demand_for_peak_h_t = make_keyed_array(model_data.index_z, model_data.index_h, model_data.index_d, model_data.index_t)
    for z in model_data.index_z, h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
        net_demand_for_peak_h_t(z, h, d, t, :) .=
            sum(
                customers.gamma(z, h) * customers.d(h, z, d, t)
            ) -
            # DG
            sum(
                (
                    customers.rho_DG(h, m, z, d, t) *
                    customers.x_DG_E_my(first(model_data.index_y), h, z, m) + sum(
                        customers.rho_DG(h, m, z, d, t) * customers.x_DG_new_my(Symbol(Int(y)), h, z, m) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    )
                ) for m in customers.index_m
            ) + 
            # remove aggregated behind-the-meter storage/pv generation/consumption since they're front-of-the-meter now
            customers.rho_DG(h, :BTMStorage, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_stor_capacity(z, h) + 
            customers.rho_DG(h, :BTMPV, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_pv_capacity(z, h) -
            # green technology subscription
            # if customer h in zone z build green tech, that should be able to offest their capacity purchase obligation
            sum(
                utility.rho_C_my(j, z, d, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, z, h) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):reg_year)
                for j in model_data.index_j
            )
    end

    # excluding green tech generation offset when it comes to sharing T&D costs
    net_demand_for_peak_wo_green_tech_h_t = make_keyed_array(model_data.index_z, model_data.index_h, model_data.index_d, model_data.index_t)
    for z in model_data.index_z, h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
        net_demand_for_peak_wo_green_tech_h_t(z, h, d, t, :) .=
            sum(
                customers.gamma(z, h) * customers.d(h, z, d, t)
            ) -
            # DG
            sum(
                (
                    customers.rho_DG(h, m, z, d, t) *
                    customers.x_DG_E_my(first(model_data.index_y), h, z, m) + sum(
                        customers.rho_DG(h, m, z, d, t) * customers.x_DG_new_my(Symbol(Int(y)), h, z, m) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    )
                ) for m in customers.index_m
            ) + 
            # remove aggregated behind-the-meter storage/pv generation/consumption since they're front-of-the-meter now
            customers.rho_DG(h, :BTMStorage, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_stor_capacity(z, h) + 
            customers.rho_DG(h, :BTMPV, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_pv_capacity(z, h)
    end

    # non-coincident peak load
    net_peak_load_h = make_keyed_array(model_data.index_z, model_data.index_h)
    for z in model_data.index_z, h in model_data.index_h
        net_peak_load_h(z, h, :) .=
            findmax(Dict((d, t) => net_demand_for_peak_h_t(z, h, d, t) for d in model_data.index_d, t in model_data.index_t))[1]
    end

    net_peak_load_wo_green_tech_h = make_keyed_array(model_data.index_z, model_data.index_h)
    for z in model_data.index_z, h in model_data.index_h
        net_peak_load_wo_green_tech_h(z, h, :) .=
            findmax(Dict((d, t) => net_demand_for_peak_wo_green_tech_h_t(z, h, d, t) for d in model_data.index_d, t in model_data.index_t))[1]
    end


    # Cost Classification/Allocation
    # for cost allocation under wholesale markets, make sure the costs charged match the load being charged
    local_net_load = make_keyed_array(model_data.index_z, model_data.index_h, model_data.index_d, model_data.index_t)
    for z in model_data.index_z, h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
        local_net_load(z, h, d, t, :) .= 
            customers.gamma(z, h) * customers.d(h, z, d, t) - 
                sum(
                    customers.rho_DG(h, m, z, d, t) *
                    customers.x_DG_E_my(first(model_data.index_y), h, z, m) + sum(
                        customers.rho_DG(h, m, z, d, t) * customers.x_DG_new_my(Symbol(Int(y)), h, z, m) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    )
                    for m in customers.index_m
                ) + 
                # remove aggregated behind-the-meter storage/pv generation/consumption since they're front-of-the-meter now
                customers.rho_DG(h, :BTMStorage, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_stor_capacity(z, h) + 
                customers.rho_DG(h, :BTMPV, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_pv_capacity(z, h)
    end

    energy_purchase_quantity_detailed = make_keyed_array(model_data.index_z, model_data.index_d, model_data.index_t)
    for z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
        energy_purchase_quantity_detailed(z, d, t, :) .=
            # positive net load
            # from the perspective of load purchasing from wholesale markets
            # aggregate across consumer class, then take the maximum of net load and zero.
            max(0.0, 
                sum(
                    local_net_load(z, h, d, t) for h in model_data.index_h
                )
            ) + 
            # exogenous export
            max(0.0, ipp.eximport_my(reg_year_index, z, d, t)) + 
            # endogenous export 
            sum(max(0.0, ipp.trans_topology(l, z) * ipp.flow_my(reg_year_index, l, d, t)) for l in utility.index_l) + 
            # storage charge
            sum(ipp.charge_E_my(reg_year_index, p, s, z, d, t) for s in ipp.index_stor_existing, p in ipp.index_p) +
            sum(ipp.charge_C_my(reg_year_index, p, s, z, d, t) for s in ipp.index_stor_new, p in ipp.index_p)
    end

    local_net_load_tou = make_keyed_array(model_data.index_z, model_data.index_h, regulator.index_rate_tou)
    for z in model_data.index_z, h in model_data.index_h, tou in regulator.index_rate_tou
        local_net_load_tou_temp = 0.0
        for i in 1:size(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :])[1]
            d = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_d][i])
            t = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_t][i])

            local_net_load_tou_temp = local_net_load_tou_temp + 
                model_data.omega(d) * delta_t * 
                (customers.gamma(z, h) * customers.d(h, z, d, t) - 
                sum(
                    customers.rho_DG(h, m, z, d, t) *
                    customers.x_DG_E_my(first(model_data.index_y), h, z, m) + sum(
                        customers.rho_DG(h, m, z, d, t) * customers.x_DG_new_my(Symbol(Int(y)), h, z, m) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    )
                    for m in customers.index_m
                ) + 
                # remove aggregated behind-the-meter storage/pv generation/consumption since they're front-of-the-meter now
                customers.rho_DG(h, :BTMStorage, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_stor_capacity(z, h) + 
                customers.rho_DG(h, :BTMPV, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_pv_capacity(z, h)
                )
        end
        local_net_load_tou(z, h, tou, :) .= local_net_load_tou_temp
    end

    energy_purchase_quantity_detailed_tou = make_keyed_array(model_data.index_z, regulator.index_rate_tou)
    for z in model_data.index_z, tou in regulator.index_rate_tou
        energy_purchase_quantity_detailed_tou_temp = 0.0
        for i in 1:size(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :])[1]
            d = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_d][i])
            t = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_t][i])

            energy_purchase_quantity_detailed_tou_temp = energy_purchase_quantity_detailed_tou_temp +
                # positive net load
                model_data.omega(d) * delta_t * 
                max(0.0, 
                    sum( 
                        (customers.gamma(z, h) * customers.d(h, z, d, t) - 
                        sum(
                            customers.rho_DG(h, m, z, d, t) *
                            customers.x_DG_E_my(first(model_data.index_y), h, z, m) + sum(
                                customers.rho_DG(h, m, z, d, t) * customers.x_DG_new_my(Symbol(Int(y)), h, z, m) for
                                y in model_data.year(first(model_data.index_y_fix)):reg_year
                            )
                            for m in customers.index_m
                        ) + 
                        # remove aggregated behind-the-meter storage/pv generation/consumption since they're front-of-the-meter now
                        customers.rho_DG(h, :BTMStorage, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_stor_capacity(z, h) + 
                        customers.rho_DG(h, :BTMPV, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * total_der_pv_capacity(z, h)
                        ) for h in model_data.index_h
                    )
                ) + 
                # exogenous export
                model_data.omega(d) * delta_t * max(0.0, ipp.eximport_my(reg_year_index, z, d, t)) + 
                # endogenous export 
                model_data.omega(d) * delta_t * sum(max(0.0, ipp.trans_topology(l, z) * ipp.flow_my(reg_year_index, l, d, t)) for l in utility.index_l) + 
                # storage charge
                model_data.omega(d) * delta_t * sum(ipp.charge_E_my(reg_year_index, p, s, z, d, t) for s in ipp.index_stor_existing, p in ipp.index_p) +
                model_data.omega(d) * delta_t * sum(ipp.charge_C_my(reg_year_index, p, s, z, d, t) for s in ipp.index_stor_new, p in ipp.index_p)
        end
        energy_purchase_quantity_detailed_tou(z, tou, :) .= energy_purchase_quantity_detailed_tou_temp
    end

    energy_purchase_cost = make_keyed_array(model_data.index_z, model_data.index_h)
    for z in model_data.index_z, h in model_data.index_h
        energy_purchase_cost(z, h, :) .=
            sum(
                model_data.omega(d) * delta_t * ipp.LMP_my(reg_year_index, z, d, t) *
                local_net_load(z, h, d, t) for d in model_data.index_d, t in model_data.index_t
            )
    end

    rec_purchase_cost = make_keyed_array(model_data.index_z)
    for z in model_data.index_z
        rec_purchase_cost(z, :) .=
            # incorporate REC cost into energy purchase cost
            sum(
                regulator.REC(z, reg_year_index) *
                model_data.omega(d) * delta_t *
                (
                    sum(
                        ipp.y_E_my(reg_year_index, p, rps, z, d, t) for rps in ipp.index_rps,
                        p in ipp.index_p
                    ) + sum(
                        ipp.y_C_my(reg_year_index, p, rps, z, d, t) for rps in ipp.index_rps,
                        p in ipp.index_p
                    )
                ) for d in model_data.index_d, t in model_data.index_t
            )
    end

    energy_purchase_cost_t = make_keyed_array(model_data.index_z, model_data.index_h, regulator.index_rate_tou)
    for z in model_data.index_z, h in model_data.index_h, tou in regulator.index_rate_tou
        energy_purchase_cost_t_temp = 0.0
        for i in 1:size(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :])[1]
            d = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_d][i])
            t = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_t][i])

            energy_purchase_cost_t_temp = energy_purchase_cost_t_temp + 
                model_data.omega(d) * delta_t * ipp.LMP_my(reg_year_index, z, d, t) *
                local_net_load(z, h, d, t)
        end
        energy_purchase_cost_t(z, h, tou, :) .= energy_purchase_cost_t_temp
    end

    rec_purchase_cost_t = make_keyed_array(model_data.index_z, regulator.index_rate_tou)
    for z in model_data.index_z, tou in regulator.index_rate_tou
        rec_purchase_cost_t_temp = 0.0
        for i in 1:size(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :])[1]
            d = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_d][i])
            t = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_rate_tou.==String(tou)), :index_t][i])

            rec_purchase_cost_t_temp = rec_purchase_cost_t_temp + 
                # incorporate REC cost into energy purchase cost
                regulator.REC(z, reg_year_index) *
                model_data.omega(d) * delta_t *
                (
                    sum(
                        ipp.y_E_my(reg_year_index, p, rps, z, d, t) for rps in ipp.index_rps,
                        p in ipp.index_p
                    ) + sum(
                        ipp.y_C_my(reg_year_index, p, rps, z, d, t) for rps in ipp.index_rps,
                        p in ipp.index_p
                    )
                )
        end
        rec_purchase_cost_t(z, tou, :) .= rec_purchase_cost_t_temp
    end

    capacity_purchase_cost = ipp.capacity_price(reg_year_index) * ipp.ucap_total(reg_year_index)

    # rate-making and settltment related to consumer contracted green technologies:
    # consumers enter contracts with green developers, consumers pay for contract price (which we assume to be a rate-of-return model)
    # green-tech's contribution to energy and capacity markets will be settled at market prices (since they're modeled on the supply-side)
    # but those revenues shall be passed down to consumers since they paid off green-tech on contract prices
    # the question is should regulator consider those pass-down in rate-making?
    # current assumption is regulator does not consider those pass-down in rate-making, but those pass-down shall be considered in consumers' welfare calculation
    # this is a reasonale assumption because contracts with green developers are made on individual basis
    # it's unlikely that those contracts will be factored into rate-making for a certain consumer class
    # note that this treatment is different from BTM resources

    energy_cost_allocation_h = make_keyed_array(model_data.index_z, model_data.index_h)
    for z in model_data.index_z, h in model_data.index_h
        energy_cost_allocation_h(z, h, :) .= 
            energy_purchase_cost(z, h) + 
            rec_purchase_cost(z) * 
            sum(
                model_data.omega(d) * delta_t * local_net_load(z, h, d, t)
                for d in model_data.index_d, t in model_data.index_t
            ) / 
            sum(
                model_data.omega(d) * delta_t * 
                energy_purchase_quantity_detailed(z, d, t) for d in model_data.index_d, t in model_data.index_t
            ) +
            der_excess_cost_h(z, h)
    end

    #=
        energy_cost_allocation_sector = make_keyed_array(model_data.index_z, model_data.index_sector)
        energy_cost_allocation_sector should be summed over h that belongs to a sector
    =#

    # allocate capacity_purchase_cost and regulator.othercost (T&D costs) by net peak load without green tech generation offset;
    # The reason only net peak load of a customer type h of zone z is considered here is because export, green tech, and storage
    # are all considered in the supply-side of capacity markets.
    demand_cost_allocation_h = make_keyed_array(model_data.index_z, model_data.index_h)
    for z in model_data.index_z, h in model_data.index_h
        demand_cost_allocation_h(z, h, :) .= 
        capacity_purchase_cost * net_peak_load_wo_green_tech_h(z, h) / sum(net_peak_load_wo_green_tech_h(z, h) for z in model_data.index_z, h in model_data.index_h) + 
        regulator.othercost(z, reg_year_index) * net_peak_load_wo_green_tech_h(z, h) / (
            sum(net_peak_load_wo_green_tech_h(z, h) for h in model_data.index_h)
        )
    end
    replace!(demand_cost_allocation_h, NaN => 0.0)

    demand_cost_allocation_capacity_h = make_keyed_array(model_data.index_z, model_data.index_h)
    for z in model_data.index_z, h in model_data.index_h
        demand_cost_allocation_capacity_h(z, h, :) .=
        capacity_purchase_cost * net_peak_load_wo_green_tech_h(z, h) / sum(net_peak_load_wo_green_tech_h(z, h) for z in model_data.index_z, h in model_data.index_h)
    end
    replace!(demand_cost_allocation_capacity_h, NaN => 0.0)

    demand_cost_allocation_othercost_h = make_keyed_array(model_data.index_z, model_data.index_h)
    for z in model_data.index_z, h in model_data.index_h
        demand_cost_allocation_othercost_h(z, h, :) .=
            regulator.othercost(z, reg_year_index) * net_peak_load_wo_green_tech_h(z, h) / (
                sum(net_peak_load_wo_green_tech_h(z, h) for h in model_data.index_h)
            )
    end
    replace!(demand_cost_allocation_othercost_h, NaN => 0.0)


    energy_cost_allocation_h_t = make_keyed_array(model_data.index_z, model_data.index_h, regulator.index_rate_tou)
    for z in model_data.index_z, h in model_data.index_h, tou in regulator.index_rate_tou
        energy_cost_allocation_h_t(z, h, tou, :) .=
            energy_purchase_cost_t(z, h, tou) + 
            rec_purchase_cost_t(z, tou) * local_net_load_tou(z, h, tou) / energy_purchase_quantity_detailed_tou(z, tou) +
            der_excess_cost_h_t(z, h, tou)
    end
    replace!(energy_cost_allocation_h_t, NaN => 0.0)


    p_before = ParamArray(regulator.p, "p_before")
    fill!(p_before, NaN)
    for z in model_data.index_z, h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
        p_before(z, h, d, t, :) .= regulator.p_my(reg_year_index, z, h, d, t)
    end

    p_before_wavg = ParamArray(regulator.p_td, "p_before_wavg")
    fill!(p_before_wavg, NaN)  # TODO DT: debug only
    for z in model_data.index_z, h in model_data.index_h
        p_before_wavg(z, h, :) .= 
            sum(regulator.p_my(reg_year_index, z, h, d, t) * model_data.omega(d) * delta_t * customers.d(h, z, d, t) for d in model_data.index_d, t in model_data.index_t) / 
            sum(model_data.omega(d) * delta_t * customers.d(h, z, d, t) for d in model_data.index_d, t in model_data.index_t)
    end

    # TODO: Call a function instead of using if-then
    # TODO: the demonimator need to be further thought through (in the case without green-tech, it's the same)
    if regulator_opts.rate_design isa FlatRate
        fill!(regulator.p, NaN)
        for h in model_data.index_h, t in model_data.index_t
            for z in model_data.index_z, h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
                regulator.p(z, h, d, t, :) .=
                    (energy_cost_allocation_h(z, h) + demand_cost_allocation_capacity_h(z, h)) /
                    net_demand_wo_green_tech_h_wo_loss(z, h) +
                    demand_cost_allocation_othercost_h(z, h) / net_demand_wo_green_tech_h_wo_loss(z, h)
            end
            replace!(regulator.p.values, NaN => 0.0)
        end
    elseif regulator_opts.rate_design isa TOU
        fill!(regulator.p, NaN)
        for z in model_data.index_z, h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
            tou = Symbol(regulator.rep_day_time_tou_mapping[(regulator.rep_day_time_tou_mapping.index_d.==String(d)) .& (regulator.rep_day_time_tou_mapping.index_t.==String(t)), :index_rate_tou][1])
            regulator.p(z, h, d, t, :) .=
                energy_cost_allocation_h_t(z, h, tou) / net_demand_wo_green_tech_h_t_wo_loss(z, h, tou) +
                demand_cost_allocation_capacity_h(z, h) / net_demand_wo_green_tech_h_wo_loss(z, h) +
                demand_cost_allocation_othercost_h(z, h) / net_demand_wo_green_tech_h_wo_loss(z, h)
        end
        replace!(regulator.p.values, NaN => 0.0)
    end

    # fill!(regulator.p_regression, NaN)
    # for h in model_data.index_h
    #     regulator.p_regression(h, :) .=
    #         (energy_cost_allocation_h(h) + demand_cost_allocation_capacity_h(h)) /
    #         net_demand_h_wo_loss(h)
    # end

    # fill!(regulator.p_td, NaN)
    # for h in model_data.index_h
    #     regulator.p_td(h, :) .=
    #         demand_cost_allocation_othercost_h(h) / net_demand_wo_green_tech_h_wo_loss(h)
    # end

    # TODO: Call a function instead of using if-then
    if regulator_opts.net_metering_policy isa ExcessRetailRate
        regulator.p_ex = ParamArray(regulator.p)
    elseif regulator_opts.net_metering_policy isa ExcessMarginalCost
        fill!(regulator.p_ex, NaN)
        for z in model_data.index_z, h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
            regulator.p_ex(z, h, d, t, :) .= ipp.LMP_my(reg_year_index, z, d, t)
        end
    elseif regulator_opts.net_metering_policy isa ExcessZero
        fill!(regulator.p_ex, 0.0)
    end

    for z in model_data.index_z, h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
        regulator.p_my(reg_year_index, z, h, d, t, :) .= regulator.p(z, h, d, t)
        regulator.p_ex_my(reg_year_index, z, h, d, t, :) .= regulator.p_ex(z, h, d, t)
    end

    # for h in model_data.index_h
    #     regulator.p_my_regression(reg_year_index, h, :) .= regulator.p_regression(h)
    #     regulator.p_my_td(reg_year_index, h, :) .= regulator.p_td(h)
    # end

    p_after_wavg = ParamArray(regulator.p_td, "p_after_wavg")
    fill!(p_after_wavg, NaN)  # TODO DT: debug only
    for z in model_data.index_z, h in model_data.index_h
        p_after_wavg(z, h, :) .= 
            sum(regulator.p_my(reg_year_index, z, h, d, t) * model_data.omega(d) * delta_t * customers.d(h, z, d, t) for d in model_data.index_d, t in model_data.index_t) / 
            sum(model_data.omega(d) * delta_t * customers.d(h, z, d, t) for d in model_data.index_d, t in model_data.index_t)
    end

    # @info "Original retail price" p_before
    # @info "Original DER excess rate" p_ex_before
    @info "New retail price" regulator.p
    @info "New DER excess rate" regulator.p_ex

    # return compute_difference_percentage_one_norm([
    #     (p_before_wavg.values, p_after_wavg.values),
    # ])
    return compute_difference_percentage_maximum_one_norm([
        (p_before.values, regulator.p.values),
    ])
end

function save_results(
    regulator::Regulator,
    regulator_opts::RegulatorOptions,
    hem_opts::HEMOptions{<:MarketStructure},
    export_file_path::AbstractString,
)
    # Primal Variables
    save_param(
        regulator.p_my.values,
        [:Year, :Zone, :CustomerType, :Day, :Time],
        :Price,
        joinpath(export_file_path, "p.csv"),
    )
    save_param(
        regulator.p_ex_my.values,
        [:Year, :Zone, :CustomerType, :Day, :Time],
        :Price,
        joinpath(export_file_path, "p_ex.csv"),
    )
end
