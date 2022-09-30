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
    # Parameters
    "planning reserve (fraction)"
    r::ParamScalar
    "allowed return on investment (fraction)"
    z::ParamScalar
    distribution_cost::ParamArray
    administration_cost::ParamArray
    transmission_cost::ParamArray
    interconnection_cost::ParamArray
    system_cost::ParamArray
    "other cost not related to the optimization problem"
    othercost::ParamArray
    "Renewable Energy Credits"
    REC::ParamScalar

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

    distribution_cost = read_param(
        "distribution_cost",
        input_filename,
        "distribution_cost",
        model_data.index_y,
        description = "distribution cost (dollar)",
    )

    administration_cost = read_param(
        "administration_cost",
        input_filename,
        "administration_cost",
        model_data.index_y,
        description = "administration cost (dollar)",
    )

    transmission_cost = read_param(
        "transmission_cost",
        input_filename,
        "transmission_cost",
        model_data.index_y,
        description = "transmission cost (dollar)",
    )

    interconnection_cost = read_param(
        "interconnection_cost",
        input_filename,
        "interconnection_cost",
        model_data.index_y,
        description = "interconnection cost (dollar)",
    )

    system_cost = read_param(
        "system_cost",
        input_filename,
        "system_cost",
        model_data.index_y,
        description = "system cost (dollar)",
    )

    return Regulator(
        id,
        ParamScalar("r", 0.12, description = "planning reserve (fraction)"),
        ParamScalar("z", 0.112, description = "allowed return on investment (fraction)"),
        distribution_cost,
        administration_cost,
        transmission_cost,
        interconnection_cost,
        system_cost,
        initialize_param(
            "othercost",
            model_data.index_y;
            description = "other cost not related to the optimization problem",
        ),
        ParamScalar("REC", 30.0, description = "Renewable Energy Credits"),
        initialize_param(
            "p",
            model_data.index_h,
            model_data.index_t;
            value = 10.0,
            description = "retail price",
        ),
        initialize_param(
            "p_ex",
            model_data.index_h,
            model_data.index_t;
            value = 10.0,
            description = "DER excess generation rate",
        ),
        initialize_param(
            "p_eximport",
            model_data.index_t;
            value = 10.0,
            description = "retail price of import/export",
        ),
        initialize_param(
            "p_green",
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
            model_data.index_h,
            model_data.index_t;
            value = 10.0,
            description = "multi-year retail price",
        ),
        initialize_param(
            "p_ex_my",
            model_data.index_y,
            model_data.index_h,
            model_data.index_t;
            value = 10.0,
            description = "multi-year DER excess generation rate",
        ),
        initialize_param(
            "p_eximport_my",
            model_data.index_y,
            model_data.index_t;
            value = 10.0,
            description = "multi-year retail price of import/export",
        ),
        initialize_param(
            "revenue_req_my",
            model_data.index_y;
            description = "revenue (requirement) of utility company by year",
        ),
        initialize_param(
            "cost_my",
            model_data.index_y;
            description = "cost of utility company (without return on equity) by year",
        ),
        initialize_param(
            "debt_interest_my",
            model_data.index_y;
            description = "debt interest by year",
        ),
        initialize_param(
            "income_tax_my",
            model_data.index_y;
            description = "income tax by year",
        ),
        initialize_param(
            "operational_cost_my",
            model_data.index_y;
            description = "operational cost by year",
        ),
        initialize_param(
            "depreciation_my",
            model_data.index_y;
            description = "accounting depreciation by year",
        ),
        initialize_param(
            "depreciation_tax_my",
            model_data.index_y;
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
            model_data.index_h;
            value = 0.0,
            description = "T&D component charge",
        ),
        initialize_param(
            "p_my_td",
            model_data.index_y,
            model_data.index_h;
            value = 0.0,
            description = "multi-year T&D component charge",
        ),
    )
end

get_id(x::Regulator) = x.id

function get_file_prefix(agent::Regulator)
    return "REC$(agent.REC.value)"
end

# although Customer is subtype of Agent, 
# Vector{Customer} is not subtype of Vector{Agent}
# But if a vector of customers c1, c2, c3 is defined 
# using the syntax Agent[c1, c2, c3], calling 
# this function will work. Can also:
# Vector{Agent}([c1, c2, c3])

function solve_agent_problem!(
    regulator::Regulator,
    regulator_opts::RegulatorOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{VerticallyIntegratedUtility},
    agent_store::AgentStore,
    w_iter,
)

    for y in model_data.index_y_fix
        regulator.othercost(y, :) .= regulator.distribution_cost(y) + regulator.administration_cost(y) + regulator.transmission_cost(y) + regulator.interconnection_cost(y) + regulator.system_cost(y)
    end

    utility = get_agent(Utility, agent_store)
    customers = get_agent(CustomerGroup, agent_store)
    green_developer = get_agent(GreenDeveloper, agent_store)

    # the year regulator is making a rate case
    reg_year = model_data.year(first(model_data.index_y))
    reg_year_index = Symbol(Int(reg_year))
    # retirement up to the rate case year
    reg_retirement = KeyedArray(
        [
            sum(
                utility.x_R_my(Symbol(Int(y_symbol)), k) for
                y_symbol in model_data.year(first(model_data.index_y_fix)):reg_year
            ) for k in utility.index_k_existing
        ];
        [get_pair(utility.index_k_existing)]...
    )

    # pure volumetric rate
    energy_cost =
        sum(
            model_data.omega(t) *
            (utility.v_E_my(reg_year_index, k, t) * utility.y_E_my(reg_year_index, k, t))
            for k in utility.index_k_existing, t in model_data.index_t
        ) + sum(
            model_data.omega(t) *
            (utility.v_C_my(reg_year_index, k, t) * utility.y_C_my(reg_year_index, k, t))
            for k in utility.index_k_new, t in model_data.index_t
        )
    fixed_om =
        sum(
            utility.fom_E_my(reg_year_index, k) * (utility.x_E_my(k) - reg_retirement(k))
            for k in utility.index_k_existing
        ) + sum(
            utility.fom_C_my(Symbol(Int(y_symbol)), k) *
            utility.x_C_my(Symbol(Int(y_symbol)), k) for
            y_symbol in model_data.year(first(model_data.index_y_fix)):reg_year,
            k in utility.index_k_new
        )
    operational_cost = energy_cost + fixed_om
    working_capital = utility.DaysofWC / 365 * operational_cost

    # calculate ADIT and rate base (no working capital) for new builds, "reg_year-y+1" represents the number of years since the new investment is made
    ADITNew = KeyedArray(
        [
            sum(
                utility.CapEx_my(Symbol(Int(y)), k) *
                utility.x_C_my(Symbol(Int(y)), k) *
                (
                    utility.CumuTaxDepre_new_my(Symbol(Int(reg_year - y + 1)), k) -
                    utility.CumuAccoutDepre_new_my(Symbol(Int(reg_year - y + 1)), k)
                ) *
                utility.Tax +
                utility.ITC_new_my(Symbol(Int(y)), k) *
                utility.CapEx_my(Symbol(Int(y)), k) *
                utility.x_C_my(Symbol(Int(y)), k) *
                (1 - utility.CumuITCAmort_new_my(Symbol(Int(reg_year - y + 1)), k)) for
                y in model_data.year(first(model_data.index_y_fix)):reg_year
            ) for k in utility.index_k_new
        ];
        [get_pair(utility.index_k_new)]...
    )
    RateBaseNoWC_new = KeyedArray(
        [
            sum(
                utility.CapEx_my(Symbol(Int(y)), k) *
                utility.x_C_my(Symbol(Int(y)), k) *
                (1 - utility.CumuAccoutDepre_new_my(Symbol(Int(reg_year - y + 1)), k))
                for y in model_data.year(first(model_data.index_y_fix)):reg_year
            ) - ADITNew(k) for k in utility.index_k_new
        ];
        [get_pair(utility.index_k_new)]...
    )

    # calculate total rate base for the year of rate making
    rate_base =
        sum(
            utility.RateBaseNoWC_existing_my(reg_year_index, k) *
            (utility.x_E_my(k) - reg_retirement(k)) for k in utility.index_k_existing
        ) +
        sum(RateBaseNoWC_new(k) for k in utility.index_k_new) +
        working_capital
    debt_interest = rate_base * utility.DebtRatio * utility.COD
    # calculate total depreciation
    depreciation =
    # annual depreciation on existing units that have not retired
        sum(
            utility.CapEx_existing_my(k) *
            (utility.x_E_my(k) - reg_retirement(k)) *
            utility.AnnualAccoutDepre_existing_my(reg_year_index, k) +
            # existing units that are retired this year will incur their regular annual depreciation, as well as the remaining un-depreciated asset
            utility.CapEx_existing_my(k) *
            utility.x_R_my(reg_year_index, k) *
            (
                utility.AnnualAccoutDepre_existing_my(reg_year_index, k) + 1 -
                utility.CumuAccoutDepre_existing_my(reg_year_index, k)
            ) for k in utility.index_k_existing
        ) +
        # annual depreciation on new units
        sum(
            utility.CapEx_my(Symbol(Int(y)), k) *
            utility.x_C_my(Symbol(Int(y)), k) *
            utility.AnnualAccoutDepre_new_my(Symbol(Int(reg_year - y + 1)), k) for
            y in model_data.year(first(model_data.index_y_fix)):reg_year,
            k in utility.index_k_new
        )
    # calculate total tax depreciation
    depreciation_tax =
    # annual depreciation on existing units that have not retired
        sum(
            utility.CapEx_existing_my(k) *
            (utility.x_E_my(k) - reg_retirement(k)) *
            utility.AnnualTaxDepre_existing_my(reg_year_index, k) +
            # existing units that are retired this year will incur their regular annual depreciation, as well as the remaining un-depreciated asset
            utility.CapEx_existing_my(k) *
            utility.x_R_my(reg_year_index, k) *
            (
                utility.AnnualTaxDepre_existing_my(reg_year_index, k) + 1 -
                utility.CumuTaxDepre_existing_my(reg_year_index, k)
            ) for k in utility.index_k_existing
        ) +
        # annual depreciation on new units
        sum(
            utility.CapEx_my(Symbol(Int(y)), k) *
            utility.x_C_my(Symbol(Int(y)), k) *
            utility.AnnualTaxDepre_new_my(Symbol(Int(reg_year - y + 1)), k) for
            y in model_data.year(first(model_data.index_y_fix)):reg_year,
            k in utility.index_k_new
        )
    return_to_equity = rate_base * (1 - utility.DebtRatio) * utility.COE
    #=
    income_tax = (return_to_equity -
        sum(utility.CapEx_my(reg_year_index,k)*utility.x_C_my(reg_year_index,k)*utility.ITC_new_my(reg_year_index,k) for k in utility.index_k_new)) /
        (1-utility.Tax) - return_to_equity
    =#
    income_tax =
        (
            return_to_equity * utility.Tax +
            (depreciation - depreciation_tax) * utility.Tax - 
            sum(
                utility.ITC_existing_my(k) *
                utility.CapEx_existing_my(k) *
                (utility.x_E_my(k) - reg_retirement(k)) *
                utility.AnnualITCAmort_existing_my(reg_year_index, k) +
                # existing units that are retired this year will incur their regular annual depreciation, as well as the remaining un-depreciated asset
                utility.ITC_existing_my(k) *
                utility.CapEx_existing_my(k) *
                utility.x_R_my(reg_year_index, k) *
                (
                    utility.AnnualITCAmort_existing_my(reg_year_index, k) + 1 -
                    utility.CumuITCAmort_existing_my(reg_year_index, k)
                ) for k in utility.index_k_existing
            ) -
            sum(
                utility.ITC_new_my(Symbol(Int(y)), k) *
                utility.CapEx_my(Symbol(Int(y)), k) *
                utility.x_C_my(Symbol(Int(y)), k) *
                utility.AnnualITCAmort_new_my(Symbol(Int(reg_year - y + 1)), k) for
                y in model_data.year(first(model_data.index_y_fix)):reg_year, k in utility.index_k_new
            )
        ) / (1 - utility.Tax)
    # calculate revenue requirement
    revenue_requirement =
        debt_interest + return_to_equity + income_tax + operational_cost + depreciation

    regulator.revenue_req_my(reg_year_index, :) .= revenue_requirement
    regulator.cost_my(reg_year_index, :) .=
        debt_interest + income_tax + operational_cost + depreciation

    regulator.debt_interest_my(reg_year_index, :) .= debt_interest
    regulator.income_tax_my(reg_year_index, :) .= income_tax
    regulator.operational_cost_my(reg_year_index, :) .= operational_cost
    regulator.depreciation_my(reg_year_index, :) .= depreciation
    regulator.depreciation_tax_my(reg_year_index, :) .= depreciation_tax

    # when it comes to net demand, calculate two values: the one with loss is used for cost allocation;
    # the one without loss is used for rate calculation.
    net_demand_w_loss = (
        # demand
        # when it comes to sharing the revenue requirement (cost), use load including distribution loss
        # e.g. utility generation is 100 MW, 50 MW to serve load (including distribution loss),
        #      50 MW for export. It makes sense to allocate the same cost for internal load and export.
        sum(
            customers.gamma(h) * model_data.omega(t) * customers.d(h, t) for
            h in model_data.index_h, t in model_data.index_t
        ) +
        # export
        sum(
            model_data.omega(t) * utility.eximport_my(reg_year_index, t) for
            t in model_data.index_t
        ) -
        # DG
        # when it comes to cost allocation, simply use total DER generation to offset total load;
        # this does not consider enery offset at the household level, which is only considered when 
        # calculating retail rates.
        sum(
            model_data.omega(t) * (
                customers.rho_DG(h, m, t) *
                customers.x_DG_E_my(first(model_data.index_y), h, m) + sum(
                    customers.rho_DG(h, m, t) * customers.x_DG_new_my(Symbol(Int(y)), h, m) for
                    y in model_data.year(first(model_data.index_y_fix)):reg_year
                )
            ) for t in model_data.index_t, h in model_data.index_h,
            m in customers.index_m
        ) -
        # green technology subscription
        sum(
            model_data.omega(t) * utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
            model_data.year(first(model_data.index_y_fix)):reg_year)
            for t in model_data.index_t, j in model_data.index_j, h in model_data.index_h
        )
    )

    # In the presence of distribution loss, multiply load (including distribution loss) by (1 - loss factor),
    # DER generation above this value shall be compansated for DER excess credits
    # e.g. DER generation is 100 MW, load (without loss) is 95 MW, receive 5 MW excess credits
    der_excess_cost_h = KeyedArray(
        [
            sum(
                model_data.omega(t) *
                regulator.p_ex(h, t) *
                (
                    max(
                        0,
                        customers.rho_DG(h, m, t) * customers.Opti_DG_E(h, m) -
                        customers.d(h, t) * (1 - utility.loss_dist),
                    ) * customers.x_DG_E_my(first(model_data.index_y), h, m) /
                    customers.Opti_DG_E(h, m) + sum(
                        max(
                            0,
                            customers.rho_DG(h, m, t) *
                            customers.Opti_DG_my(Symbol(Int(y)), h, m) -
                            customers.d(h, t) * (1 - utility.loss_dist),
                        ) * customers.x_DG_new_my(Symbol(Int(y)), h, m) /
                        customers.Opti_DG_my(Symbol(Int(y)), h, m) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    )
                ) for t in model_data.index_t, m in customers.index_m
            ) for h in model_data.index_h
        ];
        [get_pair(model_data.index_h)]...
    )

    net_demand_h_w_loss = KeyedArray(
        # demand
        [
            sum(
                customers.gamma(h) * model_data.omega(t) * customers.d(h, t) for
                t in model_data.index_t
            ) -
            # DG
            sum(
                model_data.omega(t) * (
                    customers.rho_DG(h, m, t) *
                    customers.x_DG_E_my(first(model_data.index_y), h, m) + sum(
                        customers.rho_DG(h, m, t) *
                        customers.x_DG_new_my(Symbol(Int(y)), h, m) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    )
                ) for t in model_data.index_t, m in customers.index_m
            ) -
            # green technology subscription
            sum(
                model_data.omega(t) * utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):reg_year)
                for t in model_data.index_t, j in model_data.index_j
            )
            for h in model_data.index_h
        ];
        [get_pair(model_data.index_h)]...
    )

    net_demand_h_wo_loss = KeyedArray(
        # demand
        [
            sum(
                customers.gamma(h) *
                model_data.omega(t) *
                customers.d(h, t) *
                (1 - utility.loss_dist) for t in model_data.index_t
            ) -
            # DG
            # since this net demand is for rate calculation, we need to consider energy offset at the household level.
            # e.g. two household, with 100 MW load each (without loss), if one of them has DER and generated 
            # 120 MW of energy, he/she does not need to pay for energy, but the other one still have to pay 
            # 100 MW instead of 80 MW.
            sum(
                model_data.omega(t) * (
                    min(
                        customers.rho_DG(h, m, t) * customers.Opti_DG_E(h, m),
                        customers.d(h, t) * (1 - utility.loss_dist),
                    ) * customers.x_DG_E_my(first(model_data.index_y), h, m) /
                    customers.Opti_DG_E(h, m) + sum(
                        min(
                            customers.rho_DG(h, m, t) *
                            customers.Opti_DG_my(Symbol(Int(y)), h, m),
                            customers.d(h, t) * (1 - utility.loss_dist),
                        ) * customers.x_DG_new_my(Symbol(Int(y)), h, m) /
                        customers.Opti_DG_my(Symbol(Int(y)), h, m) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    )
                ) for t in model_data.index_t, m in customers.index_m
            ) -
            # green technology subscription
            sum(
                model_data.omega(t) * utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):reg_year)
                for t in model_data.index_t, j in model_data.index_j
            )
            for h in model_data.index_h
        ];
        [get_pair(model_data.index_h)]...
    )

    net_demand_wo_green_tech_h_wo_loss = KeyedArray(
        [
            # demand
            sum(
                customers.gamma(h) *
                model_data.omega(t) *
                customers.d(h, t) *
                (1 - utility.loss_dist) for t in model_data.index_t
            ) -
            # DG
            # since this net demand is for rate calculation, we need to consider energy offset at the household level.
            # e.g. two household, with 100 MW load each (without loss), if one of them has DER and generated 
            # 120 MW of energy, he/she does not need to pay for energy, but the other one still have to pay 
            # 100 MW instead of 80 MW.
            sum(
                model_data.omega(t) * (
                    min(
                        customers.rho_DG(h, m, t) * customers.Opti_DG_E(h, m),
                        customers.d(h, t) * (1 - utility.loss_dist),
                    ) * customers.x_DG_E_my(first(model_data.index_y), h, m) /
                    customers.Opti_DG_E(h, m) + sum(
                        min(
                            customers.rho_DG(h, m, t) *
                            customers.Opti_DG_my(Symbol(Int(y)), h, m),
                            customers.d(h, t) * (1 - utility.loss_dist),
                        ) * customers.x_DG_new_my(Symbol(Int(y)), h, m) /
                        customers.Opti_DG_my(Symbol(Int(y)), h, m) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    )
                ) for t in model_data.index_t, m in customers.index_m
            )
            for h in model_data.index_h
        ];
        [get_pair(model_data.index_h)]...
    )

    #=
    energy_cost_t = Dict(t => sum(utility.v_E(k, t)*utility.y_E(k, t) for k in utility.index_k_existing) +
        sum(utility.v_C(k, t)*utility.y_C(k, t) for k in utility.index_k_new) for t in model_data.index_t)
    =#
    energy_cost_t = KeyedArray(
        [
            sum(
                utility.v_E_my(reg_year_index, k, t) * utility.y_E_my(reg_year_index, k, t) for k in utility.index_k_existing
            ) + sum(
                utility.v_C_my(reg_year_index, k, t) * utility.y_C_my(reg_year_index, k, t) for k in utility.index_k_new
            ) for t in model_data.index_t
        ];
        [get_pair(model_data.index_t)]...
    )

    net_demand_t_w_loss = KeyedArray(
        # demand
        [
            sum(customers.gamma(h) * customers.d(h, t) for h in model_data.index_h) +
            # export
            utility.eximport_my(reg_year_index, t) -
            # DG
            sum(
                customers.rho_DG(h, m, t) *
                customers.x_DG_E_my(first(model_data.index_y), h, m) + sum(
                    customers.rho_DG(h, m, t) * customers.x_DG_new_my(Symbol(Int(y)), h, m) for
                    y in model_data.year(first(model_data.index_y_fix)):reg_year
                ) for h in model_data.index_h, m in customers.index_m
            ) -
            # green technology subscription
            sum(
                utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):reg_year)
                for h in model_data.index_h, j in model_data.index_j
            )
            for t in model_data.index_t
        ];
        [get_pair(model_data.index_t)]...
    )

    der_excess_cost_h_t = make_keyed_array(model_data.index_h, model_data.index_t)
    for h in model_data.index_h, t in model_data.index_t
        der_excess_cost_h_t(h, t, :) .= sum(
            regulator.p_ex(h, t) * (
                max(
                    0,
                    customers.rho_DG(h, m, t) * customers.Opti_DG_E(h, m) -
                    customers.d(h, t) * (1 - utility.loss_dist),
                ) * customers.x_DG_E_my(first(model_data.index_y), h, m) /
                customers.Opti_DG_E(h, m) + sum(
                    max(
                        0,
                        customers.rho_DG(h, m, t) *
                        customers.Opti_DG_my(Symbol(Int(y)), h, m) -
                        customers.d(h, t) * (1 - utility.loss_dist),
                    ) * customers.x_DG_new_my(Symbol(Int(y)), h, m) /
                    customers.Opti_DG_my(Symbol(Int(y)), h, m) for
                    y in model_data.year(first(model_data.index_y_fix)):reg_year
                )
            ) for m in customers.index_m
        )
    end

    net_demand_h_t_w_loss = make_keyed_array(model_data.index_h, model_data.index_t)
    for h in model_data.index_h, t in model_data.index_t
        net_demand_h_t_w_loss(h, t, :) .=
        # demand
            customers.gamma(h) * customers.d(h, t) -
            # DG
            sum(
                customers.rho_DG(h, m, t) *
                customers.x_DG_E_my(first(model_data.index_y), h, m) + sum(
                    customers.rho_DG(h, m, t) * customers.x_DG_new_my(Symbol(Int(y)), h, m)
                    for y in model_data.year(first(model_data.index_y_fix)):reg_year
                ) for m in customers.index_m
            ) -
            sum(
                utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):reg_year)
                for j in model_data.index_j
            )
    end

    net_demand_h_t_wo_loss = make_keyed_array(model_data.index_h, model_data.index_t)
    for h in model_data.index_h, t in model_data.index_t
        net_demand_h_t_wo_loss(h, t, :) .=
        # demand
            customers.gamma(h) * customers.d(h, t) * (1 - utility.loss_dist) -
            # DG
            sum(
                min(
                    customers.rho_DG(h, m, t) * customers.Opti_DG_E(h, m),
                    customers.d(h, t) * (1 - utility.loss_dist),
                ) * customers.x_DG_E_my(first(model_data.index_y), h, m) /
                customers.Opti_DG_E(h, m) + sum(
                    min(
                        customers.rho_DG(h, m, t) *
                        customers.Opti_DG_my(Symbol(Int(y)), h, m),
                        customers.d(h, t) * (1 - utility.loss_dist),
                    ) * customers.x_DG_new_my(Symbol(Int(y)), h, m) /
                    customers.Opti_DG_my(Symbol(Int(y)), h, m) for
                    y in model_data.year(first(model_data.index_y_fix)):reg_year
                ) for m in customers.index_m
            ) -
            sum(
                utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):reg_year)
                for j in model_data.index_j
            )
    end

    # for the purpose of calculating net peak load, use load including distribution loss
    net_demand_for_peak_h_t = make_keyed_array(model_data.index_h, model_data.index_t)
    for h in model_data.index_h, t in model_data.index_t
        net_demand_for_peak_h_t(h, t, :) .=
            customers.gamma(h) * customers.d(h, t) - sum(
                customers.rho_DG(h, m, t) *
                customers.x_DG_E_my(first(model_data.index_y), h, m) + sum(
                    customers.rho_DG(h, m, t) * customers.x_DG_new_my(Symbol(Int(y)), h, m)
                    for y in model_data.year(first(model_data.index_y_fix)):reg_year
                ) for m in customers.index_m
            ) -
            sum(
                utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):reg_year)
                for j in model_data.index_j
            )
    end

    # excluding green tech generation offset when it comes to sharing T&D costs
    net_demand_for_peak_wo_green_tech_h_t = make_keyed_array(model_data.index_h, model_data.index_t)
    for h in model_data.index_h, t in model_data.index_t
        net_demand_for_peak_wo_green_tech_h_t(h, t, :) .=
            customers.gamma(h) * customers.d(h, t) - sum(
                customers.rho_DG(h, m, t) *
                customers.x_DG_E_my(first(model_data.index_y), h, m) + sum(
                    customers.rho_DG(h, m, t) * customers.x_DG_new_my(Symbol(Int(y)), h, m) for
                    y in model_data.year(first(model_data.index_y_fix)):reg_year
                ) for m in customers.index_m
            )
    end

    net_peak_load_h = KeyedArray(
        [
            findmax(Dict(t => net_demand_for_peak_h_t(h, t) for t in model_data.index_t))[1] for h in model_data.index_h
        ];
        [get_pair(model_data.index_h)]...
    )

    net_peak_load_wo_green_tech_h = KeyedArray(
        [ 
            findmax(Dict(t => net_demand_for_peak_wo_green_tech_h_t(h, t) for t in model_data.index_t))[1] for h in model_data.index_h
        ];
        [get_pair(model_data.index_h)]...
    )

    # Cost Classification/Allocation
    energy_cost_allocation_h = KeyedArray(
        [
            energy_cost * net_demand_h_w_loss(h) / net_demand_w_loss + der_excess_cost_h(h) for h in model_data.index_h
        ];
        [get_pair(model_data.index_h)]...
    )
    energy_cost_allocation_eximport =
        energy_cost * sum(
            model_data.omega(t) * utility.eximport_my(reg_year_index, t) for
            t in model_data.index_t
        ) / net_demand_w_loss

    # allocate (revenue_requirement - energy_cost) by net peak load with green tech generation offset;
    # allocate regulator.othercost (T&D costs) by net peak load without green tech generation offset;
    demand_cost_allocation_h = KeyedArray(
        [
            (revenue_requirement - energy_cost) * net_peak_load_h(h) / (
                sum(net_peak_load_h(h) for h in model_data.index_h) +
                utility.Peak_eximport_my(reg_year_index)
            ) + 
            regulator.othercost(reg_year_index) * net_peak_load_wo_green_tech_h(h) / (
                sum(net_peak_load_wo_green_tech_h(h) for h in model_data.index_h) +
                utility.Peak_eximport_my(reg_year_index)
            )
            for h in model_data.index_h
        ];
        [get_pair(model_data.index_h)]...
    )
    demand_cost_allocation_capacity_h = KeyedArray(
        [
            (revenue_requirement - energy_cost) * net_peak_load_h(h) / (
                sum(net_peak_load_h(h) for h in model_data.index_h) +
                utility.Peak_eximport_my(reg_year_index)
            )
            for h in model_data.index_h
        ];
        [get_pair(model_data.index_h)]...
    )
    demand_cost_allocation_othercost_h = KeyedArray(
        [
            regulator.othercost(reg_year_index) * net_peak_load_wo_green_tech_h(h) / (
                sum(net_peak_load_wo_green_tech_h(h) for h in model_data.index_h) +
                utility.Peak_eximport_my(reg_year_index)
            )
            for h in model_data.index_h
        ];
        [get_pair(model_data.index_h)]...
    )

    demand_cost_allocation_eximport =
        (revenue_requirement - energy_cost) *
        utility.Peak_eximport_my(reg_year_index) / (
            sum(net_peak_load_h(h) for h in model_data.index_h) +
            utility.Peak_eximport_my(reg_year_index)
        ) +
        regulator.othercost(reg_year_index) *
        utility.Peak_eximport_my(reg_year_index) / (
            sum(net_peak_load_wo_green_tech_h(h) for h in model_data.index_h) +
            utility.Peak_eximport_my(reg_year_index)
        )

    energy_cost_allocation_h_t = make_keyed_array(model_data.index_h, model_data.index_t)
    for h in model_data.index_h, t in model_data.index_t
        energy_cost_allocation_h_t(h, t, :) .=
            energy_cost_t(t) * net_demand_h_t_w_loss(h, t) / net_demand_t_w_loss(t) +
            der_excess_cost_h_t(h, t)
    end
    energy_cost_allocation_eximport_t = KeyedArray(
        [
            energy_cost_t(t) * utility.eximport_my(reg_year_index, t) /
            net_demand_t_w_loss(t) for t in model_data.index_t
        ];
        [get_pair(model_data.index_t)]...
    )

    # compute the retail price
    p_before = ParamArray(regulator.p, "p_before")
    fill!(p_before, NaN)  # TODO DT: debug only
    for h in model_data.index_h, t in model_data.index_t
        p_before(h, t, :) .= regulator.p_my(reg_year_index, h, t)
    end

    p_before_wavg = ParamArray(regulator.p_td, "p_before_wavg")
    fill!(p_before_wavg, NaN)  # TODO DT: debug only
    for h in model_data.index_h
        p_before_wavg(h, :) .= 
            sum(regulator.p_my(reg_year_index, h, t) * model_data.omega(t) * customers.d(h, t) for t in model_data.index_t) / 
            sum(model_data.omega(t) * customers.d(h, t) for t in model_data.index_t)
    end

    p_ex_before = ParamArray(regulator.p_ex, "p_ex_before")
    fill!(p_ex_before, NaN)
    for h in model_data.index_h, t in model_data.index_t
        p_ex_before(h, t, :) .= regulator.p_ex_my(reg_year_index, h, t)
    end

    # TODO: Call a function instead of using if-then
    if regulator_opts.rate_design isa FlatRate
        fill!(regulator.p, NaN)
        for h in model_data.index_h, t in model_data.index_t
            regulator.p(h, t, :) .=
                (energy_cost_allocation_h(h) + demand_cost_allocation_capacity_h(h)) /
                net_demand_h_wo_loss(h) +
                demand_cost_allocation_othercost_h(h) / net_demand_wo_green_tech_h_wo_loss(h)
        end
        fill!(regulator.p_eximport, NaN)
        for t in model_data.index_t
            regulator.p_eximport(t, :) .=
                (energy_cost_allocation_eximport + demand_cost_allocation_eximport) /
                sum(
                    model_data.omega(t) * utility.eximport_my(reg_year_index, t) for
                    t in model_data.index_t
                )
        end
    elseif regulator_opts.rate_design isa TOU
        fill!(regulator.p, NaN)
        for h in model_data.index_h, t in model_data.index_t
            regulator.p(h, t, :) .=
                energy_cost_allocation_h_t(h, t) / net_demand_h_t_wo_loss(h, t) +
                demand_cost_allocation_capacity_h(h) / net_demand_h_wo_loss(h) +
                demand_cost_allocation_othercost_h(h) / net_demand_wo_green_tech_h_wo_loss(h)
        end
        fill!(regulator.p_eximport, NaN)
        for t in model_data.index_t
            regulator.p_eximport(t, :) .=
                energy_cost_allocation_eximport_t(t) /
                utility.eximport_my(reg_year_index, t) +
                demand_cost_allocation_eximport / sum(
                    model_data.omega(t) * utility.eximport_my(reg_year_index, t) for
                    t in model_data.index_t
                )
        end
    end

    fill!(regulator.p_regression, NaN)
    for h in model_data.index_h
        regulator.p_regression(h, :) .=
            (energy_cost_allocation_h(h) + demand_cost_allocation_capacity_h(h)) /
            net_demand_h_wo_loss(h)
    end

    fill!(regulator.p_td, NaN)
    for h in model_data.index_h
        regulator.p_td(h, :) .=
            demand_cost_allocation_othercost_h(h) / net_demand_wo_green_tech_h_wo_loss(h)
    end

    # TODO: Call a function instead of using if-then
    if regulator_opts.net_metering_policy isa ExcessRetailRate
        regulator.p_ex = ParamArray(regulator.p)
    elseif regulator_opts.net_metering_policy isa ExcessMarginalCost
        fill!(regulator.p_ex, NaN)
        for h in model_data.index_h, t in model_data.index_t
            regulator.p_ex(h, t, :) .=
                utility.miu_my(reg_year_index, t) /
                (model_data.omega(t) * utility.pvf_onm(reg_year_index))
        end
    elseif regulator_opts.net_metering_policy isa ExcessZero
        fill!(regulator.p_ex, 0.0)
    end

    for h in model_data.index_h, t in model_data.index_t
        regulator.p_my(reg_year_index, h, t, :) .= regulator.p(h, t)
        regulator.p_ex_my(reg_year_index, h, t, :) .= regulator.p_ex(h, t)
        regulator.p_eximport_my(reg_year_index, t, :) .= regulator.p_eximport(t)
    end

    for h in model_data.index_h
        regulator.p_my_regression(reg_year_index, h, :) .= regulator.p_regression(h)
        regulator.p_my_td(reg_year_index, h, :) .= regulator.p_td(h)
    end

    p_after_wavg = ParamArray(regulator.p_td, "p_after_wavg")
    fill!(p_after_wavg, NaN)  # TODO DT: debug only
    for h in model_data.index_h
        p_after_wavg(h, :) .= 
            sum(regulator.p_my(reg_year_index, h, t) * model_data.omega(t) * customers.d(h, t) for t in model_data.index_t) / 
            sum(model_data.omega(t) * customers.d(h, t) for t in model_data.index_t)
    end

    @info "Original retail price" p_before
    @info "Original DER excess rate" p_ex_before
    @info "New retail price" regulator.p
    @info "New DER excess rate" regulator.p_ex

    return compute_difference_percentage_one_norm([
        (p_before_wavg.values, p_after_wavg.values),
    ])
end

function solve_agent_problem!(
    regulator::Regulator,
    regulator_opts::RegulatorOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{WholesaleMarket},
    agent_store::AgentStore,
    w_iter,
)

    for y in model_data.index_y_fix
        regulator.othercost(y, :) .= (
            regulator.distribution_cost(y) + 
            regulator.administration_cost(y) + 
            regulator.transmission_cost(y) + 
            regulator.interconnection_cost(y) + 
            regulator.system_cost(y)
        )
    end

    customers = get_agent(CustomerGroup, agent_store)
    ipp = get_agent(IPPGroup, agent_store)
    utility = get_agent(Utility, agent_store)
    green_developer = get_agent(GreenDeveloper, agent_store)

    # the year regulator is making a rate case
    reg_year = model_data.year(first(model_data.index_y))
    reg_year_index = Symbol(Int(reg_year))

    net_demand_w_loss = (
        # demand
        # when it comes to sharing the revenue requirement (cost), use load including distribution loss
        # e.g. utility generation is 100 MW, 50 MW to serve load (including distribution loss),
        #      50 MW for export. It makes sense to allocate the same cost for internal load and export.
        sum(
            customers.gamma(h) * model_data.omega(t) * customers.d(h, t) for
            h in model_data.index_h, t in model_data.index_t
        ) +
        # export
        sum(
            model_data.omega(t) * ipp.eximport_my(reg_year_index, t) for
            t in model_data.index_t
        ) -
        # DG
        # when it comes to cost allocation, simply use total DER generation to offset total load;
        # this does not consider enery offset at the household level, which is only considered when 
        # calculating retail rates.
        sum(
            model_data.omega(t) * (
                customers.rho_DG(h, m, t) *
                customers.x_DG_E_my(first(model_data.index_y), h, m) + sum(
                    customers.rho_DG(h, m, t) * customers.x_DG_new_my(Symbol(Int(y)), h, m) for
                    y in model_data.year(first(model_data.index_y_fix)):reg_year
                )
            ) for t in model_data.index_t, h in model_data.index_h,
            m in customers.index_m
        ) -
        # green technology subscription
        sum(
            model_data.omega(t) * utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
            model_data.year(first(model_data.index_y_fix)):reg_year)
            for t in model_data.index_t, j in model_data.index_j, h in model_data.index_h
        )
    )

    der_excess_cost_h = KeyedArray(
        [
            sum(
                model_data.omega(t) *
                regulator.p_ex(h, t) *
                (
                    max(
                        0,
                        customers.rho_DG(h, m, t) * customers.Opti_DG_E(h, m) -
                        customers.d(h, t) * (1 - utility.loss_dist),
                    ) * customers.x_DG_E_my(first(model_data.index_y), h, m) /
                    customers.Opti_DG_E(h, m) + sum(
                        max(
                            0,
                            customers.rho_DG(h, m, t) *
                            customers.Opti_DG_my(Symbol(Int(y)), h, m) -
                            customers.d(h, t) * (1 - utility.loss_dist),
                        ) * customers.x_DG_new_my(Symbol(Int(y)), h, m) /
                        customers.Opti_DG_my(Symbol(Int(y)), h, m) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    )
                ) for t in model_data.index_t, m in customers.index_m
            ) for h in model_data.index_h
        ];
        [get_pair(model_data.index_h)]...
    )

    net_demand_h_w_loss = KeyedArray(
        # demand
        [
            sum(
                customers.gamma(h) * model_data.omega(t) * customers.d(h, t) for
                t in model_data.index_t
            ) -
            # DG
            sum(
                model_data.omega(t) * (
                    customers.rho_DG(h, m, t) *
                    customers.x_DG_E_my(first(model_data.index_y), h, m) + sum(
                        customers.rho_DG(h, m, t) *
                        customers.x_DG_new_my(Symbol(Int(y)), h, m) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    )
                ) for t in model_data.index_t, m in customers.index_m
            ) -
            # green technology subscription
            sum(
                model_data.omega(t) * utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):reg_year)
                for t in model_data.index_t, j in model_data.index_j
            )
            for h in model_data.index_h
        ];
        [get_pair(model_data.index_h)]...
    )

    net_demand_h_wo_loss = KeyedArray(
        # demand
        [
            sum(
                customers.gamma(h) *
                model_data.omega(t) *
                customers.d(h, t) *
                (1 - utility.loss_dist) for t in model_data.index_t
            ) -
            # DG
            # since this net demand is for rate calculation, we need to consider energy offset at the household level.
            # e.g. two household, with 100 MW load each (without loss), if one of them has DER and generated 
            # 120 MW of energy, he/she does not need to pay for energy, but the other one still have to pay 
            # 100 MW instead of 80 MW.
            sum(
                model_data.omega(t) * (
                    min(
                        customers.rho_DG(h, m, t) * customers.Opti_DG_E(h, m),
                        customers.d(h, t) * (1 - utility.loss_dist),
                    ) * customers.x_DG_E_my(first(model_data.index_y), h, m) /
                    customers.Opti_DG_E(h, m) + sum(
                        min(
                            customers.rho_DG(h, m, t) *
                            customers.Opti_DG_my(Symbol(Int(y)), h, m),
                            customers.d(h, t) * (1 - utility.loss_dist),
                        ) * customers.x_DG_new_my(Symbol(Int(y)), h, m) /
                        customers.Opti_DG_my(Symbol(Int(y)), h, m) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    )
                ) for t in model_data.index_t, m in customers.index_m
            ) -
            # green technology subscription
            sum(
                model_data.omega(t) * utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):reg_year)
                for t in model_data.index_t, j in model_data.index_j
            )
            for h in model_data.index_h
        ];
        [get_pair(model_data.index_h)]...
    )

    net_demand_wo_green_tech_h_wo_loss = KeyedArray(
        [
        # demand
            sum(
                customers.gamma(h) *
                model_data.omega(t) *
                customers.d(h, t) *
                (1 - utility.loss_dist) for t in model_data.index_t
            ) -
            # DG
            # since this net demand is for rate calculation, we need to consider energy offset at the household level.
            # e.g. two household, with 100 MW load each (without loss), if one of them has DER and generated 
            # 120 MW of energy, he/she does not need to pay for energy, but the other one still have to pay 
            # 100 MW instead of 80 MW.
            sum(
                model_data.omega(t) * (
                    min(
                        customers.rho_DG(h, m, t) * customers.Opti_DG_E(h, m),
                        customers.d(h, t) * (1 - utility.loss_dist),
                    ) * customers.x_DG_E_my(first(model_data.index_y), h, m) /
                    customers.Opti_DG_E(h, m) + sum(
                        min(
                            customers.rho_DG(h, m, t) *
                            customers.Opti_DG_my(Symbol(Int(y)), h, m),
                            customers.d(h, t) * (1 - utility.loss_dist),
                        ) * customers.x_DG_new_my(Symbol(Int(y)), h, m) /
                        customers.Opti_DG_my(Symbol(Int(y)), h, m) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    )
                ) for t in model_data.index_t, m in customers.index_m
            )
            for h in model_data.index_h
        ];
        [get_pair(model_data.index_h)]...
    )

    net_demand_t_w_loss = KeyedArray(
        # demand
        [
            sum(customers.gamma(h) * customers.d(h, t) for h in model_data.index_h) +
            # export
            utility.eximport_my(reg_year_index, t) -
            # DG
            sum(
                customers.rho_DG(h, m, t) *
                customers.x_DG_E_my(first(model_data.index_y), h, m) + sum(
                    customers.rho_DG(h, m, t) * customers.x_DG_new_my(Symbol(Int(y)), h, m) for
                    y in model_data.year(first(model_data.index_y_fix)):reg_year
                ) for h in model_data.index_h, m in customers.index_m
            ) -
            # green technology subscription
            sum(
                utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):reg_year)
                for h in model_data.index_h, j in model_data.index_j
            )
            for t in model_data.index_t
        ];
        [get_pair(model_data.index_t)]...
    )

    der_excess_cost_h_t = make_keyed_array(model_data.index_h, model_data.index_t)
    for h in model_data.index_h, t in model_data.index_t
        der_excess_cost_h_t(h, t, :) .= sum(
            regulator.p_ex(h, t) * (
                max(
                    0,
                    customers.rho_DG(h, m, t) * customers.Opti_DG_E(h, m) -
                    customers.d(h, t) * (1 - utility.loss_dist),
                ) * customers.x_DG_E_my(first(model_data.index_y), h, m) /
                customers.Opti_DG_E(h, m) + sum(
                    max(
                        0,
                        customers.rho_DG(h, m, t) *
                        customers.Opti_DG_my(Symbol(Int(y)), h, m) -
                        customers.d(h, t) * (1 - utility.loss_dist),
                    ) * customers.x_DG_new_my(Symbol(Int(y)), h, m) /
                    customers.Opti_DG_my(Symbol(Int(y)), h, m) for
                    y in model_data.year(first(model_data.index_y_fix)):reg_year
                )
            ) for m in customers.index_m
        )
    end

    net_demand_h_t_w_loss = make_keyed_array(model_data.index_h, model_data.index_t)
    for h in model_data.index_h, t in model_data.index_t
        net_demand_h_t_w_loss(h, t, :) .=
        # demand
            customers.gamma(h) * customers.d(h, t) -
            # DG
            sum(
                customers.rho_DG(h, m, t) *
                customers.x_DG_E_my(first(model_data.index_y), h, m) + sum(
                    customers.rho_DG(h, m, t) * customers.x_DG_new_my(Symbol(Int(y)), h, m)
                    for y in model_data.year(first(model_data.index_y_fix)):reg_year
                ) for m in customers.index_m
            ) -
            sum(
                utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):reg_year)
                for j in model_data.index_j
            )
    end

    net_demand_h_t_wo_loss = make_keyed_array(model_data.index_h, model_data.index_t)
    for h in model_data.index_h, t in model_data.index_t
        net_demand_h_t_wo_loss(h, t, :) .=
        # demand
            customers.gamma(h) * customers.d(h, t) * (1 - utility.loss_dist) -
            # DG
            sum(
                min(
                    customers.rho_DG(h, m, t) * customers.Opti_DG_E(h, m),
                    customers.d(h, t) * (1 - utility.loss_dist),
                ) * customers.x_DG_E_my(first(model_data.index_y), h, m) /
                customers.Opti_DG_E(h, m) + sum(
                    min(
                        customers.rho_DG(h, m, t) *
                        customers.Opti_DG_my(Symbol(Int(y)), h, m),
                        customers.d(h, t) * (1 - utility.loss_dist),
                    ) * customers.x_DG_new_my(Symbol(Int(y)), h, m) /
                    customers.Opti_DG_my(Symbol(Int(y)), h, m) for
                    y in model_data.year(first(model_data.index_y_fix)):reg_year
                ) for m in customers.index_m
            ) -
            sum(
                utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):reg_year)
                for j in model_data.index_j
            )
    end

    # for the purpose of calculating net peak load, use load including distribution loss
    net_demand_for_peak_h_t = make_keyed_array(model_data.index_h, model_data.index_t)
    for h in model_data.index_h, t in model_data.index_t
        net_demand_for_peak_h_t(h, t, :) .=
            customers.gamma(h) * customers.d(h, t) - sum(
                customers.rho_DG(h, m, t) *
                customers.x_DG_E_my(first(model_data.index_y), h, m) + sum(
                    customers.rho_DG(h, m, t) * customers.x_DG_new_my(Symbol(Int(y)), h, m)
                    for y in model_data.year(first(model_data.index_y_fix)):reg_year
                ) for m in customers.index_m
            ) -
            sum(
                utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):reg_year)
                for j in model_data.index_j
            )
    end

    # excluding green tech generation offset when it comes to sharing T&D costs
    net_demand_for_peak_wo_green_tech_h_t = make_keyed_array(model_data.index_h, model_data.index_t)
    for h in model_data.index_h, t in model_data.index_t
        net_demand_for_peak_wo_green_tech_h_t(h, t, :) .=
                customers.gamma(h) * customers.d(h, t) - sum(
                    customers.rho_DG(h, m, t) *
                    customers.x_DG_E_my(first(model_data.index_y), h, m) + sum(
                        customers.rho_DG(h, m, t) * customers.x_DG_new_my(Symbol(Int(y)), h, m) for
                        y in model_data.year(first(model_data.index_y_fix)):reg_year
                    ) for m in customers.index_m
                )
    end

    net_peak_load_h = KeyedArray(
        [
            findmax(Dict(t => net_demand_for_peak_h_t(h, t) for t in model_data.index_t))[1] for h in model_data.index_h
        ];
        [get_pair(model_data.index_h)]...
    )

    net_peak_load_wo_green_tech_h = KeyedArray(
        [ 
            findmax(Dict(t => net_demand_for_peak_wo_green_tech_h_t(h, t) for t in model_data.index_t))[1] for h in model_data.index_h
        ];
        [get_pair(model_data.index_h)]...
    )

    # Cost Classification/Allocation
    energy_purchase_cost =
        sum(
            ipp.miu_my(reg_year_index, t) * (
                sum(
                    ipp.y_E_my(reg_year_index, p, k, t) for k in ipp.index_k_existing,
                    p in ipp.index_p
                ) + sum(
                    ipp.y_C_my(reg_year_index, p, k, t) for k in ipp.index_k_new,
                    p in ipp.index_p
                )
            ) for t in model_data.index_t
        ) +
        # incorporate REC cost into energy purchase cost
        sum(
            regulator.REC *
            model_data.omega(t) *
            (
                sum(
                    ipp.y_E_my(reg_year_index, p, rps, t) for rps in ipp.index_rps,
                    p in ipp.index_p
                ) + sum(
                    ipp.y_C_my(reg_year_index, p, rps, t) for rps in ipp.index_rps,
                    p in ipp.index_p
                )
            ) for t in model_data.index_t
        )
    energy_purchase_cost_t = KeyedArray(
        [
            ipp.miu_my(reg_year_index, t) / model_data.omega(t) * (
                sum(
                    ipp.y_E_my(reg_year_index, p, k, t) for k in ipp.index_k_existing,
                    p in ipp.index_p
                ) + sum(
                    ipp.y_C_my(reg_year_index, p, k, t) for k in ipp.index_k_new,
                    p in ipp.index_p
                )
            ) +
            # incorporate REC cost into energy purchase cost
            regulator.REC * (
                sum(
                    ipp.y_E_my(reg_year_index, p, rps, t) for rps in ipp.index_rps,
                    p in ipp.index_p
                ) + sum(
                    ipp.y_C_my(reg_year_index, p, rps, t) for rps in ipp.index_rps,
                    p in ipp.index_p
                )
            ) for t in model_data.index_t
        ];
        [get_pair(model_data.index_t)]...
    )

    capacity_purchase_cost = sum(
        ipp.ucap(reg_year_index, p) * (
            ipp.Capacity_intercept_my(reg_year_index) +
            ipp.Capacity_slope_my(reg_year_index) *
            sum(ipp.ucap(reg_year_index, p) for p in ipp.index_p)
        ) for p in ipp.index_p
    )

    energy_cost_allocation_h = KeyedArray(
        [
            energy_purchase_cost * net_demand_h_w_loss(h) / net_demand_w_loss +
            der_excess_cost_h(h) for h in model_data.index_h
        ];
        [get_pair(model_data.index_h)]...
    )
    energy_cost_allocation_eximport =
        energy_purchase_cost * sum(
            model_data.omega(t) * utility.eximport_my(reg_year_index, t) for
            t in model_data.index_t
        ) / net_demand_w_loss

    # allocate capacity_purchase_cost by net peak load with green tech generation offset;
    # allocate regulator.othercost (T&D costs) by net peak load without green tech generation offset;
    demand_cost_allocation_h = KeyedArray(
        [
            capacity_purchase_cost * net_peak_load_h(h) / (
                sum(net_peak_load_h(h) for h in model_data.index_h) +
                utility.Peak_eximport_my(reg_year_index)
            ) +
            regulator.othercost(reg_year_index) * net_peak_load_wo_green_tech_h(h) / (
                sum(net_peak_load_wo_green_tech_h(h) for h in model_data.index_h) +
                utility.Peak_eximport_my(reg_year_index)
            )
            for h in model_data.index_h
        ];
        [get_pair(model_data.index_h)]...
    )
    demand_cost_allocation_capacity_h = KeyedArray(
        [
            capacity_purchase_cost * net_peak_load_h(h) / (
                sum(net_peak_load_h(h) for h in model_data.index_h) +
                utility.Peak_eximport_my(reg_year_index)
            )
            for h in model_data.index_h
        ];
        [get_pair(model_data.index_h)]...
    )
    demand_cost_allocation_othercost_h = KeyedArray(
        [
            regulator.othercost(reg_year_index) * net_peak_load_wo_green_tech_h(h) / (
                sum(net_peak_load_wo_green_tech_h(h) for h in model_data.index_h) +
                utility.Peak_eximport_my(reg_year_index)
            )
            for h in model_data.index_h
        ];
        [get_pair(model_data.index_h)]...
    )

    demand_cost_allocation_eximport =
        capacity_purchase_cost *
        utility.Peak_eximport_my(reg_year_index) / (
            sum(net_peak_load_h(h) for h in model_data.index_h) +
            utility.Peak_eximport_my(reg_year_index)
        ) +
        regulator.othercost(reg_year_index) *
        utility.Peak_eximport_my(reg_year_index) / (
            sum(net_peak_load_wo_green_tech_h(h) for h in model_data.index_h) +
            utility.Peak_eximport_my(reg_year_index)
        )
    energy_cost_allocation_h_t = make_keyed_array(model_data.index_h, model_data.index_t)
    for h in model_data.index_h, t in model_data.index_t
        energy_cost_allocation_h_t(h, t, :) .=
            energy_purchase_cost_t(t) * net_demand_h_t_w_loss(h, t) /
            net_demand_t_w_loss(t) + der_excess_cost_h_t(h, t)
    end

    energy_cost_allocation_eximport_t = KeyedArray(
        [
            energy_purchase_cost_t(t) * utility.eximport_my(reg_year_index, t) /
            net_demand_t_w_loss(t) for t in model_data.index_t
        ];
        [get_pair(model_data.index_t)]...
    )

    p_before = ParamArray(regulator.p, "p_before")
    fill!(p_before, NaN)
    for h in model_data.index_h, t in model_data.index_t
        p_before(h, t, :) .= regulator.p_my(reg_year_index, h, t)
    end

    p_before_wavg = ParamArray(regulator.p_td, "p_before_wavg")
    fill!(p_before_wavg, NaN)  # TODO DT: debug only
    for h in model_data.index_h
        p_before_wavg(h, :) .= 
            sum(regulator.p_my(reg_year_index, h, t) * model_data.omega(t) * customers.d(h, t) for t in model_data.index_t) / 
            sum(model_data.omega(t) * customers.d(h, t) for t in model_data.index_t)
    end

    p_ex_before = ParamArray(regulator.p_ex, "p_ex_before")
    fill!(p_ex_before, NaN)
    for h in model_data.index_h, t in model_data.index_t
        p_ex_before(h, t, :) .= regulator.p_ex_my(reg_year_index, h, t)
    end

    # TODO: Call a function instead of using if-then
    if regulator_opts.rate_design isa FlatRate
        fill!(regulator.p, NaN)
        for h in model_data.index_h, t in model_data.index_t
            regulator.p(h, t, :) .=
                (energy_cost_allocation_h(h) + demand_cost_allocation_capacity_h(h)) /
                net_demand_h_wo_loss(h) +
                demand_cost_allocation_othercost_h(h) / net_demand_wo_green_tech_h_wo_loss(h)
        end
        fill!(regulator.p_eximport, NaN)
        for t in model_data.index_t
            regulator.p_eximport(t, :) .=
                (energy_cost_allocation_eximport + demand_cost_allocation_eximport) /
                sum(
                    model_data.omega(t) * utility.eximport_my(reg_year_index, t) for
                    t in model_data.index_t
                )
        end
    elseif regulator_opts.rate_design isa TOU
        fill!(regulator.p, NaN)
        for h in model_data.index_h, t in model_data.index_t
            regulator.p(h, t, :) .=
                energy_cost_allocation_h_t(h, t) / net_demand_h_t_wo_loss(h, t) +
                demand_cost_allocation_capacity_h(h) / net_demand_h_wo_loss(h) +
                demand_cost_allocation_othercost_h(h) / net_demand_wo_green_tech_h_wo_loss(h)
        end
        fill!(regulator.p_eximport, NaN)
        for t in model_data.index_t
            regulator.p_eximport(t, :) .=
                energy_cost_allocation_eximport_t(t) /
                utility.eximport_my(reg_year_index, t) +
                demand_cost_allocation_eximport / sum(
                    model_data.omega(t) * utility.eximport_my(reg_year_index, t) for
                    t in model_data.index_t
                )
        end
    end

    fill!(regulator.p_regression, NaN)
    for h in model_data.index_h
        regulator.p_regression(h, :) .=
            (energy_cost_allocation_h(h) + demand_cost_allocation_capacity_h(h)) /
            net_demand_h_wo_loss(h)
    end

    fill!(regulator.p_td, NaN)
    for h in model_data.index_h
        regulator.p_td(h, :) .=
            demand_cost_allocation_othercost_h(h) / net_demand_wo_green_tech_h_wo_loss(h)
    end

    # TODO: Call a function instead of using if-then
    if regulator_opts.net_metering_policy isa ExcessRetailRate
        regulator.p_ex = ParamArray(regulator.p)
    elseif regulator_opts.net_metering_policy isa ExcessMarginalCost
        fill!(regulator.p_ex, NaN)
        for h in model_data.index_h, t in model_data.index_t
            regulator.p_ex(h, t, :) .= ipp.miu_my(reg_year_index, t) / model_data.omega(t)
        end
    elseif regulator_opts.net_metering_policy isa ExcessZero
        fill!(regulator.p_ex, 0.0)
    end

    for h in model_data.index_h, t in model_data.index_t
        regulator.p_my(reg_year_index, h, t, :) .= regulator.p(h, t)
        regulator.p_ex_my(reg_year_index, h, t, :) .= regulator.p_ex(h, t)
        regulator.p_eximport_my(reg_year_index, t, :) .= regulator.p_eximport(t)
    end

    for h in model_data.index_h
        regulator.p_my_regression(reg_year_index, h, :) .= regulator.p_regression(h)
        regulator.p_my_td(reg_year_index, h, :) .= regulator.p_td(h)
    end

    p_after_wavg = ParamArray(regulator.p_td, "p_after_wavg")
    fill!(p_after_wavg, NaN)  # TODO DT: debug only
    for h in model_data.index_h
        p_after_wavg(h, :) .= 
            sum(regulator.p_my(reg_year_index, h, t) * model_data.omega(t) * customers.d(h, t) for t in model_data.index_t) / 
            sum(model_data.omega(t) * customers.d(h, t) for t in model_data.index_t)
    end

    @info "Original retail price" p_before
    @info "Original DER excess rate" p_ex_before
    @info "New retail price" regulator.p
    @info "New DER excess rate" regulator.p_ex

    return compute_difference_percentage_one_norm([
        (p_before_wavg.values, p_after_wavg.values),
    ])
end

function save_results(
    regulator::Regulator,
    regulator_opts::RegulatorOptions,
    hem_opts::HEMOptions{<:MarketStructure},
    export_file_path::AbstractString,
    fileprefix::AbstractString,
)
    # Primal Variables
    save_param(
        regulator.p_my.values,
        [:Year, :CustomerType, :Time],
        :Price,
        joinpath(export_file_path, "$(fileprefix)_p.csv"),
    )
    save_param(
        regulator.p_ex_my.values,
        [:Year, :CustomerType, :Time],
        :Price,
        joinpath(export_file_path, "$(fileprefix)_p_ex.csv"),
    )
end
