# This file defines data and functions associated with the customer.

DER_factor = 1.0    # A scaling factor applied to existing DER penetration level

mutable struct PVAdoptionModel
    Shape::ParamAxisArray
    MeanPayback::ParamAxisArray
    Bass_p::ParamAxisArray
    Bass_q::ParamAxisArray
    Rate::ParamAxisArray
end

mutable struct GreenSubModel
    Constant::ParamAxisArray
    GreenPowerPrice_coefficient::ParamAxisArray
    EnergyRate_coefficient::ParamAxisArray
    WholesaleMarket_coefficient::ParamAxisArray
    RetailCompetition_coefficient::ParamAxisArray
    RPS_coefficient::ParamAxisArray
    WTP_coefficient::ParamAxisArray
end

"""
Constructs a PVAdoptionModel by computing Rate from the other parameters.
"""
function PVAdoptionModel(Shape, MeanPayback, Bass_p, Bass_q)
    @assert_op length(Shape) == length(MeanPayback)
    vals = deepcopy(Shape.values)
    for i in 1:length(vals)
        vals[i] = Shape[i] / MeanPayback[i]
    end

    return PVAdoptionModel(
        Shape,
        MeanPayback,
        Bass_p,
        Bass_q,
        ParamAxisArray("Rate", Shape.dims, vals), # Rate
    )
end

# function GreenSubModel(Constant, GreenPowerPrice_coefficient, EnergyRate_coefficient, WholesaleMarket_coefficient, RetailCompetition_coefficient, RPS_coefficient, WTP_coefficient)
#     return GreenSubModel(
#         Constant,
#         GreenPowerPrice_coefficient,
#         EnergyRate_coefficient,
#         WholesaleMarket_coefficient,
#         RetailCompetition_coefficient,
#         RPS_coefficient,
#         WTP_coefficient,
#     )
# end

# # declare customer decision
# abstract type ConsumerModel end
# struct DERAdoption <: ConsumerModel end
# struct SupplyChoice <: ConsumerModel end

# abstract type AbstractCustomerOptions <: AgentOptions end

# struct CustomerOptions{T <: ConsumerModel} <: AbstractCustomerOptions
#     customer_model::T
# end

abstract type AbstractCustomerGroup <: AgentGroup end

mutable struct CustomerGroup <: AbstractCustomerGroup
    id::String
    # Sets
    index_m::Dimension # behind-the-meter technologies

    # Parameters
    "number of customers of type h"
    gamma::ParamAxisArray
    "demand (MWh per representative agent per hour)"
    d::ParamAxisArray
    "multi-year demand (MWh per representative agent per hour)"
    d_my::ParamAxisArray
    x_DG_E::ParamAxisArray
    "Existing DER at year y. This is a cumulative number but without x_DG_new_my built by this module"
    x_DG_E_my::ParamAxisArray
    Opti_DG::ParamAxisArray
    Opti_DG_E::ParamAxisArray
    Opti_DG_my::ParamAxisArray
    # "DER generation by a representative customer h and DER technology m"
    # DERGen::ParamAxisArray
    CapEx_DG::ParamAxisArray
    CapEx_DG_my::ParamAxisArray
    FOM_DG::ParamAxisArray
    FOM_DG_my::ParamAxisArray
    rho_DG::ParamAxisArray
    "Annualization factor for net consumer surplus of PV installation"
    delta::ParamScalar
    PeakLoad::ParamAxisArray
    PeakLoad_my::ParamAxisArray

    # Primal Variables
    x_DG_new::ParamAxisArray
    x_DG_new_my::ParamAxisArray    # Annual new DER build (not cumulative)
    x_green_sub::ParamAxisArray
    x_green_sub_my::ParamAxisArray
    x_green_sub_incremental_my::ParamAxisArray

    # Auxiliary Variables
    Payback::ParamAxisArray
    MarketShare::ParamAxisArray
    MaxDG::ParamAxisArray
    F::ParamAxisArray
    year::ParamAxisArray
    A::ParamAxisArray
    ConPVNetSurplus::ParamAxisArray
    ConPVNetSurplus_my::ParamAxisArray

    GreenTechIntercept::ParamAxisArray
    GreenTechSlope::ParamAxisArray
    pv_adoption_model::PVAdoptionModel
    green_sub_model::GreenSubModel

    pvf::Any
    rooftop::Any
    MaxDG_my::ParamAxisArray

    RetailCompetition::ParamAxisArray
    WTP_green_power::ParamAxisArray

    ConGreenPowerNetSurplus_pre_proportion_my::ParamAxisArray
    ConGreenPowerNetSurplus_post_proportion_my::ParamAxisArray
    ConGreenPowerNetSurplus_cumu_my::ParamAxisArray
end

function CustomerGroup(input_filename::AbstractString, model_data::HEMData; id = DEFAULT_ID)
    index_m = read_set(
        input_filename,
        "index_m",
        "index_m",
        prose_name = "behind-the-meter technologies m",
    )

    gamma = read_param(
        "gamma",
        input_filename,
        "Gamma",
        model_data.index_h,
        description = "number of customers of type h",
    )
    demand =
        read_param("d", input_filename, "Demand", model_data.index_t, [model_data.index_h])
    demand_my = read_param(
        "d_my",
        input_filename,
        "Demandmy",
        model_data.index_t,
        [model_data.index_y, model_data.index_h],
    )
    x_DG_E = read_param(
        "x_DG_E",
        input_filename,
        "ExistingDER",
        index_m,
        [model_data.index_h],
        description = "existing DG capacity",
    )
    for h in model_data.index_h, m in index_m
        x_DG_E[h, m] = x_DG_E[h, m] * DER_factor
    end
    x_DG_E_my = read_param(
        "x_DG_E_my",
        input_filename,
        "ExistingDERmy",
        index_m,
        [model_data.index_y, model_data.index_h],
    )
    Opti_DG =
        read_param("Opti_DG", input_filename, "OptimalDER", index_m, [model_data.index_h])
    Opti_DG_my = read_param(
        "Opti_DG_my",
        input_filename,
        "OptimalDERmy",
        index_m,
        [model_data.index_y, model_data.index_h],
    )
    rho_DG = read_param(
        "rho_DG",
        input_filename,
        "AvailabilityDER",
        model_data.index_t,
        [model_data.index_h, index_m],
    )
    # # Define total DER generation per individual customer per hour
    # DERGen = initialize_param("DERGen", model_data.index_h, model_data.index_t, value = 1.0)
    # for h in model_data.index_h, t in model_data.index_t
    #     if sum(rho_DG[h, m, t] * Opti_DG[h, m] for m in index_m) != 0.0
    #         DERGen[h, t] = sum(rho_DG[h, m, t] * Opti_DG[h, m] for m in index_m)
    #     else
    #         DERGen[h, t] = 1.0
    #     end
    # end
    # Calculate maximum demand for each customer type
    MaxLoad = AxisArray(
        [
            gamma[h] * findmax(Dict(t => demand[h, t] for t in model_data.index_t))[1]
            for h in model_data.index_h
        ],
        model_data.index_h.elements,
    )
    MaxLoad_my = make_axis_array(model_data.index_y, model_data.index_h)
    for y in model_data.index_y, h in model_data.index_h
        MaxLoad_my[y, h] =
            gamma[h] * findmax(Dict(t => demand_my[y, h, t] for t in model_data.index_t))[1]
    end

    pv_adoption_model = PVAdoptionModel(
        initialize_param("Shape", model_data.index_h, index_m, value = 1.7), # Shape
        initialize_param("MeanPayback", model_data.index_h, index_m, value = 8.8), # MeanPayback
        ParamAxisArray(
            "Bass_p",
            (model_data.index_h,),
            AxisArray(
                [7.7E-07, 6.0E-04, 6.0E-04],
                [:Residential, :Commercial, :Industrial],
            ),
        ),
        ParamAxisArray(
            "Bass_q",
            (model_data.index_h,),
            AxisArray([0.663, 0.133, 0.133], [:Residential, :Commercial, :Industrial]),
        ),
    )

    green_sub_model = GreenSubModel(
        ParamAxisArray(
            "Constant",
            (model_data.index_h,),
            AxisArray([0.0, 0.0, 0.0], [:Residential, :Commercial, :Industrial]),
            description = "Constant in green power uptake function (regression parameter)",
        ),
        ParamAxisArray(
            "GreenPowerPrice_coefficient",
            (model_data.index_h,),
            AxisArray([0.0, -0.55, -0.55], [:Residential, :Commercial, :Industrial]),
            description = "Sum of PPA and REC prices (regression parameter)",
        ),
        ParamAxisArray(
            "EnergyRate_coefficient",
            (model_data.index_h,),
            AxisArray([0.0, 0.0, 0.0], [:Residential, :Commercial, :Industrial]),
            description = "Weighted mean C&I volumetric (\$/MWh) rate (regression parameter)",
        ),
        ParamAxisArray(
            "WholesaleMarket_coefficient",
            (model_data.index_h,),
            AxisArray([0.0, 0.14, 0.14], [:Residential, :Commercial, :Industrial]),
            description = "% of load served by an ISO (regression parameter)",
        ),
        ParamAxisArray(
            "RetailCompetition_coefficient",
            (model_data.index_h,),
            AxisArray([0.0, 0.16, 0.16], [:Residential, :Commercial, :Industrial]),
            description = "% of C&I customers that are eligible for retail choice (regression parameter)",
        ),
        ParamAxisArray(
            "RPS_coefficient",
            (model_data.index_h,),
            AxisArray([0.0, 0.42, 0.42], [:Residential, :Commercial, :Industrial]),
            description = "RPS percentage requirement in 2019 (regression parameter)",
        ),
        ParamAxisArray(
            "WTP_coefficient",
            (model_data.index_h,),
            AxisArray([0.0, 0.0, 0.0], [:Residential, :Commercial, :Industrial]),
            description = "% of customers willing to pay for renewable energy at the state level (regression parameter)",
        ),
    )

    # Customer financing
    debt_ratio =
        read_param("debt_ratio", input_filename, "CustomerDebtRatio", model_data.index_h)
    cost_of_debt =
        read_param("cost_of_debt", input_filename, "CustomerCOD", model_data.index_h)
    cost_of_equity =
        read_param("cost_of_equity", input_filename, "CustomerCOE", model_data.index_h)
    tax_rate = read_param("tax_rate", input_filename, "CustomerTax", model_data.index_h)
    atwacc = AxisArray(
        [
            debt_ratio[h] * cost_of_debt[h] * (1 - tax_rate[h]) +
            (1 - debt_ratio[h]) * cost_of_equity[h] for h in model_data.index_h
        ],
        model_data.index_h.elements,
    )
    CRF = AxisArray(
        [
            atwacc[h] * (1 + atwacc[h])^20 / ((1 + atwacc[h])^20 - 1) for
            h in model_data.index_h
        ],
        model_data.index_h.elements,
    )
    pvf = AxisArray([1 / CRF[h] for h in model_data.index_h], model_data.index_h.elements)

    return CustomerGroup(
        id,
        index_m,
        gamma,
        demand,
        demand_my,
        x_DG_E,
        x_DG_E_my,
        Opti_DG,
        Opti_DG,
        Opti_DG_my,
        # DERGen,
        read_param("CapEx_DG", input_filename, "CapExDER", index_m, [model_data.index_h]),
        read_param(
            "CapEx_DG_my",
            input_filename,
            "CapExDERmy",
            index_m,
            [model_data.index_y, model_data.index_h],
        ),
        read_param("FOM_DG", input_filename, "FOMDER", index_m, [model_data.index_h]),
        read_param(
            "FOM_DG_my",
            input_filename,
            "FOMDERmy",
            index_m,
            [model_data.index_y, model_data.index_h],
        ),
        rho_DG,
        ParamScalar("delta", 0.05),
        ParamAxisArray("PeakLoad", (model_data.index_h,), MaxLoad),
        ParamAxisArray(
            "PeakLoad_my",
            Tuple(push!(copy([model_data.index_y]), model_data.index_h)),
            MaxLoad_my,
        ),
        initialize_param("x_DG_new", model_data.index_h, index_m),
        initialize_param("x_DG_new_my", model_data.index_y, model_data.index_h, index_m),
        initialize_param("x_green_sub", model_data.index_h, value = 10.0),
        initialize_param("x_green_sub_my", model_data.index_y, model_data.index_h, value = 100.0),
        initialize_param("x_green_sub_incremental_my", model_data.index_y, model_data.index_h, value = 0.0),
        initialize_param("Payback", model_data.index_h, index_m),
        initialize_param("MarketShare", model_data.index_h, index_m),
        initialize_param("MaxDG", model_data.index_h, index_m),
        initialize_param("F", model_data.index_h, index_m),
        initialize_param("year", model_data.index_h, index_m),
        initialize_param("A", model_data.index_h, index_m),
        initialize_param("ConPVNetSurplus", model_data.index_h, index_m),
        initialize_param(
            "ConPVNetSurplus_my",
            model_data.index_y,
            model_data.index_h,
            index_m,
        ),
        initialize_param(
            "GreenTechIntercept",
            model_data.index_h,
            model_data.index_j,
            value = 3.5,
        ), # Intercept of green tech demand curve
        initialize_param(
            "GreenTechSlope",
            model_data.index_h,
            model_data.index_j,
            value = -0.07,
        ),
        pv_adoption_model,
        green_sub_model,
        pvf,
        read_param("rooftop", input_filename, "RooftopDER", index_m, [model_data.index_h]),
        initialize_param("MaxDG_my", model_data.index_y, model_data.index_h, index_m),
        read_param("RetailCompetition", input_filename, "RetailCompetition", model_data.index_y),
        read_param("WTP_green_power", input_filename, "WTP", model_data.index_y),
        initialize_param(
            "ConGreenPowerNetSurplus_pre_proportion_my",
            model_data.index_y,
            model_data.index_h,
        ),
        initialize_param(
            "ConGreenPowerNetSurplus_post_proportion_my",
            model_data.index_y,
            model_data.index_h,
        ),
        initialize_param(
            "ConGreenPowerNetSurplus_cumu_my",
            model_data.index_y,
            model_data.index_h,
        ),
    )
end

get_id(x::CustomerGroup) = x.id

function solve_agent_problem!(
    customers::CustomerGroup,
    customers_opts::AgentOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{<:MarketStructure, DERUseCase},
    agent_store::AgentStore,
    w_iter,
)
    regulator = get_agent(Regulator, agent_store)
    utility = get_agent(Utility, agent_store)

    # the year consumer is making DER investment decision
    reg_year = model_data.year[first(model_data.index_y)]
    reg_year_index = Symbol(Int(reg_year))

    x_DG_before = ParamAxisArray(customers.x_DG_new, "x_DG_before")
    fill!(x_DG_before, NaN)
    for h in model_data.index_h, m in customers.index_m
        x_DG_before[h, m] = customers.x_DG_new_my[reg_year_index, h, m]
    end

    adopt_model = customers.pv_adoption_model

    # update all the annual parameters to the solve year (so we don't have to change the majority of the functions)
    for h in model_data.index_h
        customers.PeakLoad[h] = customers.PeakLoad_my[reg_year_index, h]
    end
    for h in model_data.index_h, t in model_data.index_t
        customers.d[h, t] = customers.d_my[reg_year_index, h, t]
        # customers.DERGen[h,t] = customers.DERGen_my[reg_year_index,h,t]
    end
    for h in model_data.index_h, m in customers.index_m
        customers.Opti_DG[h, m] = customers.Opti_DG_my[reg_year_index, h, m]
        customers.FOM_DG[h, m] = customers.FOM_DG_my[reg_year_index, h, m]
        customers.CapEx_DG[h, m] = customers.CapEx_DG_my[reg_year_index, h, m]
        if w_iter >= 2
            customers.x_DG_E[h, m] =
                customers.x_DG_E_my[reg_year_index, h, m] + sum(
                    customers.x_DG_new_my[Symbol(Int(y)), h, m] for
                    y in model_data.year[first(model_data.index_y_fix)]:(reg_year - 1)
                )
        else
            customers.x_DG_E[h, m] = customers.x_DG_E_my[reg_year_index, h, m]
        end
    end

    # Calculate payback period of DER
    # The NetProfit represents the energy saving/credit per representative agent per DER technology, assuming the optimal DER technology size
    NetProfit = make_axis_array(model_data.index_h, customers.index_m)
    for h in model_data.index_h, m in customers.index_m
        # value of distributed generation (offset load)
        NetProfit[h, m] =
            sum(
                model_data.omega[t] *
                regulator.p[h, t] *
                min(
                    customers.d[h, t] * (1 - utility.loss_dist),
                    customers.rho_DG[h, m, t] * customers.Opti_DG[h, m],
                ) for t in model_data.index_t
            ) +
            # value of distributed generation (excess generation)
            sum(
                model_data.omega[t] *
                regulator.p_ex[h, t] *
                max(
                    0,
                    customers.rho_DG[h, m, t] * customers.Opti_DG[h, m] -
                    customers.d[h, t] * (1 - utility.loss_dist),
                ) for t in model_data.index_t
            ) -
            # cost of distributed generation 
            customers.FOM_DG[h, m] * customers.Opti_DG[h, m]
    end

    for h in model_data.index_h, m in customers.index_m
        if NetProfit[h, m] >= 0.0
            customers.Payback[h, m] =
                customers.CapEx_DG[h, m] * customers.Opti_DG[h, m] / NetProfit[h, m]
            # Calculate maximum market share and maximum DG potential (based on WTP curve)
            customers.MarketShare[h, m] =
                1.0 - Distributions.cdf(
                    Distributions.Gamma(
                        adopt_model.Shape[h, m],
                        1 / adopt_model.Rate[h, m],
                    ),
                    customers.Payback[h, m],
                )
            customers.MaxDG[h, m] =
                customers.MarketShare[h, m] * customers.gamma[h] * customers.Opti_DG[h, m]
            # Calculate the percentage of existing DER (per agent type per DER technology) as a fraction of maximum DG potential
            customers.F[h, m] = min(customers.x_DG_E[h, m] / customers.MaxDG[h, m], 1.0)
            # Back out the reference year of DER based on the percentage of existing DER
            customers.year[h, m] =
                -log(
                    (1 - customers.F[h, m]) /
                    (customers.F[h, m] * adopt_model.Bass_q[h] / adopt_model.Bass_p[h] + 1),
                ) / (adopt_model.Bass_p[h] + adopt_model.Bass_q[h])
            # Calculate incremental DG build
            customers.A[h, m] =
                (
                    1.0 - exp(
                        -(adopt_model.Bass_p[h] + adopt_model.Bass_q[h]) *
                        (customers.year[h, m] + 1),
                    )
                ) / (
                    1.0 +
                    (adopt_model.Bass_q[h] / adopt_model.Bass_p[h]) * exp(
                        -(adopt_model.Bass_p[h] + adopt_model.Bass_q[h]) *
                        (customers.year[h, m] + 1),
                    )
                )
            customers.x_DG_new[h, m] =
                max(0.0, customers.A[h, m] * customers.MaxDG[h, m] - customers.x_DG_E[h, m])
        else
            customers.x_DG_new[h, m] = 0.0
        end
    end

    for h in model_data.index_h, m in customers.index_m
        customers.x_DG_new_my[reg_year_index, h, m] = customers.x_DG_new[h, m]
        customers.MaxDG_my[reg_year_index, h, m] = customers.MaxDG[h, m]
    end

    # @info "Original new DG" x_DG_before
    # @info "New new DG" customers.x_DG_new

    return compute_difference_one_norm([(x_DG_before, customers.x_DG_new)])
end

function solve_agent_problem!(
    customers::CustomerGroup,
    customers_opts::AgentOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{<:MarketStructure, SupplyChoiceUseCase},
    agent_store::AgentStore,
    w_iter,
)
    regulator = get_agent(Regulator, agent_store)
    utility = get_agent(Utility, agent_store)
    green_developer = get_agent(GreenDeveloper, agent_store)

    # the year consumer is making green tariff subscription decision
    reg_year = model_data.year[first(model_data.index_y)]
    reg_year_index = Symbol(Int(reg_year))

    x_green_sub_before = ParamAxisArray(customers.x_green_sub, "x_green_sub_before")
    fill!(x_green_sub_before, NaN)
    for h in model_data.index_h
        x_green_sub_before[h] = customers.x_green_sub_my[reg_year_index, h]
    end

    green_sub_model = customers.green_sub_model

    # update all the annual parameters to the solve year (so we don't have to change the majority of the functions)
    for h in model_data.index_h, t in model_data.index_t
        customers.d[h, t] = customers.d_my[reg_year_index, h, t]
    end

    if typeof(hem_opts) == HEMOptions{VerticallyIntegratedUtility, SupplyChoiceUseCase}
        WholesaleMarketPerc = 0.01
    else
        WholesaleMarketPerc = 1.0
    end

    # calculate green tariff subscription (% MWh)
    GreenSubPerc = AxisArray(
        [ 
            exp(
            green_sub_model.Constant[h] + 
            green_sub_model.GreenPowerPrice_coefficient[h] * log(green_developer.ppa_my[reg_year_index, h]) + 
            green_sub_model.EnergyRate_coefficient[h] * log(regulator.p_my_regression[reg_year_index, h]) + 
            green_sub_model.WholesaleMarket_coefficient[h] * log(WholesaleMarketPerc) + 
            green_sub_model.RetailCompetition_coefficient[h] * log(customers.RetailCompetition[reg_year_index]) + 
            green_sub_model.RPS_coefficient[h] * log(utility.RPS[reg_year_index]) + 
            green_sub_model.WTP_coefficient[h] * log(customers.WTP_green_power[reg_year_index])
            ) for h in model_data.index_h
        ],
        model_data.index_h.elements,
    )

    GreenSubPerc[:Residential] = 0.0

    # is GreenSubPerc a percentage of net load? total load? shall we account for distribution loss or not?
    GreenSubMWh = AxisArray(
        [
            sum(GreenSubPerc[h] * 
            (
                customers.d[h, t] * (1 - utility.loss_dist) * model_data.omega[t] * customers.gamma[h] -
                sum(
                    customers.rho_DG[h, m, t] * customers.x_DG_E_my[reg_year_index, h, m] * model_data.omega[t] for
                    m in customers.index_m
                ) -
                sum(
                    customers.rho_DG[h, m, t] * model_data.omega[t] * sum(
                        customers.x_DG_new_my[Symbol(Int(y_symbol)), h, m] for y_symbol in
                        model_data.year[first(model_data.index_y_fix)]:model_data.year[reg_year_index]
                    ) for m in customers.index_m
                )
            ) for t in model_data.index_t)
            for h in model_data.index_h
        ],
        model_data.index_h.elementsm,
    )

    # customers.x_green_sub_my is an annual number (per the regression), however, this number cannot decrease.
    # this is to make sure the subsribed green techs (in previous years) are always paid for.

    for h in model_data.index_h
        if reg_year > model_data.year[first(model_data.index_y_fix)]
            customers.x_green_sub_my[reg_year_index, h] = max(GreenSubMWh[h], customers.x_green_sub_my[Symbol(Int(reg_year-1)), h])
            customers.x_green_sub_incremental_my[reg_year_index, h] = customers.x_green_sub_my[reg_year_index, h] - customers.x_green_sub_my[Symbol(Int(reg_year-1)), h]
        else
            customers.x_green_sub_my[reg_year_index, h] = GreenSubMWh[h]
            customers.x_green_sub_incremental_my[reg_year_index, h] = GreenSubMWh[h]
        end
    end

    return compute_difference_one_norm([(x_green_sub_before, GreenSubMWh)])

end

function solve_agent_problem!(
    customers::CustomerGroup,
    customers_opts::AgentOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{<:MarketStructure, DERSupplyChoiceUseCase},
    agent_store::AgentStore,
    w_iter,
)
    regulator = get_agent(Regulator, agent_store)
    utility = get_agent(Utility, agent_store)
    green_developer = get_agent(GreenDeveloper, agent_store)

    # the year consumer is making green tariff subscription decision
    reg_year = model_data.year[first(model_data.index_y)]
    reg_year_index = Symbol(Int(reg_year))

    x_DG_before = ParamAxisArray(customers.x_DG_new, "x_DG_before")
    fill!(x_DG_before, NaN)
    for h in model_data.index_h, m in customers.index_m
        x_DG_before[h, m] = customers.x_DG_new_my[reg_year_index, h, m]
    end

    adopt_model = customers.pv_adoption_model

    # update all the annual parameters to the solve year (so we don't have to change the majority of the functions)
    for h in model_data.index_h
        customers.PeakLoad[h] = customers.PeakLoad_my[reg_year_index, h]
    end
    for h in model_data.index_h, t in model_data.index_t
        customers.d[h, t] = customers.d_my[reg_year_index, h, t]
        # customers.DERGen[h,t] = customers.DERGen_my[reg_year_index,h,t]
    end
    for h in model_data.index_h, m in customers.index_m
        customers.Opti_DG[h, m] = customers.Opti_DG_my[reg_year_index, h, m]
        customers.FOM_DG[h, m] = customers.FOM_DG_my[reg_year_index, h, m]
        customers.CapEx_DG[h, m] = customers.CapEx_DG_my[reg_year_index, h, m]
        if w_iter >= 2
            customers.x_DG_E[h, m] =
                customers.x_DG_E_my[reg_year_index, h, m] + sum(
                    customers.x_DG_new_my[Symbol(Int(y)), h, m] for
                    y in model_data.year[first(model_data.index_y_fix)]:(reg_year - 1)
                )
        else
            customers.x_DG_E[h, m] = customers.x_DG_E_my[reg_year_index, h, m]
        end
    end

    # Calculate payback period of DER
    # The NetProfit represents the energy saving/credit per representative agent per DER technology, assuming the optimal DER technology size
    NetProfit = make_axis_array(model_data.index_h, customers.index_m)
    for h in model_data.index_h, m in customers.index_m
        NetProfit[h, m] = 
            # value of distributed generation (offset load)
            sum(
                model_data.omega[t] *
                regulator.p[h, t] *
                min(
                    customers.d[h, t] * (1 - utility.loss_dist),
                    customers.rho_DG[h, m, t] * customers.Opti_DG[h, m],
                ) for t in model_data.index_t
            ) +
            # value of distributed generation (excess generation)
            sum(
                model_data.omega[t] *
                regulator.p_ex[h, t] *
                max(
                    0,
                    customers.rho_DG[h, m, t] * customers.Opti_DG[h, m] -
                    customers.d[h, t] * (1 - utility.loss_dist),
                ) for t in model_data.index_t
            ) -
            # cost of distributed generation 
            customers.FOM_DG[h, m] * customers.Opti_DG[h, m] 
    end

    for h in model_data.index_h, m in customers.index_m
        if NetProfit[h, m] >= 0.0
            customers.Payback[h, m] =
                customers.CapEx_DG[h, m] * customers.Opti_DG[h, m] / NetProfit[h, m]
            # Calculate maximum market share and maximum DG potential (based on WTP curve)
            customers.MarketShare[h, m] =
                1.0 - Distributions.cdf(
                    Distributions.Gamma(
                        adopt_model.Shape[h, m],
                        1 / adopt_model.Rate[h, m],
                    ),
                    customers.Payback[h, m],
                )
            customers.MaxDG[h, m] =
                customers.MarketShare[h, m] * customers.gamma[h] * customers.Opti_DG[h, m]
            # Calculate the percentage of existing DER (per agent type per DER technology) as a fraction of maximum DG potential
            customers.F[h, m] = min(customers.x_DG_E[h, m] / customers.MaxDG[h, m], 1.0)
            # Back out the reference year of DER based on the percentage of existing DER
            customers.year[h, m] =
                -log(
                    (1 - customers.F[h, m]) /
                    (customers.F[h, m] * adopt_model.Bass_q[h] / adopt_model.Bass_p[h] + 1),
                ) / (adopt_model.Bass_p[h] + adopt_model.Bass_q[h])
            # Calculate incremental DG build
            customers.A[h, m] =
                (
                    1.0 - exp(
                        -(adopt_model.Bass_p[h] + adopt_model.Bass_q[h]) *
                        (customers.year[h, m] + 1),
                    )
                ) / (
                    1.0 +
                    (adopt_model.Bass_q[h] / adopt_model.Bass_p[h]) * exp(
                        -(adopt_model.Bass_p[h] + adopt_model.Bass_q[h]) *
                        (customers.year[h, m] + 1),
                    )
                )
            customers.x_DG_new[h, m] =
                max(0.0, customers.A[h, m] * customers.MaxDG[h, m] - customers.x_DG_E[h, m])
        else
            customers.x_DG_new[h, m] = 0.0
        end
    end

    for h in model_data.index_h, m in customers.index_m
        customers.x_DG_new_my[reg_year_index, h, m] = customers.x_DG_new[h, m]
        customers.MaxDG_my[reg_year_index, h, m] = customers.MaxDG[h, m]
    end

    # @info "Original new DG" x_DG_before
    # @info "New new DG" customers.x_DG_new

    x_green_sub_before = ParamAxisArray(customers.x_green_sub, "x_green_sub_before")
    fill!(x_green_sub_before, NaN)
    for h in model_data.index_h
        x_green_sub_before[h] = customers.x_green_sub_my[reg_year_index, h]
    end

    green_sub_model = customers.green_sub_model

    if typeof(hem_opts) == HEMOptions{VerticallyIntegratedUtility, DERSupplyChoiceUseCase}
        WholesaleMarketPerc = 0.01
    else
        WholesaleMarketPerc = 1.0
    end

    # calculate green tariff subscription (% MWh)
    GreenSubPerc = AxisArray(
        [ 
            exp(
            green_sub_model.Constant[h] + 
            green_sub_model.GreenPowerPrice_coefficient[h] * log(green_developer.ppa_my[reg_year_index, h]) + 
            green_sub_model.EnergyRate_coefficient[h] * log(regulator.p_my_regression[reg_year_index, h]) + 
            green_sub_model.WholesaleMarket_coefficient[h] * log(WholesaleMarketPerc) + 
            green_sub_model.RetailCompetition_coefficient[h] * log(customers.RetailCompetition[reg_year_index]) + 
            green_sub_model.RPS_coefficient[h] * log(utility.RPS[reg_year_index]) + 
            green_sub_model.WTP_coefficient[h] * log(customers.WTP_green_power[reg_year_index])
            ) for h in model_data.index_h
        ],
        model_data.index_h.elements,
    )

    GreenSubPerc[:Residential] = 0.0

    # shall we use net load here?
    GreenSubMWh = AxisArray(
        [
            sum(GreenSubPerc[h] * 
            (
                customers.d[h, t] * (1 - utility.loss_dist) * model_data.omega[t] * customers.gamma[h] -
                sum(
                    customers.rho_DG[h, m, t] * customers.x_DG_E_my[reg_year_index, h, m] * model_data.omega[t] for
                    m in customers.index_m
                ) -
                sum(
                    customers.rho_DG[h, m, t] * model_data.omega[t] * sum(
                        customers.x_DG_new_my[Symbol(Int(y_symbol)), h, m] for y_symbol in
                        model_data.year[first(model_data.index_y_fix)]:model_data.year[reg_year_index]
                    ) for m in customers.index_m
                )
            ) for t in model_data.index_t)
            for h in model_data.index_h
        ],
        model_data.index_h.elements,
    )

    # customers.x_green_sub_my is an annual number (per the regression), however, this number cannot decrease.
    # this is to make sure the subsribed green techs (in previous years) are always paid for.

    for h in model_data.index_h
        if reg_year > model_data.year[first(model_data.index_y_fix)]
            customers.x_green_sub_my[reg_year_index, h] = max(GreenSubMWh[h], customers.x_green_sub_my[Symbol(Int(reg_year-1)), h])
            customers.x_green_sub_incremental_my[reg_year_index, h] = customers.x_green_sub_my[reg_year_index, h] - customers.x_green_sub_my[Symbol(Int(reg_year-1)), h]
        else
            customers.x_green_sub_my[reg_year_index, h] = GreenSubMWh[h]
            customers.x_green_sub_incremental_my[reg_year_index, h] = GreenSubMWh[h]
        end
    end

    return compute_difference_one_norm([
        (x_green_sub_before, GreenSubMWh), 
        (x_DG_before, customers.x_DG_new)
    ])

end



function save_results(
    customers::CustomerGroup,
    customers_opts::AgentOptions,
    hem_opts::HEMOptions{<:MarketStructure, DERUseCase},
    export_file_path::AbstractString,
    fileprefix::AbstractString,
)

    # Primal Variables
    save_param(
        customers.x_DG_new_my.values,
        [:Year, :CustomerType, :DERTech],
        :Capacity_MW,
        joinpath(export_file_path, "$(fileprefix)_x_DG.csv"),
    )
end

function save_results(
    customers::CustomerGroup,
    customers_opts::AgentOptions,
    hem_opts::HEMOptions{<:MarketStructure, SupplyChoiceUseCase},
    export_file_path::AbstractString,
    fileprefix::AbstractString,
)

    # Primal Variables
    save_param(
        customers.x_green_sub_my.values,
        [:Year, :CustomerType],
        :Subscription_MWh,
        joinpath(export_file_path, "$(fileprefix)_x_green_sub.csv"),
    )
end

function save_results(
    customers::CustomerGroup,
    customers_opts::AgentOptions,
    hem_opts::HEMOptions{<:MarketStructure, DERSupplyChoiceUseCase},
    export_file_path::AbstractString,
    fileprefix::AbstractString,
)

    # Primal Variables
    save_param(
        customers.x_DG_new_my.values,
        [:Year, :CustomerType, :DERTech],
        :Capacity_MW,
        joinpath(export_file_path, "$(fileprefix)_x_DG.csv"),
    )

    save_param(
        customers.x_green_sub_my.values,
        [:Year, :CustomerType],
        :Subscription_MWh,
        joinpath(export_file_path, "$(fileprefix)_x_green_sub.csv"),
    )
end

function welfare_calculation!(
    customers::CustomerGroup,
    customers_opts::AgentOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{<:MarketStructure, DERUseCase},
    agent_store::AgentStore,
)
    adopt_model = customers.pv_adoption_model
    regulator = get_agent(Regulator, agent_store)
    utility = get_agent(Utility, agent_store)

    """
    Net Consumer Surplus Calculation:

    + Annualized Net Consumer Surplus of New DER installation (including Energy Savings and DER Excess Credits associated with New DERs)

    + Annualized Net Consumer Surplus of Existing PV installation (including Energy Savings and DER Excess Credits associated with Existing DERs)
    (Note: this term is assumed to be a constant carried over from previous years and not quantified)

    + Gross Surplus and energy consumption
    (Note: this term is assumed to be a constant and not quantified (demand is inelastic))

    - Cost of Energy Purchase
    (Note: this is the out-of-pocket payment for purchasing energy from the utility company, therefore, this term double-counted the energy savings already accounted for in the Annualized Net Consumer Surplus)

    - Energy Savings associaed with both new and existing DERs
    (Note: this term is to remove the double-counted energy savings from the terms above)

    Also note that DER Excess Credits are not listed here because they're implictly accounted for in the Annualized Net Consumer Surplus.

    """

    # The NetProfit represents the energy saving/credit per representative agent per DER technology, assuming the optimal DER technology size
    NetProfit =
        make_axis_array(model_data.index_y_fix, model_data.index_h, customers.index_m)
    for y in model_data.index_y_fix, h in model_data.index_h, m in customers.index_m
        NetProfit[y, h, m] =
        # value of distributed generation (offset load)
            sum(
                model_data.omega[t] *
                regulator.p_my[y, h, t] *
                min(
                    customers.d_my[y, h, t] * (1 - utility.loss_dist),
                    customers.rho_DG[h, m, t] * customers.Opti_DG_my[y, h, m],
                ) for t in model_data.index_t
            ) +
            # value of distributed generation (excess generation)
            sum(
                model_data.omega[t] *
                regulator.p_ex_my[y, h, t] *
                max(
                    0,
                    customers.rho_DG[h, m, t] * customers.Opti_DG_my[y, h, m] -
                    customers.d_my[y, h, t] * (1 - utility.loss_dist),
                ) for t in model_data.index_t
            ) -
            # cost of distributed generation 
            customers.FOM_DG_my[y, h, m] * customers.Opti_DG_my[y, h, m]
    end

    ######## Note that this Consumer PV Net Surplus (ConPVNetSurplus_my) only calculates the surplus for year y's new PV installer (annualized)
    for y in model_data.index_y_fix, h in model_data.index_h, m in customers.index_m
        if NetProfit[y, h, m] >= 0.0
            # Calculate total Net Consumer Surplus of PV installation
            Integral = Dict(
                (y, h, m) => QuadGK.quadgk(
                    x ->
                        customers.gamma[h] *
                        customers.Opti_DG_my[y, h, m] *
                        (
                            1 - Distributions.cdf(
                                Distributions.Gamma(
                                    adopt_model.Shape[h, m],
                                    1 / adopt_model.Rate[h, m] * NetProfit[y, h, m] /
                                    customers.Opti_DG_my[y, h, m],
                                ),
                                x,
                            )
                        ),
                    customers.CapEx_DG_my[y, h, m],
                    100 * customers.CapEx_DG_my[y, h, m],
                    rtol = 1e-8,
                ),
            )
            # Calculate annualized Net Consumer Surplus of PV installation
            if customers.MaxDG_my[y, h, m] == 0.0
                customers.ConPVNetSurplus_my[y, h, m] = 0.0
            else
                customers.ConPVNetSurplus_my[y, h, m] =
                    customers.delta * customers.x_DG_new_my[y, h, m] /
                    customers.MaxDG_my[y, h, m] * Integral[y, h, m][1]
            end
        else
            customers.ConPVNetSurplus_my[y, h, m] = 0.0
        end
    end

    # Calculate energy savings associated with new DER (including previously installed new DER) for a certain year
    #=
    EnergySaving = Dict((y,h,m) =>
        sum(model_data.omega[t] * regulator.p_my[y,h,t] * min(customers.d_my[y,h,t], customers.rho_DG[h,m,t] * customers.Opti_DG_my[y,h,m]) 
            for t in model_data.index_t) * sum(customers.x_DG_new_my[Symbol(y_star),h,m] for y_star = model_data.year[first(model_data.index_y_fix)]:model_data.year[y]) / customers.Opti_DG_my[y,h,m]
            for y in model_data.index_y_fix, h in model_data.index_h, m in customers.index_m
    )
    =#
    EnergySaving =
        make_axis_array(model_data.index_y_fix, model_data.index_h, customers.index_m)
    for y in model_data.index_y_fix, h in model_data.index_h, m in customers.index_m
        EnergySaving[y, h, m] = sum(
            sum(
                model_data.omega[t] *
                regulator.p_my[y, h, t] *
                min(
                    customers.d_my[y, h, t] * (1 - utility.loss_dist),
                    customers.rho_DG[h, m, t] *
                    customers.Opti_DG_my[Symbol(Int(y_star)), h, m],
                ) for t in model_data.index_t
            ) * customers.x_DG_new_my[Symbol(Int(y_star)), h, m] /
            customers.Opti_DG_my[Symbol(Int(y_star)), h, m] for
            y_star in model_data.year[first(model_data.index_y_fix)]:model_data.year[y]
        )
    end
    # Calculate out-of-pocket energy costs associated with new and existing DER (including previously installed new DER) for a certain year (assume Opti_DG_my is the same across years)         
    EnergyCost =
        make_axis_array(model_data.index_y_fix, model_data.index_h, customers.index_m)
    for y in model_data.index_y_fix, h in model_data.index_h, m in customers.index_m
        EnergyCost[y, h, m] = sum(
            model_data.omega[t] *
            regulator.p_my[y, h, t] *
            (
                customers.gamma[h] * customers.d_my[y, h, t] * (1 - utility.loss_dist) -
                # savings from new DERs
                sum(
                    min(
                        customers.d_my[y, h, t] * (1 - utility.loss_dist),
                        customers.rho_DG[h, m, t] *
                        customers.Opti_DG_my[Symbol(Int(y_star)), h, m],
                    ) * customers.x_DG_new_my[Symbol(Int(y_star)), h, m] /
                    customers.Opti_DG_my[Symbol(Int(y_star)), h, m] for y_star in
                    model_data.year[first(model_data.index_y_fix)]:model_data.year[y]
                )
                #= We may not need this part for existing units, because we did not remove double-counting in "EnergySaving" calculation.
                # also minus savings from existing DER here (note that surplus associated with existing DER is not available)
                - min(customers.rho_DG[h,m,t]*customers.Opti_DG_E[h,m], customers.d_my[y,h,t]) * 
                customers.x_DG_E_my[y,h,m] / customers.Opti_DG_E[h,m]
                =#
            ) for t in model_data.index_t
        )
    end

    # Calculate energy costs related to export
    EnergyCost_eximport = AxisArray(
        [
            sum(
                model_data.omega[t] *
                regulator.p_eximport_my[y, t] *
                utility.eximport_my[y, t] for t in model_data.index_t
            ) for y in model_data.index_y_fix
        ],
        model_data.index_y_fix.elements,
    )

    # Finally, calculate Net Consumer Surplus
    ConNetSurplus =
        make_axis_array(model_data.index_y_fix, model_data.index_h, customers.index_m)
    for y in model_data.index_y_fix, h in model_data.index_h, m in customers.index_m
        ConNetSurplus[y, h, m] =
            sum(
                customers.ConPVNetSurplus_my[Symbol(Int(y_star)), h, m] for
                y_star in model_data.year[first(model_data.index_y_fix)]:model_data.year[y]
            ) - EnergySaving[y, h, m] - EnergyCost[y, h, m]
    end
    # Sum of Net Consumer Surplus across customer tpye and DER technology
    TotalConNetSurplus = AxisArray(
        [
            sum(
                ConNetSurplus[y, h, m] for h in model_data.index_h, m in customers.index_m
            ) - EnergyCost_eximport[y] for y in model_data.index_y_fix
        ],
        model_data.index_y_fix.elements,
    )

    ConPVNetSurplus_PerCustomer_my =
        make_axis_array(model_data.index_y_fix, model_data.index_h, customers.index_m)
    for y in model_data.index_y_fix, h in model_data.index_h, m in customers.index_m
        ConPVNetSurplus_PerCustomer_my[y, h, m] =
            customers.ConPVNetSurplus_my[y, h, m] /
            (customers.x_DG_new_my[y, h, m] / customers.Opti_DG_my[y, h, m])
    end
    AnnualBill_PerCustomer_my = make_axis_array(model_data.index_y_fix, model_data.index_h)
    for y in model_data.index_y_fix, h in model_data.index_h
        AnnualBill_PerCustomer_my[y, h] = sum(
            model_data.omega[t] *
            regulator.p_my[y, h, t] *
            customers.d_my[y, h, t] *
            (1 - utility.loss_dist) for t in model_data.index_t
        )
    end
    AverageBill_PerCustomer_my = make_axis_array(model_data.index_y_fix, model_data.index_h)
    for y in model_data.index_y_fix, h in model_data.index_h
        AverageBill_PerCustomer_my[y, h] =
            AnnualBill_PerCustomer_my[y, h] / sum(
                model_data.omega[t] * customers.d_my[y, h, t] * (1 - utility.loss_dist) for
                t in model_data.index_t
            )
    end

    return customers.ConPVNetSurplus_my,
    customers.ConGreenPowerNetSurplus_cumu_my,
    # EnergySaving,
    # EnergyCost,
    # ConNetSurplus,
    TotalConNetSurplus
    # ConPVNetSurplus_PerCustomer_my,
    # AnnualBill_PerCustomer_my,
    # AverageBill_PerCustomer_my
end





# TODO: welfare for consumer's green tech subscription
function welfare_calculation!(
    customers::CustomerGroup,
    customers_opts::AgentOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{<:MarketStructure, SupplyChoiceUseCase},
    agent_store::AgentStore,
)
    adopt_model = customers.pv_adoption_model
    green_sub_model = customers.green_sub_model

    regulator = get_agent(Regulator, agent_store)
    utility = get_agent(Utility, agent_store)
    green_developer = get_agent(GreenDeveloper, agent_store)

    """
    Net Consumer Surplus Calculation:

    + Annualized Net Consumer Surplus of Green Power Subscription

    - Cost of Energy Purchase of all other customers

    - T&D cost shared by Green Power Subscriber

    """

    # note: may need to consider distribution loss and DPV installation?
    max_sub = make_axis_array(model_data.index_y_fix, model_data.index_h)
    for y in model_data.index_y_fix, h in model_data.index_h
        max_sub[y, h] = 
        sum(customers.d_my[y, h, t] * model_data.omega[t] * customers.gamma[h] for t in model_data.index_t)
    end

    price_at_max_sub = make_axis_array(model_data.index_y_fix, model_data.index_h)
    for y in model_data.index_y_fix, h in model_data.index_h
        price_at_max_sub[y, h] = 0.0
    end

    if typeof(hem_opts) == HEMOptions{VerticallyIntegratedUtility, SupplyChoiceUseCase}
        WholesaleMarketPerc = 0.01
    else
        WholesaleMarketPerc = 1.0
    end

    GreenSubConstant = make_axis_array(model_data.index_y_fix, model_data.index_h)
    for y in model_data.index_y_fix, h in model_data.index_h
        GreenSubConstant[y, h] = 
            green_sub_model.Constant[h] + 
            green_sub_model.EnergyRate_coefficient[h] * log(regulator.p_my_regression[y, h]) + 
            green_sub_model.WholesaleMarket_coefficient[h] * log(WholesaleMarketPerc) + 
            green_sub_model.RetailCompetition_coefficient[h] * log(customers.RetailCompetition[y]) + 
            green_sub_model.RPS_coefficient[h] * log(utility.RPS[y]) + 
            green_sub_model.WTP_coefficient[h] * log(customers.WTP_green_power[y])
    end
    
    gross_surplus_integral = make_axis_array(model_data.index_y_fix, model_data.index_h)
    for y in model_data.index_y_fix, h in model_data.index_h
        gross_surplus_integral[y, h] = 0.0
    end
    
    net_surplus = make_axis_array(model_data.index_y_fix, model_data.index_h)
    for y in model_data.index_y_fix, h in model_data.index_h
        net_surplus[y, h] = 0.0
    end

    ######## Calculate annual incremental consumer surplus associated with green power subscription. This value needs to be accumulated for all previous years.
    for y in model_data.index_y_fix, h in model_data.index_h
        if customers.x_green_sub_incremental_my[y, h] > 0.0
            # price_at_max_sub[y, h] = 
            #     exp(-GreenSubConstant[y, h]/green_sub_model.GreenPowerPrice_coefficient[h])
            # gross_surplus_rectangle[y, h] = price_at_max_sub[y, h] * max_sub[y, h]
            # check to see if price_at_max_sub[y, h] is less than green_developer.ppa_my[y, h]
            gross_surplus_integral[y, h] = QuadGK.quadgk(
                x ->
                    exp(GreenSubConstant[y, h] + log(max_sub[y, h]) + 
                    green_sub_model.GreenPowerPrice_coefficient[h] * log(x)),
                green_developer.ppa_my[y, h],
                100 * green_developer.ppa_my[y, h],
                rtol = 1e-8,
            )[1]
            net_surplus[y, h] = gross_surplus_integral[y, h]
            customers.ConGreenPowerNetSurplus_pre_proportion_my[y, h] = net_surplus[y, h]
            customers.ConGreenPowerNetSurplus_post_proportion_my[y, h] = 
                net_surplus[y, h] * customers.x_green_sub_incremental_my[y, h] / customers.x_green_sub_my[y, h]
        else
            customers.ConGreenPowerNetSurplus_pre_proportion_my[y, h] = 0.0
            customers.ConGreenPowerNetSurplus_post_proportion_my[y, h] = 0.0
        end
    end

    # calculate actual annual consumer surplus associated with green power subscription by accumulating the net CS from previous years
    for y in model_data.index_y_fix, h in model_data.index_h
        customers.ConGreenPowerNetSurplus_cumu_my[y, h] = 
            sum(customers.ConGreenPowerNetSurplus_post_proportion_my[Symbol(Int(y_star)), h] for 
            y_star in model_data.year[first(model_data.index_y_fix)]:model_data.year[y])
    end

    # Calculate energy costs of all other customers, as well as green subscribers' share of T&D cost
    # here, we do not reduce the load by the DPV generation to avoid double-counting of DPV's saving.       
    EnergyCost = make_axis_array(model_data.index_y_fix, model_data.index_h)
    for y in model_data.index_y_fix, h in model_data.index_h
        EnergyCost[y, h] = sum(
            model_data.omega[t] *
            regulator.p_my[y, h, t] *
            (
                customers.gamma[h] * customers.d_my[y, h, t] * (1 - utility.loss_dist) -
                # green power subscribers are not paying the retail rates
                sum(
                    utility.rho_C_my[j, t] * sum(green_developer.green_tech_buildout_my[Symbol(Int(y_symbol)), j, h] for y_symbol in
                    model_data.year[first(model_data.index_y_fix)]:model_data.year[y])
                    for j in model_data.index_j
                )
            ) for t in model_data.index_t
        )
    end

    # Calculate energy costs related to export
    EnergyCost_eximport = AxisArray(
        [
            sum(
            model_data.omega[t] * regulator.p_eximport_my[y, t] * utility.eximport_my[y, t] for t in model_data.index_t
            ) for y in model_data.index_y_fix
        ],
        model_data.index_y_fix.elements,
    )

    # Calculate green power subscribers' T&D cost
    Green_sub_TD_charge = make_axis_array(model_data.index_y_fix, model_data.index_h)
    for y in model_data.index_y_fix, h in model_data.index_h
        Green_sub_TD_charge[y, h] = 
            regulator.p_my_td[y, h] * 
            sum(
                model_data.omega[t] *utility.rho_C_my[j, t] * sum(green_developer.green_tech_buildout_my[Symbol(Int(y_symbol)), j, h] for y_symbol in
                model_data.year[first(model_data.index_y_fix)]:model_data.year[y])
                for j in model_data.index_j, t in model_data.index_t
            )
    end

    # Finally, calculate Net Consumer Surplus
    ConNetSurplus = make_axis_array(model_data.index_y_fix, model_data.index_h)
    for y in model_data.index_y_fix, h in model_data.index_h
        ConNetSurplus[y, h] =
            customers.ConGreenPowerNetSurplus_cumu_my[y, h] - EnergyCost[y, h] - Green_sub_TD_charge[y, h]
    end
    # Sum of Net Consumer Surplus across customer tpye and DER technology
    TotalConNetSurplus = AxisArray(
        [
            sum(
                ConNetSurplus[y, h] for h in model_data.index_h
            ) - EnergyCost_eximport[y] for y in model_data.index_y_fix
        ],
        model_data.index_y_fix.elements,
    )

    # ConPVNetSurplus_PerCustomer_my = Dict(
    #     (y, h, m) =>
    #         customers.ConPVNetSurplus_my[y, h, m] /
    #         (customers.x_DG_new_my[y, h, m] / customers.Opti_DG_my[y, h, m]) for
    #     y in model_data.index_y_fix, h in model_data.index_h, m in customers.index_m
    # )
    # AnnualBill_PerCustomer_my = Dict(
    #     (y, h) => sum(
    #         model_data.omega[t] *
    #         regulator.p_my[y, h, t] *
    #         customers.d_my[y, h, t] *
    #         (1 - utility.loss_dist) for t in model_data.index_t
    #     ) for y in model_data.index_y_fix, h in model_data.index_h
    # )
    # AverageBill_PerCustomer_my = Dict(
    #     (y, h) =>
    #         AnnualBill_PerCustomer_my[y, h] / sum(
    #             model_data.omega[t] * customers.d_my[y, h, t] * (1 - utility.loss_dist)
    #             for t in model_data.index_t
    #         ) for y in model_data.index_y_fix, h in model_data.index_h
    # )

    return customers.ConPVNetSurplus_my
    customers.ConGreenPowerNetSurplus_cumu_my,
    # EnergyCost,
    # ConNetSurplus,
    TotalConNetSurplus
    # ConPVNetSurplus_PerCustomer_my,
    # AnnualBill_PerCustomer_my,
    # AverageBill_PerCustomer_my
end


function welfare_calculation!(
    customers::CustomerGroup,
    customers_opts::AgentOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{<:MarketStructure, DERSupplyChoiceUseCase},
    agent_store::AgentStore,
)
    adopt_model = customers.pv_adoption_model
    green_sub_model = customers.green_sub_model

    regulator = get_agent(Regulator, agent_store)
    utility = get_agent(Utility, agent_store)
    green_developer = get_agent(GreenDeveloper, agent_store)

    """
    Net Consumer Surplus Calculation:

    + Annualized Net Consumer Surplus of Green Power Subscription

    - Cost of Energy Purchase of all other customers

    - T&D cost shared by Green Power Subscriber

    + Annualized Net Consumer Surplus of New DER installation (including Energy Savings and DER Excess Credits associated with New DERs)

    Also note that DER Excess Credits are not listed here because they're implictly accounted for in the Annualized Net Consumer Surplus.

    """

    # The NetProfit represents the energy saving/credit per representative agent per DER technology, assuming the optimal DER technology size
    NetProfit = make_axis_array(model_data.index_y_fix, model_data.index_h, customers.index_m)
    for y in model_data.index_y_fix, h in model_data.index_h, m in customers.index_m
        NetProfit[y, h, m] =
        # value of distributed generation (offset load)
            sum(
                model_data.omega[t] *
                regulator.p_my[y, h, t] *
                min(
                    customers.d_my[y, h, t] * (1 - utility.loss_dist),
                    customers.rho_DG[h, m, t] * customers.Opti_DG_my[y, h, m],
                ) for t in model_data.index_t
            ) +
            # value of distributed generation (excess generation)
            sum(
                model_data.omega[t] *
                regulator.p_ex_my[y, h, t] *
                max(
                    0,
                    customers.rho_DG[h, m, t] * customers.Opti_DG_my[y, h, m] -
                    customers.d_my[y, h, t] * (1 - utility.loss_dist),
                ) for t in model_data.index_t
            ) -
            # cost of distributed generation 
            customers.FOM_DG_my[y, h, m] * customers.Opti_DG_my[y, h, m]
    end

    ######## Note that this Consumer PV Net Surplus (ConPVNetSurplus_my) only calculates the surplus for year y's new PV installer (annualized)
    for y in model_data.index_y_fix, h in model_data.index_h, m in customers.index_m
        if NetProfit[y, h, m] >= 0.0
            # Calculate total Net Consumer Surplus of PV installation
            Integral = Dict(
                (y, h, m) => QuadGK.quadgk(
                    x ->
                        customers.gamma[h] *
                        customers.Opti_DG_my[y, h, m] *
                        (
                            1 - Distributions.cdf(
                                Distributions.Gamma(
                                    adopt_model.Shape[h, m],
                                    1 / adopt_model.Rate[h, m] * NetProfit[y, h, m] /
                                    customers.Opti_DG_my[y, h, m],
                                ),
                                x,
                            )
                        ),
                    customers.CapEx_DG_my[y, h, m],
                    100 * customers.CapEx_DG_my[y, h, m],
                    rtol = 1e-8,
                ),
            )
            # Calculate annualized Net Consumer Surplus of PV installation
            if customers.MaxDG_my[y, h, m] == 0.0
                customers.ConPVNetSurplus_my[y, h, m] = 0.0
            else
                customers.ConPVNetSurplus_my[y, h, m] =
                    customers.delta * customers.x_DG_new_my[y, h, m] /
                    customers.MaxDG_my[y, h, m] * Integral[y, h, m][1]
            end
        else
            customers.ConPVNetSurplus_my[y, h, m] = 0.0
        end
    end

    # note: may need to consider distribution loss and DPV installation?
    max_sub = make_axis_array(model_data.index_y_fix, model_data.index_h)
    for y in model_data.index_y_fix, h in model_data.index_h
        max_sub[y, h] = 
        sum(customers.d_my[y, h, t] * model_data.omega[t] * customers.gamma[h] for t in model_data.index_t)
    end

    price_at_max_sub = make_axis_array(model_data.index_y_fix, model_data.index_h)
    for y in model_data.index_y_fix, h in model_data.index_h
        price_at_max_sub[y, h] = 0.0
    end

    if typeof(hem_opts) == HEMOptions{VerticallyIntegratedUtility, SupplyChoiceUseCase}
        WholesaleMarketPerc = 0.01
    else
        WholesaleMarketPerc = 1.0
    end

    GreenSubConstant = make_axis_array(model_data.index_y_fix, model_data.index_h)
    for y in model_data.index_y_fix, h in model_data.index_h
        GreenSubConstant[y, h] = 
            green_sub_model.Constant[h] + 
            green_sub_model.EnergyRate_coefficient[h] * log(regulator.p_my_regression[y, h]) + 
            green_sub_model.WholesaleMarket_coefficient[h] * log(WholesaleMarketPerc) + 
            green_sub_model.RetailCompetition_coefficient[h] * log(customers.RetailCompetition[y]) + 
            green_sub_model.RPS_coefficient[h] * log(utility.RPS[y]) + 
            green_sub_model.WTP_coefficient[h] * log(customers.WTP_green_power[y])
    end
    
    gross_surplus_integral = make_axis_array(model_data.index_y_fix, model_data.index_h)
    for y in model_data.index_y_fix, h in model_data.index_h
        gross_surplus_integral[y, h] = 0.0
    end
    
    net_surplus = make_axis_array(model_data.index_y_fix, model_data.index_h)
    for y in model_data.index_y_fix, h in model_data.index_h
        net_surplus[y, h] = 0.0
    end

    ######## Calculate annual incremental consumer surplus associated with green power subscription. This value needs to be accumulated for all previous years.
    for y in model_data.index_y_fix, h in model_data.index_h
        if customers.x_green_sub_incremental_my[y, h] > 0.0
            # price_at_max_sub[y, h] = 
            #     exp(-GreenSubConstant[y, h]/green_sub_model.GreenPowerPrice_coefficient[h])
            # gross_surplus_rectangle[y, h] = price_at_max_sub[y, h] * max_sub[y, h]
            # check to see if price_at_max_sub[y, h] is less than green_developer.ppa_my[y, h]
            gross_surplus_integral[y, h] = QuadGK.quadgk(
                x ->
                    exp(GreenSubConstant[y, h] + log(max_sub[y, h]) + 
                    green_sub_model.GreenPowerPrice_coefficient[h] * log(x)),
                green_developer.ppa_my[y, h],
                100 * green_developer.ppa_my[y, h],
                rtol = 1e-8,
            )[1]
            net_surplus[y, h] = gross_surplus_integral[y, h]
            customers.ConGreenPowerNetSurplus_pre_proportion_my[y, h] = net_surplus[y, h]
            customers.ConGreenPowerNetSurplus_post_proportion_my[y, h] = 
                net_surplus[y, h] * customers.x_green_sub_incremental_my[y, h] / customers.x_green_sub_my[y, h]
        else
            customers.ConGreenPowerNetSurplus_pre_proportion_my[y, h] = 0.0
            customers.ConGreenPowerNetSurplus_post_proportion_my[y, h] = 0.0
        end
    end

    # calculate actual annual consumer surplus associated with green power subscription by accumulating the net CS from previous years
    for y in model_data.index_y_fix, h in model_data.index_h
        customers.ConGreenPowerNetSurplus_cumu_my[y, h] = 
            sum(customers.ConGreenPowerNetSurplus_post_proportion_my[Symbol(Int(y_star)), h] for 
            y_star in model_data.year[first(model_data.index_y_fix)]:model_data.year[y])
    end

    # Calculate energy costs of all other customers, as well as green subscribers' share of T&D cost
    # here, we do not reduce the load by the DPV generation to avoid double-counting of DPV's saving.
    EnergyCost = make_axis_array(model_data.index_y_fix, model_data.index_h)
    for y in model_data.index_y_fix, h in model_data.index_h
        EnergyCost[y, h] = sum(
            model_data.omega[t] *
            regulator.p_my[y, h, t] *
            (
                customers.gamma[h] * customers.d_my[y, h, t] * (1 - utility.loss_dist) -
                # green power subscribers are not paying the retail rates
                sum(
                    utility.rho_C_my[j, t] * sum(green_developer.green_tech_buildout_my[Symbol(Int(y_symbol)), j, h] for y_symbol in
                    model_data.year[first(model_data.index_y_fix)]:model_data.year[y])
                    for j in model_data.index_j
                )
            ) for t in model_data.index_t
        )
    end

    # Calculate energy costs related to export
    EnergyCost_eximport = AxisArray(
        [
            sum(
            model_data.omega[t] * regulator.p_eximport_my[y, t] * utility.eximport_my[y, t] for t in model_data.index_t
            ) for y in model_data.index_y_fix
        ],
        model_data.index_y_fix.elements,
    )

    # Calculate green power subscribers' T&D cost
    Green_sub_TD_charge = make_axis_array(model_data.index_y_fix, model_data.index_h)
    for y in model_data.index_y_fix, h in model_data.index_h
        Green_sub_TD_charge[y, h] = 
            regulator.p_my_td[y, h] * 
            sum(
                model_data.omega[t] *utility.rho_C_my[j, t] * sum(green_developer.green_tech_buildout_my[Symbol(Int(y_symbol)), j, h] for y_symbol in
                model_data.year[first(model_data.index_y_fix)]:model_data.year[y])
                for j in model_data.index_j, t in model_data.index_t
            )
    end

    # Finally, calculate Net Consumer Surplus
    ConNetSurplus = make_axis_array(model_data.index_y_fix, model_data.index_h)
    for y in model_data.index_y_fix, h in model_data.index_h
        ConNetSurplus[y, h] =
            customers.ConGreenPowerNetSurplus_cumu_my[y, h] - EnergyCost[y, h] - Green_sub_TD_charge[y, h] +
            sum(customers.ConPVNetSurplus_my[y, h, m] for m in customers.index_m)
    end
    # Sum of Net Consumer Surplus across customer tpye and DER technology
    TotalConNetSurplus = AxisArray(
        [
            sum(
                ConNetSurplus[y, h] for h in model_data.index_h
            ) - EnergyCost_eximport[y] for y in model_data.index_y_fix
        ],
        model_data.index_y_fix.elements,
    )

    # ConPVNetSurplus_PerCustomer_my = Dict(
    #     (y, h, m) =>
    #         customers.ConPVNetSurplus_my[y, h, m] /
    #         (customers.x_DG_new_my[y, h, m] / customers.Opti_DG_my[y, h, m]) for
    #     y in model_data.index_y_fix, h in model_data.index_h, m in customers.index_m
    # )
    # AnnualBill_PerCustomer_my = Dict(
    #     (y, h) => sum(
    #         model_data.omega[t] *
    #         regulator.p_my[y, h, t] *
    #         customers.d_my[y, h, t] *
    #         (1 - utility.loss_dist) for t in model_data.index_t
    #     ) for y in model_data.index_y_fix, h in model_data.index_h
    # )
    # AverageBill_PerCustomer_my = Dict(
    #     (y, h) =>
    #         AnnualBill_PerCustomer_my[y, h] / sum(
    #             model_data.omega[t] * customers.d_my[y, h, t] * (1 - utility.loss_dist)
    #             for t in model_data.index_t
    #         ) for y in model_data.index_y_fix, h in model_data.index_h
    # )

    return customers.ConPVNetSurplus_my,
    customers.ConGreenPowerNetSurplus_cumu_my,
    # EnergyCost,
    # ConNetSurplus,
    TotalConNetSurplus
    # ConPVNetSurplus_PerCustomer_my,
    # AnnualBill_PerCustomer_my,
    # AverageBill_PerCustomer_my
end