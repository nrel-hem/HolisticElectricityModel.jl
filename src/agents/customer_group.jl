# This file defines data and functions associated with the customer.

DER_factor = 1.0    # A scaling factor applied to existing DER penetration level

mutable struct PVAdoptionModel
    Shape::ParamArray
    MeanPayback::ParamArray
    Bass_p::ParamVector
    Bass_q::ParamVector

    Rate::ParamArray
end

"""
Constructs a PVAdoptionModel by computing Rate from the other parameters.
"""
function PVAdoptionModel(Shape, MeanPayback, Bass_p, Bass_q)
    return PVAdoptionModel(
        Shape,
        MeanPayback,
        Bass_p,
        Bass_q,
        ParamArray(
            "Rate",
            Shape.dims,
            Dict(key => Shape[key] / MeanPayback[key] for key in keys(Shape)),
        ), # Rate
    )
end

abstract type AbstractCustomerGroup <: AgentGroup end

mutable struct CustomerGroup <: AbstractCustomerGroup
    id::String
    # Sets
    index_m::Dimension # behind-the-meter technologies

    # Parameters
    "number of customers of type h"
    gamma::ParamVector
    "demand (MWh per representative agent per hour)"
    d::ParamArray
    "multi-year demand (MWh per representative agent per hour)"
    d_my::ParamArray
    x_DG_E::ParamArray
    "Existing DER at year y. This is a cumulative number but without x_DG_new_my built by this module"
    x_DG_E_my::ParamArray
    Opti_DG::ParamArray
    Opti_DG_E::ParamArray
    Opti_DG_my::ParamArray
    # "DER generation by a representative customer h and DER technology m"
    # DERGen::ParamArray
    CapEx_DG::ParamArray
    CapEx_DG_my::ParamArray
    FOM_DG::ParamArray
    FOM_DG_my::ParamArray
    rho_DG::ParamArray
    "Annualization factor for net consumer surplus of PV installation"
    delta::ParamScalar
    PeakLoad::ParamVector
    PeakLoad_my::ParamArray

    # Primal Variables
    x_DG_new::ParamArray
    x_DG_new_my::ParamArray    # Annual new DER build (not cumulative)
    x_green_sub::ParamArray

    # Auxiliary Variables
    Payback::ParamArray
    MarketShare::ParamArray
    MaxDG::ParamArray
    F::ParamArray
    year::ParamArray
    A::ParamArray
    ConPVNetSurplus::ParamArray
    ConPVNetSurplus_my::ParamArray

    GreenTechIntercept::ParamArray
    GreenTechSlope::ParamArray
    pv_adoption_model::PVAdoptionModel

    pvf::Any
    rooftop::Any
    MaxDG_my::ParamArray
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
    MaxLoad = Dict(
        h => gamma[h] * findmax(Dict(t => demand[h, t] for t in model_data.index_t))[1]
        for h in model_data.index_h
    )
    MaxLoad_my = Dict(
        (y, h) =>
            gamma[h] *
            findmax(Dict(t => demand_my[y, h, t] for t in model_data.index_t))[1] for
        y in model_data.index_y, h in model_data.index_h
    )

    pv_adoption_model = PVAdoptionModel(
        initialize_param("Shape", model_data.index_h, index_m, value = 1.7), # Shape
        initialize_param("MeanPayback", model_data.index_h, index_m, value = 8.8), # MeanPayback
        ParamVector(
            "Bass_p",
            model_data.index_h,
            Dict(:Residential => 7.7E-07, :Commercial => 6.0E-04, :Industrial => 6.0E-04),
        ), # Bass_p
        ParamVector(
            "Bass_q",
            model_data.index_h,
            Dict(:Residential => 0.663, :Commercial => 0.133, :Industrial => 0.133),
        ), # Bass_q
    )

    # Customer financing
    debt_ratio =
        read_param("debt_ratio", input_filename, "CustomerDebtRatio", model_data.index_h)
    cost_of_debt =
        read_param("cost_of_debt", input_filename, "CustomerCOD", model_data.index_h)
    cost_of_equity =
        read_param("cost_of_equity", input_filename, "CustomerCOE", model_data.index_h)
    tax_rate = read_param("tax_rate", input_filename, "CustomerTax", model_data.index_h)
    atwacc = Dict(
        h =>
            debt_ratio[h] * cost_of_debt[h] * (1 - tax_rate[h]) +
            (1 - debt_ratio[h]) * cost_of_equity[h] for h in model_data.index_h
    )
    CRF = Dict(
        h => atwacc[h] * (1 + atwacc[h])^20 / ((1 + atwacc[h])^20 - 1) for
        h in model_data.index_h
    )
    pvf = Dict(h => 1 / CRF[h] for h in model_data.index_h)

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
        ParamVector("PeakLoad", model_data.index_h, MaxLoad),
        ParamArray(
            "PeakLoad_my",
            Tuple(push!(copy([model_data.index_y]), model_data.index_h)),
            MaxLoad_my,
        ),
        initialize_param("x_DG_new", model_data.index_h, index_m),
        initialize_param("x_DG_new_my", model_data.index_y, model_data.index_h, index_m),
        initialize_param("x_green_sub", model_data.index_h, model_data.index_j),
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
        pvf,
        read_param("rooftop", input_filename, "RooftopDER", index_m, [model_data.index_h]),
        initialize_param("MaxDG_my", model_data.index_y, model_data.index_h, index_m),
    )
end

get_id(x::CustomerGroup) = x.id

function solve_agent_problem!(
    customers::CustomerGroup,
    customers_opts::AgentOptions,
    model_data::HEMData,
    hem_opts::HEMOptions,
    agent_store::AgentStore,
    w_iter,
)
    regulator = get_agent(Regulator, agent_store)
    utility = get_agent(Utility, agent_store)

    # the year consumer is making DER investment decision
    reg_year = model_data.year[first(model_data.index_y)]
    reg_year_index = Symbol(Int(reg_year))

    x_DG_before = ParamArray(customers.x_DG_new, "x_DG_before")
    empty!(x_DG_before)
    merge!(
        x_DG_before,
        Dict(
            (h, m) => customers.x_DG_new_my[reg_year_index, h, m] for
            h in model_data.index_h, m in customers.index_m
        ),
    )

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
    NetProfit = Dict(
        (h, m) =>
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
            customers.FOM_DG[h, m] * customers.Opti_DG[h, m] for
        h in model_data.index_h, m in customers.index_m
    )

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

function save_results(
    customers::CustomerGroup,
    customers_opts::AgentOptions,
    hem_opts::HEMOptions,
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

function welfare_calculation(
    customers::CustomerGroup,
    model_data::HEMData,
    regulator::Agent,
)
    adopt_model = customers.pv_adoption_model

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
    NetProfit = Dict(
        (h, m) =>
        # value of distributed generation (offset load)
            sum(
                model_data.omega[t] *
                regulator.p[h, t] *
                min(
                    customers.d[h, t],
                    sum(
                        customers.rho_DG[h, m, t] * customers.Opti_DG[h, m] for
                        m in customers.index_m
                    ),
                ) *
                (customers.rho_DG[h, m, t] * customers.Opti_DG[h, m]) /
                customers.DERGen[h, t] for t in model_data.index_t
            ) +
            # value of distributed generation (excess generation)
            sum(
                model_data.omega[t] *
                regulator.p_ex[h, t] *
                max(
                    0,
                    sum(
                        customers.rho_DG[h, m, t] * customers.Opti_DG[h, m] for
                        m in customers.index_m
                    ) - customers.d[h, t],
                ) *
                (customers.rho_DG[h, m, t] * customers.Opti_DG[h, m]) /
                customers.DERGen[h, t] for t in model_data.index_t
            ) -
            # cost of distributed generation 
            customers.FOM_DG[h, m] * customers.Opti_DG[h, m] for
        h in model_data.index_h, m in customers.index_m
    )

    for h in model_data.index_h, m in customers.index_m
        if NetProfit[h, m] >= 0.0
            # Calculate total Net Consumer Surplus of PV installation
            Integral = Dict(
                (h, m) => quadgk(
                    x ->
                        customers.gamma[h] *
                        customers.Opti_DG[h, m] *
                        (
                            1 - Distributions.cdf(
                                Distributions.Gamma(
                                    adopt_model.Shape[h, m],
                                    1 / adopt_model.Rate[h, m] * NetProfit[h, m] /
                                    customers.Opti_DG[h, m],
                                ),
                                x,
                            )
                        ),
                    customers.CapEx_DG[h, m],
                    100 * customers.CapEx_DG[h, m],
                    rtol = 1e-8,
                ),
            )
            # Calculate annualized Net Consumer Surplus of PV installation
            if customers.MaxDG[h, m] == 0.0
                customers.ConPVNetSurplus[h, m] = 0.0
            else
                customers.ConPVNetSurplus[h, m] =
                    customers.delta * customers.x_DG_new[h, m] / customers.MaxDG[h, m] *
                    Integral[h, m][1]
            end
        else
            customers.ConPVNetSurplus[h, m] = 0.0
        end
    end

    # Calculate energy savings associated with both new and existing DER
    EnergySaving = Dict(
        (h, m) =>
            sum(
                model_data.omega[t] *
                regulator.p[h, t] *
                min(
                    customers.d[h, t],
                    sum(
                        customers.rho_DG[h, m, t] * customers.Opti_DG[h, m] for
                        m in customers.index_m
                    ),
                ) *
                (customers.rho_DG[h, m, t] * customers.Opti_DG[h, m]) /
                customers.DERGen[h, t] for t in model_data.index_t
            ) * (customers.x_DG_E[h, m] + customers.x_DG_new[h, m]) /
            customers.Opti_DG[h, m] for h in model_data.index_h, m in customers.index_m
    )
    # Calculate out-of-pocket energy costs           
    EnergyCost = Dict(
        (h, m) => sum(
            model_data.omega[t] *
            regulator.p[h, t] *
            (
                customers.gamma[h] * customers.d[h, t] -
                min(
                    sum(
                        customers.rho_DG[h, m, t] * customers.Opti_DG[h, m] for
                        m in customers.index_m
                    ),
                    customers.d[h, t],
                ) * (customers.rho_DG[h, m, t] * customers.Opti_DG[h, m]) /
                customers.DERGen[h, t] *
                (customers.x_DG_E[h, m] + customers.x_DG_new[h, m]) /
                customers.Opti_DG[h, m]
            ) for t in model_data.index_t
        ) for h in model_data.index_h, m in customers.index_m
    )
    # Finally, calculate Net Consumer Surplus
    ConNetSurplus = Dict(
        (h, m) =>
            customers.ConPVNetSurplus[h, m] - EnergySaving[h, m] - EnergyCost[h, m] for
        h in model_data.index_h, m in customers.index_m
    )
    # Sum of Net Consumer Surplus across customer tpye and DER technology
    TotalConNetSurplus =
        sum(ConNetSurplus[h, m] for h in model_data.index_h, m in customers.index_m)

    return customers.ConPVNetSurplus,
    EnergySaving,
    EnergyCost,
    ConNetSurplus,
    TotalConNetSurplus
end

function welfare_calculation!(
    customers::CustomerGroup,
    customers_opts::AgentOptions,
    model_data::HEMData,
    hem_opts::HEMOptions,
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
    NetProfit = Dict(
        (y, h, m) =>
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
            customers.FOM_DG_my[y, h, m] * customers.Opti_DG_my[y, h, m] for
        y in model_data.index_y_fix, h in model_data.index_h, m in customers.index_m
    )

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
    EnergySaving = Dict(
        (y, h, m) => sum(
            sum(
                model_data.omega[t] *
                regulator.p_my[y, h, t] *
                min(
                    customers.d_my[y, h, t] * (1 - utility.loss_dist),
                    customers.rho_DG[h, m, t] *
                    customers.Opti_DG_my[Symbol(Int(y_star)), h, m],
                ) for t in model_data.index_t
            ) * customers.x_DG_new_my[Symbol(Int(y_star)), h, m] /
            customers.Opti_DG_my[Symbol(Int(y_star)), h, m] for y_star in
                model_data.year[first(model_data.index_y_fix)]:model_data.year[y]
        ) for y in model_data.index_y_fix, h in model_data.index_h,
        m in customers.index_m
    )
    # Calculate out-of-pocket energy costs associated with new and existing DER (including previously installed new DER) for a certain year (assume Opti_DG_my is the same across years)         
    EnergyCost = Dict(
        (y, h, m) => sum(
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
        ) for y in model_data.index_y_fix, h in model_data.index_h,
        m in customers.index_m
    )

    # Calculate energy costs related to export
    EnergyCost_eximport = Dict(
        y => sum(
            model_data.omega[t] * regulator.p_eximport_my[y, t] * utility.eximport_my[y, t] for t in model_data.index_t
        ) for y in model_data.index_y_fix
    )

    # Finally, calculate Net Consumer Surplus
    ConNetSurplus = Dict(
        (y, h, m) =>
            sum(
                customers.ConPVNetSurplus_my[Symbol(Int(y_star)), h, m] for y_star in
                    model_data.year[first(model_data.index_y_fix)]:model_data.year[y]
            ) - EnergySaving[y, h, m] - EnergyCost[y, h, m] for
        y in model_data.index_y_fix, h in model_data.index_h, m in customers.index_m
    )
    # Sum of Net Consumer Surplus across customer tpye and DER technology
    TotalConNetSurplus = Dict(
        y =>
            sum(
                ConNetSurplus[y, h, m] for h in model_data.index_h, m in customers.index_m
            ) - EnergyCost_eximport[y] for y in model_data.index_y_fix
    )

    ConPVNetSurplus_PerCustomer_my = Dict(
        (y, h, m) =>
            customers.ConPVNetSurplus_my[y, h, m] /
            (customers.x_DG_new_my[y, h, m] / customers.Opti_DG_my[y, h, m]) for
        y in model_data.index_y_fix, h in model_data.index_h, m in customers.index_m
    )
    AnnualBill_PerCustomer_my = Dict(
        (y, h) => sum(
            model_data.omega[t] *
            regulator.p_my[y, h, t] *
            customers.d_my[y, h, t] *
            (1 - utility.loss_dist) for t in model_data.index_t
        ) for y in model_data.index_y_fix, h in model_data.index_h
    )
    AverageBill_PerCustomer_my = Dict(
        (y, h) =>
            AnnualBill_PerCustomer_my[y, h] / sum(
                model_data.omega[t] * customers.d_my[y, h, t] * (1 - utility.loss_dist)
                for t in model_data.index_t
            ) for y in model_data.index_y_fix, h in model_data.index_h
    )

    return customers.ConPVNetSurplus_my,
    EnergySaving,
    EnergyCost,
    ConNetSurplus,
    TotalConNetSurplus,
    ConPVNetSurplus_PerCustomer_my,
    AnnualBill_PerCustomer_my,
    AverageBill_PerCustomer_my
end
