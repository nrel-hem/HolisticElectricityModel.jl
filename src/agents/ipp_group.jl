# This module defines the data and functions associated with the Independent Power Producer

# declare customer decision
abstract type IPPAlgorithm end
struct LagrangeDecomposition <: IPPAlgorithm end
struct MIQP <: IPPAlgorithm end

abstract type AbstractIPPOptions <: AgentOptions end

struct IPPOptions{T <: IPPAlgorithm} <: AbstractIPPOptions
    ipp_algorithm::T
end

abstract type AbstractIPPGroup <: AgentGroup end

mutable struct IPPGroup <: AbstractIPPGroup
    id::String
    # Sets
    index_k_existing::Dimension # existing bulk generation technologies
    index_k_new::Dimension # potential bulk generation technologies
    index_p::Dimension # individual independent power producers
    index_rps::Dimension # RPS-qualified technologies

    # Parameters
    "existing capacity (MW)"
    x_E::ParamAxisArray
    "fixed cost of existing capacity (\$/MW-yr)"
    f_E::ParamAxisArray
    "fixed cost of new capacity (\$/MW-yr)"
    f_C::ParamAxisArray
    "variable cost of existing capacity (\$/MWh)"
    v_E::ParamAxisArray
    "variable cost of new capacity (\$/MWh)"
    v_C::ParamAxisArray
    "availability of existing capacity (fraction)"
    rho_E::ParamAxisArray
    "availability of new capacity (fraction)"
    rho_C::ParamAxisArray
    eximport::ParamAxisArray # net export (MWh)
    Peak_eximport::ParamScalar
    "Big M Parameter"
    B1GM::ParamScalar{<:Integer}
    zeta::ParamScalar # offer price factor cap

    # Primal Variables
    y_E::ParamAxisArray
    y_C::ParamAxisArray
    x_R::ParamAxisArray
    x_C::ParamAxisArray
    miu::ParamAxisArray
    o_E::ParamAxisArray # offer price of existing capacity ($/MWh)
    o_C::ParamAxisArray # offer price of new capacity ($/MWh)
    LMP::ParamAxisArray

    # Parameters (multi-year)
    x_E_my::ParamAxisArray # existing capacity (MW)
    fom_E_my::ParamAxisArray # fixed O&M of existing capacity ($/MW-yr)
    fom_C_my::ParamAxisArray # fixed O&M of new capacity ($/MW-yr)
    CapEx_my::ParamAxisArray # capital expense of new capacity ($/MW)
    v_E_my::ParamAxisArray # variable cost of existing capacity ($/MWh)
    v_C_my::ParamAxisArray # variable cost of new capacity ($/MWh)
    rho_E_my::ParamAxisArray # availability of existing capacity (fraction)
    rho_C_my::ParamAxisArray # availability of new capacity (fraction)
    eximport_my::ParamAxisArray # net export (MWh)
    pvf_cap::ParamAxisArray # present value factor of capital expenses
    pvf_onm::ParamAxisArray # present value factor of o&m expenses
    CRF_default::ParamAxisArray
    Tax::ParamAxisArray
    DebtRatio::ParamAxisArray
    COD::ParamAxisArray
    COE::ParamAxisArray

    # Primal Variables (multi-year)
    y_E_my::ParamAxisArray
    y_C_my::ParamAxisArray
    x_R_my::ParamAxisArray
    x_C_my::ParamAxisArray
    o_E_my::ParamAxisArray # offer price of existing capacity ($/MWh)
    o_C_my::ParamAxisArray # offer price of new capacity ($/MWh)
    # Dual Variables (multi-year)
    miu_my::ParamAxisArray
    LMP_my::ParamAxisArray
    eta_my::ParamAxisArray
    lambda_my::ParamAxisArray
    u_y_E_my::ParamAxisArray
    u_y_C_my::ParamAxisArray
    u_miu_my::ParamAxisArray
    u_eta_my::ParamAxisArray
    u_lambda_my::ParamAxisArray
    x_R_cumu::ParamAxisArray
    x_C_cumu::ParamAxisArray

    # Lagrange multiplier
    L_R::ParamAxisArray
    L_C::ParamAxisArray
    x_R_my_st1::ParamAxisArray
    x_C_my_st1::ParamAxisArray
    x_R_my_st2::ParamAxisArray
    x_C_my_st2::ParamAxisArray
    # value of objective function
    obj_st1::ParamScalar
    obj_st2::ParamScalar
    obj_lower_bound::ParamScalar
    obj_upper_bound::ParamScalar
    obj_feasible::ParamScalar

    #temporary save_results
    y_E_my_temp::ParamAxisArray
    y_C_my_temp::ParamAxisArray
    x_R_my_temp::ParamAxisArray
    x_C_my_temp::ParamAxisArray
    miu_my_temp::ParamAxisArray
    LMP_my_temp::ParamAxisArray

    # capacity market parameters
    NetCONE::ParamAxisArray
    DC_length::ParamAxisArray
    capacity_credit_E_my::ParamAxisArray # capacity credit of existing resources
    capacity_credit_C_my::ParamAxisArray # capacity credit of new resources
    capacity_price::ParamAxisArray # $/MW-yr
    capacity_price_my_temp::ParamAxisArray

    Net_Load_my::ParamAxisArray
    Max_Net_Load_my::ParamAxisArray
    Reserve_req_my::ParamAxisArray
    Capacity_slope_my::ParamAxisArray
    Capacity_intercept_my::ParamAxisArray
    ucap_temp::ParamAxisArray
    ucap::ParamAxisArray

    # RPS
    RPS::ParamAxisArray

    # emission rate
    emission_rate_E_my::ParamAxisArray
    emission_rate_C_my::ParamAxisArray
end

function IPPGroup(input_filename::String, model_data::HEMData, id = DEFAULT_ID)
    index_k_existing = read_set(input_filename, "index_k_existing", "index_k_existing")
    index_k_new = read_set(input_filename, "index_k_new", "index_k_new")
    index_p = read_set(input_filename, "index_p", "index_p")
    index_rps = read_set(input_filename, "index_rps", "index_rps")

    FOMNew = read_param("FOM_new", input_filename, "FOMNewIPP", index_k_new, [index_p])
    CapExNew =
        read_param("CapEx_new", input_filename, "CapExNewIPP", index_k_new, [index_p])
    LifetimeNew =
        read_param("Lifetime_new", input_filename, "LifetimeNewIPP", index_k_new, [index_p])
    debt_ratio = read_param("DebtRatio", input_filename, "DebtRatio", index_p)
    cost_of_debt = read_param("COD", input_filename, "COD", index_p)
    cost_of_equity = read_param("COE", input_filename, "COE", index_p)
    tax_rate = read_param("Tax", input_filename, "Tax", index_p)
    atwacc = Dict(
        p =>
            debt_ratio(p) * cost_of_debt(p) * (1 - tax_rate(p)) +
            (1 - debt_ratio(p)) * cost_of_equity(p) for p in index_p
    )
    CRF = Dict(
        (p, k) =>
            atwacc[p] * (1 + atwacc[p])^LifetimeNew(p, k) /
            ((1 + atwacc[p])^LifetimeNew(p, k) - 1) for p in index_p, k in index_k_new
    )
    FixedCostNew = make_keyed_array(index_p, index_k_new)
    for p in index_p, k in index_k_new
        FixedCostNew(p, k, :) .= FOMNew(p, k) + CapExNew(p, k) * CRF[p, k]
    end

    eximport = read_param("eximport", input_filename, "Export", model_data.index_t)
    peak_eximport =
        ParamScalar("Peak_eximport", findmax(eximport)[1], description = "peak export")

    CRF_default = KeyedArray(
        [atwacc[p] * (1 + atwacc[p])^20 / ((1 + atwacc[p])^20 - 1) for p in index_p];
        [get_pair(index_p)]...
    )
    pvf_cap = make_keyed_array(model_data.index_y, index_p)
    for y in model_data.index_y, p in index_p
        pvf_cap(y, p, :) .= 1 / (1 + atwacc[p])^(model_data.year(y) - model_data.year_start)
    end
    pvf_onm = make_keyed_array(model_data.index_y, index_p)
    for y in model_data.index_y, p in index_p
        pvf_onm(y, p, :) .= 1 / (1 + atwacc[p])^(model_data.year(y) - model_data.year_start)
    end

    NetCONE = read_param("NetCONE", input_filename, "NetCONE", model_data.index_y)     # $/MW-yr
    DC_length = read_param("DC_length", input_filename, "DC_length", model_data.index_y)

    return IPPGroup(
        id,
        index_k_existing,
        index_k_new,
        index_p,
        index_rps,
        read_param(
            "x_E",
            input_filename,
            "ExistingCapacityIPP",
            index_k_existing,
            [index_p],
        ),
        read_param("f_E", input_filename, "FixedCostOldIPP", index_k_existing, [index_p]),
        ParamAxisArray("f_C", Tuple(push!(copy([index_p]), index_k_new)), FixedCostNew),
        read_param(
            "v_E",
            input_filename,
            "VariableCostOldIPP",
            model_data.index_t,
            [index_p, index_k_existing],
        ),
        read_param(
            "v_C",
            input_filename,
            "VariableCostNewIPP",
            model_data.index_t,
            [index_p, index_k_new],
        ),
        read_param(
            "rho_E",
            input_filename,
            "AvailabilityOldIPP",
            model_data.index_t,
            [index_p, index_k_existing],
        ),
        read_param(
            "rho_C",
            input_filename,
            "AvailabilityNewIPP",
            model_data.index_t,
            [index_p, index_k_new],
        ),
        eximport,
        peak_eximport,
        ParamScalar("B1GM", 1000000),
        ParamScalar("zeta", 3.0),
        initialize_param("y_E", index_p, index_k_existing, model_data.index_t),
        initialize_param("y_C", index_p, index_k_new, model_data.index_t),
        initialize_param("x_R", index_p, index_k_existing),
        initialize_param("x_C", index_p, index_k_new),
        initialize_param("miu", model_data.index_t),
        read_param(
            "o_E",
            input_filename,
            "VariableCostOldIPP",
            model_data.index_t,
            [index_p, index_k_existing],
        ),
        read_param(
            "o_C",
            input_filename,
            "VariableCostNewIPP",
            model_data.index_t,
            [index_p, index_k_new],
        ),
        initialize_param("LMP", model_data.index_t),
        read_param(
            "x_E_my",
            input_filename,
            "ExistingCapacityIPP",
            index_k_existing,
            [index_p],
        ),
        read_param(
            "fom_E_my",
            input_filename,
            "FixedCostOldIPPmy",
            index_k_existing,
            [model_data.index_y, index_p],
        ),
        read_param(
            "fom_C_my",
            input_filename,
            "FOMNewIPPmy",
            index_k_new,
            [model_data.index_y, index_p],
        ),
        read_param(
            "CapEx_my",
            input_filename,
            "CapExNewIPPmy",
            index_k_new,
            [model_data.index_y, index_p],
        ),
        read_param(
            "v_E_my",
            input_filename,
            "VariableCostOldIPPmy",
            model_data.index_t,
            [model_data.index_y, index_p, index_k_existing],
        ),
        read_param(
            "v_C_my",
            input_filename,
            "VariableCostNewIPPmy",
            model_data.index_t,
            [model_data.index_y, index_p, index_k_new],
        ),
        read_param(
            "rho_E_my",
            input_filename,
            "AvailabilityOldIPP",
            model_data.index_t,
            [index_p, index_k_existing],
        ),
        read_param(
            "rho_C_my",
            input_filename,
            "AvailabilityNewIPP",
            model_data.index_t,
            [index_p, index_k_new],
        ),
        read_param(
            "eximport_my",
            input_filename,
            "Exportmy",
            model_data.index_t,
            [model_data.index_y],
        ),
        ParamAxisArray(
            "pvf_cap",
            Tuple(push!(copy([model_data.index_y]), index_p)),
            pvf_cap,
        ),
        ParamAxisArray(
            "pvf_onm",
            Tuple(push!(copy([model_data.index_y]), index_p)),
            pvf_onm,
        ),
        ParamAxisArray("CRF_default", (index_p,), CRF_default),
        tax_rate,
        debt_ratio,
        cost_of_debt,
        cost_of_equity,
        initialize_param(
            "y_E_my",
            model_data.index_y,
            index_p,
            index_k_existing,
            model_data.index_t,
        ),
        initialize_param(
            "y_C_my",
            model_data.index_y,
            index_p,
            index_k_new,
            model_data.index_t,
        ),
        initialize_param("x_R_my", model_data.index_y, index_p, index_k_existing),
        initialize_param("x_C_my", model_data.index_y, index_p, index_k_new),
        read_param(
            "o_E_my",
            input_filename,
            "VariableCostOldIPPmy",
            model_data.index_t,
            [model_data.index_y, index_p, index_k_existing],
        ),
        read_param(
            "o_C_my",
            input_filename,
            "VariableCostNewIPPmy",
            model_data.index_t,
            [model_data.index_y, index_p, index_k_new],
        ),
        initialize_param("miu_my", model_data.index_y, model_data.index_t),
        initialize_param("LMP_my", model_data.index_y, model_data.index_t),
        initialize_param(
            "eta_my",
            model_data.index_y,
            index_p,
            index_k_existing,
            model_data.index_t,
        ),
        initialize_param(
            "lambda_my",
            model_data.index_y,
            index_p,
            index_k_new,
            model_data.index_t,
        ),
        initialize_param(
            "u_y_E_my",
            model_data.index_y,
            index_p,
            index_k_existing,
            model_data.index_t,
        ),
        initialize_param(
            "u_y_C_my",
            model_data.index_y,
            index_p,
            index_k_new,
            model_data.index_t,
        ),
        initialize_param("u_miu_my", model_data.index_y, model_data.index_t),
        initialize_param(
            "u_eta_my",
            model_data.index_y,
            index_p,
            index_k_existing,
            model_data.index_t,
        ),
        initialize_param(
            "u_lambda_my",
            model_data.index_y,
            index_p,
            index_k_new,
            model_data.index_t,
        ),
        initialize_param("x_R_cumu", index_p, index_k_existing),
        initialize_param("x_C_cumu", index_p, index_k_new),
        initialize_param("L_R", model_data.index_y, index_k_existing, value = 0.0),
        initialize_param("L_C", model_data.index_y, index_k_new, value = 0.0),
        initialize_param("x_R_my_st1", model_data.index_y, index_k_existing),
        initialize_param("x_C_my_st1", model_data.index_y, index_k_new),
        initialize_param("x_R_my_st2", model_data.index_y, index_k_existing),
        initialize_param("x_C_my_st2", model_data.index_y, index_k_new),
        ParamScalar("obj_st1", 1.0),
        ParamScalar("obj_st2", 1.0),
        ParamScalar("obj_lower_bound", -Inf),
        ParamScalar("obj_upper_bound", Inf),
        ParamScalar("obj_feasible", -Inf),
        initialize_param(
            "y_E_my_temp",
            model_data.index_y,
            index_p,
            index_k_existing,
            model_data.index_t,
        ),
        initialize_param(
            "y_C_my_temp",
            model_data.index_y,
            index_p,
            index_k_new,
            model_data.index_t,
        ),
        initialize_param("x_R_my_temp", model_data.index_y, index_p, index_k_existing),
        initialize_param("x_C_my_temp", model_data.index_y, index_p, index_k_new),
        initialize_param("miu_my_temp", model_data.index_y, model_data.index_t),
        initialize_param("LMP_my_temp", model_data.index_y, model_data.index_t),
        NetCONE,
        DC_length,
        read_param(
            "capacity_credit_E_my",
            input_filename,
            "CapacityCredit_old",
            index_k_existing,
            [model_data.index_y],
        ),
        read_param(
            "capacity_credit_C_my",
            input_filename,
            "CapacityCredit_new",
            index_k_new,
            [model_data.index_y],
        ),
        initialize_param("capacity_price", model_data.index_y),
        initialize_param("capacity_price_my_temp", model_data.index_y),
        initialize_param("Net_Load_my", model_data.index_y, model_data.index_t),
        initialize_param("Max_Net_Load_my", model_data.index_y),
        initialize_param("Reserve_req_my", model_data.index_y),
        initialize_param("Capacity_slope_my", model_data.index_y),
        initialize_param("Capacity_intercept_my", model_data.index_y),
        initialize_param("ucap_temp", model_data.index_y, index_p),
        initialize_param("ucap", model_data.index_y, index_p),
        read_param("RPS", input_filename, "RPS", model_data.index_y),
        read_param(
            "emission_rate_E_my",
            input_filename,
            "EmissionRateOldIPPmy",
            index_k_existing,
            [model_data.index_y, index_p],
        ),
        read_param(
            "emission_rate_C_my",
            input_filename,
            "EmissionRateNewIPPmy",
            index_k_new,
            [model_data.index_y, index_p],
        ),
    )
end

get_id(x::IPPGroup) = x.id

function solve_agent_problem!(
    ipps::IPPGroup,
    ipp_opts::AgentOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{VerticallyIntegratedUtility},
    agent_store::AgentStore,
    w_iter,
)
    return 0.0
end

# Lagrange decomposition of the IPP's problem
# TODO: Need to add Green Power Subscription to the Lagrange Decomposition approach!
# Subproblem 1 (Annual Investment & Retirement)
function Lagrange_Sub_Investment_Retirement_Cap(
    ipp::IPPGroup,
    ipp_opts::IPPOptions{LagrangeDecomposition},
    p_star,
    model_data::HEMData,
    hem_opts::HEMOptions{WholesaleMarket},
    agent_store::AgentStore,
    w_iter,
)
    regulator = get_agent(Regulator, agent_store)
    customers = get_agent(CustomerGroup, agent_store)

    WMDER_IPP = get_new_jump_model(hem_opts.NLP_solver)
    set_optimizer_attribute(WMDER_IPP, "print_level", 0)
    # get_optimizer_attribute(WMDER_IPP, "MIPTOL")
    # set_optimizer_attributes(WMDER_IPP, "tol" => 1e-6, "max_iter" => 500)

    # Define positive variables
    @variable(WMDER_IPP, 0 <= x_C[model_data.index_y, ipp.index_k_new] <= 5000)
    @variable(WMDER_IPP, x_R[model_data.index_y, ipp.index_k_existing] >= 0)

    for p in ipp.index_p
        for y in model_data.index_y
            if y == last(model_data.index_y.elements)
                ipp.pvf_onm(y, p, :) .= ipp.pvf_cap(y, p) / ipp.CRF_default(p)
            else
                ipp.pvf_onm(y, p, :) .= ipp.pvf_cap(y, p)
            end
        end
    end

    # adding capacity market parameters
    fill!(ipp.Net_Load_my, NaN)
    for y in model_data.index_y, t in model_data.index_t
        ipp.Net_Load_my(y, t, :) .=
            sum(customers.gamma(h) * customers.d_my(y, h, t) for h in model_data.index_h) +
            ipp.eximport_my(y, t) - sum(
                customers.rho_DG(h, m, t) * customers.x_DG_E_my(y, h, m) for
                h in model_data.index_h, m in customers.index_m
            ) - sum(
                customers.rho_DG(h, m, t) * sum(
                    customers.x_DG_new_my(Symbol(Int(y_symbol)), h, m) for y_symbol in
                    model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                ) for h in model_data.index_h, m in customers.index_m
            )
    end
    fill!(ipp.Max_Net_Load_my, NaN)
    for y in model_data.index_y
        ipp.Max_Net_Load_my(y, :) .=
            findmax(Dict(t => ipp.Net_Load_my(y, t) for t in model_data.index_t))[1]
    end

    Max_Net_Load_my_index = KeyedArray(
        [
            findmax(Dict(t => ipp.Net_Load_my(y, t) for t in model_data.index_t))[2] for
            y in model_data.index_y
        ];
        [get_pair(model_data.index_y)]...
    )

    fill!(ipp.capacity_credit_E_my, NaN)
    for y in model_data.index_y, k in ipp.index_k_existing
        ipp.capacity_credit_E_my(y, k, :) .= ipp.rho_E_my[p_star, k, Max_Net_Load_my_index(y)]
    end
    fill!(ipp.capacity_credit_C_my, NaN)
    for y in model_data.index_y, k in ipp.index_k_new
        ipp.capacity_credit_C_my(y, k, :) .= ipp.rho_C_my[p_star, k, Max_Net_Load_my_index(y)]
    end
    fill!(ipp.Reserve_req_my, NaN)
    for y in model_data.index_y
        ipp.Reserve_req_my(y, :) .= (1 + regulator.r) * ipp.Max_Net_Load_my(y)
    end

    for y in model_data.index_y
        ipp.Capacity_slope_my(y, :) .=
            -ipp.NetCONE(y) / (ipp.DC_length(y) * ipp.Reserve_req_my(y))
        ipp.Capacity_intercept_my(y, :) .=
            -ipp.Capacity_slope_my(y) * ipp.Reserve_req_my(y) + ipp.NetCONE(y)
    end

    #####################################################################

    UCAP_p_star =
        y -> begin
            sum(
                ipp.capacity_credit_E_my(y, k) * (
                    ipp.x_E_my(p_star, k) - sum(
                        x_R[Symbol(Int(y_symbol)), k] for y_symbol in
                        model_data.year(first(model_data.index_y)):model_data.year(y)
                    ) - ipp.x_R_cumu(p_star, k)
                ) for k in ipp.index_k_existing
            ) + sum(
                ipp.capacity_credit_C_my(y, k) * (
                    sum(
                        x_C[Symbol(Int(y_symbol)), k] for y_symbol in
                        model_data.year(first(model_data.index_y)):model_data.year(y)
                    ) + ipp.x_C_cumu(p_star, k)
                ) for k in ipp.index_k_new
            )
        end

    if length(ipp.index_p) >= 2
        UCAP_total =
            y -> begin
                UCAP_p_star(y) +
                sum(
                    ipp.capacity_credit_E_my(y, k) * (
                        ipp.x_E_my(p, k) - sum(
                            ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_R_cumu(p, k)
                    ) for k in ipp.index_k_existing,
                    p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))]
                ) +
                sum(
                    ipp.capacity_credit_C_my(y, k) * (
                        sum(
                            ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p, k)
                    ) for k in ipp.index_k_new,
                    p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))]
                )
            end
    else
        UCAP_total = y -> begin
            UCAP_p_star(y)
        end
    end

    objective_function = begin
        sum(
            # capacity revenue
            ipp.pvf_onm(y, p_star) * (
                UCAP_p_star(y) * (
                    ipp.Capacity_intercept_my(y) + ipp.Capacity_slope_my(y) * UCAP_total(y)
                )
            ) -
            # fixed costs
            #   fom * (cap exist - cap retiring) for gen type
            ipp.pvf_onm(y, p_star) * sum(
                ipp.fom_E_my(y, p_star, k) * (
                    ipp.x_E_my(p_star, k) - sum(
                        x_R[Symbol(Int(y_symbol)), k] for y_symbol in
                        model_data.year(first(model_data.index_y)):model_data.year(y)
                    )
                ) for k in ipp.index_k_existing
            ) -
            # fixed costs
            #   fom * cap new for gen type
            sum(
                ipp.fom_C_my(y, p_star, k) *
                x_C[y, k] *
                sum(
                    ipp.pvf_onm[Symbol(Int(y_symbol)), p_star] for y_symbol in
                    model_data.year(y):model_data.year(last(model_data.index_y))
                ) for k in ipp.index_k_new
            ) -
            # fixed costs
            #   capex * cap new for gen type
            ipp.pvf_cap(y, p_star) *
            sum(ipp.CapEx_my(y, p_star, k) * x_C[y, k] for k in ipp.index_k_new) +
            # lagrange multiplier
            sum(ipp.L_R(y, k) * x_R[y, k] for k in ipp.index_k_existing) +
            sum(ipp.L_C(y, k) * x_C[y, k] for k in ipp.index_k_new)

            for y in model_data.index_y
        )
    end

    @objective(WMDER_IPP, Max, objective_function)

    @constraint(
        WMDER_IPP,
        Eq_sigma[y in model_data.index_y, k in ipp.index_k_existing],
        ipp.x_E_my(p_star, k) - sum(
            x_R[Symbol(Int(y_symbol)), k] for
            y_symbol in model_data.year(first(model_data.index_y)):model_data.year(y)
        ) - ipp.x_R_cumu(p_star, k) >= 0
    )

    if length(ipp.index_p) >= 2
        planning_reserves_cap =
            y -> begin
                # bulk generation available capacity at time t
                sum(
                    ipp.capacity_credit_E_my(y, k) * (
                        ipp.x_E_my(p_star, k) - sum(
                            x_R[Symbol(Int(y_symbol)), k] for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_R_cumu(p_star, k)
                    ) for k in ipp.index_k_existing
                ) +
                sum(
                    ipp.capacity_credit_C_my(y, k) * (
                        sum(
                            x_C[Symbol(Int(y_symbol)), k] for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p_star, k)
                    ) for k in ipp.index_k_new
                ) +
                sum(
                    ipp.capacity_credit_E_my(y, k) * (
                        ipp.x_E_my(p, k) - sum(
                            ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_R_cumu(p, k)
                    ) for k in ipp.index_k_existing,
                    p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))]
                ) +
                sum(
                    ipp.capacity_credit_C_my(y, k) * (
                        sum(
                            ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p, k)
                    ) for k in ipp.index_k_new,
                    p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))]
                ) -
                # net_load plus planning reserve
                ipp.Reserve_req_my(y)
            end
    else
        planning_reserves_cap =
            y -> begin
                # bulk generation available capacity at time t
                sum(
                    ipp.capacity_credit_E_my(y, k) * (
                        ipp.x_E_my(p_star, k) - sum(
                            x_R[Symbol(Int(y_symbol)), k] for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_R_cumu(p_star, k)
                    ) for k in ipp.index_k_existing
                ) + sum(
                    ipp.capacity_credit_C_my(y, k) * (
                        sum(
                            x_C[Symbol(Int(y_symbol)), k] for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p_star, k)
                    ) for k in ipp.index_k_new
                ) -
                # net_load plus planning reserve
                ipp.Reserve_req_my(y)
            end
    end
    @constraint(
        WMDER_IPP,
        Eq_xi_cap[y in model_data.index_y],
        planning_reserves_cap(y) >= 0
    )

    if length(ipp.index_p) >= 2
        planning_reserves =
            (y, t) -> begin
                # bulk generation available capacity at time t
                sum(
                    ipp.rho_E_my(p_star, k, t) * (
                        ipp.x_E_my(p_star, k) - sum(
                            x_R[Symbol(Int(y_symbol)), k] for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_R_cumu(p_star, k)
                    ) for k in ipp.index_k_existing
                ) +
                sum(
                    ipp.rho_C_my(p_star, k, t) * (
                        sum(
                            x_C[Symbol(Int(y_symbol)), k] for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p_star, k)
                    ) for k in ipp.index_k_new
                ) +
                sum(
                    ipp.rho_E_my[p, k, t] * (
                        ipp.x_E_my(p, k) - sum(
                            ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_R_cumu(p, k)
                    ) for k in ipp.index_k_existing,
                    p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))]
                ) +
                sum(
                    ipp.rho_C_my[p, k, t] * (
                        sum(
                            ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p, k)
                    ) for k in ipp.index_k_new,
                    p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))]
                ) -
                # net_load plus planning reserve
                (1 + regulator.r) * (
                    sum(
                        customers.gamma(h) * customers.d_my(y, h, t) for
                        h in model_data.index_h
                    ) + ipp.eximport_my(y, t) - sum(
                        customers.rho_DG(h, m, t) * customers.x_DG_E_my(y, h, m) for
                        h in model_data.index_h, m in customers.index_m
                    ) - sum(
                        customers.rho_DG(h, m, t) * sum(
                            customers.x_DG_new_my(Symbol(Int(y_symbol)), h, m) for
                            y_symbol in
                            model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                        ) for h in model_data.index_h, m in customers.index_m
                    )
                )
            end
    else
        planning_reserves =
            (y, t) -> begin
                # bulk generation available capacity at time t
                sum(
                    ipp.rho_E_my(p_star, k, t) * (
                        ipp.x_E_my(p_star, k) - sum(
                            x_R[Symbol(Int(y_symbol)), k] for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_R_cumu(p_star, k)
                    ) for k in ipp.index_k_existing
                ) + sum(
                    ipp.rho_C_my(p_star, k, t) * (
                        sum(
                            x_C[Symbol(Int(y_symbol)), k] for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p_star, k)
                    ) for k in ipp.index_k_new
                ) -
                # net_load plus planning reserve
                (1 + regulator.r) * (
                    sum(
                        customers.gamma(h) * customers.d_my(y, h, t) for
                        h in model_data.index_h
                    ) + ipp.eximport_my(y, t) - sum(
                        customers.rho_DG(h, m, t) * customers.x_DG_E_my(y, h, m) for
                        h in model_data.index_h, m in customers.index_m
                    ) - sum(
                        customers.rho_DG(h, m, t) * sum(
                            customers.x_DG_new_my(Symbol(Int(y_symbol)), h, m) for
                            y_symbol in
                            model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                        ) for h in model_data.index_h, m in customers.index_m
                    )
                )
            end
    end
    @constraint(
        WMDER_IPP,
        Eq_xi[y in model_data.index_y, t in model_data.index_t],
        planning_reserves(y, t) >= 0
    )

    TimerOutputs.@timeit HEM_TIMER "optimize! Lagrange_Sub_Investment_Retirement_Cap" begin
        optimize!(WMDER_IPP)
    end

    for y in model_data.index_y, k in ipp.index_k_existing
        ipp.x_R_my_st1(y, k, :) .= value.(x_R[y, k])
    end
    for y in model_data.index_y, k in ipp.index_k_new
        ipp.x_C_my_st1(y, k, :) .= value.(x_C[y, k])
    end

    update!(ipp.obj_st1, objective_value(WMDER_IPP))

    # return stage 1 investment & retirement solution
    return ipp.obj_st1
end

# Lagrange decomposition of the IPP's problem
# Subproblem 2 (Hourly dispatch)
function Lagrange_Sub_Dispatch_Cap(
    ipp::IPPGroup,
    ipp_opts::IPPOptions{LagrangeDecomposition},
    p_star,
    model_data::HEMData,
    hem_opts::HEMOptions{WholesaleMarket},
    agent_store::AgentStore,
    w_iter,
)
    regulator = get_agent(Regulator, agent_store)
    customers = get_agent(CustomerGroup, agent_store)

    WMDER_IPP = get_new_jump_model(hem_opts.MIP_solver)
    # set_optimizer_attribute(WMDER_IPP, "OUTPUTLOG", 0)

    # WMDER_IPP = Model(()->Solver.Optimizer(OUTPUTLOG = 0))

    # Define positive variables
    @variable(WMDER_IPP, 0 <= x_C[model_data.index_y, ipp.index_k_new] <= 5000)
    @variable(WMDER_IPP, x_R[model_data.index_y, ipp.index_k_existing] >= 0)
    @variable(
        WMDER_IPP,
        y_E[model_data.index_y, ipp.index_p, ipp.index_k_existing, model_data.index_t] >= 0
    )
    @variable(
        WMDER_IPP,
        y_C[model_data.index_y, ipp.index_p, ipp.index_k_new, model_data.index_t] >= 0
    )
    @variable(WMDER_IPP, miu[model_data.index_y, model_data.index_t] >= 0)
    @variable(
        WMDER_IPP,
        eta[model_data.index_y, ipp.index_p, ipp.index_k_existing, model_data.index_t] >= 0
    )
    @variable(
        WMDER_IPP,
        lambda[model_data.index_y, ipp.index_p, ipp.index_k_new, model_data.index_t] >= 0
    )

    @variable(
        WMDER_IPP,
        u_y_E[model_data.index_y, ipp.index_p, ipp.index_k_existing, model_data.index_t],
        Bin
    )
    @variable(
        WMDER_IPP,
        u_y_C[model_data.index_y, ipp.index_p, ipp.index_k_new, model_data.index_t],
        Bin
    )
    @variable(WMDER_IPP, u_miu[model_data.index_y, model_data.index_t], Bin)
    @variable(
        WMDER_IPP,
        u_eta[model_data.index_y, ipp.index_p, ipp.index_k_existing, model_data.index_t],
        Bin
    )
    @variable(
        WMDER_IPP,
        u_lambda[model_data.index_y, ipp.index_p, ipp.index_k_new, model_data.index_t],
        Bin
    )

    for p in ipp.index_p
        for y in model_data.index_y
            if y == last(model_data.index_y.elements)
                ipp.pvf_onm(y, p, :) .= ipp.pvf_cap(y, p) / ipp.CRF_default(p)
            else
                ipp.pvf_onm(y, p, :) .= ipp.pvf_cap(y, p)
            end
        end
    end

    objective_function = begin
        sum(
            # Linearized revenue term
            ipp.pvf_onm(y, p_star) * (
                sum(
                    miu[y, t] * (
                        sum(
                            customers.gamma(h) * customers.d_my(y, h, t) for
                            h in model_data.index_h
                        ) + ipp.eximport_my(y, t) - sum(
                            customers.rho_DG(h, m, t) * customers.x_DG_E_my(y, h, m) for
                            h in model_data.index_h, m in customers.index_m
                        ) - sum(
                            customers.rho_DG(h, m, t) * sum(
                                customers.x_DG_new_my(Symbol(Int(y_symbol)), h, m) for
                                y_symbol in
                                model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                            ) for h in model_data.index_h, m in customers.index_m
                        )
                    ) for t in model_data.index_t
                ) - (
                    sum(
                        eta[y, p, k, t] *
                        ipp.rho_E_my[p, k, t] *
                        (
                            ipp.x_E_my(p, k) - sum(
                                ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) - ipp.x_R_cumu(p, k)
                        ) for t in model_data.index_t, k in ipp.index_k_existing,
                        p in ipp.index_p
                    ) - sum(
                        eta[y, p_star, k, t] *
                        ipp.rho_E_my(p_star, k, t) *
                        (
                            ipp.x_E_my(p_star, k) - sum(
                                ipp.x_R_my[Symbol(Int(y_symbol)), p_star, k] for
                                y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) - ipp.x_R_cumu(p_star, k)
                        ) for t in model_data.index_t, k in ipp.index_k_existing
                    )
                ) - (
                    sum(
                        lambda[y, p, k, t] *
                        ipp.rho_C_my[p, k, t] *
                        (
                            sum(
                                ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) + ipp.x_C_cumu(p, k)
                        ) for t in model_data.index_t, k in ipp.index_k_new,
                        p in ipp.index_p
                    ) - sum(
                        lambda[y, p_star, k, t] *
                        ipp.rho_C_my(p_star, k, t) *
                        (
                            sum(
                                ipp.x_C_my[Symbol(Int(y_symbol)), p_star, k] for
                                y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) + ipp.x_C_cumu(p_star, k)
                        ) for t in model_data.index_t, k in ipp.index_k_new
                    )
                ) - (
                    sum(
                        model_data.omega(t) * (ipp.v_E_my(y, p, k, t) * y_E[y, p, k, t]) for
                        t in model_data.index_t, k in ipp.index_k_existing,
                        p in ipp.index_p
                    ) - sum(
                        model_data.omega(t) *
                        (ipp.v_E_my[y, p_star, k, t] * y_E[y, p_star, k, t]) for
                        t in model_data.index_t, k in ipp.index_k_existing
                    )
                ) - (
                    sum(
                        model_data.omega(t) * (ipp.v_C_my(y, p, k, t) * y_C[y, p, k, t]) for
                        t in model_data.index_t, k in ipp.index_k_new, p in ipp.index_p
                    ) - sum(
                        model_data.omega(t) *
                        (ipp.v_C_my[y, p_star, k, t] * y_C[y, p_star, k, t]) for
                        t in model_data.index_t, k in ipp.index_k_new
                    )
                ) -

                # generation costs (this includes terms due to revenue linearization)
                #   num hrs * ((fuel + vom) * gen existing + (fuel + vom) * gen new) for t and gen type
                sum(
                    model_data.omega(t) *
                    (ipp.v_E_my[y, p_star, k, t] * y_E[y, p_star, k, t]) for
                    t in model_data.index_t, k in ipp.index_k_existing
                ) - sum(
                    model_data.omega(t) *
                    (ipp.v_C_my[y, p_star, k, t] * y_C[y, p_star, k, t]) for
                    t in model_data.index_t, k in ipp.index_k_new
                )
            ) -

            # lagrange multiplier
            sum(ipp.L_R(y, k) * x_R[y, k] for k in ipp.index_k_existing) -
            sum(ipp.L_C(y, k) * x_C[y, k] for k in ipp.index_k_new)

            for y in model_data.index_y
        )
    end

    @objective(WMDER_IPP, Max, objective_function)

    @constraint(
        WMDER_IPP,
        Eq_sigma[y in model_data.index_y, k in ipp.index_k_existing],
        ipp.x_E_my(p_star, k) - sum(
            x_R[Symbol(Int(y_symbol)), k] for
            y_symbol in model_data.year(first(model_data.index_y)):model_data.year(y)
        ) - ipp.x_R_cumu(p_star, k) >= 0
    )

    @constraint(
        WMDER_IPP,
        Eq_y_E_p_1[
            y in model_data.index_y,
            p in ipp.index_p,
            k in ipp.index_k_existing,
            t in model_data.index_t,
        ],
        y_E[y, p, k, t] <= ipp.B1GM.value * u_y_E[y, p, k, t]
    )
    @constraint(
        WMDER_IPP,
        Eq_y_E_p_2[
            y in model_data.index_y,
            p in ipp.index_p,
            k in ipp.index_k_existing,
            t in model_data.index_t,
        ],
        model_data.omega(t) * ipp.v_E_my(y, p, k, t) - miu[y, t] + eta[y, p, k, t] <=
        ipp.B1GM.value * (1 - u_y_E[y, p, k, t])
    )
    @constraint(
        WMDER_IPP,
        Eq_y_E_p_3[
            y in model_data.index_y,
            p in ipp.index_p,
            k in ipp.index_k_existing,
            t in model_data.index_t,
        ],
        model_data.omega(t) * ipp.v_E_my(y, p, k, t) - miu[y, t] + eta[y, p, k, t] >= 0
    )

    @constraint(
        WMDER_IPP,
        Eq_y_C_p_1[
            y in model_data.index_y,
            p in ipp.index_p,
            k in ipp.index_k_new,
            t in model_data.index_t,
        ],
        y_C[y, p, k, t] <= ipp.B1GM.value * u_y_C[y, p, k, t]
    )
    @constraint(
        WMDER_IPP,
        Eq_y_C_p_2[
            y in model_data.index_y,
            p in ipp.index_p,
            k in ipp.index_k_new,
            t in model_data.index_t,
        ],
        model_data.omega(t) * ipp.v_C_my(y, p, k, t) - miu[y, t] + lambda[y, p, k, t] <=
        ipp.B1GM.value * (1 - u_y_C[y, p, k, t])
    )
    @constraint(
        WMDER_IPP,
        Eq_y_C_p_3[
            y in model_data.index_y,
            p in ipp.index_p,
            k in ipp.index_k_new,
            t in model_data.index_t,
        ],
        model_data.omega(t) * ipp.v_C_my(y, p, k, t) - miu[y, t] + lambda[y, p, k, t] >= 0
    )

    supply_demand_balance =
        (y, t) -> begin
            # bulk generation at time t
            sum(y_E[y, p, k, t] for k in ipp.index_k_existing, p in ipp.index_p) +
            sum(y_C[y, p, k, t] for k in ipp.index_k_new, p in ipp.index_p) -
            # demand at time t
            sum(customers.gamma(h) * customers.d_my(y, h, t) for h in model_data.index_h) - ipp.eximport_my(y, t) +
            # existing DG generation at time t
            sum(
                customers.rho_DG(h, m, t) * customers.x_DG_E_my(y, h, m) for
                h in model_data.index_h, m in customers.index_m
            ) +
            # new DG generation at time t
            sum(
                customers.rho_DG(h, m, t) * sum(
                    customers.x_DG_new_my(Symbol(Int(y_symbol)), h, m) for y_symbol in
                    model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                ) for h in model_data.index_h, m in customers.index_m
            )
        end

    @constraint(
        WMDER_IPP,
        Eq_miu_1[y in model_data.index_y, t in model_data.index_t],
        miu[y, t] <= ipp.B1GM.value * u_miu[y, t]
    )
    @constraint(
        WMDER_IPP,
        Eq_miu_2[y in model_data.index_y, t in model_data.index_t],
        supply_demand_balance(y, t) <= ipp.B1GM.value * (1 - u_miu[y, t])
    )
    @constraint(
        WMDER_IPP,
        Eq_miu_3[y in model_data.index_y, t in model_data.index_t],
        supply_demand_balance(y, t) >= 0
    )

    @constraint(
        WMDER_IPP,
        Eq_eta_p_star_1[
            y in model_data.index_y,
            k in ipp.index_k_existing,
            t in model_data.index_t,
        ],
        eta[y, p_star, k, t] <= ipp.B1GM.value * u_eta[y, p_star, k, t]
    )
    @constraint(
        WMDER_IPP,
        Eq_eta_p_star_2[
            y in model_data.index_y,
            k in ipp.index_k_existing,
            t in model_data.index_t,
        ],
        ipp.rho_E_my(p_star, k, t) * (
            ipp.x_E_my(p_star, k) - sum(
                x_R[Symbol(Int(y_symbol)), k] for
                y_symbol in model_data.year(first(model_data.index_y)):model_data.year(y)
            ) - ipp.x_R_cumu(p_star, k)
        ) - y_E[y, p_star, k, t] <= ipp.B1GM.value * (1 - u_eta[y, p_star, k, t])
    )
    @constraint(
        WMDER_IPP,
        Eq_eta_p_star_3[
            y in model_data.index_y,
            k in ipp.index_k_existing,
            t in model_data.index_t,
        ],
        ipp.rho_E_my(p_star, k, t) * (
            ipp.x_E_my(p_star, k) - sum(
                x_R[Symbol(Int(y_symbol)), k] for
                y_symbol in model_data.year(first(model_data.index_y)):model_data.year(y)
            ) - ipp.x_R_cumu(p_star, k)
        ) - y_E[y, p_star, k, t] >= 0
    )
    if length(ipp.index_p) >= 2
        @constraint(
            WMDER_IPP,
            Eq_eta_p_other_1[
                y in model_data.index_y,
                p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))],
                k in ipp.index_k_existing,
                t in model_data.index_t,
            ],
            eta[y, p, k, t] <= ipp.B1GM.value * u_eta[y, p, k, t]
        )
        @constraint(
            WMDER_IPP,
            Eq_eta_p_other_2[
                y in model_data.index_y,
                p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))],
                k in ipp.index_k_existing,
                t in model_data.index_t,
            ],
            ipp.rho_E_my[p, k, t] * (
                ipp.x_E_my(p, k) - sum(
                    ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) - ipp.x_R_cumu(p, k)
            ) - y_E[y, p, k, t] <= ipp.B1GM.value * (1 - u_eta[y, p, k, t])
        )
        @constraint(
            WMDER_IPP,
            Eq_eta_p_other_3[
                y in model_data.index_y,
                p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))],
                k in ipp.index_k_existing,
                t in model_data.index_t,
            ],
            ipp.rho_E_my[p, k, t] * (
                ipp.x_E_my(p, k) - sum(
                    ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) - ipp.x_R_cumu(p, k)
            ) - y_E[y, p, k, t] >= 0
        )
    end

    @constraint(
        WMDER_IPP,
        Eq_lambda_p_star_1[
            y in model_data.index_y,
            k in ipp.index_k_new,
            t in model_data.index_t,
        ],
        lambda[y, p_star, k, t] <= ipp.B1GM.value * u_lambda[y, p_star, k, t]
    )
    @constraint(
        WMDER_IPP,
        Eq_lambda_p_star_2[
            y in model_data.index_y,
            k in ipp.index_k_new,
            t in model_data.index_t,
        ],
        ipp.rho_C_my(p_star, k, t) * (
            sum(
                x_C[Symbol(Int(y_symbol)), k] for
                y_symbol in model_data.year(first(model_data.index_y)):model_data.year(y)
            ) + ipp.x_C_cumu(p_star, k)
        ) - y_C[y, p_star, k, t] <= ipp.B1GM.value * (1 - u_lambda[y, p_star, k, t])
    )
    @constraint(
        WMDER_IPP,
        Eq_lambda_p_star_3[
            y in model_data.index_y,
            k in ipp.index_k_new,
            t in model_data.index_t,
        ],
        ipp.rho_C_my(p_star, k, t) * (
            sum(
                x_C[Symbol(Int(y_symbol)), k] for
                y_symbol in model_data.year(first(model_data.index_y)):model_data.year(y)
            ) + ipp.x_C_cumu(p_star, k)
        ) - y_C[y, p_star, k, t] >= 0
    )
    if length(ipp.index_p) >= 2
        @constraint(
            WMDER_IPP,
            Eq_lambda_p_other_1[
                y in model_data.index_y,
                p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))],
                k in ipp.index_k_new,
                t in model_data.index_t,
            ],
            lambda[y, p, k, t] <= ipp.B1GM.value * u_lambda[y, p, k, t]
        )
        @constraint(
            WMDER_IPP,
            Eq_lambda_p_other_2[
                y in model_data.index_y,
                p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))],
                k in ipp.index_k_new,
                t in model_data.index_t,
            ],
            ipp.rho_C_my[p, k, t] * (
                sum(
                    ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) + ipp.x_C_cumu(p, k)
            ) - y_C[y, p, k, t] <= ipp.B1GM.value * (1 - u_lambda[y, p, k, t])
        )
        @constraint(
            WMDER_IPP,
            Eq_lambda_p_other_3[
                y in model_data.index_y,
                p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))],
                k in ipp.index_k_new,
                t in model_data.index_t,
            ],
            ipp.rho_C_my[p, k, t] * (
                sum(
                    ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) + ipp.x_C_cumu(p, k)
            ) - y_C[y, p, k, t] >= 0
        )
    end

    # without this constraint, miu will essentially be capped at B1GM
    if length(ipp.index_p) >= 2
        planning_reserves =
            (y, t) -> begin
                # bulk generation available capacity at time t
                sum(
                    ipp.rho_E_my(p_star, k, t) * (
                        ipp.x_E_my(p_star, k) - sum(
                            x_R[Symbol(Int(y_symbol)), k] for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_R_cumu(p_star, k)
                    ) for k in ipp.index_k_existing
                ) +
                sum(
                    ipp.rho_C_my(p_star, k, t) * (
                        sum(
                            x_C[Symbol(Int(y_symbol)), k] for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p_star, k)
                    ) for k in ipp.index_k_new
                ) +
                sum(
                    ipp.rho_E_my[p, k, t] * (
                        ipp.x_E_my(p, k) - sum(
                            ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_R_cumu(p, k)
                    ) for k in ipp.index_k_existing,
                    p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))]
                ) +
                sum(
                    ipp.rho_C_my[p, k, t] * (
                        sum(
                            ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p, k)
                    ) for k in ipp.index_k_new,
                    p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))]
                ) -
                # net_load plus planning reserve
                (1 + regulator.r) * (
                    sum(
                        customers.gamma(h) * customers.d_my(y, h, t) for
                        h in model_data.index_h
                    ) + ipp.eximport_my(y, t) - sum(
                        customers.rho_DG(h, m, t) * customers.x_DG_E_my(y, h, m) for
                        h in model_data.index_h, m in customers.index_m
                    ) - sum(
                        customers.rho_DG(h, m, t) * sum(
                            customers.x_DG_new_my(Symbol(Int(y_symbol)), h, m) for
                            y_symbol in
                            model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                        ) for h in model_data.index_h, m in customers.index_m
                    )
                )
            end
    else
        planning_reserves =
            (y, t) -> begin
                # bulk generation available capacity at time t
                sum(
                    ipp.rho_E_my(p_star, k, t) * (
                        ipp.x_E_my(p_star, k) - sum(
                            x_R[Symbol(Int(y_symbol)), k] for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_R_cumu(p_star, k)
                    ) for k in ipp.index_k_existing
                ) + sum(
                    ipp.rho_C_my(p_star, k, t) * (
                        sum(
                            x_C[Symbol(Int(y_symbol)), k] for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p_star, k)
                    ) for k in ipp.index_k_new
                ) -
                # net_load plus planning reserve
                (1 + regulator.r) * (
                    sum(
                        customers.gamma(h) * customers.d_my(y, h, t) for
                        h in model_data.index_h
                    ) + ipp.eximport_my(y, t) - sum(
                        customers.rho_DG(h, m, t) * customers.x_DG_E_my(y, h, m) for
                        h in model_data.index_h, m in customers.index_m
                    ) - sum(
                        customers.rho_DG(h, m, t) * sum(
                            customers.x_DG_new_my(Symbol(Int(y_symbol)), h, m) for
                            y_symbol in
                            model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                        ) for h in model_data.index_h, m in customers.index_m
                    )
                )
            end
    end
    @constraint(
        WMDER_IPP,
        Eq_xi[y in model_data.index_y, t in model_data.index_t],
        planning_reserves(y, t) >= 0
    )

    if length(ipp.index_p) >= 2
        planning_reserves_cap =
            y -> begin
                # bulk generation available capacity at time t
                sum(
                    ipp.capacity_credit_E_my(y, k) * (
                        ipp.x_E_my(p_star, k) - sum(
                            x_R[Symbol(Int(y_symbol)), k] for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_R_cumu(p_star, k)
                    ) for k in ipp.index_k_existing
                ) +
                sum(
                    ipp.capacity_credit_C_my(y, k) * (
                        sum(
                            x_C[Symbol(Int(y_symbol)), k] for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p_star, k)
                    ) for k in ipp.index_k_new
                ) +
                sum(
                    ipp.capacity_credit_E_my(y, k) * (
                        ipp.x_E_my(p, k) - sum(
                            ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_R_cumu(p, k)
                    ) for k in ipp.index_k_existing,
                    p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))]
                ) +
                sum(
                    ipp.capacity_credit_C_my(y, k) * (
                        sum(
                            ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p, k)
                    ) for k in ipp.index_k_new,
                    p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))]
                ) -
                # net_load plus planning reserve
                ipp.Reserve_req_my(y)
            end
    else
        planning_reserves_cap =
            y -> begin
                # bulk generation available capacity at time t
                sum(
                    ipp.capacity_credit_E_my(y, k) * (
                        ipp.x_E_my(p_star, k) - sum(
                            x_R[Symbol(Int(y_symbol)), k] for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_R_cumu(p_star, k)
                    ) for k in ipp.index_k_existing
                ) + sum(
                    ipp.capacity_credit_C_my(y, k) * (
                        sum(
                            x_C[Symbol(Int(y_symbol)), k] for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p_star, k)
                    ) for k in ipp.index_k_new
                ) -
                # net_load plus planning reserve
                ipp.Reserve_req_my(y)
            end
    end
    @constraint(
        WMDER_IPP,
        Eq_xi_cap[y in model_data.index_y],
        planning_reserves_cap(y) >= 0
    )

    # RPS constraint
    @constraint(
        WMDER_IPP,
        Eq_rps[y in model_data.index_y],
        sum(
            model_data.omega(t) * (y_E[y, p, rps, t] + y_C[y, p, rps, t]) for
            rps in ipp.index_rps, t in model_data.index_t, p in ipp.index_p
        ) -
        ipp.RPS(y) *
        sum(model_data.omega(t) * ipp.Net_Load_my(y, t) for t in model_data.index_t) >= 0
    )

    TimerOutputs.@timeit HEM_TIMER "optimize! Lagrange_Sub_Investment_Retirement_Cap" begin
        optimize!(WMDER_IPP)
    end

    # return MOI.get(WMDER_IPP, MOI.TerminationStatus())

    for y in model_data.index_y, k in ipp.index_k_existing
        ipp.x_R_my_st2(y, k, :) .= value.(x_R[y, k])
    end
    for y in model_data.index_y, k in ipp.index_k_new
        ipp.x_C_my_st2(y, k, :) .= value.(x_C[y, k])
    end

    for y in model_data.index_y, k in ipp.index_k_existing
        ipp.x_R_my_temp(y, p_star, k, :) .= value.(x_R[y, k])
    end
    for y in model_data.index_y, k in ipp.index_k_new
        ipp.x_C_my_temp(y, p_star, k, :) .= value.(x_C[y, k])
    end
    for y in model_data.index_y,
        p in ipp.index_p,
        k in ipp.index_k_existing,
        t in model_data.index_t

        ipp.y_E_my_temp(y, p, k, t, :) .= value.(y_E[y, p, k, t])
    end
    for y in model_data.index_y,
        p in ipp.index_p,
        k in ipp.index_k_new,
        t in model_data.index_t

        ipp.y_C_my_temp(y, p, k, t, :) .= value.(y_C[y, p, k, t])
    end
    for y in model_data.index_y, t in model_data.index_t
        ipp.miu_my_temp(y, t, :) .= value.(miu[y, t])
        ipp.LMP_my_temp(y, t, :) .= ipp.miu_my_temp(y, t) / model_data.omega(t)
    end

    update!(ipp.obj_st2, objective_value(WMDER_IPP))

    # return stage 2 investment & retirement solution
    return ipp.obj_st2
end

# Lagrange decomposition of the IPP's problem
# Feasible solution
function Lagrange_Feasible_Cap(
    ipp::IPPGroup,
    ipp_opts::IPPOptions{LagrangeDecomposition},
    p_star,
    model_data::HEMData,
    hem_opts::HEMOptions{WholesaleMarket},
    agent_store::AgentStore,
    w_iter,
)
    regulator = get_agent(Regulator, agent_store)
    customers = get_agent(CustomerGroup, agent_store)

    WMDER_IPP = get_new_jump_model(hem_opts.MIP_solver)
    # set_optimizer_attribute(WMDER_IPP, "OUTPUTLOG", 0)

    # Define positive variables
    @variable(
        WMDER_IPP,
        y_E[model_data.index_y, ipp.index_p, ipp.index_k_existing, model_data.index_t] >= 0
    )
    @variable(
        WMDER_IPP,
        y_C[model_data.index_y, ipp.index_p, ipp.index_k_new, model_data.index_t] >= 0
    )
    @variable(WMDER_IPP, miu[model_data.index_y, model_data.index_t] >= 0)
    @variable(
        WMDER_IPP,
        eta[model_data.index_y, ipp.index_p, ipp.index_k_existing, model_data.index_t] >= 0
    )
    @variable(
        WMDER_IPP,
        lambda[model_data.index_y, ipp.index_p, ipp.index_k_new, model_data.index_t] >= 0
    )

    @variable(
        WMDER_IPP,
        u_y_E[model_data.index_y, ipp.index_p, ipp.index_k_existing, model_data.index_t],
        Bin
    )
    @variable(
        WMDER_IPP,
        u_y_C[model_data.index_y, ipp.index_p, ipp.index_k_new, model_data.index_t],
        Bin
    )
    @variable(WMDER_IPP, u_miu[model_data.index_y, model_data.index_t], Bin)
    @variable(
        WMDER_IPP,
        u_eta[model_data.index_y, ipp.index_p, ipp.index_k_existing, model_data.index_t],
        Bin
    )
    @variable(
        WMDER_IPP,
        u_lambda[model_data.index_y, ipp.index_p, ipp.index_k_new, model_data.index_t],
        Bin
    )

    for p in ipp.index_p
        for y in model_data.index_y
            if y == last(model_data.index_y.elements)
                ipp.pvf_onm(y, p, :) .= ipp.pvf_cap(y, p) / ipp.CRF_default(p)
            else
                ipp.pvf_onm(y, p, :) .= ipp.pvf_cap(y, p)
            end
        end
    end

    objective_function = begin
        sum(
            # Linearized revenue term
            ipp.pvf_onm(y, p_star) * (
                sum(
                    miu[y, t] * (
                        sum(
                            customers.gamma(h) * customers.d_my(y, h, t) for
                            h in model_data.index_h
                        ) + ipp.eximport_my(y, t) - sum(
                            customers.rho_DG(h, m, t) * customers.x_DG_E_my(y, h, m) for
                            h in model_data.index_h, m in customers.index_m
                        ) - sum(
                            customers.rho_DG(h, m, t) * sum(
                                customers.x_DG_new_my(Symbol(Int(y_symbol)), h, m) for
                                y_symbol in
                                model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                            ) for h in model_data.index_h, m in customers.index_m
                        )
                    ) for t in model_data.index_t
                ) - (
                    sum(
                        eta[y, p, k, t] *
                        ipp.rho_E_my[p, k, t] *
                        (
                            ipp.x_E_my(p, k) - sum(
                                ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) - ipp.x_R_cumu(p, k)
                        ) for t in model_data.index_t, k in ipp.index_k_existing,
                        p in ipp.index_p
                    ) - sum(
                        eta[y, p_star, k, t] *
                        ipp.rho_E_my(p_star, k, t) *
                        (
                            ipp.x_E_my(p_star, k) - sum(
                                ipp.x_R_my[Symbol(Int(y_symbol)), p_star, k] for
                                y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) - ipp.x_R_cumu(p_star, k)
                        ) for t in model_data.index_t, k in ipp.index_k_existing
                    )
                ) - (
                    sum(
                        lambda[y, p, k, t] *
                        ipp.rho_C_my[p, k, t] *
                        (
                            sum(
                                ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) + ipp.x_C_cumu(p, k)
                        ) for t in model_data.index_t, k in ipp.index_k_new,
                        p in ipp.index_p
                    ) - sum(
                        lambda[y, p_star, k, t] *
                        ipp.rho_C_my(p_star, k, t) *
                        (
                            sum(
                                ipp.x_C_my[Symbol(Int(y_symbol)), p_star, k] for
                                y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) + ipp.x_C_cumu(p_star, k)
                        ) for t in model_data.index_t, k in ipp.index_k_new
                    )
                ) - (
                    sum(
                        model_data.omega(t) * (ipp.v_E_my(y, p, k, t) * y_E[y, p, k, t]) for
                        t in model_data.index_t, k in ipp.index_k_existing,
                        p in ipp.index_p
                    ) - sum(
                        model_data.omega(t) *
                        (ipp.v_E_my[y, p_star, k, t] * y_E[y, p_star, k, t]) for
                        t in model_data.index_t, k in ipp.index_k_existing
                    )
                ) - (
                    sum(
                        model_data.omega(t) * (ipp.v_C_my(y, p, k, t) * y_C[y, p, k, t]) for
                        t in model_data.index_t, k in ipp.index_k_new, p in ipp.index_p
                    ) - sum(
                        model_data.omega(t) *
                        (ipp.v_C_my[y, p_star, k, t] * y_C[y, p_star, k, t]) for
                        t in model_data.index_t, k in ipp.index_k_new
                    )
                ) -

                # generation costs (this includes terms due to revenue linearization)
                #   num hrs * ((fuel + vom) * gen existing + (fuel + vom) * gen new) for t and gen type
                sum(
                    model_data.omega(t) *
                    (ipp.v_E_my[y, p_star, k, t] * y_E[y, p_star, k, t]) for
                    t in model_data.index_t, k in ipp.index_k_existing
                ) - sum(
                    model_data.omega(t) *
                    (ipp.v_C_my[y, p_star, k, t] * y_C[y, p_star, k, t]) for
                    t in model_data.index_t, k in ipp.index_k_new
                )
            )

            for y in model_data.index_y
        )
    end

    @objective(WMDER_IPP, Max, objective_function)

    @constraint(
        WMDER_IPP,
        Eq_y_E_p_1[
            y in model_data.index_y,
            p in ipp.index_p,
            k in ipp.index_k_existing,
            t in model_data.index_t,
        ],
        y_E[y, p, k, t] <= ipp.B1GM.value * u_y_E[y, p, k, t]
    )
    @constraint(
        WMDER_IPP,
        Eq_y_E_p_2[
            y in model_data.index_y,
            p in ipp.index_p,
            k in ipp.index_k_existing,
            t in model_data.index_t,
        ],
        model_data.omega(t) * ipp.v_E_my(y, p, k, t) - miu[y, t] + eta[y, p, k, t] <=
        ipp.B1GM.value * (1 - u_y_E[y, p, k, t])
    )
    @constraint(
        WMDER_IPP,
        Eq_y_E_p_3[
            y in model_data.index_y,
            p in ipp.index_p,
            k in ipp.index_k_existing,
            t in model_data.index_t,
        ],
        model_data.omega(t) * ipp.v_E_my(y, p, k, t) - miu[y, t] + eta[y, p, k, t] >= 0
    )

    @constraint(
        WMDER_IPP,
        Eq_y_C_p_1[
            y in model_data.index_y,
            p in ipp.index_p,
            k in ipp.index_k_new,
            t in model_data.index_t,
        ],
        y_C[y, p, k, t] <= ipp.B1GM.value * u_y_C[y, p, k, t]
    )
    @constraint(
        WMDER_IPP,
        Eq_y_C_p_2[
            y in model_data.index_y,
            p in ipp.index_p,
            k in ipp.index_k_new,
            t in model_data.index_t,
        ],
        model_data.omega(t) * ipp.v_C_my(y, p, k, t) - miu[y, t] + lambda(y, p, k, t) <=
        ipp.B1GM.value * (1 - u_y_C(y, p, k, t))
    )
    @constraint(
        WMDER_IPP,
        Eq_y_C_p_3[
            y in model_data.index_y,
            p in ipp.index_p,
            k in ipp.index_k_new,
            t in model_data.index_t,
        ],
        model_data.omega(t) * ipp.v_C_my(y, p, k, t) - miu[y, t] + lambda[y, p, k, t] >= 0
    )

    supply_demand_balance =
        (y, t) -> begin
            # bulk generation at time t
            sum(y_E[y, p, k, t] for k in ipp.index_k_existing, p in ipp.index_p) +
            sum(y_C[y, p, k, t] for k in ipp.index_k_new, p in ipp.index_p) -
            # demand at time t
            sum(customers.gamma(h) * customers.d_my(y, h, t) for h in model_data.index_h) - ipp.eximport_my(y, t) +
            # existing DG generation at time t
            sum(
                customers.rho_DG(h, m, t) * customers.x_DG_E_my(y, h, m) for
                h in model_data.index_h, m in customers.index_m
            ) +
            # new DG generation at time t
            sum(
                customers.rho_DG(h, m, t) * sum(
                    customers.x_DG_new_my(Symbol(Int(y_symbol)), h, m) for y_symbol in
                    model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                ) for h in model_data.index_h, m in customers.index_m
            )
        end

    @constraint(
        WMDER_IPP,
        Eq_miu_1[y in model_data.index_y, t in model_data.index_t],
        miu[y, t] <= ipp.B1GM.value * u_miu[y, t]
    )
    @constraint(
        WMDER_IPP,
        Eq_miu_2[y in model_data.index_y, t in model_data.index_t],
        supply_demand_balance(y, t) <= ipp.B1GM.value * (1 - u_miu[y, t])
    )
    @constraint(
        WMDER_IPP,
        Eq_miu_3[y in model_data.index_y, t in model_data.index_t],
        supply_demand_balance(y, t) >= 0
    )

    @constraint(
        WMDER_IPP,
        Eq_eta_p_star_1[
            y in model_data.index_y,
            k in ipp.index_k_existing,
            t in model_data.index_t,
        ],
        eta[y, p_star, k, t] <= ipp.B1GM.value * u_eta[y, p_star, k, t]
    )
    @constraint(
        WMDER_IPP,
        Eq_eta_p_star_2[
            y in model_data.index_y,
            k in ipp.index_k_existing,
            t in model_data.index_t,
        ],
        ipp.rho_E_my(p_star, k, t) * (
            ipp.x_E_my(p_star, k) - sum(
                ipp.x_R_my_st2[Symbol(Int(y_symbol)), k] for
                y_symbol in model_data.year(first(model_data.index_y)):model_data.year(y)
            ) - ipp.x_R_cumu(p_star, k)
        ) - y_E[y, p_star, k, t] <= ipp.B1GM.value * (1 - u_eta[y, p_star, k, t])
    )
    @constraint(
        WMDER_IPP,
        Eq_eta_p_star_3[
            y in model_data.index_y,
            k in ipp.index_k_existing,
            t in model_data.index_t,
        ],
        ipp.rho_E_my(p_star, k, t) * (
            ipp.x_E_my(p_star, k) - sum(
                ipp.x_R_my_st2[Symbol(Int(y_symbol)), k] for
                y_symbol in model_data.year(first(model_data.index_y)):model_data.year(y)
            ) - ipp.x_R_cumu(p_star, k)
        ) - y_E[y, p_star, k, t] >= 0
    )
    if length(ipp.index_p) >= 2
        @constraint(
            WMDER_IPP,
            Eq_eta_p_other_1[
                y in model_data.index_y,
                p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))],
                k in ipp.index_k_existing,
                t in model_data.index_t,
            ],
            eta[y, p, k, t] <= ipp.B1GM.value * u_eta[y, p, k, t]
        )
        @constraint(
            WMDER_IPP,
            Eq_eta_p_other_2[
                y in model_data.index_y,
                p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))],
                k in ipp.index_k_existing,
                t in model_data.index_t,
            ],
            ipp.rho_E_my[p, k, t] * (
                ipp.x_E_my(p, k) - sum(
                    ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) - ipp.x_R_cumu(p, k)
            ) - y_E[y, p, k, t] <= ipp.B1GM.value * (1 - u_eta[y, p, k, t])
        )
        @constraint(
            WMDER_IPP,
            Eq_eta_p_other_3[
                y in model_data.index_y,
                p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))],
                k in ipp.index_k_existing,
                t in model_data.index_t,
            ],
            ipp.rho_E_my[p, k, t] * (
                ipp.x_E_my(p, k) - sum(
                    ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) - ipp.x_R_cumu(p, k)
            ) - y_E[y, p, k, t] >= 0
        )
    end

    @constraint(
        WMDER_IPP,
        Eq_lambda_p_star_1[
            y in model_data.index_y,
            k in ipp.index_k_new,
            t in model_data.index_t,
        ],
        lambda[y, p_star, k, t] <= ipp.B1GM.value * u_lambda[y, p_star, k, t]
    )
    @constraint(
        WMDER_IPP,
        Eq_lambda_p_star_2[
            y in model_data.index_y,
            k in ipp.index_k_new,
            t in model_data.index_t,
        ],
        ipp.rho_C_my(p_star, k, t) * (
            sum(
                ipp.x_C_my_st2[Symbol(Int(y_symbol)), k] for
                y_symbol in model_data.year(first(model_data.index_y)):model_data.year(y)
            ) + ipp.x_C_cumu(p_star, k)
        ) - y_C[y, p_star, k, t] <= ipp.B1GM.value * (1 - u_lambda[y, p_star, k, t])
    )
    @constraint(
        WMDER_IPP,
        Eq_lambda_p_star_3[
            y in model_data.index_y,
            k in ipp.index_k_new,
            t in model_data.index_t,
        ],
        ipp.rho_C_my(p_star, k, t) * (
            sum(
                ipp.x_C_my_st2[Symbol(Int(y_symbol)), k] for
                y_symbol in model_data.year(first(model_data.index_y)):model_data.year(y)
            ) + ipp.x_C_cumu(p_star, k)
        ) - y_C[y, p_star, k, t] >= 0
    )
    if length(ipp.index_p) >= 2
        @constraint(
            WMDER_IPP,
            Eq_lambda_p_other_1[
                y in model_data.index_y,
                p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))],
                k in ipp.index_k_new,
                t in model_data.index_t,
            ],
            lambda[y, p, k, t] <= ipp.B1GM.value * u_lambda[y, p, k, t]
        )
        @constraint(
            WMDER_IPP,
            Eq_lambda_p_other_2[
                y in model_data.index_y,
                p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))],
                k in ipp.index_k_new,
                t in model_data.index_t,
            ],
            ipp.rho_C_my[p, k, t] * (
                sum(
                    ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) + ipp.x_C_cumu(p, k)
            ) - y_C[y, p, k, t] <= ipp.B1GM.value * (1 - u_lambda[y, p, k, t])
        )
        @constraint(
            WMDER_IPP,
            Eq_lambda_p_other_3[
                y in model_data.index_y,
                p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))],
                k in ipp.index_k_new,
                t in model_data.index_t,
            ],
            ipp.rho_C_my[p, k, t] * (
                sum(
                    ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) + ipp.x_C_cumu(p, k)
            ) - y_C[y, p, k, t] >= 0
        )
    end

    # RPS constraint
    @constraint(
        WMDER_IPP,
        Eq_rps[y in model_data.index_y],
        sum(
            model_data.omega(t) * (y_E[y, p, rps, t] + y_C[y, p, rps, t]) for
            rps in ipp.index_rps, t in model_data.index_t, p in ipp.index_p
        ) -
        ipp.RPS(y) *
        sum(model_data.omega(t) * ipp.Net_Load_my(y, t) for t in model_data.index_t) >= 0
    )

    TimerOutputs.@timeit HEM_TIMER "optimize! Lagrange_Feasible_Cap" begin
        optimize!(WMDER_IPP)
    end

    UCAP_p_star = KeyedArray(
        [
            sum(
                ipp.capacity_credit_E_my(y, k) * (
                    ipp.x_E_my(p_star, k) - sum(
                        ipp.x_R_my_st2[Symbol(Int(y_symbol)), k] for y_symbol in
                        model_data.year(first(model_data.index_y)):model_data.year(y)
                    ) - ipp.x_R_cumu(p_star, k)
                ) for k in ipp.index_k_existing
            ) + sum(
                ipp.capacity_credit_C_my(y, k) * (
                    sum(
                        ipp.x_C_my_st2[Symbol(Int(y_symbol)), k] for y_symbol in
                        model_data.year(first(model_data.index_y)):model_data.year(y)
                    ) + ipp.x_C_cumu(p_star, k)
                ) for k in ipp.index_k_new
            ) for y in model_data.index_y
        ];
        [get_pair(model_data.index_y)]...
    )

    if length(ipp.index_p) >= 2
        UCAP_total = KeyedArray(
            [
                UCAP_p_star[y] +
                sum(
                    ipp.capacity_credit_E_my(y, k) * (
                        ipp.x_E_my(p, k) - sum(
                            ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_R_cumu(p, k)
                    ) for k in ipp.index_k_existing,
                    p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))]
                ) +
                sum(
                    ipp.capacity_credit_C_my(y, k) * (
                        sum(
                            ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p, k)
                    ) for k in ipp.index_k_new,
                    p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))]
                ) for y in model_data.index_y
            ];
            [get_pair(model_data.index_y)]...
        )
    else
        UCAP_total = KeyedArray(
            [UCAP_p_star[y] for y in model_data.index_y];
            [get_pair(model_data.index_y)]...
        )
    end

    update!(
        ipp.obj_feasible,
        objective_value(WMDER_IPP) + sum(
            # capacity revenue
            ipp.pvf_onm(y, p_star) * (
                UCAP_p_star[y] *
                (ipp.Capacity_intercept_my(y) + ipp.Capacity_slope_my(y) * UCAP_total(y))
            ) -
            # fixed costs
            #   fom * (cap exist - cap retiring) for gen type
            ipp.pvf_onm(y, p_star) * sum(
                ipp.fom_E_my(y, p_star, k) * (
                    ipp.x_E_my(p_star, k) - sum(
                        ipp.x_R_my_st2[Symbol(Int(y_symbol)), k] for y_symbol in
                        model_data.year(first(model_data.index_y)):model_data.year(y)
                    )
                ) for k in ipp.index_k_existing
            ) -
            # fixed costs
            #   fom * cap new for gen type
            sum(
                ipp.fom_C_my(y, p_star, k) *
                ipp.x_C_my_st2(y, k) *
                sum(
                    ipp.pvf_onm[Symbol(Int(y_symbol)), p_star] for y_symbol in
                    model_data.year(y):model_data.year(last(model_data.index_y))
                ) for k in ipp.index_k_new
            ) -
            # fixed costs
            #   capex * cap new for gen type
            ipp.pvf_cap(y, p_star) *
            sum(ipp.CapEx_my(y, p_star, k) * ipp.x_C_my_st2(y, k) for k in ipp.index_k_new)

            for y in model_data.index_y
        ),
    )

    for y in model_data.index_y
        ipp.capacity_price_my_temp(y, :) .=
            ipp.Capacity_intercept_my(y) + ipp.Capacity_slope_my(y) * UCAP_total(y)
        ipp.ucap_temp(y, p_star, :) .= UCAP_p_star[y]
    end

    return ipp.obj_feasible
end

function solve_agent_problem_ipp_cap(
    ipp::IPPGroup,
    ipp_opts::IPPOptions{LagrangeDecomposition},
    p_star,
    model_data::HEMData,
    hem_opts::HEMOptions{WholesaleMarket},
    agent_store::AgentStore,
    w_iter,
)
    global optimality_gap = 100.0
    global tol = 0.01  #0.001
    global rho = 1.0
    global alpha = 1.0
    global lagrange_iter = 1
    global max_iter_ipp = 100   #100

    update!(ipp.obj_upper_bound, Inf)
    update!(ipp.obj_lower_bound, -Inf)

    x_R_before = ParamAxisArray(ipp.x_R_my)
    x_C_before = ParamAxisArray(ipp.x_C_my)

    while optimality_gap >= tol && lagrange_iter <= max_iter_ipp
        ipp.obj_st1 = Lagrange_Sub_Investment_Retirement_Cap(
            ipp,
            ipp_opts,
            p_star,
            model_data,
            hem_opts,
            agent_store,
            w_iter,
        )
        ipp.obj_st2 = Lagrange_Sub_Dispatch_Cap(
            ipp,
            ipp_opts,
            p_star,
            model_data,
            hem_opts,
            agent_store,
            w_iter,
        )
        ipp.obj_feasible = Lagrange_Feasible_Cap(
            ipp,
            ipp_opts,
            p_star,
            model_data,
            hem_opts,
            agent_store,
            w_iter,
        )

        #=
        if ipp.obj_st1 + ipp.obj_st2 < ipp.obj_lower_bound
            save_results_debug(ipp, exportfilepath, "Results", temprep, lagrange_iter)
            @info "Glitch: Upper Bound is smaller than Lower Bound"
            break
        end
        =#

        # if the supposedly upper bound is less than lower bound, modify lagrange multiplier to move away from this solution
        while ipp.obj_st1 + ipp.obj_st2 < ipp.obj_feasible
            for y in model_data.index_y, k in ipp.index_k_existing
                ipp.L_R(y, k, :) .= ipp.L_R(y, k) + 1e-9
            end
            for y in model_data.index_y, k in ipp.index_k_new
                ipp.L_C(y, k, :) .= ipp.L_C(y, k) + 1e-9
            end
            ipp.obj_st1 = Lagrange_Sub_Investment_Retirement_Cap(
                ipp,
                ipp_opts,
                p_star,
                model_data,
                hem_opts,
                agent_store,
                w_iter,
            )
            ipp.obj_st2 = Lagrange_Sub_Dispatch_Cap(
                ipp,
                ipp_opts,
                p_star,
                model_data,
                hem_opts,
                agent_store,
                w_iter,
            )
        end

        if ipp.obj_st1 + ipp.obj_st2 < ipp.obj_upper_bound &&
           ipp.obj_st1 + ipp.obj_st2 > ipp.obj_lower_bound
            update!(ipp.obj_upper_bound, ipp.obj_st1 + ipp.obj_st2)
        end

        # feasible = ipp.obj_feasible
        # objective_st1 = ipp.obj_st1
        # objective_st2 = ipp.obj_st2

        if ipp.obj_feasible > ipp.obj_lower_bound
            # @info "Iteration $lagrange_iter feasibility $feasible"
            update!(ipp.obj_lower_bound, ipp.obj_feasible.value)
            for y in model_data.index_y, p in ipp.index_p, k in ipp.index_k_existing
                ipp.x_R_my(y, p, k, :) .= ipp.x_R_my_temp[y, p, k]
            end
            for y in model_data.index_y, p in ipp.index_p, k in ipp.index_k_new
                ipp.x_C_my(y, p, k, :) .= ipp.x_C_my_temp[y, p, k]
            end
            for y in model_data.index_y,
                p in ipp.index_p,
                k in ipp.index_k_existing,
                t in model_data.index_t

                ipp.y_E_my(y, p, k, t, :) .= ipp.y_E_my_temp(y, p, k, t)
            end
            for y in model_data.index_y,
                p in ipp.index_p,
                k in ipp.index_k_new,
                t in model_data.index_t

                ipp.y_C_my(y, p, k, t, :) .= ipp.y_C_my_temp(y, p, k, t)
            end
            for y in model_data.index_y, t in model_data.index_t
                ipp.miu_my(y, t, :) .= ipp.miu_my_temp(y, t)
            end
            for y in model_data.index_y, t in model_data.index_t
                ipp.LMP_my(y, t, :) .= ipp.LMP_my_temp(y, t)
            end
            for y in model_data.index_y
                ipp.capacity_price(y, :) .= ipp.capacity_price_my_temp(y)
                ipp.ucap(y, p_star, :) .= ipp.ucap_temp(y, p_star)
            end
        end

        optimality_gap = (ipp.obj_upper_bound - ipp.obj_lower_bound) / ipp.obj_upper_bound
        UB = ipp.obj_upper_bound.value
        LB = ipp.obj_lower_bound.value

        alpha =
            rho * (ipp.obj_upper_bound - ipp.obj_lower_bound) / (
                sum(
                    (ipp.x_C_my_st1(y, k) - ipp.x_C_my_st2(y, k))^2 for
                    y in model_data.index_y, k in ipp.index_k_new
                ) + sum(
                    (ipp.x_R_my_st1(y, k) - ipp.x_R_my_st2(y, k))^2 for
                    y in model_data.index_y, k in ipp.index_k_existing
                )
            )

        for y in model_data.index_y, k in ipp.index_k_existing
            ipp.L_R(y, k, :) .=
                ipp.L_R(y, k) - alpha * (ipp.x_R_my_st1(y, k) - ipp.x_R_my_st2(y, k))
        end
        for y in model_data.index_y, k in ipp.index_k_new
            ipp.L_C(y, k, :) .=
                ipp.L_C(y, k) - alpha * (ipp.x_C_my_st1(y, k) - ipp.x_C_my_st2(y, k))
        end

        lagrange_iter += 1
        @info "Lagrange iteration $lagrange_iter optimality gap: $optimality_gap"
        @info "Iteration $lagrange_iter upper bound: $UB"
        @info "Iteration $lagrange_iter lower bound: $LB"

        #=
        @info "Iteration $lagrange_iter stage 1 objective: $objective_st1"
        @info "Iteration $lagrange_iter stage 2 objective: $objective_st2"
        =#

    end

    return compute_difference_percentage_one_norm([(x_R_before, ipp.x_R_my), (x_C_before, ipp.x_C_my)])
end

function solve_agent_problem_ipp_cap(
    ipp::IPPGroup,
    ipp_opts::IPPOptions{MIQP},
    p_star,
    model_data::HEMData,
    hem_opts::HEMOptions{WholesaleMarket},
    agent_store::AgentStore,
    w_iter,
)
    x_R_before = ParamAxisArray(ipp.x_R_my)
    x_C_before = ParamAxisArray(ipp.x_C_my)

    utility = get_agent(Utility, agent_store)
    regulator = get_agent(Regulator, agent_store)
    customers = get_agent(CustomerGroup, agent_store)
    green_developer = get_agent(GreenDeveloper, agent_store)

    WMDER_IPP = get_new_jump_model(hem_opts.MIP_solver)
    # set_optimizer_attribute(WMDER_IPP, "OUTPUTLOG", 0)

    # get_optimizer_attribute(WMDER_IPP, "MIPTOL")
    # set_optimizer_attributes(WMDER_IPP, "tol" => 1e-6, "max_iter" => 500)

    # Define positive variables
    @variable(WMDER_IPP, x_C[model_data.index_y, ipp.index_k_new] >= 0)
    @variable(WMDER_IPP, x_R[model_data.index_y, ipp.index_k_existing] >= 0)
    @variable(
        WMDER_IPP,
        y_E[model_data.index_y, ipp.index_p, ipp.index_k_existing, model_data.index_t] >= 0
    )
    @variable(
        WMDER_IPP,
        y_C[model_data.index_y, ipp.index_p, ipp.index_k_new, model_data.index_t] >= 0
    )
    @variable(WMDER_IPP, miu[model_data.index_y, model_data.index_t] >= 0)
    @variable(
        WMDER_IPP,
        eta[model_data.index_y, ipp.index_p, ipp.index_k_existing, model_data.index_t] >= 0
    )
    @variable(
        WMDER_IPP,
        lambda[model_data.index_y, ipp.index_p, ipp.index_k_new, model_data.index_t] >= 0
    )

    @variable(
        WMDER_IPP,
        u_y_E[model_data.index_y, ipp.index_p, ipp.index_k_existing, model_data.index_t],
        Bin
    )
    @variable(
        WMDER_IPP,
        u_y_C[model_data.index_y, ipp.index_p, ipp.index_k_new, model_data.index_t],
        Bin
    )
    @variable(WMDER_IPP, u_miu[model_data.index_y, model_data.index_t], Bin)
    @variable(
        WMDER_IPP,
        u_eta[model_data.index_y, ipp.index_p, ipp.index_k_existing, model_data.index_t],
        Bin
    )
    @variable(
        WMDER_IPP,
        u_lambda[model_data.index_y, ipp.index_p, ipp.index_k_new, model_data.index_t],
        Bin
    )

    # setvalue(x_C[Symbol("2019"), :GasCC], 0.0)
    # setvalue(x_C[Symbol("2019"), :GasCT], 0.0)
    # setvalue(x_C[Symbol("2019"), :UPV], 762.961280674661)
    # setvalue(x_R[Symbol("2019"), :GasCC], 437.665817196685)
    # setvalue(x_R[Symbol("2019"), :GasCT], 0.0)
    # setvalue(x_R[Symbol("2019"), :UPV], 0.0)

    for p in ipp.index_p
        for y in model_data.index_y
            if y == last(model_data.index_y.elements)
                ipp.pvf_onm(y, p, :) .= ipp.pvf_cap(y, p) / ipp.CRF_default(p)
            else
                ipp.pvf_onm(y, p, :) .= ipp.pvf_cap(y, p)
            end
        end
    end

    # adding capacity market parameters
    fill!(ipp.Net_Load_my, NaN)
    for y in model_data.index_y, t in model_data.index_t
        ipp.Net_Load_my(y, t, :) .=
            sum(customers.gamma(h) * customers.d_my(y, h, t) for h in model_data.index_h) +
            ipp.eximport_my(y, t) - sum(
                customers.rho_DG(h, m, t) * customers.x_DG_E_my(y, h, m) for
                h in model_data.index_h, m in customers.index_m
            ) - sum(
                customers.rho_DG(h, m, t) * sum(
                    customers.x_DG_new_my(Symbol(Int(y_symbol)), h, m) for y_symbol in
                    model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                ) for h in model_data.index_h, m in customers.index_m
            )
    end
    fill!(ipp.Max_Net_Load_my, NaN)
    for y in model_data.index_y
        ipp.Max_Net_Load_my(y, :) .=
            findmax(Dict(t => ipp.Net_Load_my(y, t) for t in model_data.index_t))[1]
    end

    Max_Net_Load_my_index = Dict(
        y => findmax(Dict(t => ipp.Net_Load_my(y, t) for t in model_data.index_t))[2]
        for y in model_data.index_y
    )
    fill!(ipp.capacity_credit_E_my, NaN)
    for y in model_data.index_y, k in ipp.index_k_existing
        ipp.capacity_credit_E_my(y, k, :) .= ipp.rho_E_my[p_star, k, Max_Net_Load_my_index(y)]
    end
    fill!(ipp.capacity_credit_C_my, NaN)
    for y in model_data.index_y, k in ipp.index_k_new
        ipp.capacity_credit_C_my(y, k, :) .= ipp.rho_C_my[p_star, k, Max_Net_Load_my_index(y)]
    end
    fill!(ipp.Reserve_req_my, NaN)
    for y in model_data.index_y
        ipp.Reserve_req_my(y, :) .= (1 + regulator.r) * ipp.Max_Net_Load_my(y)
    end

    for y in model_data.index_y
        ipp.Capacity_slope_my(y, :) .=
            -ipp.NetCONE(y) / (ipp.DC_length(y) * ipp.Reserve_req_my(y))
        ipp.Capacity_intercept_my(y, :) .=
            -ipp.Capacity_slope_my(y) * ipp.Reserve_req_my(y) + ipp.NetCONE(y)
    end
    #=
    ipp.Capacity_slope_my = Dict(y => - ipp.NetCONE(y) / (ipp.DC_length(y)*ipp.Reserve_req_my(y))
        for y in model_data.index_y)
    ipp.Capacity_intercept_my = Dict(y => - ipp.Capacity_slope_my(y) * ipp.Reserve_req_my(y) + ipp.NetCONE(y)
        for y in model_data.index_y)
    =#
    #####################################################################

    UCAP_p_star =
        y -> begin
            sum(
                ipp.capacity_credit_E_my(y, k) * (
                    ipp.x_E_my(p_star, k) - sum(
                        x_R[Symbol(Int(y_symbol)), k] for y_symbol in
                        model_data.year(first(model_data.index_y)):model_data.year(y)
                    ) - ipp.x_R_cumu(p_star, k)
                ) for k in ipp.index_k_existing
            ) + sum(
                ipp.capacity_credit_C_my(y, k) * (
                    sum(
                        x_C[Symbol(Int(y_symbol)), k] for y_symbol in
                        model_data.year(first(model_data.index_y)):model_data.year(y)
                    ) + ipp.x_C_cumu(p_star, k)
                ) for k in ipp.index_k_new
            )
        end

    if length(ipp.index_p) >= 2
        UCAP_total =
            y -> begin
                UCAP_p_star(y) +
                sum(
                    ipp.capacity_credit_E_my(y, k) * (
                        ipp.x_E_my(p, k) - sum(
                            ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_R_cumu(p, k)
                    ) for k in ipp.index_k_existing,
                    p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))]
                ) +
                sum(
                    ipp.capacity_credit_C_my(y, k) * (
                        sum(
                            ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p, k)
                    ) for k in ipp.index_k_new,
                    p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))]
                ) +
                # green technology subscription
                sum(
                    ipp.capacity_credit_C_my[y, j] * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
                    model_data.year(first(model_data.index_y_fix)):model_data.year(y))
                    for j in model_data.index_j, h in model_data.index_h
                )
            end
    else
        UCAP_total = y -> begin
            UCAP_p_star(y) +
            # green technology subscription
            sum(
                ipp.capacity_credit_C_my[y, j] * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):model_data.year(y))
                for j in model_data.index_j, h in model_data.index_h
            )
        end
    end

    objective_function = begin
        sum(

            # capacity revenue
            ipp.pvf_onm(y, p_star) * (
                UCAP_p_star(y) * (
                    ipp.Capacity_intercept_my(y) + ipp.Capacity_slope_my(y) * UCAP_total(y)
                )
            ) -

            # fixed costs
            #   fom * (cap exist - cap retiring) for gen type
            ipp.pvf_onm(y, p_star) * sum(
                ipp.fom_E_my(y, p_star, k) * (
                    ipp.x_E_my(p_star, k) - sum(
                        x_R[Symbol(Int(y_symbol)), k] for y_symbol in
                        model_data.year(first(model_data.index_y)):model_data.year(y)
                    )
                ) for k in ipp.index_k_existing
            ) -
            # fixed costs
            #   fom * cap new for gen type
            sum(
                ipp.fom_C_my(y, p_star, k) *
                x_C[y, k] *
                sum(
                    ipp.pvf_onm[Symbol(Int(y_symbol)), p_star] for y_symbol in
                    model_data.year(y):model_data.year(last(model_data.index_y))
                ) for k in ipp.index_k_new
            ) -
            # fixed costs
            #   capex * cap new for gen type
            ipp.pvf_cap(y, p_star) *
            sum(ipp.CapEx_my(y, p_star, k) * x_C[y, k] for k in ipp.index_k_new) +

            # Linearized profit term
            ipp.pvf_onm(y, p_star) * (
                sum(
                    miu[y, t] * (
                        sum(
                            customers.gamma(h) * customers.d_my(y, h, t) for
                            h in model_data.index_h
                        ) + ipp.eximport_my(y, t) - sum(
                            customers.rho_DG(h, m, t) * customers.x_DG_E_my(y, h, m) for
                            h in model_data.index_h, m in customers.index_m
                        ) - sum(
                            customers.rho_DG(h, m, t) * sum(
                                customers.x_DG_new_my(Symbol(Int(y_symbol)), h, m) for
                                y_symbol in
                                model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                            ) for h in model_data.index_h, m in customers.index_m
                        ) -
                        # green technology subscription at time t
                        sum(
                            utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
                            model_data.year(first(model_data.index_y_fix)):model_data.year(y))
                            for j in model_data.index_j, h in model_data.index_h
                        )
                    ) for t in model_data.index_t
                ) - (
                    sum(
                        eta[y, p, k, t] *
                        ipp.rho_E_my[p, k, t] *
                        (
                            ipp.x_E_my(p, k) - sum(
                                ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) - ipp.x_R_cumu(p, k)
                        ) for t in model_data.index_t, k in ipp.index_k_existing,
                        p in ipp.index_p
                    ) - sum(
                        eta[y, p_star, k, t] *
                        ipp.rho_E_my(p_star, k, t) *
                        (
                            ipp.x_E_my(p_star, k) - sum(
                                ipp.x_R_my[Symbol(Int(y_symbol)), p_star, k] for
                                y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) - ipp.x_R_cumu(p_star, k)
                        ) for t in model_data.index_t, k in ipp.index_k_existing
                    )
                ) - (
                    sum(
                        lambda[y, p, k, t] *
                        ipp.rho_C_my[p, k, t] *
                        (
                            sum(
                                ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) + ipp.x_C_cumu(p, k)
                        ) for t in model_data.index_t, k in ipp.index_k_new,
                        p in ipp.index_p
                    ) - sum(
                        lambda[y, p_star, k, t] *
                        ipp.rho_C_my(p_star, k, t) *
                        (
                            sum(
                                ipp.x_C_my[Symbol(Int(y_symbol)), p_star, k] for
                                y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) + ipp.x_C_cumu(p_star, k)
                        ) for t in model_data.index_t, k in ipp.index_k_new
                    )
                ) - (
                    sum(
                        model_data.omega(t) * (ipp.v_E_my(y, p, k, t) * y_E[y, p, k, t]) for
                        t in model_data.index_t, k in ipp.index_k_existing,
                        p in ipp.index_p
                    ) - sum(
                        model_data.omega(t) *
                        (ipp.v_E_my[y, p_star, k, t] * y_E[y, p_star, k, t]) for
                        t in model_data.index_t, k in ipp.index_k_existing
                    )
                ) - (
                    sum(
                        model_data.omega(t) * (ipp.v_C_my(y, p, k, t) * y_C[y, p, k, t]) for
                        t in model_data.index_t, k in ipp.index_k_new, p in ipp.index_p
                    ) - sum(
                        model_data.omega(t) *
                        (ipp.v_C_my[y, p_star, k, t] * y_C[y, p_star, k, t]) for
                        t in model_data.index_t, k in ipp.index_k_new
                    )
                ) -

                # generation costs (this includes terms due to revenue linearization)
                #   num hrs * ((fuel + vom) * gen existing + (fuel + vom) * gen new) for t and gen type
                sum(
                    model_data.omega(t) *
                    (ipp.v_E_my[y, p_star, k, t] * y_E[y, p_star, k, t]) for
                    t in model_data.index_t, k in ipp.index_k_existing
                ) - sum(
                    model_data.omega(t) *
                    (ipp.v_C_my[y, p_star, k, t] * y_C[y, p_star, k, t]) for
                    t in model_data.index_t, k in ipp.index_k_new
                )
            )

            for y in model_data.index_y
        )
    end

    @objective(WMDER_IPP, Max, objective_function)

    @constraint(
        WMDER_IPP,
        Eq_sigma[y in model_data.index_y, k in ipp.index_k_existing],
        ipp.x_E_my(p_star, k) - sum(
            x_R[Symbol(Int(y_symbol)), k] for
            y_symbol in model_data.year(first(model_data.index_y)):model_data.year(y)
        ) - ipp.x_R_cumu(p_star, k) >= 0
    )

    @constraint(
        WMDER_IPP,
        Eq_y_E_p_1[
            y in model_data.index_y,
            p in ipp.index_p,
            k in ipp.index_k_existing,
            t in model_data.index_t,
        ],
        y_E[y, p, k, t] <= ipp.B1GM.value * u_y_E[y, p, k, t]
    )
    @constraint(
        WMDER_IPP,
        Eq_y_E_p_2[
            y in model_data.index_y,
            p in ipp.index_p,
            k in ipp.index_k_existing,
            t in model_data.index_t,
        ],
        model_data.omega(t) * ipp.v_E_my(y, p, k, t) - miu[y, t] + eta[y, p, k, t] <=
        ipp.B1GM.value * (1 - u_y_E[y, p, k, t])
    )
    @constraint(
        WMDER_IPP,
        Eq_y_E_p_3[
            y in model_data.index_y,
            p in ipp.index_p,
            k in ipp.index_k_existing,
            t in model_data.index_t,
        ],
        model_data.omega(t) * ipp.v_E_my(y, p, k, t) - miu[y, t] + eta[y, p, k, t] >= 0
    )

    @constraint(
        WMDER_IPP,
        Eq_y_C_p_1[
            y in model_data.index_y,
            p in ipp.index_p,
            k in ipp.index_k_new,
            t in model_data.index_t,
        ],
        y_C[y, p, k, t] <= ipp.B1GM.value * u_y_C[y, p, k, t]
    )
    @constraint(
        WMDER_IPP,
        Eq_y_C_p_2[
            y in model_data.index_y,
            p in ipp.index_p,
            k in ipp.index_k_new,
            t in model_data.index_t,
        ],
        model_data.omega(t) * ipp.v_C_my(y, p, k, t) - miu[y, t] + lambda[y, p, k, t] <=
        ipp.B1GM.value * (1 - u_y_C[y, p, k, t])
    )
    @constraint(
        WMDER_IPP,
        Eq_y_C_p_3[
            y in model_data.index_y,
            p in ipp.index_p,
            k in ipp.index_k_new,
            t in model_data.index_t,
        ],
        model_data.omega(t) * ipp.v_C_my(y, p, k, t) - miu[y, t] + lambda[y, p, k, t] >= 0
    )

    supply_demand_balance =
        (y, t) -> begin
            # bulk generation at time t
            sum(y_E[y, p, k, t] for k in ipp.index_k_existing, p in ipp.index_p) +
            sum(y_C[y, p, k, t] for k in ipp.index_k_new, p in ipp.index_p) -
            # demand at time t
            sum(customers.gamma(h) * customers.d_my(y, h, t) for h in model_data.index_h) - ipp.eximport_my(y, t) +
            # existing DG generation at time t
            sum(
                customers.rho_DG(h, m, t) * customers.x_DG_E_my(y, h, m) for
                h in model_data.index_h, m in customers.index_m
            ) +
            # new DG generation at time t
            sum(
                customers.rho_DG(h, m, t) * sum(
                    customers.x_DG_new_my(Symbol(Int(y_symbol)), h, m) for y_symbol in
                    model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                ) for h in model_data.index_h, m in customers.index_m
            ) +
            # green technology subscription at time t
            sum(
                utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):model_data.year(y))
                for j in model_data.index_j, h in model_data.index_h
            )
        end

    @constraint(
        WMDER_IPP,
        Eq_miu_1[y in model_data.index_y, t in model_data.index_t],
        miu[y, t] <= ipp.B1GM.value * u_miu[y, t]
    )
    @constraint(
        WMDER_IPP,
        Eq_miu_2[y in model_data.index_y, t in model_data.index_t],
        supply_demand_balance(y, t) <= ipp.B1GM.value * (1 - u_miu[y, t])
    )
    @constraint(
        WMDER_IPP,
        Eq_miu_3[y in model_data.index_y, t in model_data.index_t],
        supply_demand_balance(y, t) >= 0
    )

    @constraint(
        WMDER_IPP,
        Eq_eta_p_star_1[
            y in model_data.index_y,
            k in ipp.index_k_existing,
            t in model_data.index_t,
        ],
        eta[y, p_star, k, t] <= ipp.B1GM.value * u_eta[y, p_star, k, t]
    )
    @constraint(
        WMDER_IPP,
        Eq_eta_p_star_2[
            y in model_data.index_y,
            k in ipp.index_k_existing,
            t in model_data.index_t,
        ],
        ipp.rho_E_my(p_star, k, t) * (
            ipp.x_E_my(p_star, k) - sum(
                x_R[Symbol(Int(y_symbol)), k] for
                y_symbol in model_data.year(first(model_data.index_y)):model_data.year(y)
            ) - ipp.x_R_cumu(p_star, k)
        ) - y_E[y, p_star, k, t] <= ipp.B1GM.value * (1 - u_eta[y, p_star, k, t])
    )
    @constraint(
        WMDER_IPP,
        Eq_eta_p_star_3[
            y in model_data.index_y,
            k in ipp.index_k_existing,
            t in model_data.index_t,
        ],
        ipp.rho_E_my(p_star, k, t) * (
            ipp.x_E_my(p_star, k) - sum(
                x_R[Symbol(Int(y_symbol)), k] for
                y_symbol in model_data.year(first(model_data.index_y)):model_data.year(y)
            ) - ipp.x_R_cumu(p_star, k)
        ) - y_E[y, p_star, k, t] >= 0
    )
    if length(ipp.index_p) >= 2
        @constraint(
            WMDER_IPP,
            Eq_eta_p_other_1[
                y in model_data.index_y,
                p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))],
                k in ipp.index_k_existing,
                t in model_data.index_t,
            ],
            eta[y, p, k, t] <= ipp.B1GM.value * u_eta[y, p, k, t]
        )
        @constraint(
            WMDER_IPP,
            Eq_eta_p_other_2[
                y in model_data.index_y,
                p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))],
                k in ipp.index_k_existing,
                t in model_data.index_t,
            ],
            ipp.rho_E_my[p, k, t] * (
                ipp.x_E_my(p, k) - sum(
                    ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) - ipp.x_R_cumu(p, k)
            ) - y_E[y, p, k, t] <= ipp.B1GM.value * (1 - u_eta[y, p, k, t])
        )
        @constraint(
            WMDER_IPP,
            Eq_eta_p_other_3[
                y in model_data.index_y,
                p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))],
                k in ipp.index_k_existing,
                t in model_data.index_t,
            ],
            ipp.rho_E_my[p, k, t] * (
                ipp.x_E_my(p, k) - sum(
                    ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) - ipp.x_R_cumu(p, k)
            ) - y_E[y, p, k, t] >= 0
        )
    end

    @constraint(
        WMDER_IPP,
        Eq_lambda_p_star_1[
            y in model_data.index_y,
            k in ipp.index_k_new,
            t in model_data.index_t,
        ],
        lambda[y, p_star, k, t] <= ipp.B1GM.value * u_lambda[y, p_star, k, t]
    )
    @constraint(
        WMDER_IPP,
        Eq_lambda_p_star_2[
            y in model_data.index_y,
            k in ipp.index_k_new,
            t in model_data.index_t,
        ],
        ipp.rho_C_my(p_star, k, t) * (
            sum(
                x_C[Symbol(Int(y_symbol)), k] for
                y_symbol in model_data.year(first(model_data.index_y)):model_data.year(y)
            ) + ipp.x_C_cumu(p_star, k)
        ) - y_C[y, p_star, k, t] <= ipp.B1GM.value * (1 - u_lambda[y, p_star, k, t])
    )
    @constraint(
        WMDER_IPP,
        Eq_lambda_p_star_3[
            y in model_data.index_y,
            k in ipp.index_k_new,
            t in model_data.index_t,
        ],
        ipp.rho_C_my(p_star, k, t) * (
            sum(
                x_C[Symbol(Int(y_symbol)), k] for
                y_symbol in model_data.year(first(model_data.index_y)):model_data.year(y)
            ) + ipp.x_C_cumu(p_star, k)
        ) - y_C[y, p_star, k, t] >= 0
    )
    if length(ipp.index_p) >= 2
        @constraint(
            WMDER_IPP,
            Eq_lambda_p_other_1[
                y in model_data.index_y,
                p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))],
                k in ipp.index_k_new,
                t in model_data.index_t,
            ],
            lambda[y, p, k, t] <= ipp.B1GM.value * u_lambda[y, p, k, t]
        )
        @constraint(
            WMDER_IPP,
            Eq_lambda_p_other_2[
                y in model_data.index_y,
                p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))],
                k in ipp.index_k_new,
                t in model_data.index_t,
            ],
            ipp.rho_C_my[p, k, t] * (
                sum(
                    ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) + ipp.x_C_cumu(p, k)
            ) - y_C[y, p, k, t] <= ipp.B1GM.value * (1 - u_lambda[y, p, k, t])
        )
        @constraint(
            WMDER_IPP,
            Eq_lambda_p_other_3[
                y in model_data.index_y,
                p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))],
                k in ipp.index_k_new,
                t in model_data.index_t,
            ],
            ipp.rho_C_my[p, k, t] * (
                sum(
                    ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) + ipp.x_C_cumu(p, k)
            ) - y_C[y, p, k, t] >= 0
        )
    end

    if length(ipp.index_p) >= 2
        planning_reserves =
            (y, t) -> begin
                # bulk generation available capacity at time t
                sum(
                    ipp.rho_E_my(p_star, k, t) * (
                        ipp.x_E_my(p_star, k) - sum(
                            x_R[Symbol(Int(y_symbol)), k] for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_R_cumu(p_star, k)
                    ) for k in ipp.index_k_existing
                ) +
                sum(
                    ipp.rho_C_my(p_star, k, t) * (
                        sum(
                            x_C[Symbol(Int(y_symbol)), k] for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p_star, k)
                    ) for k in ipp.index_k_new
                ) +
                sum(
                    ipp.rho_E_my[p, k, t] * (
                        ipp.x_E_my(p, k) - sum(
                            ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_R_cumu(p, k)
                    ) for k in ipp.index_k_existing,
                    p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))]
                ) +
                sum(
                    ipp.rho_C_my[p, k, t] * (
                        sum(
                            ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p, k)
                    ) for k in ipp.index_k_new,
                    p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))]
                ) +
                # green technology subscription
                sum(
                    utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
                    model_data.year(first(model_data.index_y_fix)):model_data.year(y))
                    for j in model_data.index_j, h in model_data.index_h
                ) -
                # net_load plus planning reserve
                (1 + regulator.r) * (
                    sum(
                        customers.gamma(h) * customers.d_my(y, h, t) for
                        h in model_data.index_h
                    ) + ipp.eximport_my(y, t) - sum(
                        customers.rho_DG(h, m, t) * customers.x_DG_E_my(y, h, m) for
                        h in model_data.index_h, m in customers.index_m
                    ) - sum(
                        customers.rho_DG(h, m, t) * sum(
                            customers.x_DG_new_my(Symbol(Int(y_symbol)), h, m) for
                            y_symbol in
                            model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                        ) for h in model_data.index_h, m in customers.index_m
                    )
                )
            end
    else
        planning_reserves =
            (y, t) -> begin
                # bulk generation available capacity at time t
                sum(
                    ipp.rho_E_my(p_star, k, t) * (
                        ipp.x_E_my(p_star, k) - sum(
                            x_R[Symbol(Int(y_symbol)), k] for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_R_cumu(p_star, k)
                    ) for k in ipp.index_k_existing
                ) + sum(
                    ipp.rho_C_my(p_star, k, t) * (
                        sum(
                            x_C[Symbol(Int(y_symbol)), k] for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p_star, k)
                    ) for k in ipp.index_k_new
                ) +
                # green technology subscription
                sum(
                    utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
                    model_data.year(first(model_data.index_y_fix)):model_data.year(y))
                    for j in model_data.index_j, h in model_data.index_h
                ) -
                # net_load plus planning reserve
                (1 + regulator.r) * (
                    sum(
                        customers.gamma(h) * customers.d_my(y, h, t) for
                        h in model_data.index_h
                    ) + ipp.eximport_my(y, t) - sum(
                        customers.rho_DG(h, m, t) * customers.x_DG_E_my(y, h, m) for
                        h in model_data.index_h, m in customers.index_m
                    ) - sum(
                        customers.rho_DG(h, m, t) * sum(
                            customers.x_DG_new_my(Symbol(Int(y_symbol)), h, m) for
                            y_symbol in
                            model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                        ) for h in model_data.index_h, m in customers.index_m
                    )
                )
            end
    end
    @constraint(
        WMDER_IPP,
        Eq_xi[y in model_data.index_y, t in model_data.index_t],
        planning_reserves(y, t) >= 0
    )

    if length(ipp.index_p) >= 2
        planning_reserves_cap =
            y -> begin
                # bulk generation available capacity at time t
                sum(
                    ipp.capacity_credit_E_my(y, k) * (
                        ipp.x_E_my(p_star, k) - sum(
                            x_R[Symbol(Int(y_symbol)), k] for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_R_cumu(p_star, k)
                    ) for k in ipp.index_k_existing
                ) +
                sum(
                    ipp.capacity_credit_C_my(y, k) * (
                        sum(
                            x_C[Symbol(Int(y_symbol)), k] for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p_star, k)
                    ) for k in ipp.index_k_new
                ) +
                sum(
                    ipp.capacity_credit_E_my(y, k) * (
                        ipp.x_E_my(p, k) - sum(
                            ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_R_cumu(p, k)
                    ) for k in ipp.index_k_existing,
                    p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))]
                ) +
                sum(
                    ipp.capacity_credit_C_my(y, k) * (
                        sum(
                            ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p, k)
                    ) for k in ipp.index_k_new,
                    p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))]
                ) +
                # green technology subscription
                sum(
                    ipp.capacity_credit_C_my[y, j] * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
                    model_data.year(first(model_data.index_y_fix)):model_data.year(y))
                    for j in model_data.index_j, h in model_data.index_h
                ) -
                # net_load plus planning reserve
                ipp.Reserve_req_my(y)
            end
    else
        planning_reserves_cap =
            y -> begin
                # bulk generation available capacity at time t
                sum(
                    ipp.capacity_credit_E_my(y, k) * (
                        ipp.x_E_my(p_star, k) - sum(
                            x_R[Symbol(Int(y_symbol)), k] for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_R_cumu(p_star, k)
                    ) for k in ipp.index_k_existing
                ) + sum(
                    ipp.capacity_credit_C_my(y, k) * (
                        sum(
                            x_C[Symbol(Int(y_symbol)), k] for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p_star, k)
                    ) for k in ipp.index_k_new
                ) +
                # green technology subscription
                sum(
                    ipp.capacity_credit_C_my[y, j] * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
                    model_data.year(first(model_data.index_y_fix)):model_data.year(y))
                    for j in model_data.index_j, h in model_data.index_h
                ) -
                # net_load plus planning reserve
                ipp.Reserve_req_my(y)
            end
    end
    @constraint(
        WMDER_IPP,
        Eq_xi_cap[y in model_data.index_y],
        planning_reserves_cap(y) >= 0
    )

    # RPS constraint
    @constraint(
        WMDER_IPP,
        Eq_rps[y in model_data.index_y],
        sum(
            model_data.omega(t) * (y_E[y, p, rps, t] + y_C[y, p, rps, t]) for
            rps in ipp.index_rps, t in model_data.index_t, p in ipp.index_p
        ) -
        ipp.RPS(y) *
        sum(model_data.omega(t) * ipp.Net_Load_my(y, t) for t in model_data.index_t) >= 0
    )

    TimerOutputs.@timeit HEM_TIMER "optimize! Lagrange_Sub_Dispatch_Cap" begin
        optimize!(WMDER_IPP)
    end

    for y in model_data.index_y, k in ipp.index_k_existing
        ipp.x_R_my(y, p_star, k, :) .= value.(x_R[y, k])
    end
    for y in model_data.index_y, p in ipp.index_p, k in ipp.index_k_new
        ipp.x_C_my(y, p_star, k, :) .= value.(x_C[y, k])
    end
    for y in model_data.index_y,
        p in ipp.index_p,
        k in ipp.index_k_existing,
        t in model_data.index_t

        ipp.y_E_my(y, p, k, t, :) .= value.(y_E[y, p, k, t])
    end
    for y in model_data.index_y,
        p in ipp.index_p,
        k in ipp.index_k_new,
        t in model_data.index_t

        ipp.y_C_my(y, p, k, t, :) .= value.(y_C[y, p, k, t])
    end
    for y in model_data.index_y, t in model_data.index_t
        ipp.miu_my(y, t, :) .= value.(miu[y, t])
    end
    for y in model_data.index_y, t in model_data.index_t
        ipp.LMP_my(y, t, :) .= ipp.miu_my(y, t) / model_data.omega(t)
    end

    ###### running MILP only ######
    # UCAP_p_star_solution = Dict(y =>
    #     sum(ipp.capacity_credit_E_my[y,k]*(ipp.x_E_my[p_star,k]-sum(ipp.x_R_my[Symbol(y_symbol),p_star,k] for y_symbol = model_data.year(first(model_data.index_y)):model_data.year(y)) - ipp.x_R_cumu[p_star,k]) for k in ipp.index_k_existing) + 
    #     sum(ipp.capacity_credit_C_my[y,k]*(sum(ipp.x_C_my[Symbol(y_symbol),p_star,k] for y_symbol = model_data.year(first(model_data.index_y)):model_data.year(y)) + ipp.x_C_cumu[p_star,k]) for k in ipp.index_k_new)
    #     for y in model_data.index_y
    # )
    # capacity_profit = Dict(y =>
    #     ipp.pvf_onm[y,p_star] * UCAP_p_star_solution(y) * (ipp.Capacity_intercept_my(y) + ipp.Capacity_slope_my(y) * UCAP_p_star_solution(y))
    #     for y in model_data.index_y
    # )
    # total_profit = Dict(y => capacity_profit(y) + objective_value(WMDER_IPP) for y in model_data.index_y)

    ###### running MIQP ######
    UCAP_p_star = Dict(
        y =>
            sum(
                ipp.capacity_credit_E_my(y, k) * (
                    ipp.x_E_my(p_star, k) - sum(
                        ipp.x_R_my[Symbol(Int(y_symbol)), p_star, k] for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                    ) - ipp.x_R_cumu(p_star, k)
                ) for k in ipp.index_k_existing
            ) + sum(
                ipp.capacity_credit_C_my(y, k) * (
                    sum(
                        ipp.x_C_my[Symbol(Int(y_symbol)), p_star, k] for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                    ) + ipp.x_C_cumu(p_star, k)
                ) for k in ipp.index_k_new
            ) for y in model_data.index_y
    )

    if length(ipp.index_p) >= 2
        UCAP_total = Dict(
            y =>
                UCAP_p_star[y] +
                sum(
                    ipp.capacity_credit_E_my(y, k) * (
                        ipp.x_E_my(p, k) - sum(
                            ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_R_cumu(p, k)
                    ) for k in ipp.index_k_existing,
                    p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))]
                ) +
                sum(
                    ipp.capacity_credit_C_my(y, k) * (
                        sum(
                            ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p, k)
                    ) for k in ipp.index_k_new,
                    p in ipp.index_p[Not(findall(x -> x == p_star, ipp.index_p))]
                ) +
                # green technology subscription
                sum(
                    ipp.capacity_credit_C_my[y, j] * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
                    model_data.year(first(model_data.index_y_fix)):model_data.year(y))
                    for j in model_data.index_j, h in model_data.index_h
                )
                for y in model_data.index_y
        )
    else
        UCAP_total = Dict(y => UCAP_p_star[y] +
        # green technology subscription
        sum(
            ipp.capacity_credit_C_my[y, j] * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
            model_data.year(first(model_data.index_y_fix)):model_data.year(y))
            for j in model_data.index_j, h in model_data.index_h
        ) for y in model_data.index_y)
    end

    for y in model_data.index_y
        ipp.capacity_price(y, :) .=
            ipp.Capacity_intercept_my(y) + ipp.Capacity_slope_my(y) * UCAP_total[y]
        ipp.ucap(y, p_star, :) .= UCAP_p_star[y]
    end

    return compute_difference_percentage_one_norm([
        (x_R_before, ipp.x_R_my),
        (x_C_before, ipp.x_C_my),
    ])
end

function solve_agent_problem!(
    ipp::IPPGroup,
    ipp_opts::AgentOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{WholesaleMarket},
    agent_store::AgentStore,
    w_iter,
)
    diff = 0.0

    for p in ipp.index_p
        TimerOutputs.@timeit HEM_TIMER "solve_agent_problem_ipp_cap" begin
            diff += solve_agent_problem_ipp_cap(
                ipp,
                ipp_opts,
                p,
                model_data,
                hem_opts,
                agent_store,
                w_iter,
            )
        end
    end

    # report change in key variables from previous iteration to this one
    return diff
end

function save_results(
    ipps::IPPGroup,
    ipp_opts::AgentOptions,
    hem_opts::HEMOptions{WholesaleMarket},
    export_file_path::AbstractString,
    fileprefix::AbstractString,
)
    # Primal Variables
    save_param(
        ipps.y_E_my.values,
        [:Year, :IPP, :GenTech, :Time],
        :Generation_MWh,
        joinpath(export_file_path, "$(fileprefix)_y_E.csv"),
    )
    save_param(
        ipps.y_C_my.values,
        [:Year, :IPP, :GenTech, :Time],
        :Generation_MWh,
        joinpath(export_file_path, "$(fileprefix)_y_C.csv"),
    )
    save_param(
        ipps.x_R_my.values,
        [:Year, :IPP, :GenTech],
        :Capacity_MW,
        joinpath(export_file_path, "$(fileprefix)_x_R.csv"),
    )
    save_param(
        ipps.x_C_my.values,
        [:Year, :IPP, :GenTech],
        :Capacity_MW,
        joinpath(export_file_path, "$(fileprefix)_x_C.csv"),
    )
    save_param(
        ipps.LMP_my.values,
        [:Year, :Time],
        :MarginalCost,
        joinpath(export_file_path, "$(fileprefix)_LMP.csv"),
    )
end

function welfare_calculation!(
    ipp::IPPGroup,
    ipp_opts::AgentOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{WholesaleMarket},
    agent_store::AgentStore,
)
    regulator = get_agent(Regulator, agent_store)
    utility = get_agent(Utility, agent_store)

    IPP_Revenue_p = make_keyed_array(model_data.index_y_fix, ipp.index_p)
    for y in model_data.index_y_fix, p in ipp.index_p
        IPP_Revenue_p(y, p, :) .=
        # Linearized revenue term 
            sum(
                ipp.miu_my(y, t) * (
                    sum(ipp.y_E_my(y, p, k, t) for k in ipp.index_k_existing) +
                    sum(ipp.y_C_my(y, p, k, t) for k in ipp.index_k_new)
                ) for t in model_data.index_t
            ) +
            ipp.ucap(y, p) * (
                ipp.Capacity_intercept_my(y) +
                ipp.Capacity_slope_my(y) * sum(ipp.ucap(y, p) for p in ipp.index_p)
            ) +
            # REC revenue
            sum(
                regulator.REC *
                model_data.omega(t) *
                (
                    sum(ipp.y_E_my[y, p, rps, t] for rps in ipp.index_rps) +
                    sum(ipp.y_C_my[y, p, rps, t] for rps in ipp.index_rps)
                ) for t in model_data.index_t
            )
    end

    IPP_Revenue_total = KeyedArray(
        [
            sum(IPP_Revenue_p(y, p) for p in ipp.index_p) + regulator.othercost(y) for
            y in model_data.index_y_fix
        ];
        [get_pair(model_data.index_y_fix)]...
    )

    energy_cost = make_keyed_array(model_data.index_y_fix, ipp.index_p)
    for y in model_data.index_y_fix, p in ipp.index_p
        energy_cost(y, p, :) .=
            sum(
                model_data.omega(t) * (ipp.v_E_my(y, p, k, t) * ipp.y_E_my(y, p, k, t)) for
                t in model_data.index_t, k in ipp.index_k_existing
            ) + sum(
                model_data.omega(t) * (ipp.v_C_my(y, p, k, t) * ipp.y_C_my(y, p, k, t)) for
                t in model_data.index_t, k in ipp.index_k_new
            )
    end
    fixed_om = make_keyed_array(model_data.index_y_fix, ipp.index_p)
    for y in model_data.index_y_fix, p in ipp.index_p
        fixed_om(y, p, :) .=
            sum(
                ipp.fom_E_my[y, p, k] * (
                    ipp.x_E_my(p, k) - sum(
                        ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                        model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                    )
                ) for k in ipp.index_k_existing
            ) + sum(
                ipp.fom_C_my(Symbol(Int(y_symbol)), p, k) *
                ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for k in ipp.index_k_new,
                y_symbol in
                model_data.year(first(model_data.index_y_fix)):model_data.year(y)
            )
    end
    operational_cost = make_keyed_array(model_data.index_y_fix, ipp.index_p)
    for y in model_data.index_y_fix, p in ipp.index_p
        operational_cost(y, p, :) .= energy_cost(y, p) + fixed_om(y, p)
    end
    working_capital = make_keyed_array(model_data.index_y_fix, ipp.index_p)
    for y in model_data.index_y_fix, p in ipp.index_p
        working_capital(y, p, :) .= utility.DaysofWC / 365 * operational_cost(y, p)
    end

    # assume the ipps' new depreciation schedule is the same as utility's
    ADITNew = make_keyed_array(model_data.index_y_fix, ipp.index_p, ipp.index_k_new)
    for y in model_data.index_y_fix, p in ipp.index_p, k in ipp.index_k_new
        ADITNew[y, p, k] = sum(
            ipp.CapEx_my(Symbol(Int(y_symbol)), p, k) *
            ipp.x_C_my(Symbol(Int(y_symbol)), p, k) *
            (
                utility.CumuTaxDepre_new_my[
                    Symbol(Int(model_data.year(y) - y_symbol + 1)),
                    k,
                ] - utility.CumuAccoutDepre_new_my[
                    Symbol(Int(model_data.year(y) - y_symbol + 1)),
                    k,
                ]
            ) *
            ipp.Tax(p) +
            utility.ITC_new_my[Symbol(Int(y_symbol)), k] *
            ipp.CapEx_my(Symbol(Int(y_symbol)), p, k) *
            ipp.x_C_my(Symbol(Int(y_symbol)), p, k) *
            (
                1 - utility.CumuITCAmort_new_my[
                    Symbol(Int(model_data.year(y) - y_symbol + 1)),
                    k,
                ]
            ) for y_symbol in
            model_data.year(first(model_data.index_y_fix)):model_data.year(y)
        )
    end
    RateBaseNoWC_new = make_keyed_array(model_data.index_y_fix, ipp.index_p, ipp.index_k_new)
    for y in model_data.index_y_fix, p in ipp.index_p, k in ipp.index_k_new
        RateBaseNoWC_new[y, p, k] =
            sum(
                ipp.CapEx_my(Symbol(Int(y_symbol)), p, k) *
                ipp.x_C_my(Symbol(Int(y_symbol)), p, k) *
                (
                    1 - utility.CumuAccoutDepre_new_my[
                        Symbol(Int(model_data.year(y) - y_symbol + 1)),
                        k,
                    ]
                ) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):model_data.year(y)
            ) - ADITNew[y, p, k]
    end

    rate_base = make_keyed_array(model_data.index_y_fix, ipp.index_p)
    for y in model_data.index_y_fix, p in ipp.index_p
        rate_base(y, p, :) .=
            sum(
                utility.RateBaseNoWC_existing_my(y, k) * (
                    ipp.x_E_my(p, k) - sum(
                        ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                        model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                    )
                ) for k in ipp.index_k_existing
            ) +
            sum(RateBaseNoWC_new[y, p, k] for k in ipp.index_k_new) +
            working_capital(y, p)
    end
    debt_interest = make_keyed_array(model_data.index_y_fix, ipp.index_p)
    for y in model_data.index_y_fix, p in ipp.index_p
        debt_interest(y, p, :) .= rate_base(y, p) * ipp.DebtRatio(p) * ipp.COD(p)
    end

    depreciation = make_keyed_array(model_data.index_y_fix, ipp.index_p)
    for y in model_data.index_y_fix, p in ipp.index_p
        depreciation(y, p, :) .=
            sum(
                utility.CapEx_existing_my[k] *
                (
                    ipp.x_E_my(p, k) - sum(
                        ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                        model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                    )
                ) *
                utility.AnnualAccoutDepre_existing_my(y, k) +
                utility.CapEx_existing_my[k] *
                ipp.x_R_my[y, p, k] *
                (
                    utility.AnnualAccoutDepre_existing_my(y, k) + 1 -
                    utility.CumuAccoutDepre_existing_my(y, k)
                ) for k in ipp.index_k_existing
            ) + sum(
                ipp.CapEx_my(Symbol(Int(y_symbol)), p, k) *
                ipp.x_C_my(Symbol(Int(y_symbol)), p, k) *
                utility.AnnualAccoutDepre_new_my[
                    Symbol(Int(model_data.year(y) - y_symbol + 1)),
                    k,
                ] for y_symbol in
                model_data.year(first(model_data.index_y_fix)):model_data.year(y),
                k in ipp.index_k_new
            )
    end

    depreciation_tax = make_keyed_array(model_data.index_y_fix, ipp.index_p)
    for y in model_data.index_y_fix, p in ipp.index_p
        depreciation_tax(y, p, :) .=
            sum(
                utility.CapEx_existing_my[k] *
                (
                    ipp.x_E_my(p, k) - sum(
                        ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                        model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                    )
                ) *
                utility.AnnualTaxDepre_existing_my(y, k) +
                utility.CapEx_existing_my[k] *
                ipp.x_R_my[y, p, k] *
                (
                    utility.AnnualTaxDepre_existing_my(y, k) + 1 -
                    utility.CumuTaxDepre_existing_my(y, k)
                ) for k in ipp.index_k_existing
            ) + sum(
                ipp.CapEx_my(Symbol(Int(y_symbol)), p, k) *
                ipp.x_C_my(Symbol(Int(y_symbol)), p, k) *
                utility.AnnualTaxDepre_new_my[
                    Symbol(Int(model_data.year(y) - y_symbol + 1)),
                    k,
                ] for y_symbol in
                model_data.year(first(model_data.index_y_fix)):model_data.year(y),
                k in ipp.index_k_new
            )
    end

    income_tax = make_keyed_array(model_data.index_y_fix, ipp.index_p)
    for y in model_data.index_y_fix, p in ipp.index_p
        income_tax(y, p, :) .=
            (
                IPP_Revenue_p(y, p) - debt_interest(y, p) - operational_cost(y, p) -
                depreciation_tax(y, p)
            ) * ipp.Tax(p) - 
            sum(
                utility.ITC_existing_my[k] *
                utility.CapEx_existing_my[k] *
                (
                    ipp.x_E_my(p, k) - sum(
                        ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                        model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                    )
                ) *
                utility.AnnualITCAmort_existing_my(y, k) +
                # existing units that are retired this year will incur their regular annual depreciation, as well as the remaining un-depreciated asset
                utility.ITC_existing_my[k] *
                utility.CapEx_existing_my[k] *
                ipp.x_R_my[y, p, k] *
                (
                    utility.AnnualITCAmort_existing_my(y, k) + 1 -
                    utility.CumuITCAmort_existing_my(y, k)
                ) for k in ipp.index_k_existing
            ) -
            sum(
                utility.ITC_new_my[Symbol(Int(y_symbol)), k] *
                ipp.CapEx_my(Symbol(Int(y_symbol)), p, k) *
                ipp.x_C_my(Symbol(Int(y_symbol)), p, k) *
                utility.AnnualITCAmort_new_my[Symbol(Int(model_data.year(y) - y_symbol + 1)), k] for
                y_symbol in model_data.year(first(model_data.index_y_fix)):model_data.year(y), k in ipp.index_k_new
            )
    end

    IPP_Cost_p = make_keyed_array(model_data.index_y_fix, ipp.index_p)
    for y in model_data.index_y_fix, p in ipp.index_p
        IPP_Cost_p(y, p, :) .=
            debt_interest(y, p) +
            income_tax(y, p) +
            operational_cost(y, p) +
            depreciation(y, p)
    end

    IPP_Cost_total = KeyedArray(
        [
            sum(IPP_Cost_p(y, p) for p in ipp.index_p) + regulator.othercost(y) for
            y in model_data.index_y_fix
        ];
        [get_pair(model_data.index_y_fix)]...
    )

    IPP_debt_interest_my = KeyedArray(
        [sum(debt_interest(y, p) for p in ipp.index_p) for y in model_data.index_y_fix];
        [get_pair(model_data.index_y_fix)]...
    )
    IPP_income_tax_my = KeyedArray(
        [sum(income_tax(y, p) for p in ipp.index_p) for y in model_data.index_y_fix];
        [get_pair(model_data.index_y_fix)]...
    )
    IPP_operational_cost_my = KeyedArray(
        [sum(operational_cost(y, p) for p in ipp.index_p) for y in model_data.index_y_fix];
        [get_pair(model_data.index_y_fix)]...
    )
    IPP_depreciation_my = KeyedArray(
        [sum(depreciation(y, p) for p in ipp.index_p) for y in model_data.index_y_fix];
        [get_pair(model_data.index_y_fix)]...
    )
    IPP_depreciation_tax_my = KeyedArray(
        [sum(depreciation_tax(y, p) for p in ipp.index_p) for y in model_data.index_y_fix];
        [get_pair(model_data.index_y_fix)]...
    )
    IPP_total_emission_my = KeyedArray(
        [
            sum(
                model_data.omega(t) * (
                    sum(
                        ipp.y_E_my(y, p, k, t) * ipp.emission_rate_E_my[y, p, k] for
                        k in utility.index_k_existing, p in ipp.index_p
                    ) + sum(
                        ipp.y_C_my(y, p, k, t) * ipp.emission_rate_C_my[y, p, k] for
                        k in utility.index_k_new, p in ipp.index_p
                    )
                ) for t in model_data.index_t
            ) * 0.000453592 for y in model_data.index_y_fix
        ];
        [get_pair(model_data.index_y_fix)]...
    )

    return IPP_Revenue_total,
    IPP_Cost_total,
    IPP_debt_interest_my,
    IPP_income_tax_my,
    IPP_operational_cost_my,
    IPP_depreciation_my,
    IPP_depreciation_tax_my,
    IPP_total_emission_my,
    IPP_Revenue_p,
    IPP_Cost_p
end

"""
Update IPPGroup cumulative parameters
"""
function update_cumulative!(model_data::HEMData, ipp::IPPGroup)
    for p in ipp.index_p, k in ipp.index_k_existing
        ipp.x_R_cumu(p,k, :) .= ipp.x_R_cumu(p,k) + ipp.x_R_my(first(model_data.index_y),p,k)
    end

    for p in ipp.index_p, k in ipp.index_k_new
        ipp.x_C_cumu(p,k, :) .= ipp.x_C_cumu(p,k) + ipp.x_C_my(first(model_data.index_y),p,k)
    end
end
