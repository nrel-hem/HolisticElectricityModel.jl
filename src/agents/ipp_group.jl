# This module defines the data and functions associated with the Independent Power Producer

# declare customer decision
abstract type IPPAlgorithm end
struct LagrangeDecomposition <: IPPAlgorithm end
struct MIQP <: IPPAlgorithm end
struct MPPDCMER <: IPPAlgorithm end

abstract type AbstractIPPOptions <: AgentOptions end

struct IPPOptions{T <: IPPAlgorithm} <: AbstractIPPOptions
    ipp_algorithm::T
    solvers::Dict{String, <:HEMSolver}
end

"""
Construct IPPOptions with solvers defined as MOI.OptimizerWithAttributes instances or
HEMSolver instances.

Review the functions that call `get_new_jump_model` for Dict key requirements.

# Examples
```julia
julia> ipp_opts = IPPOptions(
    LagrangeDecomposition(),
    Dict(
        "Lagrange_Sub_Investment_Retirement_Cap" => JuMP.optimizer_with_attributes(
            Ipopt.Optimizer,
            "print_level" => 0,
        ),
        "Lagrange_Sub_Dispatch_Cap" => JuMP.optimizer_with_attributes(
            Xpress.Optimizer,
        ),
        "Lagrange_Feasible_Cap" => JuMP.optimizer_with_attributes(
            Xpress.Optimizer,
        )
    )
)
```
"""
function IPPOptions(algorithm::IPPAlgorithm, solvers::Dict)
    hem_solvers = Dict{String, HEMSolver}()
    for (key, val) in solvers
        if val isa MOI.OptimizerWithAttributes
            hem_solvers[key] = AnySolver(val)
        else
            hem_solvers[key] = val
        end
    end

    return IPPOptions(algorithm, hem_solvers)
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
    x_E::ParamArray
    "fixed cost of existing capacity (\$/MW-yr)"
    f_E::ParamArray
    "fixed cost of new capacity (\$/MW-yr)"
    f_C::ParamArray
    "variable cost of existing capacity (\$/MWh)"
    v_E::ParamArray
    "variable cost of new capacity (\$/MWh)"
    v_C::ParamArray
    "availability of existing capacity (fraction)"
    rho_E::ParamArray
    "availability of new capacity (fraction)"
    rho_C::ParamArray
    eximport::ParamArray # net export (MWh)
    Peak_eximport::ParamScalar
    "Big M Parameter"
    B1GM::ParamScalar{<:Integer}
    zeta::ParamScalar # offer price factor cap

    # Primal Variables
    y_E::ParamArray
    y_C::ParamArray
    x_R::ParamArray
    x_C::ParamArray
    miu::ParamArray
    o_E::ParamArray # offer price of existing capacity ($/MWh)
    o_C::ParamArray # offer price of new capacity ($/MWh)
    LMP::ParamArray

    # Parameters (multi-year)
    x_E_my::ParamArray # existing capacity (MW)
    fom_E_my::ParamArray # fixed O&M of existing capacity ($/MW-yr)
    fom_C_my::ParamArray # fixed O&M of new capacity ($/MW-yr)
    CapEx_my::ParamArray # capital expense of new capacity ($/MW)
    v_E_my::ParamArray # variable cost of existing capacity ($/MWh)
    v_C_my::ParamArray # variable cost of new capacity ($/MWh)
    rho_E_my::ParamArray # availability of existing capacity (fraction)
    rho_C_my::ParamArray # availability of new capacity (fraction)
    eximport_my::ParamArray # net export (MWh)
    pvf_cap::ParamArray # present value factor of capital expenses
    pvf_onm::ParamArray # present value factor of o&m expenses
    CRF_default::ParamArray
    Tax::ParamArray
    DebtRatio::ParamArray
    COD::ParamArray
    COE::ParamArray

    # Primal Variables (multi-year)
    y_E_my::ParamArray
    y_C_my::ParamArray
    x_R_my::ParamArray
    x_C_my::ParamArray
    o_E_my::ParamArray # offer price of existing capacity ($/MWh)
    o_C_my::ParamArray # offer price of new capacity ($/MWh)
    # Dual Variables (multi-year)
    miu_my::ParamArray
    LMP_my::ParamArray
    eta_my::ParamArray
    lambda_my::ParamArray
    u_y_E_my::ParamArray
    u_y_C_my::ParamArray
    u_miu_my::ParamArray
    u_eta_my::ParamArray
    u_lambda_my::ParamArray
    x_R_cumu::ParamArray
    x_C_cumu::ParamArray

    # Lagrange multiplier
    L_R::ParamArray
    L_C::ParamArray
    x_R_my_st1::ParamArray
    x_C_my_st1::ParamArray
    x_R_my_st2::ParamArray
    x_C_my_st2::ParamArray
    # value of objective function
    obj_st1::ParamScalar
    obj_st2::ParamScalar
    obj_lower_bound::ParamScalar
    obj_upper_bound::ParamScalar
    obj_feasible::ParamScalar

    #temporary save_results
    y_E_my_temp::ParamArray
    y_C_my_temp::ParamArray
    x_R_my_temp::ParamArray
    x_C_my_temp::ParamArray
    miu_my_temp::ParamArray
    LMP_my_temp::ParamArray

    # capacity market parameters
    NetCONE::ParamArray
    DC_length::ParamArray
    capacity_credit_E_my::ParamArray # capacity credit of existing resources
    capacity_credit_C_my::ParamArray # capacity credit of new resources
    capacity_price::ParamArray # $/MW-yr
    capacity_price_my_temp::ParamArray

    Net_Load_my::ParamArray
    Max_Net_Load_my::ParamArray
    Reserve_req_my::ParamArray
    Capacity_slope_my::ParamArray
    Capacity_intercept_my::ParamArray
    ucap_temp::ParamArray
    ucap::ParamArray

    # RPS
    RPS::ParamArray

    # emission rate
    emission_rate_E_my::ParamArray
    emission_rate_C_my::ParamArray
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
        ParamArray("f_C", Tuple(push!(copy([index_p]), index_k_new)), FixedCostNew),
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
        ParamArray(
            "pvf_cap",
            Tuple(push!(copy([model_data.index_y]), index_p)),
            pvf_cap,
        ),
        ParamArray(
            "pvf_onm",
            Tuple(push!(copy([model_data.index_y]), index_p)),
            pvf_onm,
        ),
        ParamArray("CRF_default", (index_p,), CRF_default),
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

    WMDER_IPP = get_new_jump_model(ipp_opts.solvers["Lagrange_Sub_Investment_Retirement_Cap"])

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
        ipp.capacity_credit_E_my(y, k, :) .= ipp.rho_E_my(p_star, k, Max_Net_Load_my_index(y))
    end
    fill!(ipp.capacity_credit_C_my, NaN)
    for y in model_data.index_y, k in ipp.index_k_new
        ipp.capacity_credit_C_my(y, k, :) .= ipp.rho_C_my(p_star, k, Max_Net_Load_my_index(y))
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
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
                ) +
                sum(
                    ipp.capacity_credit_C_my(y, k) * (
                        sum(
                            ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p, k)
                    ) for k in ipp.index_k_new,
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
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
                    ipp.pvf_onm(Symbol(Int(y_symbol)), p_star) for y_symbol in
                    model_data.year(y):model_data.year(last(model_data.index_y.elements))
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
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
                ) +
                sum(
                    ipp.capacity_credit_C_my(y, k) * (
                        sum(
                            ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p, k)
                    ) for k in ipp.index_k_new,
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
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
                    ipp.rho_E_my(p, k, t) * (
                        ipp.x_E_my(p, k) - sum(
                            ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_R_cumu(p, k)
                    ) for k in ipp.index_k_existing,
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
                ) +
                sum(
                    ipp.rho_C_my(p, k, t) * (
                        sum(
                            ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p, k)
                    ) for k in ipp.index_k_new,
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
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

    WMDER_IPP = get_new_jump_model(ipp_opts.solvers["Lagrange_Sub_Dispatch_Cap"])

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
                        ipp.rho_E_my(p, k, t) *
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
                                ipp.x_R_my(Symbol(Int(y_symbol)), p_star, k) for
                                y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) - ipp.x_R_cumu(p_star, k)
                        ) for t in model_data.index_t, k in ipp.index_k_existing
                    )
                ) - (
                    sum(
                        lambda[y, p, k, t] *
                        ipp.rho_C_my(p, k, t) *
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
                                ipp.x_C_my(Symbol(Int(y_symbol)), p_star, k) for
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
                        (ipp.v_E_my(y, p_star, k, t) * y_E[y, p_star, k, t]) for
                        t in model_data.index_t, k in ipp.index_k_existing
                    )
                ) - (
                    sum(
                        model_data.omega(t) * (ipp.v_C_my(y, p, k, t) * y_C[y, p, k, t]) for
                        t in model_data.index_t, k in ipp.index_k_new, p in ipp.index_p
                    ) - sum(
                        model_data.omega(t) *
                        (ipp.v_C_my(y, p_star, k, t) * y_C[y, p_star, k, t]) for
                        t in model_data.index_t, k in ipp.index_k_new
                    )
                ) -

                # generation costs (this includes terms due to revenue linearization)
                #   num hrs * ((fuel + vom) * gen existing + (fuel + vom) * gen new) for t and gen type
                sum(
                    model_data.omega(t) *
                    (ipp.v_E_my(y, p_star, k, t) * y_E[y, p_star, k, t]) for
                    t in model_data.index_t, k in ipp.index_k_existing
                ) - sum(
                    model_data.omega(t) *
                    (ipp.v_C_my(y, p_star, k, t) * y_C[y, p_star, k, t]) for
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
                p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p))),
                k in ipp.index_k_existing,
                t in model_data.index_t,
            ],
            eta[y, p, k, t] <= ipp.B1GM.value * u_eta[y, p, k, t]
        )
        @constraint(
            WMDER_IPP,
            Eq_eta_p_other_2[
                y in model_data.index_y,
                p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p))),
                k in ipp.index_k_existing,
                t in model_data.index_t,
            ],
            ipp.rho_E_my(p, k, t) * (
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
                p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p))),
                k in ipp.index_k_existing,
                t in model_data.index_t,
            ],
            ipp.rho_E_my(p, k, t) * (
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
                p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p))),
                k in ipp.index_k_new,
                t in model_data.index_t,
            ],
            lambda[y, p, k, t] <= ipp.B1GM.value * u_lambda[y, p, k, t]
        )
        @constraint(
            WMDER_IPP,
            Eq_lambda_p_other_2[
                y in model_data.index_y,
                p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p))),
                k in ipp.index_k_new,
                t in model_data.index_t,
            ],
            ipp.rho_C_my(p, k, t) * (
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
                p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p))),
                k in ipp.index_k_new,
                t in model_data.index_t,
            ],
            ipp.rho_C_my(p, k, t) * (
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
                    ipp.rho_E_my(p, k, t) * (
                        ipp.x_E_my(p, k) - sum(
                            ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_R_cumu(p, k)
                    ) for k in ipp.index_k_existing,
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
                ) +
                sum(
                    ipp.rho_C_my(p, k, t) * (
                        sum(
                            ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p, k)
                    ) for k in ipp.index_k_new,
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
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
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
                ) +
                sum(
                    ipp.capacity_credit_C_my(y, k) * (
                        sum(
                            ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p, k)
                    ) for k in ipp.index_k_new,
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
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

    WMDER_IPP = get_new_jump_model(ipp_opts.solvers["Lagrange_Feasible_Cap"])

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
                        ipp.rho_E_my(p, k, t) *
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
                                ipp.x_R_my(Symbol(Int(y_symbol)), p_star, k) for
                                y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) - ipp.x_R_cumu(p_star, k)
                        ) for t in model_data.index_t, k in ipp.index_k_existing
                    )
                ) - (
                    sum(
                        lambda[y, p, k, t] *
                        ipp.rho_C_my(p, k, t) *
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
                                ipp.x_C_my(Symbol(Int(y_symbol)), p_star, k) for
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
                        (ipp.v_E_my(y, p_star, k, t) * y_E[y, p_star, k, t]) for
                        t in model_data.index_t, k in ipp.index_k_existing
                    )
                ) - (
                    sum(
                        model_data.omega(t) * (ipp.v_C_my(y, p, k, t) * y_C[y, p, k, t]) for
                        t in model_data.index_t, k in ipp.index_k_new, p in ipp.index_p
                    ) - sum(
                        model_data.omega(t) *
                        (ipp.v_C_my(y, p_star, k, t) * y_C[y, p_star, k, t]) for
                        t in model_data.index_t, k in ipp.index_k_new
                    )
                ) -

                # generation costs (this includes terms due to revenue linearization)
                #   num hrs * ((fuel + vom) * gen existing + (fuel + vom) * gen new) for t and gen type
                sum(
                    model_data.omega(t) *
                    (ipp.v_E_my(y, p_star, k, t) * y_E[y, p_star, k, t]) for
                    t in model_data.index_t, k in ipp.index_k_existing
                ) - sum(
                    model_data.omega(t) *
                    (ipp.v_C_my(y, p_star, k, t) * y_C[y, p_star, k, t]) for
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
                ipp.x_R_my_st2(Symbol(Int(y_symbol)), k) for
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
                ipp.x_R_my_st2(Symbol(Int(y_symbol)), k) for
                y_symbol in model_data.year(first(model_data.index_y)):model_data.year(y)
            ) - ipp.x_R_cumu(p_star, k)
        ) - y_E[y, p_star, k, t] >= 0
    )
    if length(ipp.index_p) >= 2
        @constraint(
            WMDER_IPP,
            Eq_eta_p_other_1[
                y in model_data.index_y,
                p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p))),
                k in ipp.index_k_existing,
                t in model_data.index_t,
            ],
            eta[y, p, k, t] <= ipp.B1GM.value * u_eta[y, p, k, t]
        )
        @constraint(
            WMDER_IPP,
            Eq_eta_p_other_2[
                y in model_data.index_y,
                p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p))),
                k in ipp.index_k_existing,
                t in model_data.index_t,
            ],
            ipp.rho_E_my(p, k, t) * (
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
                p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p))),
                k in ipp.index_k_existing,
                t in model_data.index_t,
            ],
            ipp.rho_E_my(p, k, t) * (
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
                ipp.x_C_my_st2(Symbol(Int(y_symbol)), k) for
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
                ipp.x_C_my_st2(Symbol(Int(y_symbol)), k) for
                y_symbol in model_data.year(first(model_data.index_y)):model_data.year(y)
            ) + ipp.x_C_cumu(p_star, k)
        ) - y_C[y, p_star, k, t] >= 0
    )
    if length(ipp.index_p) >= 2
        @constraint(
            WMDER_IPP,
            Eq_lambda_p_other_1[
                y in model_data.index_y,
                p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p))),
                k in ipp.index_k_new,
                t in model_data.index_t,
            ],
            lambda[y, p, k, t] <= ipp.B1GM.value * u_lambda[y, p, k, t]
        )
        @constraint(
            WMDER_IPP,
            Eq_lambda_p_other_2[
                y in model_data.index_y,
                p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p))),
                k in ipp.index_k_new,
                t in model_data.index_t,
            ],
            ipp.rho_C_my(p, k, t) * (
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
                p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p))),
                k in ipp.index_k_new,
                t in model_data.index_t,
            ],
            ipp.rho_C_my(p, k, t) * (
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
                        ipp.x_R_my_st2(Symbol(Int(y_symbol)), k) for y_symbol in
                        model_data.year(first(model_data.index_y)):model_data.year(y)
                    ) - ipp.x_R_cumu(p_star, k)
                ) for k in ipp.index_k_existing
            ) + sum(
                ipp.capacity_credit_C_my(y, k) * (
                    sum(
                        ipp.x_C_my_st2(Symbol(Int(y_symbol)), k) for y_symbol in
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
                UCAP_p_star(y) +
                sum(
                    ipp.capacity_credit_E_my(y, k) * (
                        ipp.x_E_my(p, k) - sum(
                            ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_R_cumu(p, k)
                    ) for k in ipp.index_k_existing,
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
                ) +
                sum(
                    ipp.capacity_credit_C_my(y, k) * (
                        sum(
                            ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p, k)
                    ) for k in ipp.index_k_new,
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
                ) for y in model_data.index_y
            ];
            [get_pair(model_data.index_y)]...
        )
    else
        UCAP_total = KeyedArray(
            [UCAP_p_star(y) for y in model_data.index_y];
            [get_pair(model_data.index_y)]...
        )
    end

    update!(
        ipp.obj_feasible,
        objective_value(WMDER_IPP) + sum(
            # capacity revenue
            ipp.pvf_onm(y, p_star) * (
                UCAP_p_star(y) *
                (ipp.Capacity_intercept_my(y) + ipp.Capacity_slope_my(y) * UCAP_total(y))
            ) -
            # fixed costs
            #   fom * (cap exist - cap retiring) for gen type
            ipp.pvf_onm(y, p_star) * sum(
                ipp.fom_E_my(y, p_star, k) * (
                    ipp.x_E_my(p_star, k) - sum(
                        ipp.x_R_my_st2(Symbol(Int(y_symbol)), k) for y_symbol in
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
                    ipp.pvf_onm(Symbol(Int(y_symbol)), p_star) for y_symbol in
                    model_data.year(y):model_data.year(last(model_data.index_y.elements))
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
        ipp.ucap_temp(y, p_star, :) .= UCAP_p_star(y)
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

    x_R_before = ParamArray(ipp.x_R_my)
    x_C_before = ParamArray(ipp.x_C_my)

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
                ipp.x_R_my(y, p, k, :) .= ipp.x_R_my_temp(y, p, k)
            end
            for y in model_data.index_y, p in ipp.index_p, k in ipp.index_k_new
                ipp.x_C_my(y, p, k, :) .= ipp.x_C_my_temp(y, p, k)
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
    x_R_before = ParamArray(ipp.x_R_my)
    x_C_before = ParamArray(ipp.x_C_my)

    utility = get_agent(Utility, agent_store)
    regulator = get_agent(Regulator, agent_store)
    customers = get_agent(CustomerGroup, agent_store)
    green_developer = get_agent(GreenDeveloper, agent_store)

    WMDER_IPP = get_new_jump_model(ipp_opts.solvers["solve_agent_problem_ipp_cap"])

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

    Max_Net_Load_my_index = KeyedArray(
        [
            findmax(Dict(t => ipp.Net_Load_my(y, t) for t in model_data.index_t))[2] for
            y in model_data.index_y
        ];
        [get_pair(model_data.index_y)]...
    )


    fill!(ipp.capacity_credit_E_my, NaN)
    for y in model_data.index_y, k in ipp.index_k_existing
        ipp.capacity_credit_E_my(y, k, :) .= ipp.rho_E_my(p_star, k, Max_Net_Load_my_index(y))
    end
    fill!(ipp.capacity_credit_C_my, NaN)
    for y in model_data.index_y, k in ipp.index_k_new
        ipp.capacity_credit_C_my(y, k, :) .= ipp.rho_C_my(p_star, k, Max_Net_Load_my_index(y))
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
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
                ) +
                sum(
                    ipp.capacity_credit_C_my(y, k) * (
                        sum(
                            ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p, k)
                    ) for k in ipp.index_k_new,
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
                ) +
                # green technology subscription
                sum(
                    ipp.capacity_credit_C_my(y, j) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
                    model_data.year(first(model_data.index_y_fix)):model_data.year(y))
                    for j in model_data.index_j, h in model_data.index_h
                )
            end
    else
        UCAP_total = y -> begin
            UCAP_p_star(y) +
            # green technology subscription
            sum(
                ipp.capacity_credit_C_my(y, j) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
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
                    ipp.pvf_onm(Symbol(Int(y_symbol)), p_star) for y_symbol in
                    model_data.year(y):model_data.year(last(model_data.index_y.elements))
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
                        ipp.rho_E_my(p, k, t) *
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
                                ipp.x_R_my(Symbol(Int(y_symbol)), p_star, k) for
                                y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) - ipp.x_R_cumu(p_star, k)
                        ) for t in model_data.index_t, k in ipp.index_k_existing
                    )
                ) - (
                    sum(
                        lambda[y, p, k, t] *
                        ipp.rho_C_my(p, k, t) *
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
                                ipp.x_C_my(Symbol(Int(y_symbol)), p_star, k) for
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
                        (ipp.v_E_my(y, p_star, k, t) * y_E[y, p_star, k, t]) for
                        t in model_data.index_t, k in ipp.index_k_existing
                    )
                ) - (
                    sum(
                        model_data.omega(t) * (ipp.v_C_my(y, p, k, t) * y_C[y, p, k, t]) for
                        t in model_data.index_t, k in ipp.index_k_new, p in ipp.index_p
                    ) - sum(
                        model_data.omega(t) *
                        (ipp.v_C_my(y, p_star, k, t) * y_C[y, p_star, k, t]) for
                        t in model_data.index_t, k in ipp.index_k_new
                    )
                ) -

                # generation costs (this includes terms due to revenue linearization)
                #   num hrs * ((fuel + vom) * gen existing + (fuel + vom) * gen new) for t and gen type
                sum(
                    model_data.omega(t) *
                    (ipp.v_E_my(y, p_star, k, t) * y_E[y, p_star, k, t]) for
                    t in model_data.index_t, k in ipp.index_k_existing
                ) - sum(
                    model_data.omega(t) *
                    (ipp.v_C_my(y, p_star, k, t) * y_C[y, p_star, k, t]) for
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
                p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p))),
                k in ipp.index_k_existing,
                t in model_data.index_t,
            ],
            eta[y, p, k, t] <= ipp.B1GM.value * u_eta[y, p, k, t]
        )
        @constraint(
            WMDER_IPP,
            Eq_eta_p_other_2[
                y in model_data.index_y,
                p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p))),
                k in ipp.index_k_existing,
                t in model_data.index_t,
            ],
            ipp.rho_E_my(p, k, t) * (
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
                p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p))),
                k in ipp.index_k_existing,
                t in model_data.index_t,
            ],
            ipp.rho_E_my(p, k, t) * (
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
                p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p))),
                k in ipp.index_k_new,
                t in model_data.index_t,
            ],
            lambda[y, p, k, t] <= ipp.B1GM.value * u_lambda[y, p, k, t]
        )
        @constraint(
            WMDER_IPP,
            Eq_lambda_p_other_2[
                y in model_data.index_y,
                p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p))),
                k in ipp.index_k_new,
                t in model_data.index_t,
            ],
            ipp.rho_C_my(p, k, t) * (
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
                p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p))),
                k in ipp.index_k_new,
                t in model_data.index_t,
            ],
            ipp.rho_C_my(p, k, t) * (
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
                    ipp.rho_E_my(p, k, t) * (
                        ipp.x_E_my(p, k) - sum(
                            ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_R_cumu(p, k)
                    ) for k in ipp.index_k_existing,
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
                ) +
                sum(
                    ipp.rho_C_my(p, k, t) * (
                        sum(
                            ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p, k)
                    ) for k in ipp.index_k_new,
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
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
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
                ) +
                sum(
                    ipp.capacity_credit_C_my(y, k) * (
                        sum(
                            ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p, k)
                    ) for k in ipp.index_k_new,
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
                ) +
                # green technology subscription
                sum(
                    ipp.capacity_credit_C_my(y, j) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
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
                    ipp.capacity_credit_C_my(y, j) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
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
    #     sum(ipp.capacity_credit_E_my(y,k)*(ipp.x_E_my(p_star,k)-sum(ipp.x_R_my(Symbol(y_symbol),p_star,k) for y_symbol = model_data.year(first(model_data.index_y)):model_data.year(y)) - ipp.x_R_cumu[p_star,k]) for k in ipp.index_k_existing) + 
    #     sum(ipp.capacity_credit_C_my(y,k)*(sum(ipp.x_C_my(Symbol(y_symbol),p_star,k) for y_symbol = model_data.year(first(model_data.index_y)):model_data.year(y)) + ipp.x_C_cumu[p_star,k]) for k in ipp.index_k_new)
    #     for y in model_data.index_y
    # )
    # capacity_profit = Dict(y =>
    #     ipp.pvf_onm(y,p_star) * UCAP_p_star_solution(y) * (ipp.Capacity_intercept_my(y) + ipp.Capacity_slope_my(y) * UCAP_p_star_solution(y))
    #     for y in model_data.index_y
    # )
    # total_profit = Dict(y => capacity_profit(y) + objective_value(WMDER_IPP) for y in model_data.index_y)

    ###### running MIQP ######
    UCAP_p_star = KeyedArray(
        [
            sum(
                ipp.capacity_credit_E_my(y, k) * (
                    ipp.x_E_my(p_star, k) - sum(
                        ipp.x_R_my(Symbol(Int(y_symbol)), p_star, k) for y_symbol in
                        model_data.year(first(model_data.index_y)):model_data.year(y)
                    ) - ipp.x_R_cumu(p_star, k)
                ) for k in ipp.index_k_existing
            ) + sum(
                ipp.capacity_credit_C_my(y, k) * (
                    sum(
                        ipp.x_C_my(Symbol(Int(y_symbol)), p_star, k) for y_symbol in
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
                UCAP_p_star(y) +
                sum(
                    ipp.capacity_credit_E_my(y, k) * (
                        ipp.x_E_my(p, k) - sum(
                            ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_R_cumu(p, k)
                    ) for k in ipp.index_k_existing,
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
                ) +
                sum(
                    ipp.capacity_credit_C_my(y, k) * (
                        sum(
                            ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p, k)
                    ) for k in ipp.index_k_new,
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
                ) +
                # green technology subscription
                sum(
                    ipp.capacity_credit_C_my(y, j) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
                    model_data.year(first(model_data.index_y_fix)):model_data.year(y))
                    for j in model_data.index_j, h in model_data.index_h
                ) 
                for y in model_data.index_y
            ];
            [get_pair(model_data.index_y)]...
        )
    else
        UCAP_total = KeyedArray(
            [UCAP_p_star(y) +
            # green technology subscription
            sum(
                ipp.capacity_credit_C_my(y, j) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):model_data.year(y))
                for j in model_data.index_j, h in model_data.index_h
            ) for y in model_data.index_y];
            [get_pair(model_data.index_y)]...
        )
    end

    for y in model_data.index_y
        ipp.capacity_price(y, :) .=
            ipp.Capacity_intercept_my(y) + ipp.Capacity_slope_my(y) * UCAP_total(y)
        ipp.ucap(y, p_star, :) .= UCAP_p_star(y)
    end

    # report lower level duality gap
    lower_level_primal_obj = DataFrame()
    lower_level_dual_obj = DataFrame()
    lower_level_duality_gap = DataFrame()
    for y in model_data.index_y
        lower_level_primal_obj[!, y] = Vector{Float64}(undef, 1)
        lower_level_dual_obj[!, y] = Vector{Float64}(undef, 1)
        lower_level_duality_gap[!, y] = Vector{Float64}(undef, 1)
    end
    for y in model_data.index_y
        lower_level_primal_obj[1, y] = 
        sum(
            model_data.omega(t) *
            (ipp.v_E_my(y, p, k, t) * value.(y_E[y, p, k, t])) for
            t in model_data.index_t, k in ipp.index_k_existing, p in ipp.index_p
        ) + 
        sum(
            model_data.omega(t) *
            (ipp.v_C_my(y, p, k, t) * value.(y_C[y, p, k, t])) for
            t in model_data.index_t, k in ipp.index_k_new, p in ipp.index_p
        )

        if length(ipp.index_p) >= 2
            lower_level_dual_obj[1, y] = 
            sum(
                value.(miu[y, t]) * (
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
            ) - 
            (
                sum(
                    value.(eta[y, p, k, t]) *
                    ipp.rho_E_my(p, k, t) *
                    (
                        ipp.x_E_my(p, k) - ipp.x_R_cumu(p, k)
                    ) for t in model_data.index_t, k in ipp.index_k_existing, p in ipp.index_p
                ) + sum(
                    value.(lambda[y, p, k, t]) *
                    ipp.rho_C_my(p, k, t) *
                    ipp.x_C_cumu(p, k) 
                    for t in model_data.index_t, k in ipp.index_k_new, p in ipp.index_p
                )
            ) + sum(
                value.(eta[y, p, k, t]) *
                ipp.rho_E_my(p, k, t) *
                sum(
                    ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) for t in model_data.index_t, k in ipp.index_k_existing,
                p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
            ) - sum(
                value.(lambda[y, p, k, t]) *
                ipp.rho_C_my(p, k, t) *
                sum(
                    ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) for t in model_data.index_t, k in ipp.index_k_new,
                p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
            ) + sum(
                ipp.rho_E_my(p_star, k, t) *
                value.(eta[y, p_star, k, t]) *
                sum(
                    value.(x_R[Symbol(Int(y_symbol)), k]) for
                    y_symbol in model_data.year(first(model_data.index_y)):model_data.year(y)
                ) for t in model_data.index_t, k in ipp.index_k_existing
            ) - sum(
                ipp.rho_C_my(p_star, k, t) *
                value.(lambda[y, p_star, k, t]) *
                sum(
                    value.(x_C[Symbol(Int(y_symbol)), k]) for
                    y_symbol in model_data.year(first(model_data.index_y)):model_data.year(y)
                ) for t in model_data.index_t, k in ipp.index_k_new
            )
        else
            lower_level_dual_obj[1, y] = 
            sum(
                value.(miu[y, t]) * (
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
            ) - 
            (
                sum(
                    value.(eta[y, p, k, t]) *
                    ipp.rho_E_my(p, k, t) *
                    (
                        ipp.x_E_my(p, k) - ipp.x_R_cumu(p, k)
                    ) for t in model_data.index_t, k in ipp.index_k_existing, p in ipp.index_p
                ) + sum(
                    value.(lambda[y, p, k, t]) *
                    ipp.rho_C_my(p, k, t) *
                    ipp.x_C_cumu(p, k) 
                    for t in model_data.index_t, k in ipp.index_k_new, p in ipp.index_p
                )
            ) + sum(
                ipp.rho_E_my(p_star, k, t) *
                value.(eta[y, p_star, k, t]) *
                sum(
                    value.(x_R[Symbol(Int(y_symbol)), k]) for
                    y_symbol in model_data.year(first(model_data.index_y)):model_data.year(y)
                ) for t in model_data.index_t, k in ipp.index_k_existing
            ) - sum(
                ipp.rho_C_my(p_star, k, t) *
                value.(lambda[y, p_star, k, t]) *
                sum(
                    value.(x_C[Symbol(Int(y_symbol)), k]) for
                    y_symbol in model_data.year(first(model_data.index_y)):model_data.year(y)
                ) for t in model_data.index_t, k in ipp.index_k_new
            )
        end
    end
    lower_level_duality_gap = (lower_level_primal_obj .- lower_level_dual_obj) ./ abs.(lower_level_primal_obj)
    @info "lower level primal obj is $(lower_level_primal_obj)"
    @info "lower level dual obj is $(lower_level_dual_obj)"
    @info "lower level duality gap is $(lower_level_duality_gap)"

    return compute_difference_percentage_one_norm([
        (x_R_before, ipp.x_R_my),
        (x_C_before, ipp.x_C_my),
    ])
end



function solve_agent_problem_ipp_cap(
    ipp::IPPGroup,
    ipp_opts::IPPOptions{MPPDCMER},
    p_star,
    model_data::HEMData,
    hem_opts::HEMOptions{WholesaleMarket},
    agent_store::AgentStore,
    w_iter,
)
    x_R_before = ParamArray(ipp.x_R_my)
    x_C_before = ParamArray(ipp.x_C_my)

    utility = get_agent(Utility, agent_store)
    regulator = get_agent(Regulator, agent_store)
    customers = get_agent(CustomerGroup, agent_store)
    green_developer = get_agent(GreenDeveloper, agent_store)

    WMDER_IPP = get_new_jump_model(ipp_opts.solvers["solve_agent_problem_ipp_mppdc"])
    MPPDCMER_lower = get_new_jump_model(ipp_opts.solvers["solve_agent_problem_ipp_mppdc_mccormic_lower"])

    # first, use lower level optimization results to set variable bounds for McCormick-envelope Relaxation
    @variable(
        MPPDCMER_lower,
        y_E_bounds[model_data.index_y, ipp.index_p, ipp.index_k_existing, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower,
        y_C_bounds[model_data.index_y, ipp.index_p, ipp.index_k_new, model_data.index_t] >= 0
    )

    objective_function_lower = begin
        sum(
            sum(
                model_data.omega(t) *
                (ipp.v_E_my(y, p, k, t) * y_E_bounds[y, p, k, t]) for
                t in model_data.index_t, k in ipp.index_k_existing, p in ipp.index_p
            ) + 
            sum(
                model_data.omega(t) *
                (ipp.v_C_my(y, p, k, t) * y_C_bounds[y, p, k, t]) for
                t in model_data.index_t, k in ipp.index_k_new, p in ipp.index_p
            )
            for y in model_data.index_y
        )
    end

    @objective(MPPDCMER_lower, Min, objective_function_lower)

    supply_demand_balance_lower =
        (y, t) -> begin
            # bulk generation at time t
            sum(y_E_bounds[y, p, k, t] for k in ipp.index_k_existing, p in ipp.index_p) +
            sum(y_C_bounds[y, p, k, t] for k in ipp.index_k_new, p in ipp.index_p) -
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
        MPPDCMER_lower,
        Eq_primal_feasible_supplydemandbalance_lower[y in model_data.index_y, t in model_data.index_t],
        supply_demand_balance_lower(y, t) >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_gen_max_E_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            k in ipp.index_k_existing,
            t in model_data.index_t,
        ],
        ipp.rho_E_my(p, k, t) * (
            ipp.x_E_my(p, k) - sum(
                ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                model_data.year(first(model_data.index_y)):model_data.year(y)
            ) - ipp.x_R_cumu(p, k)
        ) - y_E_bounds[y, p, k, t] >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_gen_max_C_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            k in ipp.index_k_new,
            t in model_data.index_t,
        ],
        ipp.rho_C_my(p, k, t) * (
            sum(
                ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                model_data.year(first(model_data.index_y)):model_data.year(y)
            ) + ipp.x_C_cumu(p, k)
        ) - y_C_bounds[y, p, k, t] >= 0
    )

    TimerOutputs.@timeit HEM_TIMER "optimize! MPPDC McCormic Envelope Relaxation lower level" begin
        optimize!(MPPDCMER_lower)
    end

    eta_param = initialize_param("eta_param", model_data.index_y, ipp.index_k_existing, model_data.index_t)

    for y in model_data.index_y, k in ipp.index_k_existing, t in model_data.index_t
        eta_param(y, k, t) = dual.(Eq_primal_feasible_gen_max_E_lower[y, p_star,k,t])
    end

    lambda_param = initialize_param("lambda_param", model_data.index_y, ipp.index_k_new, model_data.index_t)

    for y in model_data.index_y, k in ipp.index_k_new, t in model_data.index_t
        lambda_param(y, k, t) = dual.(Eq_primal_feasible_gen_max_C_lower[y, p_star,k,t])
    end

    # test bounds McCormick-envelope Relaxation
    # eta_param = CSV.read(joinpath("/home/nguo/HolisticElectricityModel-Data/outputs/ba_1_base_2018_future_1_ipps_1", "eta.csv"), DataFrame)
    # lambda_param = CSV.read(joinpath("/home/nguo/HolisticElectricityModel-Data/outputs/ba_1_base_2018_future_1_ipps_1", "lambda.csv"), DataFrame)
 
    eta_upper_bound_adj = 1.001
    eta_lower_bound_adj = 0.999
    lambda_upper_bound_adj = 1.001
    lambda_lower_bound_adj = 0.999

    # eta_upper_bound_adj = 10.0
    # eta_lower_bound_adj = 10.0
    # lambda_upper_bound_adj = 10.0
    # lambda_lower_bound_adj = 10.0

    eta_L = initialize_param("eta_L", model_data.index_y, ipp.index_k_existing, model_data.index_t)
    eta_U = initialize_param("eta_U", model_data.index_y, ipp.index_k_existing, model_data.index_t)

    for y in model_data.index_y
        for k in ipp.index_k_existing
            for t in model_data.index_t
                eta_U(y, k, t) = eta_param(y, k, t) * eta_upper_bound_adj
                eta_L(y, k, t) = eta_param(y, k, t) * eta_lower_bound_adj
            end
        end
    end

    # for y in model_data.index_y
    #     for k in ipp.index_k_existing
    #         for t in model_data.index_t
    #             eta_U[y, k, t] = eta_param[(eta_param.Year .== model_data.year[y]) .& (eta_param.IPP .== "ipp1") .& (eta_param.GenTech .== string(k)) .& (eta_param.Time .== string(t)), "eta"][1] + 10.0
    #             eta_L[y, k, t] = eta_param[(eta_param.Year .== model_data.year[y]) .& (eta_param.IPP .== "ipp1") .& (eta_param.GenTech .== string(k)) .& (eta_param.Time .== string(t)), "eta"][1] - 10.0
    #         end
    #     end
    # end

    X_R_cumu_L = initialize_param("x_R_L", model_data.index_y, ipp.index_k_existing)
    X_R_cumu_U = initialize_param("x_R_U", model_data.index_y, ipp.index_k_existing)
    for y in model_data.index_y
        for k in ipp.index_k_existing
            # need to have constraints in the optimization as well
            X_R_cumu_U(y, k) = ipp.x_E_my(p_star, k) - ipp.x_R_cumu(p_star, k)
        end
    end

    lambda_L = initialize_param("lambda_L", model_data.index_y, ipp.index_k_new, model_data.index_t)
    lambda_U = initialize_param("lambda_U", model_data.index_y, ipp.index_k_new, model_data.index_t)

    for y in model_data.index_y
        for k in ipp.index_k_new
            for t in model_data.index_t
                lambda_U(y, k, t) = lambda_param(y, k, t) * lambda_upper_bound_adj
                lambda_L(y, k, t) = lambda_param(y, k, t) * lambda_lower_bound_adj
            end
        end
    end

    # for y in model_data.index_y
    #     for k in ipp.index_k_new
    #         for t in model_data.index_t
    #             lambda_U[y, k, t] = lambda_param[(lambda_param.Year .== model_data.year[y]) .& (lambda_param.IPP .== "ipp1") .& (lambda_param.GenTech .== string(k)) .& (lambda_param.Time .== string(t)), "Generation_MWh"][1] + 10.0
    #             lambda_L[y, k, t] = lambda_param[(lambda_param.Year .== model_data.year[y]) .& (lambda_param.IPP .== "ipp1") .& (lambda_param.GenTech .== string(k)) .& (lambda_param.Time .== string(t)), "Generation_MWh"][1] - 10.0
    #         end
    #     end
    # end

    X_C_cumu_L = initialize_param("x_C_L", model_data.index_y, ipp.index_k_new)
    X_C_cumu_U = initialize_param("x_C_U", model_data.index_y, ipp.index_k_new)
    for y in model_data.index_y
        for k in ipp.index_k_new
            # need to have constraints in the optimization as well
            # this hard-coded number needs to change
            X_C_cumu_U(y, k) = 5000.0
        end
    end


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
        mce_E_p_star[model_data.index_y, ipp.index_k_existing, model_data.index_t],
    )
    @variable(
        WMDER_IPP,
        mce_C_p_star[model_data.index_y, ipp.index_k_new, model_data.index_t],
    )
    # cumulative retirement of p_star
    @variable(WMDER_IPP, x_R_mce[model_data.index_y, ipp.index_k_existing] >= 0)
    # cumulative investment of p_star
    @variable(WMDER_IPP, 0 <= x_C_mce[model_data.index_y, ipp.index_k_new] <= 5000.0)

    @variable(WMDER_IPP, UCAP_p_star[model_data.index_y])

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

    Max_Net_Load_my_index = KeyedArray(
        [
            findmax(Dict(t => ipp.Net_Load_my(y, t) for t in model_data.index_t))[2] for
            y in model_data.index_y
        ];
        [get_pair(model_data.index_y)]...
    )


    fill!(ipp.capacity_credit_E_my, NaN)
    for y in model_data.index_y, k in ipp.index_k_existing
        ipp.capacity_credit_E_my(y, k, :) .= ipp.rho_E_my(p_star, k, Max_Net_Load_my_index(y))
    end
    fill!(ipp.capacity_credit_C_my, NaN)
    for y in model_data.index_y, k in ipp.index_k_new
        ipp.capacity_credit_C_my(y, k, :) .= ipp.rho_C_my(p_star, k, Max_Net_Load_my_index(y))
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
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
                ) +
                sum(
                    ipp.capacity_credit_C_my(y, k) * (
                        sum(
                            ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p, k)
                    ) for k in ipp.index_k_new,
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
                ) +
                # green technology subscription
                sum(
                    ipp.capacity_credit_C_my(y, j) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
                    model_data.year(first(model_data.index_y_fix)):model_data.year(y))
                    for j in model_data.index_j, h in model_data.index_h
                )
            end
    else
        UCAP_total = y -> begin
            UCAP_p_star(y) +
            # green technology subscription
            sum(
                ipp.capacity_credit_C_my(y, j) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
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
                    ipp.pvf_onm(Symbol(Int(y_symbol)), p_star) for y_symbol in
                    model_data.year(y):model_data.year(last(model_data.index_y.elements))
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
                        ipp.rho_E_my(p, k, t) *
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
                                ipp.x_R_my(Symbol(Int(y_symbol)), p_star, k) for
                                y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) - ipp.x_R_cumu(p_star, k)
                        ) for t in model_data.index_t, k in ipp.index_k_existing
                    )
                ) - (
                    sum(
                        lambda[y, p, k, t] *
                        ipp.rho_C_my(p, k, t) *
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
                                ipp.x_C_my(Symbol(Int(y_symbol)), p_star, k) for
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
                        (ipp.v_E_my(y, p_star, k, t) * y_E[y, p_star, k, t]) for
                        t in model_data.index_t, k in ipp.index_k_existing
                    )
                ) - (
                    sum(
                        model_data.omega(t) * (ipp.v_C_my(y, p, k, t) * y_C[y, p, k, t]) for
                        t in model_data.index_t, k in ipp.index_k_new, p in ipp.index_p
                    ) - sum(
                        model_data.omega(t) *
                        (ipp.v_C_my(y, p_star, k, t) * y_C[y, p_star, k, t]) for
                        t in model_data.index_t, k in ipp.index_k_new
                    )
                ) -

                # generation costs (this includes terms due to revenue linearization)
                #   num hrs * ((fuel + vom) * gen existing + (fuel + vom) * gen new) for t and gen type
                sum(
                    model_data.omega(t) *
                    (ipp.v_E_my(y, p_star, k, t) * y_E[y, p_star, k, t]) for
                    t in model_data.index_t, k in ipp.index_k_existing
                ) - sum(
                    model_data.omega(t) *
                    (ipp.v_C_my(y, p_star, k, t) * y_C[y, p_star, k, t]) for
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

    # dual feasible constraints
    @constraint(
        WMDER_IPP,
        Eq_dual_feasible_E[
            y in model_data.index_y,
            p in ipp.index_p,
            k in ipp.index_k_existing,
            t in model_data.index_t,
        ],
        model_data.omega(t) * ipp.v_E_my(y, p, k, t) - miu[y, t] + eta[y, p, k, t] >= 0
    )
    @constraint(
        WMDER_IPP,
        Eq_dual_feasible_C[
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
        Eq_primal_feasible_supplydemandbalance[y in model_data.index_y, t in model_data.index_t],
        supply_demand_balance(y, t) >= 0
    )
    
    @constraint(
        WMDER_IPP,
        Eq_primal_feasible_gen_max_E[
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
            Eq_primal_feasible_other_gen_max_E[
                y in model_data.index_y,
                p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p))),
                k in ipp.index_k_existing,
                t in model_data.index_t,
            ],
            ipp.rho_E_my(p, k, t) * (
                ipp.x_E_my(p, k) - sum(
                    ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) - ipp.x_R_cumu(p, k)
            ) - y_E[y, p, k, t] >= 0
        )
    end

    @constraint(
        WMDER_IPP,
        Eq_primal_feasible_gen_max_C[
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
            Eq_primal_feasible_other_gen_max_C[
                y in model_data.index_y,
                p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p))),
                k in ipp.index_k_new,
                t in model_data.index_t,
            ],
            ipp.rho_C_my(p, k, t) * (
                sum(
                    ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) + ipp.x_C_cumu(p, k)
            ) - y_C[y, p, k, t] >= 0
        )
    end

    # strong-duality constraints
    if length(ipp.index_p) >= 2
        @constraint(
            WMDER_IPP,
            Eq_strong_duality[
                y in model_data.index_y
            ],
            sum(
                model_data.omega(t) *
                (ipp.v_E_my(y, p, k, t) * y_E[y, p, k, t]) for
                t in model_data.index_t, k in ipp.index_k_existing, p in ipp.index_p
            ) + 
            sum(
                model_data.omega(t) *
                (ipp.v_C_my(y, p, k, t) * y_C[y, p, k, t]) for
                t in model_data.index_t, k in ipp.index_k_new, p in ipp.index_p
            ) == 
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
                            model_data.year(first(model_data.index_y_fix)):model_data.year[y]
                        ) for h in model_data.index_h, m in customers.index_m
                    ) -
                    # green technology subscription at time t
                    sum(
                        utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
                        model_data.year(first(model_data.index_y_fix)):model_data.year(y))
                        for j in model_data.index_j, h in model_data.index_h
                    )
                ) for t in model_data.index_t
            ) - 
            (
                sum(
                    eta[y, p, k, t] *
                    ipp.rho_E_my(p, k, t) *
                    (
                        ipp.x_E_my(p, k) - ipp.x_R_cumu(p, k)
                    ) for t in model_data.index_t, k in ipp.index_k_existing, p in ipp.index_p
                ) + sum(
                    lambda[y, p, k, t] *
                    ipp.rho_C_my(p, k, t) *
                    ipp.x_C_cumu(p, k) 
                    for t in model_data.index_t, k in ipp.index_k_new, p in ipp.index_p
                )
            ) + sum(
                eta[y, p, k, t] *
                ipp.rho_E_my(p, k, t) *
                sum(
                    ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) for t in model_data.index_t, k in ipp.index_k_existing,
                p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
            ) - sum(
                lambda[y, p, k, t] *
                ipp.rho_C_my(p, k, t) *
                sum(
                    ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) for t in model_data.index_t, k in ipp.index_k_new,
                p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
            ) + sum(
                ipp.rho_E_my(p_star, k, t) *
                mce_E_p_star[y, k, t] for t in model_data.index_t, k in ipp.index_k_existing
            ) - sum(
                ipp.rho_C_my(p_star, k, t) *
                mce_C_p_star[y, k, t] for t in model_data.index_t, k in ipp.index_k_new
            )
        )
    else
        @constraint(
            WMDER_IPP,
            Eq_strong_duality[
                y in model_data.index_y
            ],
            sum(
                model_data.omega(t) *
                (ipp.v_E_my(y, p, k, t) * y_E[y, p, k, t]) for
                t in model_data.index_t, k in ipp.index_k_existing, p in ipp.index_p
            ) + 
            sum(
                model_data.omega(t) *
                (ipp.v_C_my(y, p, k, t) * y_C[y, p, k, t]) for
                t in model_data.index_t, k in ipp.index_k_new, p in ipp.index_p
            ) == 
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
            ) - 
            (
                sum(
                    eta[y, p, k, t] *
                    ipp.rho_E_my(p, k, t) *
                    (
                        ipp.x_E_my(p, k) - ipp.x_R_cumu(p, k)
                    ) for t in model_data.index_t, k in ipp.index_k_existing, p in ipp.index_p
                ) + sum(
                    lambda[y, p, k, t] *
                    ipp.rho_C_my(p, k, t) *
                    ipp.x_C_cumu(p, k) 
                    for t in model_data.index_t, k in ipp.index_k_new, p in ipp.index_p
                )
            ) + sum(
                ipp.rho_E_my(p_star, k, t) *
                mce_E_p_star[y, k, t] for t in model_data.index_t, k in ipp.index_k_existing
            ) - sum(
                ipp.rho_C_my(p_star, k, t) *
                mce_C_p_star[y, k, t] for t in model_data.index_t, k in ipp.index_k_new
            )
        )
    end

    # linking constraints on total retirement of p_star at year y
    @constraint(
        WMDER_IPP,
        Eq_cumu_R_mce[
            y in model_data.index_y,
            k in ipp.index_k_existing
        ],
        x_R_mce[y, k] == sum(
            x_R[Symbol(Int(y_symbol)), k] for
            y_symbol in model_data.year(first(model_data.index_y)):model_data.year(y)
        )
    )
    # McCormick-envelope relaxation of mce_E_p_star
    @constraint(
        WMDER_IPP,
        Eq_MCE_Relax_LL_R[
            y in model_data.index_y,
            k in ipp.index_k_existing,
            t in model_data.index_t
        ],
        mce_E_p_star[y, k, t] >= X_R_cumu_L(y, k) * eta[y, p_star, k, t] + 
        eta_L(y, k, t) * x_R_mce[y, k] - eta_L(y, k, t) * X_R_cumu_L(y, k)
    )
    @constraint(
        WMDER_IPP,
        Eq_MCE_Relax_UU_R[
            y in model_data.index_y,
            k in ipp.index_k_existing,
            t in model_data.index_t
        ],
        mce_E_p_star[y, k, t] >= - eta_U(y, k, t) * X_R_cumu_U(y, k) +
        eta_U(y, k, t) * x_R_mce[y, k] + eta[y, p_star, k, t] * X_R_cumu_U(y, k)
    )
    @constraint(
        WMDER_IPP,
        Eq_MCE_Relax_LU_R[
            y in model_data.index_y,
            k in ipp.index_k_existing,
            t in model_data.index_t
        ],
        mce_E_p_star[y, k, t] <= X_R_cumu_U(y, k) * eta[y, p_star, k, t] -
        eta_L(y, k, t) * X_R_cumu_U(y, k) + eta_L(y, k, t) * x_R_mce[y, k]
    )
    @constraint(
        WMDER_IPP,
        Eq_MCE_Relax_UL_R[
            y in model_data.index_y,
            k in ipp.index_k_existing,
            t in model_data.index_t
        ],
        mce_E_p_star[y, k, t] <= eta_U(y, k, t) * x_R_mce[y, k] -
        eta_U(y, k, t) * X_R_cumu_L(y, k) + X_R_cumu_L(y, k) * eta[y, p_star, k, t]
    )


    # linking constraints on total investment of p_star at year y
    @constraint(
        WMDER_IPP,
        Eq_cumu_C_mce[
            y in model_data.index_y,
            k in ipp.index_k_new,
        ],
        x_C_mce[y, k] == sum(
            x_C[Symbol(Int(y_symbol)), k] for
            y_symbol in model_data.year(first(model_data.index_y)):model_data.year(y)
        )
    )

    # McCormick-envelope relaxation of mce_C_p_star

    @constraint(
        WMDER_IPP,
        Eq_MCE_Relax_LL_C[
            y in model_data.index_y,
            k in ipp.index_k_new,
            t in model_data.index_t
        ],
        mce_C_p_star[y, k, t] >= X_C_cumu_L(y, k) * lambda[y, p_star, k, t] + 
        lambda_L(y, k, t) * x_C_mce[y, k] - lambda_L(y, k, t) * X_C_cumu_L(y, k)
    )
    @constraint(
        WMDER_IPP,
        Eq_MCE_Relax_UU_C[
            y in model_data.index_y,
            k in ipp.index_k_new,
            t in model_data.index_t
        ],
        mce_C_p_star[y, k, t] >= - lambda_U(y, k, t) * X_C_cumu_U(y, k) +
        lambda_U(y, k, t) * x_C_mce[y, k] + lambda[y, p_star, k, t] * X_C_cumu_U(y, k)
    )
    @constraint(
        WMDER_IPP,
        Eq_MCE_Relax_LU_C[
            y in model_data.index_y,
            k in ipp.index_k_new,
            t in model_data.index_t
        ],
        mce_C_p_star[y, k, t] <= X_C_cumu_U(y, k) * lambda[y, p_star, k, t] -
        lambda_L(y, k, t) * X_C_cumu_U(y, k) + lambda_L(y, k, t) * x_C_mce[y, k]
    )
    @constraint(
        WMDER_IPP,
        Eq_MCE_Relax_UL_C[
            y in model_data.index_y,
            k in ipp.index_k_new,
            t in model_data.index_t
        ],
        mce_C_p_star[y, k, t] <= lambda_U(y, k, t) * x_C_mce[y, k] -
        lambda_U(y, k, t) * X_C_cumu_L(y, k) + X_C_cumu_L(y, k) * lambda[y, p_star, k, t]
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
                    ipp.rho_E_my(p, k, t) * (
                        ipp.x_E_my(p, k) - sum(
                            ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_R_cumu(p, k)
                    ) for k in ipp.index_k_existing,
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
                ) +
                sum(
                    ipp.rho_C_my(p, k, t) * (
                        sum(
                            ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p, k)
                    ) for k in ipp.index_k_new,
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
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
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
                ) +
                sum(
                    ipp.capacity_credit_C_my(y, k) * (
                        sum(
                            ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p, k)
                    ) for k in ipp.index_k_new,
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
                ) +
                # green technology subscription
                sum(
                    ipp.capacity_credit_C_my(y, j) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
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
                    ipp.capacity_credit_C_my(y, j) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
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
    #     sum(ipp.capacity_credit_E_my(y,k)*(ipp.x_E_my(p_star,k)-sum(ipp.x_R_my(Symbol(y_symbol),p_star,k) for y_symbol = model_data.year(first(model_data.index_y)):model_data.year(y)) - ipp.x_R_cumu[p_star,k]) for k in ipp.index_k_existing) + 
    #     sum(ipp.capacity_credit_C_my(y,k)*(sum(ipp.x_C_my(Symbol(y_symbol),p_star,k) for y_symbol = model_data.year(first(model_data.index_y)):model_data.year(y)) + ipp.x_C_cumu[p_star,k]) for k in ipp.index_k_new)
    #     for y in model_data.index_y
    # )
    # capacity_profit = Dict(y =>
    #     ipp.pvf_onm(y,p_star) * UCAP_p_star_solution(y) * (ipp.Capacity_intercept_my(y) + ipp.Capacity_slope_my(y) * UCAP_p_star_solution(y))
    #     for y in model_data.index_y
    # )
    # total_profit = Dict(y => capacity_profit(y) + objective_value(WMDER_IPP) for y in model_data.index_y)

    ###### running MIQP ######
    UCAP_p_star = KeyedArray(
        [
            sum(
                ipp.capacity_credit_E_my(y, k) * (
                    ipp.x_E_my(p_star, k) - sum(
                        ipp.x_R_my(Symbol(Int(y_symbol)), p_star, k) for y_symbol in
                        model_data.year(first(model_data.index_y)):model_data.year(y)
                    ) - ipp.x_R_cumu(p_star, k)
                ) for k in ipp.index_k_existing
            ) + sum(
                ipp.capacity_credit_C_my(y, k) * (
                    sum(
                        ipp.x_C_my(Symbol(Int(y_symbol)), p_star, k) for y_symbol in
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
                UCAP_p_star(y) +
                sum(
                    ipp.capacity_credit_E_my(y, k) * (
                        ipp.x_E_my(p, k) - sum(
                            ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_R_cumu(p, k)
                    ) for k in ipp.index_k_existing,
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
                ) +
                sum(
                    ipp.capacity_credit_C_my(y, k) * (
                        sum(
                            ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p, k)
                    ) for k in ipp.index_k_new,
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
                ) +
                # green technology subscription
                sum(
                    ipp.capacity_credit_C_my(y, j) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
                    model_data.year(first(model_data.index_y_fix)):model_data.year(y))
                    for j in model_data.index_j, h in model_data.index_h
                ) 
                for y in model_data.index_y
            ];
            [get_pair(model_data.index_y)]...
        )
    else
        UCAP_total = KeyedArray(
            [UCAP_p_star(y) +
            # green technology subscription
            sum(
                ipp.capacity_credit_C_my(y, j) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):model_data.year(y))
                for j in model_data.index_j, h in model_data.index_h
            ) for y in model_data.index_y];
            [get_pair(model_data.index_y)]...
        )
    end

    for y in model_data.index_y
        ipp.capacity_price(y, :) .=
            ipp.Capacity_intercept_my(y) + ipp.Capacity_slope_my(y) * UCAP_total(y)
        ipp.ucap(y, p_star, :) .= UCAP_p_star(y)
    end


    # report lower level duality gap
    lower_level_primal_obj = DataFrame()
    lower_level_dual_obj = DataFrame()
    lower_level_duality_gap = DataFrame()
    for y in model_data.index_y
        lower_level_primal_obj[!, y] = Vector{Float64}(undef, 1)
        lower_level_dual_obj[!, y] = Vector{Float64}(undef, 1)
        lower_level_duality_gap[!, y] = Vector{Float64}(undef, 1)
    end
    for y in model_data.index_y
        lower_level_primal_obj[1, y] = 
        sum(
            model_data.omega(t) *
            (ipp.v_E_my(y, p, k, t) * value.(y_E[y, p, k, t])) for
            t in model_data.index_t, k in ipp.index_k_existing, p in ipp.index_p
        ) + 
        sum(
            model_data.omega(t) *
            (ipp.v_C_my(y, p, k, t) * value.(y_C[y, p, k, t])) for
            t in model_data.index_t, k in ipp.index_k_new, p in ipp.index_p
        )

        if length(ipp.index_p) >= 2
            lower_level_dual_obj[1, y] = 
            sum(
                value.(miu[y, t]) * (
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
            ) - 
            (
                sum(
                    value.(eta[y, p, k, t]) *
                    ipp.rho_E_my(p, k, t) *
                    (
                        ipp.x_E_my(p, k) - ipp.x_R_cumu(p, k)
                    ) for t in model_data.index_t, k in ipp.index_k_existing, p in ipp.index_p
                ) + sum(
                    value.(lambda[y, p, k, t]) *
                    ipp.rho_C_my(p, k, t) *
                    ipp.x_C_cumu(p, k) 
                    for t in model_data.index_t, k in ipp.index_k_new, p in ipp.index_p
                )
            ) + sum(
                value.(eta[y, p, k, t]) *
                ipp.rho_E_my(p, k, t) *
                sum(
                    ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) for t in model_data.index_t, k in ipp.index_k_existing,
                p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
            ) - sum(
                value.(lambda[y, p, k, t]) *
                ipp.rho_C_my(p, k, t) *
                sum(
                    ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) for t in model_data.index_t, k in ipp.index_k_new,
                p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
            ) + sum(
                ipp.rho_E_my(p_star, k, t) *
                value.(eta[y, p_star, k, t]) *
                value.(x_R_mce[y, k]) for t in model_data.index_t, k in ipp.index_k_existing
            ) - sum(
                ipp.rho_C_my(p_star, k, t) *
                value.(lambda[y, p_star, k, t]) *
                value.(x_C_mce[y, k]) for t in model_data.index_t, k in ipp.index_k_new
            )
        else
            lower_level_dual_obj[1, y] = 
            sum(
                value.(miu[y, t]) * (
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
            ) - 
            (
                sum(
                    value.(eta[y, p, k, t]) *
                    ipp.rho_E_my(p, k, t) *
                    (
                        ipp.x_E_my(p, k) - ipp.x_R_cumu(p, k)
                    ) for t in model_data.index_t, k in ipp.index_k_existing, p in ipp.index_p
                ) + sum(
                    value.(lambda[y, p, k, t]) *
                    ipp.rho_C_my(p, k, t) *
                    ipp.x_C_cumu(p, k) 
                    for t in model_data.index_t, k in ipp.index_k_new, p in ipp.index_p
                )
            ) + sum(
                ipp.rho_E_my(p_star, k, t) *
                value.(eta[y, p_star, k, t]) *
                value.(x_R_mce[y, k]) for t in model_data.index_t, k in ipp.index_k_existing
            ) - sum(
                ipp.rho_C_my(p_star, k, t) *
                value.(lambda[y, p_star, k, t]) *
                value.(x_C_mce[y, k]) for t in model_data.index_t, k in ipp.index_k_new
            )
        end
    end
    lower_level_duality_gap = (lower_level_primal_obj .- lower_level_dual_obj) ./ abs.(lower_level_primal_obj)
    @info "lower level primal obj is $(lower_level_primal_obj)"
    @info "lower level dual obj is $(lower_level_dual_obj)"
    @info "lower level duality gap is $(lower_level_duality_gap)"

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
)
    # Primal Variables
    save_param(
        ipps.y_E_my.values,
        [:Year, :IPP, :GenTech, :Time],
        :Generation_MWh,
        joinpath(export_file_path, "y_E.csv"),
    )
    save_param(
        ipps.y_C_my.values,
        [:Year, :IPP, :GenTech, :Time],
        :Generation_MWh,
        joinpath(export_file_path, "y_C.csv"),
    )
    save_param(
        ipps.x_R_my.values,
        [:Year, :IPP, :GenTech],
        :Capacity_MW,
        joinpath(export_file_path, "x_R.csv"),
    )
    save_param(
        ipps.x_C_my.values,
        [:Year, :IPP, :GenTech],
        :Capacity_MW,
        joinpath(export_file_path, "x_C.csv"),
    )
    save_param(
        ipps.LMP_my.values,
        [:Year, :Time],
        :MarginalCost,
        joinpath(export_file_path, "LMP.csv"),
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
                    sum(ipp.y_E_my(y, p, rps, t) for rps in ipp.index_rps) +
                    sum(ipp.y_C_my(y, p, rps, t) for rps in ipp.index_rps)
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
                ipp.fom_E_my(y, p, k) * (
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
        ADITNew(y, p, k) = sum(
            ipp.CapEx_my(Symbol(Int(y_symbol)), p, k) *
            ipp.x_C_my(Symbol(Int(y_symbol)), p, k) *
            (
                utility.CumuTaxDepre_new_my(
                    Symbol(Int(model_data.year(y) - y_symbol + 1)),
                    k,
                 ) - utility.CumuAccoutDepre_new_my(
                    Symbol(Int(model_data.year(y) - y_symbol + 1)),
                    k,
                 )
            ) *
            ipp.Tax(p) +
            utility.ITC_new_my(Symbol(Int(y_symbol)), k) *
            ipp.CapEx_my(Symbol(Int(y_symbol)), p, k) *
            ipp.x_C_my(Symbol(Int(y_symbol)), p, k) *
            (
                1 - utility.CumuITCAmort_new_my(
                    Symbol(Int(model_data.year(y) - y_symbol + 1)),
                    k,
                )
            ) for y_symbol in
            model_data.year(first(model_data.index_y_fix)):model_data.year(y)
        )
    end
    RateBaseNoWC_new = make_keyed_array(model_data.index_y_fix, ipp.index_p, ipp.index_k_new)
    for y in model_data.index_y_fix, p in ipp.index_p, k in ipp.index_k_new
        RateBaseNoWC_new(y, p, k) =
            sum(
                ipp.CapEx_my(Symbol(Int(y_symbol)), p, k) *
                ipp.x_C_my(Symbol(Int(y_symbol)), p, k) *
                (
                    1 - utility.CumuAccoutDepre_new_my(
                        Symbol(Int(model_data.year(y) - y_symbol + 1)),
                        k,
                    )
                ) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):model_data.year(y)
            ) - ADITNew(y, p, k)
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
            sum(RateBaseNoWC_new(y, p, k) for k in ipp.index_k_new) +
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
                utility.CapEx_existing_my(k) *
                (
                    ipp.x_E_my(p, k) - sum(
                        ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                        model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                    )
                ) *
                utility.AnnualAccoutDepre_existing_my(y, k) +
                utility.CapEx_existing_my(k) *
                ipp.x_R_my(y, p, k) *
                (
                    utility.AnnualAccoutDepre_existing_my(y, k) + 1 -
                    utility.CumuAccoutDepre_existing_my(y, k)
                ) for k in ipp.index_k_existing
            ) + sum(
                ipp.CapEx_my(Symbol(Int(y_symbol)), p, k) *
                ipp.x_C_my(Symbol(Int(y_symbol)), p, k) *
                utility.AnnualAccoutDepre_new_my(
                    Symbol(Int(model_data.year(y) - y_symbol + 1)),
                    k,
                 ) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):model_data.year(y),
                k in ipp.index_k_new
            )
    end

    depreciation_tax = make_keyed_array(model_data.index_y_fix, ipp.index_p)
    for y in model_data.index_y_fix, p in ipp.index_p
        depreciation_tax(y, p, :) .=
            sum(
                utility.CapEx_existing_my(k) *
                (
                    ipp.x_E_my(p, k) - sum(
                        ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                        model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                    )
                ) *
                utility.AnnualTaxDepre_existing_my(y, k) +
                utility.CapEx_existing_my(k) *
                ipp.x_R_my(y, p, k) *
                (
                    utility.AnnualTaxDepre_existing_my(y, k) + 1 -
                    utility.CumuTaxDepre_existing_my(y, k)
                ) for k in ipp.index_k_existing
            ) + sum(
                ipp.CapEx_my(Symbol(Int(y_symbol)), p, k) *
                ipp.x_C_my(Symbol(Int(y_symbol)), p, k) *
                utility.AnnualTaxDepre_new_my(
                    Symbol(Int(model_data.year(y) - y_symbol + 1)),
                    k,
                 ) for y_symbol in
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
                utility.ITC_existing_my(k) *
                utility.CapEx_existing_my(k) *
                (
                    ipp.x_E_my(p, k) - sum(
                        ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
                        model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                    )
                ) *
                utility.AnnualITCAmort_existing_my(y, k) +
                # existing units that are retired this year will incur their regular annual depreciation, as well as the remaining un-depreciated asset
                utility.ITC_existing_my(k) *
                utility.CapEx_existing_my(k) *
                ipp.x_R_my(y, p, k) *
                (
                    utility.AnnualITCAmort_existing_my(y, k) + 1 -
                    utility.CumuITCAmort_existing_my(y, k)
                ) for k in ipp.index_k_existing
            ) -
            sum(
                utility.ITC_new_my(Symbol(Int(y_symbol)), k) *
                ipp.CapEx_my(Symbol(Int(y_symbol)), p, k) *
                ipp.x_C_my(Symbol(Int(y_symbol)), p, k) *
                utility.AnnualITCAmort_new_my(Symbol(Int(model_data.year(y) - y_symbol + 1)), k) for
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
                        ipp.y_E_my(y, p, k, t) * ipp.emission_rate_E_my(y, p, k) for
                        k in utility.index_k_existing, p in ipp.index_p
                    ) + sum(
                        ipp.y_C_my(y, p, k, t) * ipp.emission_rate_C_my(y, p, k) for
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
