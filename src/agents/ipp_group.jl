# This module defines the data and functions associated with the Independent Power Producer

# declare customer decision
abstract type IPPAlgorithm end
struct LagrangeDecomposition <: IPPAlgorithm end
struct MIQP <: IPPAlgorithm end
struct MPPDCMER <: IPPAlgorithm end
struct MPPDCMERTransStorage <: IPPAlgorithm end

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
    index_stor_existing::Dimension # existing bulk storage technologies
    index_stor_new::Dimension # potential bulk storage technologies
    index_p::Dimension # individual independent power producers
    index_rps::Dimension # RPS-qualified technologies
    index_l::Dimension # transmission lines

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
    x_stor_E_my::ParamArray # existing storage capacity (MW)
    fom_E_my::ParamArray # fixed O&M of existing capacity ($/MW-yr)
    fom_C_my::ParamArray # fixed O&M of new capacity ($/MW-yr)
    fom_stor_E_my::ParamArray # fixed O&M of existing storage capacity ($/MW-yr)
    fom_stor_C_my::ParamArray # fixed O&M of new storage capacity ($/MW-yr)
    CapEx_my::ParamArray # capital expense of new capacity ($/MW)
    CapEx_stor_my::ParamArray # capital expense of new storage capacity ($/MW)
    rte_stor_E_my::ParamArray # round trip efficiency of existing storage ($/MW)
    rte_stor_C_my::ParamArray # round trip efficiency of new storage ($/MW)
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
    initial_energy_existing_my::ParamArray
    initial_energy_new_my::ParamArray
    stor_duration_existing::ParamArray
    stor_duration_new::ParamArray
    trans_topology::ParamArray
    trans_capacity::ParamArray

    # Primal Variables (multi-year)
    y_E_my::ParamArray
    y_C_my::ParamArray
    x_R_my::ParamArray
    x_C_my::ParamArray
    x_stor_R_my::ParamArray
    x_stor_C_my::ParamArray
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
    x_stor_R_cumu::ParamArray
    x_stor_C_cumu::ParamArray

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
    capacity_credit_stor_E_my::ParamArray # capacity credit of existing resources
    capacity_credit_stor_C_my::ParamArray # capacity credit of new resources
    capacity_price::ParamArray # $/MW-yr
    capacity_price_my_temp::ParamArray

    Net_Load_my::ParamArray
    Max_Net_Load_my::ParamArray
    Reserve_req_my::ParamArray
    Capacity_slope_my::ParamArray
    Capacity_intercept_my::ParamArray
    ucap_temp::ParamArray
    ucap::ParamArray
    ucap_total::ParamArray

    # RPS
    RPS::ParamArray

    # emission rate
    emission_rate_E_my::ParamArray
    emission_rate_C_my::ParamArray

    # storage and transmission related dispatch
    charge_E_my::ParamArray
    discharge_E_my::ParamArray
    charge_C_my::ParamArray
    discharge_C_my::ParamArray
    energy_E_my::ParamArray
    energy_C_my::ParamArray
    flow_my::ParamArray

    Max_Net_Load_my_dict::Dict

    # McCormic bounds
    eta_L_vec::Vector{}
    eta_U_vec::Vector{}
    lambda_L_vec::Vector{}
    lambda_U_vec::Vector{}
    theta_E_energy_L_vec::Vector{}
    theta_E_energy_U_vec::Vector{}
    theta_E_discharge_L_vec::Vector{}
    theta_E_discharge_U_vec::Vector{}
    theta_E_charge_L_vec::Vector{}
    theta_E_charge_U_vec::Vector{}
    pi_E_charge_L_vec::Vector{}
    pi_E_charge_U_vec::Vector{}
    kappa_E_L_vec::Vector{}
    kappa_E_U_vec::Vector{}
    theta_C_energy_L_vec::Vector{}
    theta_C_energy_U_vec::Vector{}
    theta_C_discharge_L_vec::Vector{}
    theta_C_discharge_U_vec::Vector{}
    theta_C_charge_L_vec::Vector{}
    theta_C_charge_U_vec::Vector{}
    pi_C_charge_L_vec::Vector{}
    pi_C_charge_U_vec::Vector{}
    kappa_C_L_vec::Vector{}
    kappa_C_U_vec::Vector{}
    
end

function IPPGroup(input_filename::String, model_data::HEMData, id = DEFAULT_ID)
    index_k_existing = read_set(input_filename, "index_k_existing", "index_k_existing")
    index_k_new = read_set(input_filename, "index_k_new", "index_k_new")
    index_stor_existing = read_set(input_filename, "index_stor_existing", "index_stor_existing")
    index_stor_new = read_set(input_filename, "index_stor_new", "index_stor_new")
    index_p = read_set(input_filename, "index_p", "index_p")
    index_rps = read_set(input_filename, "index_rps", "index_rps")
    index_l = read_set(input_filename, "index_l", "index_l")

    min_max = Dimension(
        "min_max",
        Symbol.(["min", "max"]),
        prose_name = "min_max",
        description = "minimum and maximum capacity of transmission lines",
    )

    FOMNew = read_param("FOM_new", input_filename, "FOMNewIPP", index_k_new, [index_p, model_data.index_z])
    CapExNew =
        read_param("CapEx_new", input_filename, "CapExNewIPP", index_k_new, [index_p, model_data.index_z])
    LifetimeNew =
        read_param("Lifetime_new", input_filename, "LifetimeNewIPP", index_k_new, [index_p])
    LifetimeStorNew =
        read_param("LifetimeStor_new", input_filename, "LifetimeStorNewIPP", index_stor_new, [index_p])
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
    FixedCostNew = make_keyed_array(index_p, model_data.index_z, index_k_new)
    for p in index_p, z in model_data.index_z, k in index_k_new
        FixedCostNew(p, z, k, :) .= FOMNew(p, z, k) + CapExNew(p, z, k) * CRF[p, k]
    end

    eximport = read_param("eximport", input_filename, "Export", model_data.index_t, [model_data.index_z, model_data.index_d])
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
        index_stor_existing,
        index_stor_new,
        index_p,
        index_rps,
        index_l,
        read_param(
            "x_E",
            input_filename,
            "ExistingCapacityIPP",
            index_k_existing,
            [index_p, model_data.index_z],
        ),
        read_param("f_E", input_filename, "FixedCostOldIPP", index_k_existing, [index_p, model_data.index_z]),
        ParamArray("f_C", Tuple(push!(copy([index_p, model_data.index_z]), index_k_new)), FixedCostNew),
        read_param(
            "v_E",
            input_filename,
            "VariableCostOldIPP",
            model_data.index_t,
            [index_p, index_k_existing, model_data.index_z, model_data.index_d],
        ),
        read_param(
            "v_C",
            input_filename,
            "VariableCostNewIPP",
            model_data.index_t,
            [index_p, index_k_new, model_data.index_z, model_data.index_d],
        ),
        read_param(
            "rho_E",
            input_filename,
            "AvailabilityOldIPP",
            model_data.index_t,
            [index_p, index_k_existing, model_data.index_z, model_data.index_d],
        ),
        read_param(
            "rho_C",
            input_filename,
            "AvailabilityNewIPP",
            model_data.index_t,
            [index_p, index_k_new, model_data.index_z, model_data.index_d],
        ),
        eximport,
        peak_eximport,
        ParamScalar("B1GM", 1000000),
        ParamScalar("zeta", 3.0),
        initialize_param("y_E", index_p, index_k_existing, model_data.index_z, model_data.index_d, model_data.index_t),
        initialize_param("y_C", index_p, index_k_new, model_data.index_z, model_data.index_d, model_data.index_t),
        initialize_param("x_R", index_p, index_k_existing, model_data.index_z),
        initialize_param("x_C", index_p, index_k_new, model_data.index_z),
        initialize_param("miu", model_data.index_z, model_data.index_d, model_data.index_t),
        read_param(
            "o_E",
            input_filename,
            "VariableCostOldIPP",
            model_data.index_t,
            [index_p, index_k_existing, model_data.index_z, model_data.index_d],
        ),
        read_param(
            "o_C",
            input_filename,
            "VariableCostNewIPP",
            model_data.index_t,
            [index_p, index_k_new, model_data.index_z, model_data.index_d],
        ),
        initialize_param("LMP", model_data.index_z, model_data.index_d, model_data.index_t),
        read_param(
            "x_E_my",
            input_filename,
            "ExistingCapacityIPP",
            index_k_existing,
            [index_p, model_data.index_z],
        ),
        read_param(
            "x_stor_E_my",
            input_filename,
            "ExistingStorCapacityIPP",
            index_stor_existing,
            [index_p, model_data.index_z],
        ),
        read_param(
            "fom_E_my",
            input_filename,
            "FixedCostOldIPPmy",
            index_k_existing,
            [model_data.index_y, index_p, model_data.index_z],
        ),
        read_param(
            "fom_C_my",
            input_filename,
            "FOMNewIPPmy",
            index_k_new,
            [model_data.index_y, index_p, model_data.index_z],
        ),
        read_param(
            "fom_stor_E_my",
            input_filename,
            "FixedCostStorOldIPPmy",
            index_stor_existing,
            [model_data.index_y, index_p, model_data.index_z],
        ),
        read_param(
            "fom_stor_C_my",
            input_filename,
            "StorFOMNewIPPmy",
            index_stor_new,
            [model_data.index_y, index_p, model_data.index_z],
        ),
        read_param(
            "CapEx_my",
            input_filename,
            "CapExNewIPPmy",
            index_k_new,
            [model_data.index_y, index_p, model_data.index_z],
        ),
        read_param(
            "CapEx_stor_my",
            input_filename,
            "StorCapExNewIPPmy",
            index_stor_new,
            [model_data.index_y, index_p, model_data.index_z],
        ),
        read_param(
            "rte_stor_E_my",
            input_filename,
            "StorRTEOldIPPmy",
            index_stor_existing,
            [model_data.index_y, index_p, model_data.index_z],
        ),
        read_param(
            "rte_stor_C_my",
            input_filename,
            "StorRTENewIPPmy",
            index_stor_new,
            [model_data.index_y, index_p, model_data.index_z],
        ),
        read_param(
            "v_E_my",
            input_filename,
            "VariableCostOldIPPmy",
            model_data.index_t,
            [model_data.index_y, index_p, index_k_existing, model_data.index_z, model_data.index_d],
        ),
        read_param(
            "v_C_my",
            input_filename,
            "VariableCostNewIPPmy",
            model_data.index_t,
            [model_data.index_y, index_p, index_k_new, model_data.index_z, model_data.index_d],
        ),
        read_param(
            "rho_E_my",
            input_filename,
            "AvailabilityOldIPP",
            model_data.index_t,
            [index_p, index_k_existing, model_data.index_z, model_data.index_d],
        ),
        read_param(
            "rho_C_my",
            input_filename,
            "AvailabilityNewIPP",
            model_data.index_t,
            [index_p, index_k_new, model_data.index_z, model_data.index_d],
        ),
        read_param(
            "eximport_my",
            input_filename,
            "Exportmy",
            model_data.index_t,
            [model_data.index_y, model_data.index_z, model_data.index_d],
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
        read_param(
            "initial_energy_existing_my",
            input_filename,
            "ExistingStorInitialEnergyIPP",
            model_data.index_d,
            [model_data.index_y, index_p, index_stor_existing, model_data.index_z],
        ),
        read_param(
            "initial_energy_new_my",
            input_filename,
            "NewStorInitialEnergyIPP",
            model_data.index_d,
            [model_data.index_y, index_p, index_stor_new, model_data.index_z],
        ),
        read_param(
            "stor_duration_existing",
            input_filename,
            "ExistingStorDuration",
            index_stor_existing,
        ),
        read_param(
            "stor_duration_new",
            input_filename,
            "NewStorDuration",
            index_stor_new,
        ),
        read_param(
            "trans_topology",
            input_filename,
            "TransmissionTopology",
            model_data.index_z,
            [index_l],
        ),
        read_param(
            "trans_capacity",
            input_filename,
            "TransmissionCapacity",
            min_max,
            [index_l],
        ),
        initialize_param(
            "y_E_my",
            model_data.index_y,
            index_p,
            index_k_existing,
            model_data.index_z,
            model_data.index_d,
            model_data.index_t,
        ),
        initialize_param(
            "y_C_my",
            model_data.index_y,
            index_p,
            index_k_new,
            model_data.index_z,
            model_data.index_d,
            model_data.index_t,
        ),
        initialize_param("x_R_my", model_data.index_y, index_p, index_k_existing, model_data.index_z),
        initialize_param("x_C_my", model_data.index_y, index_p, index_k_new, model_data.index_z),
        initialize_param("x_stor_R_my", model_data.index_y, index_p, index_stor_existing, model_data.index_z),
        initialize_param("x_stor_C_my", model_data.index_y, index_p, index_stor_new, model_data.index_z),
        read_param(
            "o_E_my",
            input_filename,
            "VariableCostOldIPPmy",
            model_data.index_t,
            [model_data.index_y, index_p, index_k_existing, model_data.index_z, model_data.index_d],
        ),
        read_param(
            "o_C_my",
            input_filename,
            "VariableCostNewIPPmy",
            model_data.index_t,
            [model_data.index_y, index_p, index_k_new, model_data.index_z, model_data.index_d],
        ),
        initialize_param("miu_my", model_data.index_y, model_data.index_z, model_data.index_d, model_data.index_t),
        initialize_param("LMP_my", model_data.index_y, model_data.index_z, model_data.index_d, model_data.index_t),
        initialize_param(
            "eta_my",
            model_data.index_y,
            index_p,
            index_k_existing,
            model_data.index_z,
            model_data.index_d,
            model_data.index_t,
        ),
        initialize_param(
            "lambda_my",
            model_data.index_y,
            index_p,
            index_k_new,
            model_data.index_z,
            model_data.index_d,
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
        initialize_param("x_R_cumu", index_p, index_k_existing, model_data.index_z),
        initialize_param("x_C_cumu", index_p, index_k_new, model_data.index_z),
        initialize_param("x_stor_R_cumu", index_p, index_stor_existing, model_data.index_z),
        initialize_param("x_stor_C_cumu", index_p, index_stor_new, model_data.index_z),
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
            model_data.index_z,
            model_data.index_d,
            model_data.index_t,
        ),
        initialize_param(
            "y_C_my_temp",
            model_data.index_y,
            index_p,
            index_k_new,
            model_data.index_z,
            model_data.index_d,
            model_data.index_t,
        ),
        initialize_param("x_R_my_temp", model_data.index_y, index_p, index_k_existing, model_data.index_z),
        initialize_param("x_C_my_temp", model_data.index_y, index_p, index_k_new, model_data.index_z),
        initialize_param("miu_my_temp", model_data.index_y, model_data.index_z, model_data.index_d, model_data.index_t),
        initialize_param("LMP_my_temp", model_data.index_y, model_data.index_z, model_data.index_d, model_data.index_t),
        NetCONE,
        DC_length,
        read_param(
            "capacity_credit_E_my",
            input_filename,
            "CapacityCredit_old",
            index_k_existing,
            [model_data.index_y, model_data.index_z],
        ),
        read_param(
            "capacity_credit_C_my",
            input_filename,
            "CapacityCredit_new",
            index_k_new,
            [model_data.index_y, model_data.index_z],
        ),
        read_param(
            "capacity_credit_stor_E_my",
            input_filename,
            "CapacityCreditStor_old",
            index_stor_existing,
            [model_data.index_y, model_data.index_z],
        ),
        read_param(
            "capacity_credit_stor_C_my",
            input_filename,
            "CapacityCreditStor_new",
            index_stor_new,
            [model_data.index_y, model_data.index_z],
        ),
        initialize_param("capacity_price", model_data.index_y),
        initialize_param("capacity_price_my_temp", model_data.index_y),
        initialize_param("Net_Load_my", model_data.index_y, model_data.index_z, model_data.index_d, model_data.index_t),
        initialize_param("Max_Net_Load_my", model_data.index_y, model_data.index_z),
        initialize_param("Reserve_req_my", model_data.index_y, model_data.index_z),
        initialize_param("Capacity_slope_my", model_data.index_y),
        initialize_param("Capacity_intercept_my", model_data.index_y),
        initialize_param("ucap_temp", model_data.index_y, index_p),
        initialize_param("ucap", model_data.index_y, index_p),
        initialize_param("ucap_total", model_data.index_y),
        read_param("RPS", input_filename, "RPS", model_data.index_y),
        read_param(
            "emission_rate_E_my",
            input_filename,
            "EmissionRateOldIPPmy",
            index_k_existing,
            [model_data.index_y, index_p, model_data.index_z],
        ),
        read_param(
            "emission_rate_C_my",
            input_filename,
            "EmissionRateNewIPPmy",
            index_k_new,
            [model_data.index_y, index_p, model_data.index_z],
        ),
        initialize_param(
            "charge_E_my",
            model_data.index_y,
            index_p,
            index_stor_existing,
            model_data.index_z,
            model_data.index_d,
            model_data.index_t,
        ),
        initialize_param(
            "discharge_E_my",
            model_data.index_y,
            index_p,
            index_stor_existing,
            model_data.index_z,
            model_data.index_d,
            model_data.index_t,
        ),
        initialize_param(
            "charge_C_my",
            model_data.index_y,
            index_p,
            index_stor_new,
            model_data.index_z,
            model_data.index_d,
            model_data.index_t,
        ),
        initialize_param(
            "discharge_C_my",
            model_data.index_y,
            index_p,
            index_stor_new,
            model_data.index_z,
            model_data.index_d,
            model_data.index_t,
        ),
        initialize_param(
            "energy_E_my",
            model_data.index_y,
            index_p,
            index_stor_existing,
            model_data.index_z,
            model_data.index_d,
            model_data.index_t,
        ),
        initialize_param(
            "energy_C_my",
            model_data.index_y,
            index_p,
            index_stor_new,
            model_data.index_z,
            model_data.index_d,
            model_data.index_t,
        ),
        initialize_param(
            "flow_my",
            model_data.index_y,
            index_l,
            model_data.index_d,
            model_data.index_t,
        ),
        Dict(),
        repeat([[]], 24)...
    )
end

get_id(x::IPPGroup) = x.id

function solve_agent_problem!(
    ipps::IPPGroup,
    ipp_opts::IPPOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{VIU},
    agent_store::AgentStore,
    w_iter,
    window_length,
    jump_model,
    export_file_path,
    update_results::Bool
)
    return 0.0
end


##################### add transmission and storage #####################
function solve_agent_problem_ipp_cap(
    ipp::IPPGroup,
    ipp_opts::IPPOptions{MPPDCMERTransStorage},
    p_star,
    model_data::HEMData,
    hem_opts::HEMOptions{WM},
    agent_store::AgentStore,
    w_iter,
    window_length,
    jump_model
)
    x_R_before = ParamArray(ipp.x_R_my)
    x_C_before = ParamArray(ipp.x_C_my)
    delta_t = get_delta_t(model_data)
    # the aggregator problem hasn't solved yet, so use last year's participation rates
    reg_year_dera, reg_year_index_dera = get_prev_reg_year(model_data, w_iter)
    reg_year_dera_pre, reg_year_index_dera_pre = get_prev_two_reg_year(model_data, w_iter)

    # utility = get_agent(Utility, agent_store)
    regulator = get_agent(Regulator, agent_store)
    customers = get_agent(CustomerGroup, agent_store)
    der_aggregator = get_agent(DERAggregator, agent_store)
    green_developer = get_agent(GreenDeveloper, agent_store)

    WMDER_IPP = get_new_jump_model(ipp_opts.solvers["solve_agent_problem_ipp_mppdc"])
    MPPDCMER_lower = get_new_jump_model(ipp_opts.solvers["solve_agent_problem_ipp_mppdc_mccormic_lower"])
    MPPDCMER_lower_dual = get_new_jump_model(ipp_opts.solvers["solve_agent_problem_ipp_mppdc_mccormic_lower"])

    if w_iter >= 2
        model_data.index_y.elements =
                model_data.index_y_fix.elements[w_iter-1:(w_iter-1 + window_length - 1)]
    end

    # first, use lower level optimization results to set variable bounds for McCormick-envelope Relaxation
    @variable(
        MPPDCMER_lower,
        y_E_bounds[model_data.index_y, ipp.index_p, ipp.index_k_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower,
        y_C_bounds[model_data.index_y, ipp.index_p, ipp.index_k_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower,
        charge_E_bounds[model_data.index_y, ipp.index_p, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower,
        discharge_E_bounds[model_data.index_y, ipp.index_p, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower,
        charge_C_bounds[model_data.index_y, ipp.index_p, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower,
        discharge_C_bounds[model_data.index_y, ipp.index_p, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower,
        energy_E_bounds[model_data.index_y, ipp.index_p, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower,
        energy_C_bounds[model_data.index_y, ipp.index_p, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower,
        flow_bounds[model_data.index_y, ipp.index_l, model_data.index_d, model_data.index_t]
    )

    objective_function_lower = begin
        sum(
            sum(
                model_data.omega(d) * delta_t *
                (ipp.v_E_my(y, p, k, z, d, t) * y_E_bounds[y, p, k, z, d, t]) for
                d in model_data.index_d, t in model_data.index_t, z in model_data.index_z, k in ipp.index_k_existing, p in ipp.index_p
            ) + 
            sum(
                model_data.omega(d) * delta_t *
                (ipp.v_C_my(y, p, k, z, d, t) * y_C_bounds[y, p, k, z, d, t]) for
                d in model_data.index_d, t in model_data.index_t, z in model_data.index_z, k in ipp.index_k_new, p in ipp.index_p
            )
            for y in model_data.index_y
        )
    end

    @objective(MPPDCMER_lower, Min, objective_function_lower)

    supply_demand_balance_lower =
        (y, z, d, t) -> begin
            # bulk generation at time t
            sum(y_E_bounds[y, p, k, z, d, t] for k in ipp.index_k_existing, p in ipp.index_p) +
            sum(y_C_bounds[y, p, k, z, d, t] for k in ipp.index_k_new, p in ipp.index_p) -
            # flow out of zone z
            sum(ipp.trans_topology(l, z) * flow_bounds[y, l, d, t] for l in ipp.index_l) +
            # battery discharge
            sum(discharge_E_bounds[y, p, s, z, d, t] for s in ipp.index_stor_existing, p in ipp.index_p) +
            sum(discharge_C_bounds[y, p, s, z, d, t] for s in ipp.index_stor_new, p in ipp.index_p) -
            # battery charge
            sum(charge_E_bounds[y, p, s, z, d, t] for s in ipp.index_stor_existing, p in ipp.index_p) -
            sum(charge_C_bounds[y, p, s, z, d, t] for s in ipp.index_stor_new, p in ipp.index_p) -
            # demand at time t
            sum(customers.gamma(z, h) * customers.d_my(y, h, z, d, t) for h in model_data.index_h) -
            ipp.eximport_my(y, z, d, t) +
            # total DG generation at time t
            sum(
                customers.rho_DG(h, m, z, d, t) * customers.total_der_capacity_my(y, z, h, m) for
                h in model_data.index_h, m in customers.index_m
            ) - 
            # remove aggregated behind-the-meter pv/storage generation/consumption since they're front-of-the-meter now
            sum(
                customers.rho_DG(h, m, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera_pre, z) *
                customers.total_pv_stor_capacity_my(y, z, h, m) for h in model_data.index_h, m in (:BTMStorage, :BTMPV)
            ) +
            # green technology subscription at time t
            sum(
                ipp.rho_C_my(Symbol("ipp1"), j, z, d, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, z, h) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):model_data.year(y))
                for j in model_data.index_j, h in model_data.index_h
            )
        end

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_supplydemandbalance_lower[y in model_data.index_y, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t],
        supply_demand_balance_lower(y, z, d, t) == 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_gen_max_E_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            k in ipp.index_k_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        ipp.rho_E_my(p, k, z, d, t) * (
            ipp.x_E_my(p, z, k) - ipp.x_R_cumu(p, k, z)
        ) - y_E_bounds[y, p, k, z, d, t] >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_gen_max_C_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            k in ipp.index_k_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        ipp.rho_C_my(p, k, z, d, t) * (
            ipp.x_C_cumu(p, k, z)
        ) - y_C_bounds[y, p, k, z, d, t] >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_flow_lower_lower[
            y in model_data.index_y,
            l in ipp.index_l,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        flow_bounds[y, l, d, t] - ipp.trans_capacity(l, :min) >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_flow_upper_lower[
            y in model_data.index_y,
            l in ipp.index_l,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        ipp.trans_capacity(l, :max) - flow_bounds[y, l, d, t] >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_energy_E_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t.elements[2:end],
        ],
        energy_E_bounds[y, p, s, z, d, t] == energy_E_bounds[y, p, s, z, d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)-delta_t), model_data.time.values)][1]] - discharge_E_bounds[y, p, s, z, d, t] / ipp.rte_stor_E_my(y, p, z, s) * delta_t +
            charge_E_bounds[y, p, s, z, d, t] * delta_t
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_energy_E_0_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in [model_data.index_t.elements[1]],
        ],
        energy_E_bounds[y, p, s, z, d, t] == ipp.initial_energy_existing_my(y, p, s, z, d) - discharge_E_bounds[y, p, s, z, d, t] / ipp.rte_stor_E_my(y, p, z, s) * delta_t +
            charge_E_bounds[y, p, s, z, d, t] * delta_t
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_energy_upper_bound_E_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        ipp.stor_duration_existing(s) * (
            ipp.x_stor_E_my(p, z, s) - ipp.x_stor_R_cumu(p, s, z)
        ) - energy_E_bounds[y, p, s, z, d, t] >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_discharge_upper_bound_E_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        ipp.rte_stor_E_my(y, p, z, s) * (
            ipp.x_stor_E_my(p, z, s) - ipp.x_stor_R_cumu(p, s, z)
        ) - discharge_E_bounds[y, p, s, z, d, t] >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_charge_upper_bound_E_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        ipp.x_stor_E_my(p, z, s) - ipp.x_stor_R_cumu(p, s, z) - charge_E_bounds[y, p, s, z, d, t] >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_discharge_energy_upper_bound_E_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t.elements[2:end],
        ],
        energy_E_bounds[y, p, s, z, d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)-delta_t), model_data.time.values)][1]] >= discharge_E_bounds[y, p, s, z, d, t] * delta_t
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_discharge_energy_upper_bound_E_0_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in [model_data.index_t.elements[1]],
        ],
        ipp.initial_energy_existing_my(y, p, s, z, d) >= discharge_E_bounds[y, p, s, z, d, t] * delta_t
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_charge_energy_upper_bound_E_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t.elements[2:end],
        ],
        ipp.stor_duration_existing(s) * (
            ipp.x_stor_E_my(p, z, s) - ipp.x_stor_R_cumu(p, s, z)
        ) -
        energy_E_bounds[y, p, s, z, d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)-delta_t), model_data.time.values)][1]] -
        charge_E_bounds[y, p, s, z, d, t] * delta_t >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_charge_energy_upper_bound_E_0_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in [model_data.index_t.elements[1]],
        ],
        ipp.stor_duration_existing(s) * (
            ipp.x_stor_E_my(p, z, s) - ipp.x_stor_R_cumu(p, s, z)
        ) -
        ipp.initial_energy_existing_my(y, p, s, z, d) - charge_E_bounds[y, p, s, z, d, t] * delta_t >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_charge_discharge_upper_bound_E_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t.elements,
        ],
            ipp.x_stor_E_my(p, z, s) - ipp.x_stor_R_cumu(p, s, z) -  
            charge_E_bounds[y, p, s, z, d, t] - discharge_E_bounds[y, p, s, z, d, t] / ipp.rte_stor_E_my(y, p, z, s) >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_energy_C_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t.elements[2:end],
        ],
        energy_C_bounds[y, p, s, z, d, t] == energy_C_bounds[y, p, s, z, d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)-delta_t), model_data.time.values)][1]] - discharge_C_bounds[y, p, s, z, d, t] / ipp.rte_stor_C_my(y, p, z, s) * delta_t +
            charge_C_bounds[y, p, s, z, d, t] * delta_t
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_energy_C_0_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in [model_data.index_t.elements[1]],
        ],
        energy_C_bounds[y, p, s, z, d, t] == ipp.initial_energy_new_my(y, p, s, z, d) - discharge_C_bounds[y, p, s, z, d, t] / ipp.rte_stor_C_my(y, p, z, s) * delta_t +
            charge_C_bounds[y, p, s, z, d, t] * delta_t
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_energy_upper_bound_C_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        ipp.stor_duration_new(s) * (
            ipp.x_stor_C_cumu(p, s, z)
        ) - energy_C_bounds[y, p, s, z, d, t] >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_discharge_upper_bound_C_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        ipp.rte_stor_C_my(y, p, z, s) * (
            ipp.x_stor_C_cumu(p, s, z)
        ) - discharge_C_bounds[y, p, s, z, d, t] >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_charge_upper_bound_C_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        ipp.x_stor_C_cumu(p, s, z) - charge_C_bounds[y, p, s, z, d, t] >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_discharge_energy_upper_bound_C_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t.elements[2:end],
        ],
        energy_C_bounds[y, p, s, z, d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)-delta_t), model_data.time.values)][1]] - discharge_C_bounds[y, p, s, z, d, t] * delta_t >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_discharge_energy_upper_bound_C_0_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in [model_data.index_t.elements[1]],
        ],
        ipp.initial_energy_new_my(y, p, s, z, d) - discharge_C_bounds[y, p, s, z, d, t] * delta_t >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_charge_energy_upper_bound_C_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t.elements[2:end],
        ],
        ipp.stor_duration_new(s) * (
            ipp.x_stor_C_cumu(p, s, z)
        ) -
        energy_C_bounds[y, p, s, z, d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)-delta_t), model_data.time.values)][1]] -
        charge_C_bounds[y, p, s, z, d, t] * delta_t >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_charge_energy_upper_bound_C_0_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in [model_data.index_t.elements[1]],
        ],
        ipp.stor_duration_new(s) * (
            ipp.x_stor_C_cumu(p, s, z)
        ) -
        ipp.initial_energy_new_my(y, p, s, z, d) - charge_C_bounds[y, p, s, z, d, t] * delta_t >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_charge_discharge_upper_bound_C_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t.elements,
        ],
        ipp.x_stor_C_cumu(p, s, z) - charge_C_bounds[y, p, s, z, d, t] - discharge_C_bounds[y, p, s, z, d, t] / ipp.rte_stor_C_my(y, p, z, s) >= 0
    )

    push!(jump_model, MPPDCMER_lower)

    @info("MPPDC lower level primal")
    TimerOutputs.@timeit HEM_TIMER "optimize! MPPDC McCormic Envelope Relaxation lower level primal" begin
        optimize!(MPPDCMER_lower)
    end

    # TEMPORARY CODE FOR TESTING
    # objective_value(MPPDCMER_lower) # calling this errors the program if the optimization failed
    # dual_model = dualize(MPPDCMER_lower; dual_names = DualNames("dual", ""))
    # f = open("lower_level_dual.txt","w"); print(f, dual_model); close(f)

    ############ lower-level dual problem ############
    @variable(MPPDCMER_lower_dual, miu_lower[model_data.index_y, model_data.index_z, model_data.index_d, model_data.index_t])
    @variable(
        MPPDCMER_lower_dual,
        eta_lower[model_data.index_y, ipp.index_p, ipp.index_k_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower_dual,
        lambda_lower[model_data.index_y, ipp.index_p, ipp.index_k_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )

    @variable(
        MPPDCMER_lower_dual,
        iota_min_lower[model_data.index_y, ipp.index_l, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower_dual,
        iota_max_lower[model_data.index_y, ipp.index_l, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower_dual,
        psi_E_lower[model_data.index_y, ipp.index_p, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t]
    )
    @variable(
        MPPDCMER_lower_dual,
        theta_E_energy_lower[model_data.index_y, ipp.index_p, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower_dual,
        theta_E_discharge_lower[model_data.index_y, ipp.index_p, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower_dual,
        theta_E_charge_lower[model_data.index_y, ipp.index_p, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower_dual,
        pi_E_discharge_lower[model_data.index_y, ipp.index_p, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower_dual,
        pi_E_charge_lower[model_data.index_y, ipp.index_p, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower_dual,
        kappa_E_lower[model_data.index_y, ipp.index_p, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower_dual,
        psi_C_lower[model_data.index_y, ipp.index_p, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t]
    )
    @variable(
        MPPDCMER_lower_dual,
        theta_C_energy_lower[model_data.index_y, ipp.index_p, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower_dual,
        theta_C_discharge_lower[model_data.index_y, ipp.index_p, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower_dual,
        theta_C_charge_lower[model_data.index_y, ipp.index_p, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower_dual,
        pi_C_discharge_lower[model_data.index_y, ipp.index_p, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower_dual,
        pi_C_charge_lower[model_data.index_y, ipp.index_p, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower_dual,
        kappa_C_lower[model_data.index_y, ipp.index_p, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )

    objective_function_lower_dual = begin
        sum(
            sum(
                miu_lower[y, z, d, t] * (
                    sum(
                        customers.gamma(z, h) * customers.d_my(y, h, z, d, t) for h in model_data.index_h
                    ) + ipp.eximport_my(y, z, d, t) -
                    # total DG generation at time t
                    sum(
                        customers.rho_DG(h, m, z, d, t) * customers.total_der_capacity_my(y, z, h, m) for
                        h in model_data.index_h, m in customers.index_m
                    ) +
                    # remove aggregated behind-the-meter pv/storage generation/consumption since they're front-of-the-meter now
                    sum(
                        customers.rho_DG(h, m, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera_pre, z) * 
                        customers.total_pv_stor_capacity_my(y, z, h, m) for h in model_data.index_h, m in (:BTMStorage, :BTMPV)
                    ) -
                    # green technology subscription at time t
                    sum(
                        ipp.rho_C_my(Symbol("ipp1"), j, z, d, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, z, h) for y_symbol in
                        model_data.year(first(model_data.index_y_fix)):model_data.year(y))
                        for j in model_data.index_j, h in model_data.index_h
                    )
                ) for z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
            ) - 
            (
                sum(
                    eta_lower[y, p, k, z, d, t] *
                    ipp.rho_E_my(p, k, z, d, t) * (
                        ipp.x_E_my(p, z, k) - ipp.x_R_cumu(p, k, z)
                    ) for p in ipp.index_p, k in ipp.index_k_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    lambda_lower[y, p, k, z, d, t] *
                    ipp.rho_C_my(p, k, z, d, t) * (
                        ipp.x_C_cumu(p, k, z)
                    ) for p in ipp.index_p, k in ipp.index_k_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) + 
            (
                sum(iota_min_lower[y, l, d, t] * ipp.trans_capacity(l, :min) for l in ipp.index_l, d in model_data.index_d, t in model_data.index_t) - 
                sum(iota_max_lower[y, l, d, t] * ipp.trans_capacity(l, :max) for l in ipp.index_l, d in model_data.index_d, t in model_data.index_t)
            ) - 
            (
                sum(
                    theta_E_energy_lower[y, p, s, z, d, t] *
                    ipp.stor_duration_existing(s) * (
                        ipp.x_stor_E_my(p, z, s) - ipp.x_stor_R_cumu(p, s, z)
                    ) for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    theta_E_discharge_lower[y, p, s, z, d, t] *
                    ipp.rte_stor_E_my(y, p, z, s) * (
                        ipp.x_stor_E_my(p, z, s) - ipp.x_stor_R_cumu(p, s, z)
                    ) for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    theta_E_charge_lower[y, p, s, z, d, t] *
                    (
                        ipp.x_stor_E_my(p, z, s) - ipp.x_stor_R_cumu(p, s, z)
                    ) for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    pi_E_charge_lower[y, p, s, z, d, t] *
                    ipp.stor_duration_existing(s) * (
                        ipp.x_stor_E_my(p, z, s) - ipp.x_stor_R_cumu(p, s, z)
                    ) for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    kappa_E_lower[y, p, s, z, d, t] *
                    (
                        ipp.x_stor_E_my(p, z, s) - ipp.x_stor_R_cumu(p, s, z)
                    ) for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    theta_C_energy_lower[y, p, s, z, d, t] *
                    ipp.stor_duration_new(s) * (
                        ipp.x_stor_C_cumu(p, s, z)
                    )
                    for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    theta_C_discharge_lower[y, p, s, z, d, t] *
                    ipp.rte_stor_C_my(y, p, z, s) * (
                        ipp.x_stor_C_cumu(p, s, z)
                    ) 
                    for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    theta_C_charge_lower[y, p, s, z, d, t] * (
                        ipp.x_stor_C_cumu(p, s, z)
                    ) 
                    for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    pi_C_charge_lower[y, p, s, z, d, t] *
                    ipp.stor_duration_new(s) * (
                        ipp.x_stor_C_cumu(p, s, z)
                    ) 
                    for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    kappa_C_lower[y, p, s, z, d, t] * (
                        ipp.x_stor_C_cumu(p, s, z)
                    ) 
                    for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) + 
            sum(psi_E_lower[y, p, s, z, d, model_data.index_t.elements[1]] * ipp.initial_energy_existing_my(y, p, s, z, d)
                for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d) - 
            sum(pi_E_discharge_lower[y, p, s, z, d, model_data.index_t.elements[1]] * ipp.initial_energy_existing_my(y, p, s, z, d)
                for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d) + 
            sum(pi_E_charge_lower[y, p, s, z, d, model_data.index_t.elements[1]] * ipp.initial_energy_existing_my(y, p, s, z, d)
                for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d) + 
            sum(psi_C_lower[y, p, s, z, d, model_data.index_t.elements[1]] * ipp.initial_energy_new_my(y, p, s, z, d)
                for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d) - 
            sum(pi_C_discharge_lower[y, p, s, z, d, model_data.index_t.elements[1]] * ipp.initial_energy_new_my(y, p, s, z, d)
                for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d) + 
            sum(pi_C_charge_lower[y, p, s, z, d, model_data.index_t.elements[1]] * ipp.initial_energy_new_my(y, p, s, z, d)
                for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d)
            for y in model_data.index_y
        )
    end

    @objective(MPPDCMER_lower_dual, Max, objective_function_lower_dual)

    # dual feasible constraints
    @constraint(
        MPPDCMER_lower_dual,
        Eq_dual_feasible_y_E_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            k in ipp.index_k_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        model_data.omega(d) * delta_t * ipp.v_E_my(y, p, k, z, d, t) - miu_lower[y, z, d, t] + eta_lower[y, p, k, z, d, t] >= 0
    )
    @constraint(
        MPPDCMER_lower_dual,
        Eq_dual_feasible_y_C_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            k in ipp.index_k_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        model_data.omega(d) * delta_t * ipp.v_C_my(y, p, k, z, d, t) - miu_lower[y, z, d, t] + lambda_lower[y, p, k, z, d, t] >= 0
    )
    @constraint(
        MPPDCMER_lower_dual,
        Eq_dual_feasible_flow_lower[
            y in model_data.index_y,
            l in ipp.index_l,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        sum(miu_lower[y, z, d, t] * ipp.trans_topology(l, z) for z in model_data.index_z) - iota_min_lower[y, l, d, t] + iota_max_lower[y, l, d, t] == 0
    )
    @constraint(
        MPPDCMER_lower_dual,
        Eq_dual_feasible_discharge_E_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        - miu_lower[y, z, d, t] - psi_E_lower[y, p, s, z, d, t] / ipp.rte_stor_E_my(y, p, z, s) * delta_t + theta_E_discharge_lower[y, p, s, z, d, t] + pi_E_discharge_lower[y, p, s, z, d, t] * delta_t + kappa_E_lower[y, p, s, z, d, t] / ipp.rte_stor_E_my(y, p, z, s) >= 0
    )
    @constraint(
        MPPDCMER_lower_dual,
        Eq_dual_feasible_charge_E_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        miu_lower[y, z, d, t] + psi_E_lower[y, p, s, z, d, t] * delta_t + theta_E_charge_lower[y, p, s, z, d, t] + pi_E_charge_lower[y, p, s, z, d, t] * delta_t + kappa_E_lower[y, p, s, z, d, t] >= 0
    )
    @constraint(
        MPPDCMER_lower_dual,
        Eq_dual_feasible_energy_E_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t.elements[1:end-1],
        ],
        psi_E_lower[y, p, s, z, d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)+delta_t), model_data.time.values)][1]] - psi_E_lower[y, p, s, z, d, t] + theta_E_energy_lower[y, p, s, z, d, t] 
        - pi_E_discharge_lower[y, p, s, z, d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)+delta_t), model_data.time.values)][1]] 
        + pi_E_charge_lower[y, p, s, z, d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)+delta_t), model_data.time.values)][1]] >= 0
    )
    @constraint(
        MPPDCMER_lower_dual,
        Eq_dual_feasible_energy_last_E_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in [model_data.index_t.elements[end]],
        ],
        - psi_E_lower[y, p, s, z, d, t] + theta_E_energy_lower[y, p, s, z, d, t] >= 0
    )
    @constraint(
        MPPDCMER_lower_dual,
        Eq_dual_feasible_discharge_C_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        - miu_lower[y, z, d, t] - psi_C_lower[y, p, s, z, d, t] / ipp.rte_stor_C_my(y, p, z, s) * delta_t + theta_C_discharge_lower[y, p, s, z, d, t] + pi_C_discharge_lower[y, p, s, z, d, t] * delta_t + kappa_C_lower[y, p, s, z, d, t] / ipp.rte_stor_C_my(y, p, z, s) >= 0
    )
    @constraint(
        MPPDCMER_lower_dual,
        Eq_dual_feasible_charge_C_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        miu_lower[y, z, d, t] + psi_C_lower[y, p, s, z, d, t] * delta_t + theta_C_charge_lower[y, p, s, z, d, t] + pi_C_charge_lower[y, p, s, z, d, t] * delta_t + kappa_C_lower[y, p, s, z, d, t] >= 0
    )
    @constraint(
        MPPDCMER_lower_dual,
        Eq_dual_feasible_energy_C_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t.elements[1:end-1],
        ],
        psi_C_lower[y, p, s, z, d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)+delta_t), model_data.time.values)][1]] - psi_C_lower[y, p, s, z, d, t] + theta_C_energy_lower[y, p, s, z, d, t] 
        - pi_C_discharge_lower[y, p, s, z, d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)+delta_t), model_data.time.values)][1]] 
        + pi_C_charge_lower[y, p, s, z, d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)+delta_t), model_data.time.values)][1]] >= 0
    )
    @constraint(
        MPPDCMER_lower_dual,
        Eq_dual_feasible_energy_last_C_lower[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in [model_data.index_t.elements[end]],
        ],
        - psi_C_lower[y, p, s, z, d, t] + theta_C_energy_lower[y, p, s, z, d, t] >= 0
    )

    jump_model[end] = MPPDCMER_lower_dual

    @info("MPPDC lower level dual")
    TimerOutputs.@timeit HEM_TIMER "optimize! MPPDC McCormic Envelope Relaxation lower level dual problem" begin
        optimize!(MPPDCMER_lower_dual)
    end

    # TEMPORARY CODE FOR TESTING
    # objective_value(MPPDCMER_lower_dual) # calling this errors the program if the optimization failed
    # f = open("lower_level_dual_my_version.txt","w"); print(f, MPPDCMER_lower_dual); close(f)
    # abs.(value.(miu_lower).data) .- abs.(dual.(Eq_primal_feasible_supplydemandbalance_lower).data)

    model_data.index_y.elements =
                model_data.index_y_fix.elements[w_iter:(w_iter + window_length - 1)]

    if termination_status(MPPDCMER_lower) == OPTIMAL
        # Lower-level problem completed successfully. Calculate and save McCormick bounds.

        eta_param = initialize_param("eta_param", model_data.index_y, ipp.index_k_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        fill!(eta_param, NaN)
        for y in model_data.index_y, k in ipp.index_k_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
            y_minus = (w_iter >= 2) ? model_data.year(y) - 1 : model_data.year(y)
            eta_param(y, k, z, d, t, :) .= abs(dual.(Eq_primal_feasible_gen_max_E_lower[Symbol(y_minus), p_star, k, z, d, t]))
        end

        lambda_param = initialize_param("lambda_param", model_data.index_y, ipp.index_k_new, model_data.index_z, model_data.index_d, model_data.index_t)
        fill!(lambda_param, NaN)
        for y in model_data.index_y, k in ipp.index_k_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
            y_minus = (w_iter >= 2) ? model_data.year(y) - 1 : model_data.year(y)
            lambda_param(y, k, z, d, t, :) .= abs(dual.(Eq_primal_feasible_gen_max_C_lower[Symbol(y_minus), p_star, k, z, d, t]))
        end

        theta_E_energy_param = initialize_param("theta_E_energy_param", model_data.index_y, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        fill!(theta_E_energy_param, NaN)
        for y in model_data.index_y, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
            y_minus = (w_iter >= 2) ? model_data.year(y) - 1 : model_data.year(y)
            theta_E_energy_param(y, s, z, d, t, :) .= abs(dual.(Eq_primal_feasible_energy_upper_bound_E_lower[Symbol(y_minus), p_star, s, z, d, t]))
        end

        theta_E_discharge_param = initialize_param("theta_E_discharge_param", model_data.index_y, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        fill!(theta_E_discharge_param, NaN)
        for y in model_data.index_y, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
            y_minus = (w_iter >= 2) ? model_data.year(y) - 1 : model_data.year(y)
            theta_E_discharge_param(y, s, z, d, t, :) .= abs(dual.(Eq_primal_feasible_discharge_upper_bound_E_lower[Symbol(y_minus), p_star, s, z, d, t]))
        end

        theta_E_charge_param = initialize_param("theta_E_charge_param", model_data.index_y, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        fill!(theta_E_charge_param, NaN)
        for y in model_data.index_y, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
            y_minus = (w_iter >= 2) ? model_data.year(y) - 1 : model_data.year(y)
            theta_E_charge_param(y, s, z, d, t, :) .= abs(dual.(Eq_primal_feasible_charge_upper_bound_E_lower[Symbol(y_minus), p_star, s, z, d, t]))
        end

        pi_E_charge_param = initialize_param("pi_E_charge_param", model_data.index_y, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        fill!(pi_E_charge_param, NaN)
        for y in model_data.index_y, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t.elements[2:end]
            y_minus = (w_iter >= 2) ? model_data.year(y) - 1 : model_data.year(y)
            pi_E_charge_param(y, s, z, d, t, :) .= abs(dual.(Eq_primal_feasible_charge_energy_upper_bound_E_lower[Symbol(y_minus), p_star, s, z, d, t]))
        end
        for y in model_data.index_y, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in [model_data.index_t.elements[1]]
            y_minus = (w_iter >= 2) ? model_data.year(y) - 1 : model_data.year(y)
            pi_E_charge_param(y, s, z, d, t, :) .= abs(dual.(Eq_primal_feasible_charge_energy_upper_bound_E_0_lower[Symbol(y_minus), p_star, s, z, d, t]))
        end

        kappa_E_param = initialize_param("kappa_E_param", model_data.index_y, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        fill!(kappa_E_param, NaN)
        for y in model_data.index_y, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
            y_minus = (w_iter >= 2) ? model_data.year(y) - 1 : model_data.year(y)
            kappa_E_param(y, s, z, d, t, :) .= abs(dual.(Eq_primal_feasible_charge_discharge_upper_bound_E_lower[Symbol(y_minus), p_star, s, z, d, t]))
        end

        theta_C_energy_param = initialize_param("theta_C_energy_param", model_data.index_y, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        fill!(theta_C_energy_param, NaN)
        for y in model_data.index_y, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
            y_minus = (w_iter >= 2) ? model_data.year(y) - 1 : model_data.year(y)
            theta_C_energy_param(y, s, z, d, t, :) .= abs(dual.(Eq_primal_feasible_energy_upper_bound_C_lower[Symbol(y_minus), p_star, s, z, d, t]))
        end

        theta_C_discharge_param = initialize_param("theta_C_discharge_param", model_data.index_y, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        fill!(theta_C_discharge_param, NaN)
        for y in model_data.index_y, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
            y_minus = (w_iter >= 2) ? model_data.year(y) - 1 : model_data.year(y)
            theta_C_discharge_param(y, s, z, d, t, :) .= abs(dual.(Eq_primal_feasible_discharge_upper_bound_C_lower[Symbol(y_minus), p_star, s, z, d, t]))
        end

        theta_C_charge_param = initialize_param("theta_C_charge_param", model_data.index_y, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        fill!(theta_C_charge_param, NaN)
        for y in model_data.index_y, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
            y_minus = (w_iter >= 2) ? model_data.year(y) - 1 : model_data.year(y)
            theta_C_charge_param(y, s, z, d, t, :) .= abs(dual.(Eq_primal_feasible_charge_upper_bound_C_lower[Symbol(y_minus), p_star, s, z, d, t]))
        end

        pi_C_charge_param = initialize_param("pi_C_charge_param", model_data.index_y, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        fill!(pi_C_charge_param, NaN)
        for y in model_data.index_y, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t.elements[2:end]
            y_minus = (w_iter >= 2) ? model_data.year(y) - 1 : model_data.year(y)
            pi_C_charge_param(y, s, z, d, t, :) .= abs(dual.(Eq_primal_feasible_charge_energy_upper_bound_C_lower[Symbol(y_minus), p_star, s, z, d, t]))
        end
        for y in model_data.index_y, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in [model_data.index_t.elements[1]]
            y_minus = (w_iter >= 2) ? model_data.year(y) - 1 : model_data.year(y)
            pi_C_charge_param(y, s, z, d, t, :) .= abs(dual.(Eq_primal_feasible_charge_energy_upper_bound_C_0_lower[Symbol(y_minus), p_star, s, z, d, t]))
        end

        kappa_C_param = initialize_param("kappa_C_param", model_data.index_y, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        fill!(kappa_C_param, NaN)
        for y in model_data.index_y, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
            y_minus = (w_iter >= 2) ? model_data.year(y) - 1 : model_data.year(y)
            kappa_C_param(y, s, z, d, t, :) .= abs(dual.(Eq_primal_feasible_charge_discharge_upper_bound_C_lower[Symbol(y_minus), p_star, s, z, d, t]))
        end

        # TEMPORARY CODE FOR TESTING
        # test bounds McCormick-envelope Relaxation
        # eta_param = CSV.read(joinpath("/home/nguo/HolisticElectricityModel-Data/outputs/ba_1_base_2018_future_1_ipps_1", "eta.csv"), DataFrame)
        # lambda_param = CSV.read(joinpath("/home/nguo/HolisticElectricityModel-Data/outputs/ba_1_base_2018_future_1_ipps_1", "lambda.csv"), DataFrame)

        # adjust these bounds may make the problem easier to solve! But it may also hurt the duality gap.
        eta_upper_bound_adj = 1.3 # 1.1
        eta_lower_bound_adj = 0.7 # 0.8
        lambda_upper_bound_adj = 1.3
        lambda_lower_bound_adj = 0.7
        battery_upper_bound_adj = 1.3
        battery_lower_bound_adj = 0.7

        eta_L = initialize_param("eta_L", deepcopy(model_data.index_y), ipp.index_k_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        eta_U = initialize_param("eta_U", deepcopy(model_data.index_y), ipp.index_k_existing, model_data.index_z, model_data.index_d, model_data.index_t)

        for y in model_data.index_y
            for k in ipp.index_k_existing
                for z in model_data.index_z
                    for d in model_data.index_d
                        for t in model_data.index_t
                            eta_U(y, k, z, d, t, :) .= eta_param(y, k, z, d, t) * eta_upper_bound_adj
                            eta_L(y, k, z, d, t, :) .= eta_param(y, k, z, d, t) * eta_lower_bound_adj
                        end
                    end
                end
            end
        end
        push!(ipp.eta_U_vec, eta_U)
        push!(ipp.eta_L_vec, eta_L)

        lambda_L = initialize_param("lambda_L", deepcopy(model_data.index_y), ipp.index_k_new, model_data.index_z, model_data.index_d, model_data.index_t)
        lambda_U = initialize_param("lambda_U", deepcopy(model_data.index_y), ipp.index_k_new, model_data.index_z, model_data.index_d, model_data.index_t)
        
        for y in model_data.index_y
            for k in ipp.index_k_new
                for z in model_data.index_z
                    for d in model_data.index_d
                        for t in model_data.index_t
                            lambda_U(y, k, z, d, t, :) .= lambda_param(y, k, z, d, t) * lambda_upper_bound_adj
                            lambda_L(y, k, z, d, t, :) .= lambda_param(y, k, z, d, t) * lambda_lower_bound_adj
                        end
                    end
                end
            end
        end

        push!(ipp.lambda_U_vec, lambda_U)
        push!(ipp.lambda_L_vec, lambda_L)

        theta_E_energy_L = initialize_param("theta_E_energy_L", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        theta_E_energy_U = initialize_param("theta_E_energy_U", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        theta_E_discharge_L = initialize_param("theta_E_discharge_L", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        theta_E_discharge_U = initialize_param("theta_E_discharge_U", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        theta_E_charge_L = initialize_param("theta_E_charge_L", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        theta_E_charge_U = initialize_param("theta_E_charge_U", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        pi_E_charge_L = initialize_param("pi_E_charge_L", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        pi_E_charge_U = initialize_param("pi_E_charge_U", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        kappa_E_L = initialize_param("kappa_E_L", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        kappa_E_U = initialize_param("kappa_E_U", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)

        for y in model_data.index_y
            for k in ipp.index_stor_existing
                for z in model_data.index_z
                    for d in model_data.index_d
                        for t in model_data.index_t
                            theta_E_energy_U(y, k, z, d, t, :) .= theta_E_energy_param(y, k, z, d, t) * battery_upper_bound_adj
                            theta_E_energy_L(y, k, z, d, t, :) .= theta_E_energy_param(y, k, z, d, t) * battery_lower_bound_adj
                            theta_E_discharge_U(y, k, z, d, t, :) .= theta_E_discharge_param(y, k, z, d, t) * battery_upper_bound_adj
                            theta_E_discharge_L(y, k, z, d, t, :) .= theta_E_discharge_param(y, k, z, d, t) * battery_lower_bound_adj
                            theta_E_charge_U(y, k, z, d, t, :) .= theta_E_charge_param(y, k, z, d, t) * battery_upper_bound_adj
                            theta_E_charge_L(y, k, z, d, t, :) .= theta_E_charge_param(y, k, z, d, t) * battery_lower_bound_adj
                            pi_E_charge_U(y, k, z, d, t, :) .= pi_E_charge_param(y, k, z, d, t) * battery_upper_bound_adj
                            pi_E_charge_L(y, k, z, d, t, :) .= pi_E_charge_param(y, k, z, d, t) * battery_lower_bound_adj
                            kappa_E_U(y, k, z, d, t, :) .= kappa_E_param(y, k, z, d, t) * battery_upper_bound_adj
                            kappa_E_L(y, k, z, d, t, :) .= kappa_E_param(y, k, z, d, t) * battery_lower_bound_adj
                        end
                    end
                end
            end
        end

        push!(ipp.theta_E_energy_L_vec, theta_E_energy_L)
        push!(ipp.theta_E_energy_U_vec, theta_E_energy_U)
        push!(ipp.theta_E_discharge_L_vec, theta_E_discharge_L)
        push!(ipp.theta_E_discharge_U_vec, theta_E_discharge_U)
        push!(ipp.theta_E_charge_L_vec, theta_E_charge_L)
        push!(ipp.theta_E_charge_U_vec, theta_E_charge_U)
        push!(ipp.pi_E_charge_L_vec, pi_E_charge_L)
        push!(ipp.pi_E_charge_U_vec, pi_E_charge_U)
        push!(ipp.kappa_E_L_vec, kappa_E_L)
        push!(ipp.kappa_E_U_vec, kappa_E_U)

        theta_C_energy_L = initialize_param("theta_C_energy_L", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        theta_C_energy_U = initialize_param("theta_C_energy_U", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        theta_C_discharge_L = initialize_param("theta_C_discharge_L", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        theta_C_discharge_U = initialize_param("theta_C_discharge_U", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        theta_C_charge_L = initialize_param("theta_C_charge_L", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        theta_C_charge_U = initialize_param("theta_C_charge_U", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        pi_C_charge_L = initialize_param("pi_C_charge_L", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        pi_C_charge_U = initialize_param("pi_C_charge_U", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        kappa_C_L = initialize_param("kappa_C_L", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        kappa_C_U = initialize_param("kappa_C_U", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)

        for y in model_data.index_y
            for k in ipp.index_stor_new
                for z in model_data.index_z
                    for d in model_data.index_d
                        for t in model_data.index_t
                            theta_C_energy_U(y, k, z, d, t, :) .= theta_C_energy_param(y, k, z, d, t) * battery_upper_bound_adj
                            theta_C_energy_L(y, k, z, d, t, :) .= theta_C_energy_param(y, k, z, d, t) * battery_lower_bound_adj
                            theta_C_discharge_U(y, k, z, d, t, :) .= theta_C_discharge_param(y, k, z, d, t) * battery_upper_bound_adj
                            theta_C_discharge_L(y, k, z, d, t, :) .= theta_C_discharge_param(y, k, z, d, t) * battery_lower_bound_adj
                            theta_C_charge_U(y, k, z, d, t, :) .= theta_C_charge_param(y, k, z, d, t) * battery_upper_bound_adj
                            theta_C_charge_L(y, k, z, d, t, :) .= theta_C_charge_param(y, k, z, d, t) * battery_lower_bound_adj
                            pi_C_charge_U(y, k, z, d, t, :) .= pi_C_charge_param(y, k, z, d, t) * battery_upper_bound_adj
                            pi_C_charge_L(y, k, z, d, t, :) .= pi_C_charge_param(y, k, z, d, t) * battery_lower_bound_adj
                            kappa_C_U(y, k, z, d, t, :) .= kappa_C_param(y, k, z, d, t) * battery_upper_bound_adj
                            kappa_C_L(y, k, z, d, t, :) .= kappa_C_param(y, k, z, d, t) * battery_lower_bound_adj
                        end
                    end
                end
            end
        end

        push!(ipp.theta_C_energy_L_vec, theta_C_energy_L)
        push!(ipp.theta_C_energy_U_vec, theta_C_energy_U)
        push!(ipp.theta_C_discharge_L_vec, theta_C_discharge_L)
        push!(ipp.theta_C_discharge_U_vec, theta_C_discharge_U)
        push!(ipp.theta_C_charge_L_vec, theta_C_charge_L)
        push!(ipp.theta_C_charge_U_vec, theta_C_charge_U)
        push!(ipp.pi_C_charge_L_vec, pi_C_charge_L)
        push!(ipp.pi_C_charge_U_vec, pi_C_charge_U)
        push!(ipp.kappa_C_L_vec, kappa_C_L)
        push!(ipp.kappa_C_U_vec, kappa_C_U)
    else
        # Lower-level problem failed. Use previous McCormick bounds.

        eta_L = initialize_param("eta_L", deepcopy(model_data.index_y), ipp.index_k_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        eta_U = initialize_param("eta_U", deepcopy(model_data.index_y), ipp.index_k_existing, model_data.index_z, model_data.index_d, model_data.index_t)

        for y in model_data.index_y
            i = 1
            y_before = ipp.eta_U_vec[end].dims[1][i]

            for k in ipp.index_k_existing
                for z in model_data.index_z
                    for d in model_data.index_d
                        for t in model_data.index_t
                            eta_U(y, k, z, d, t, :) .= ipp.eta_U_vec[end](y_before, k, z, d, t)
                            eta_L(y, k, z, d, t, :) .= ipp.eta_L_vec[end](y_before, k, z, d, t)
                        end
                    end
                end
            end

            i += 1
        end

        lambda_L = initialize_param("lambda_L", deepcopy(model_data.index_y), ipp.index_k_new, model_data.index_z, model_data.index_d, model_data.index_t)
        lambda_U = initialize_param("lambda_U", deepcopy(model_data.index_y), ipp.index_k_new, model_data.index_z, model_data.index_d, model_data.index_t)

        for y in model_data.index_y
            i = 1
            y_before = ipp.eta_U_vec[end].dims[1][i]

            for k in ipp.index_k_new
                for z in model_data.index_z
                    for d in model_data.index_d
                        for t in model_data.index_t
                            lambda_U(y, k, z, d, t, :) .= ipp.lambda_U_vec[end](y_before, k, z, d, t)
                            lambda_L(y, k, z, d, t, :) .= ipp.lambda_L_vec[end](y_before, k, z, d, t)
                        end
                    end
                end
            end

            i += 1
        end

        theta_E_energy_L = initialize_param("theta_E_energy_L", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        theta_E_energy_U = initialize_param("theta_E_energy_U", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        theta_E_discharge_L = initialize_param("theta_E_discharge_L", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        theta_E_discharge_U = initialize_param("theta_E_discharge_U", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        theta_E_charge_L = initialize_param("theta_E_charge_L", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        theta_E_charge_U = initialize_param("theta_E_charge_U", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        pi_E_charge_L = initialize_param("pi_E_charge_L", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        pi_E_charge_U = initialize_param("pi_E_charge_U", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        kappa_E_L = initialize_param("kappa_E_L", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        kappa_E_U = initialize_param("kappa_E_U", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        
        for y in model_data.index_y
            i = 1
            y_before = ipp.eta_U_vec[end].dims[1][i]

            for k in ipp.index_stor_existing
                for z in model_data.index_z
                    for d in model_data.index_d
                        for t in model_data.index_t
                            theta_E_energy_U(y, k, z, d, t, :) .= ipp.theta_E_energy_U_vec[end](y_before, k, z, d, t)
                            theta_E_energy_L(y, k, z, d, t, :) .= ipp.theta_E_energy_L_vec[end](y_before, k, z, d, t)
                            theta_E_discharge_U(y, k, z, d, t, :) .= ipp.theta_E_discharge_U_vec[end](y_before, k, z, d, t)
                            theta_E_discharge_L(y, k, z, d, t, :) .= ipp.theta_E_discharge_L_vec[end](y_before, k, z, d, t)
                            theta_E_charge_U(y, k, z, d, t, :) .= ipp.theta_E_charge_U_vec[end](y_before, k, z, d, t)
                            theta_E_charge_L(y, k, z, d, t, :) .= ipp.theta_E_charge_L_vec[end](y_before, k, z, d, t)
                            pi_E_charge_U(y, k, z, d, t, :) .= ipp.pi_E_charge_U_vec[end](y_before, k, z, d, t)
                            pi_E_charge_L(y, k, z, d, t, :) .= ipp.pi_E_charge_L_vec[end](y_before, k, z, d, t)
                            kappa_E_U(y, k, z, d, t, :) .= ipp.kappa_E_U_vec[end](y_before, k, z, d, t)
                            kappa_E_L(y, k, z, d, t, :) .= ipp.kappa_E_L_vec[end](y_before, k, z, d, t)
                        end
                    end
                end
            end

            i += 1
        end

        theta_C_energy_L = initialize_param("theta_C_energy_L", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        theta_C_energy_U = initialize_param("theta_C_energy_U", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        theta_C_discharge_L = initialize_param("theta_C_discharge_L", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        theta_C_discharge_U = initialize_param("theta_C_discharge_U", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        theta_C_charge_L = initialize_param("theta_C_charge_L", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        theta_C_charge_U = initialize_param("theta_C_charge_U", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        pi_C_charge_L = initialize_param("pi_C_charge_L", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        pi_C_charge_U = initialize_param("pi_C_charge_U", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        kappa_C_L = initialize_param("kappa_C_L", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        kappa_C_U = initialize_param("kappa_C_U", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)

        for y in model_data.index_y
            i = 1
            y_before = ipp.eta_U_vec[end].dims[1][i]

            for k in ipp.index_stor_new
                for z in model_data.index_z
                    for d in model_data.index_d
                        for t in model_data.index_t
                            theta_C_energy_U(y, k, z, d, t, :) .= ipp.theta_C_energy_U_vec[end](y_before, k, z, d, t)
                            theta_C_energy_L(y, k, z, d, t, :) .= ipp.theta_C_energy_L_vec[end](y_before, k, z, d, t)
                            theta_C_discharge_U(y, k, z, d, t, :) .= ipp.theta_C_discharge_U_vec[end](y_before, k, z, d, t)
                            theta_C_discharge_L(y, k, z, d, t, :) .= ipp.theta_C_discharge_L_vec[end](y_before, k, z, d, t)
                            theta_C_charge_U(y, k, z, d, t, :) .= ipp.theta_C_charge_U_vec[end](y_before, k, z, d, t)
                            theta_C_charge_L(y, k, z, d, t, :) .= ipp.theta_C_charge_L_vec[end](y_before, k, z, d, t)
                            pi_C_charge_U(y, k, z, d, t, :) .= ipp.pi_C_charge_U_vec[end](y_before, k, z, d, t)
                            pi_C_charge_L(y, k, z, d, t, :) .= ipp.pi_C_charge_L_vec[end](y_before, k, z, d, t)
                            kappa_C_U(y, k, z, d, t, :) .= ipp.kappa_C_U_vec[end](y_before, k, z, d, t)
                            kappa_C_L(y, k, z, d, t, :) .= ipp.kappa_C_L_vec[end](y_before, k, z, d, t)
                        end
                    end
                end
            end

            i += 1
        end

    end


    X_R_cumu_L = initialize_param("x_R_L", model_data.index_y, ipp.index_k_existing, model_data.index_z)
    X_R_cumu_U = initialize_param("x_R_U", model_data.index_y, ipp.index_k_existing, model_data.index_z)
    fill!(X_R_cumu_U, NaN)
    X_R_stor_cumu_L = initialize_param("x_R_stor_L", model_data.index_y, ipp.index_stor_existing, model_data.index_z)
    X_R_stor_cumu_U = initialize_param("x_R_stor_U", model_data.index_y, ipp.index_stor_existing, model_data.index_z)
    fill!(X_R_stor_cumu_U, NaN)
    for y in model_data.index_y
        for k in ipp.index_k_existing
            for z in model_data.index_z
                # need to have constraints in the optimization as well
                if k == Symbol("nuclear")
                    # X_R_cumu_U(y, k, z, :) .= min(ipp.x_E_my(p_star, z, k) - ipp.x_R_cumu(p_star, k, z), 1.0)
                    X_R_cumu_U(y, k, z, :) .= ipp.x_E_my(p_star, z, k) - ipp.x_R_cumu(p_star, k, z)
                else
                    X_R_cumu_U(y, k, z, :) .= ipp.x_E_my(p_star, z, k) - ipp.x_R_cumu(p_star, k, z)
                end
            end
        end
    end
    for y in model_data.index_y
        for s in ipp.index_stor_existing
            for z in model_data.index_z
                # need to have constraints in the optimization as well
                X_R_stor_cumu_U(y, s, z, :) .= ipp.x_stor_E_my(p_star, z, s) - ipp.x_stor_R_cumu(p_star, s, z)
            end
        end
    end

    X_C_cumu_L = initialize_param("x_C_L", model_data.index_y, ipp.index_k_new, model_data.index_z)
    X_C_cumu_U = initialize_param("x_C_U", model_data.index_y, ipp.index_k_new, model_data.index_z)
    X_C_stor_cumu_L = initialize_param("x_C_stor_L", model_data.index_y, ipp.index_stor_new, model_data.index_z)
    X_C_stor_cumu_U = initialize_param("x_C_stor_U", model_data.index_y, ipp.index_stor_new, model_data.index_z)
    for y in model_data.index_y
        for k in ipp.index_k_new
            for z in model_data.index_z
                # need to have constraints in the optimization as well
                # this hard-coded number needs to change
                X_C_cumu_U(y, k, z, :) .= 5000.0
            end
        end
    end
    for y in model_data.index_y
        for s in ipp.index_stor_new
            for z in model_data.index_z
                # need to have constraints in the optimization as well
                # this hard-coded number needs to change
                X_C_stor_cumu_U(y, s, z, :) .= 5000.0
            end
        end
    end

    # Define positive variables
    @variable(WMDER_IPP, x_C[model_data.index_y, ipp.index_k_new, model_data.index_z] >= 0)
    @variable(WMDER_IPP, x_R[model_data.index_y, ipp.index_k_existing, model_data.index_z] >= 0)
    @variable(WMDER_IPP, x_stor_C[model_data.index_y, ipp.index_stor_new, model_data.index_z] >= 0)
    @variable(WMDER_IPP, x_stor_R[model_data.index_y, ipp.index_stor_existing, model_data.index_z] >= 0)

    @variable(
        WMDER_IPP,
        y_E[model_data.index_y, ipp.index_p, ipp.index_k_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        WMDER_IPP,
        y_C[model_data.index_y, ipp.index_p, ipp.index_k_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(WMDER_IPP, miu[model_data.index_y, model_data.index_z, model_data.index_d, model_data.index_t])
    @variable(
        WMDER_IPP,
        eta[model_data.index_y, ipp.index_p, ipp.index_k_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        WMDER_IPP,
        lambda[model_data.index_y, ipp.index_p, ipp.index_k_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    # variables that represent bilibear terms
    @variable(
        WMDER_IPP,
        mce_eta_x_R_p_star[model_data.index_y, ipp.index_k_existing, model_data.index_z, model_data.index_d, model_data.index_t],
    )
    @variable(
        WMDER_IPP,
        mce_lambda_x_C_p_star[model_data.index_y, ipp.index_k_new, model_data.index_z, model_data.index_d, model_data.index_t],
    )
    @variable(
        WMDER_IPP,
        mce_theta_E_energy_x_stor_R_p_star[model_data.index_y, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t],
    )
    @variable(
        WMDER_IPP,
        mce_theta_E_discharge_x_stor_R_p_star[model_data.index_y, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t],
    )
    @variable(
        WMDER_IPP,
        mce_theta_E_charge_x_stor_R_p_star[model_data.index_y, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t],
    )
    @variable(
        WMDER_IPP,
        mce_pi_E_charge_x_stor_R_p_star[model_data.index_y, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t],
    )
    @variable(
        WMDER_IPP,
        mce_kappa_E_x_stor_R_p_star[model_data.index_y, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t],
    )

    @variable(
        WMDER_IPP,
        mce_theta_C_energy_x_stor_C_p_star[model_data.index_y, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t],
    )
    @variable(
        WMDER_IPP,
        mce_theta_C_discharge_x_stor_C_p_star[model_data.index_y, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t],
    )
    @variable(
        WMDER_IPP,
        mce_theta_C_charge_x_stor_C_p_star[model_data.index_y, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t],
    )
    @variable(
        WMDER_IPP,
        mce_pi_C_charge_x_stor_C_p_star[model_data.index_y, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t],
    )
    @variable(
        WMDER_IPP,
        mce_kappa_C_x_stor_C_p_star[model_data.index_y, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t],
    )
    # cumulative retirement of p_star
    @variable(WMDER_IPP, x_R_mce[model_data.index_y, ipp.index_k_existing, model_data.index_z] >= 0)
    @variable(WMDER_IPP, x_stor_R_mce[model_data.index_y, ipp.index_stor_existing, model_data.index_z] >= 0)
    # cumulative investment of p_star
    @variable(WMDER_IPP, x_C_mce[model_data.index_y, ipp.index_k_new, model_data.index_z] >= 0)
    @variable(WMDER_IPP, x_stor_C_mce[model_data.index_y, ipp.index_stor_new, model_data.index_z] >= 0)

    # adjust upper bound of X_C for specific technology
    for y in model_data.index_y, z in model_data.index_z, k in ipp.index_k_existing
        set_upper_bound(x_R_mce[y, k, z], X_R_cumu_U(y, k, z))
    end
    for y in model_data.index_y, z in model_data.index_z, s in ipp.index_stor_existing
        set_upper_bound(x_stor_R_mce[y, s, z], X_R_stor_cumu_U(y, s, z))
    end

    for y in model_data.index_y, z in model_data.index_z, k in ipp.index_k_new
        set_upper_bound(x_C_mce[y, k, z], X_C_cumu_U(y, k, z))
    end
    for y in model_data.index_y, z in model_data.index_z, s in ipp.index_stor_new
        set_upper_bound(x_stor_C_mce[y, s, z], X_C_stor_cumu_U(y, s, z))
    end

    # the constraint that makes multi-year look-ahead hard to solve
    for y in model_data.index_y, z in model_data.index_z
        fix(x_R[y, :nuclear, z], 0.0; force = true)
        fix(x_R[y, :dera_pv, z], 0.0; force = true)
        fix(x_stor_R[y, :der_aggregator, z], 0.0; force = true)
    end

    constraint_scaling = 1.0

    @variable(
        WMDER_IPP,
        charge_E[model_data.index_y, ipp.index_p, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        WMDER_IPP,
        discharge_E[model_data.index_y, ipp.index_p, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        WMDER_IPP,
        charge_C[model_data.index_y, ipp.index_p, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        WMDER_IPP,
        discharge_C[model_data.index_y, ipp.index_p, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        WMDER_IPP,
        energy_E[model_data.index_y, ipp.index_p, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        WMDER_IPP,
        energy_C[model_data.index_y, ipp.index_p, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        WMDER_IPP,
        flow[model_data.index_y, ipp.index_l, model_data.index_d, model_data.index_t]
    )

    @variable(
        WMDER_IPP,
        iota_min[model_data.index_y, ipp.index_l, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        WMDER_IPP,
        iota_max[model_data.index_y, ipp.index_l, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        WMDER_IPP,
        psi_E[model_data.index_y, ipp.index_p, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t]
    )
    @variable(
        WMDER_IPP,
        theta_E_energy[model_data.index_y, ipp.index_p, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        WMDER_IPP,
        theta_E_discharge[model_data.index_y, ipp.index_p, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        WMDER_IPP,
        theta_E_charge[model_data.index_y, ipp.index_p, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        WMDER_IPP,
        pi_E_discharge[model_data.index_y, ipp.index_p, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        WMDER_IPP,
        pi_E_charge[model_data.index_y, ipp.index_p, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        WMDER_IPP,
        kappa_E[model_data.index_y, ipp.index_p, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        WMDER_IPP,
        psi_C[model_data.index_y, ipp.index_p, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t]
    )
    @variable(
        WMDER_IPP,
        theta_C_energy[model_data.index_y, ipp.index_p, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        WMDER_IPP,
        theta_C_discharge[model_data.index_y, ipp.index_p, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        WMDER_IPP,
        theta_C_charge[model_data.index_y, ipp.index_p, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        WMDER_IPP,
        pi_C_discharge[model_data.index_y, ipp.index_p, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        WMDER_IPP,
        pi_C_charge[model_data.index_y, ipp.index_p, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        WMDER_IPP,
        kappa_C[model_data.index_y, ipp.index_p, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )

    @variable(WMDER_IPP, UCAP_p_star[model_data.index_y])

    @variable(
        WMDER_IPP,
        flow_cap[model_data.index_y, ipp.index_l]
    )

    # TEMPORARY CODE FOR TESTING
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
    for y in model_data.index_y, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
        ipp.Net_Load_my(y, z, d, t, :) .=
            sum(customers.gamma(z, h) * customers.d_my(y, h, z, d, t) for h in model_data.index_h) - 
            # total DG generation at time t
            sum(
                customers.rho_DG(h, m, z, d, t) * customers.total_der_capacity_my(y, z, h, m) for
                h in model_data.index_h, m in customers.index_m
            ) +
            # remove aggregated behind-the-meter pv/storage generation/consumption since they're front-of-the-meter now
            sum(
                customers.rho_DG(h, m, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) *
                customers.total_pv_stor_capacity_my(y, z, h, m) for h in model_data.index_h, m in (:BTMStorage, :BTMPV)
            )
    end
    fill!(ipp.Max_Net_Load_my, NaN)
    Max_Net_Load_my_dict = Dict()
    for y in model_data.index_y, z in model_data.index_z
        ipp.Max_Net_Load_my(y, z, :) .=
            findmax(Dict((d, t) => ipp.Net_Load_my(y, z, d, t) for d in model_data.index_d, t in model_data.index_t))[1]
        push!(
            Max_Net_Load_my_dict,
            (y, z) => findmax(Dict((d, t) => ipp.Net_Load_my(y, z, d, t) for d in model_data.index_d, t in model_data.index_t))[2],
        )
    end
    ipp.Max_Net_Load_my_dict = Max_Net_Load_my_dict

    # commented out for now, use default capacity credit
    # Max_Net_Load_my_index = KeyedArray(
    #     [
    #         findmax(Dict(t => ipp.Net_Load_my(y, t) for t in model_data.index_t))[2] for
    #         y in model_data.index_y
    #     ];
    #     [get_pair(model_data.index_y)]...
    # )

    # fill!(ipp.capacity_credit_E_my, NaN)
    # for y in model_data.index_y, k in ipp.index_k_existing
    #     ipp.capacity_credit_E_my(y, k, :) .= ipp.rho_E_my(p_star, k, Max_Net_Load_my_index(y))
    # end
    # fill!(ipp.capacity_credit_C_my, NaN)
    # for y in model_data.index_y, k in ipp.index_k_new
    #     ipp.capacity_credit_C_my(y, k, :) .= ipp.rho_C_my(p_star, k, Max_Net_Load_my_index(y))
    # end
    fill!(ipp.Reserve_req_my, NaN)
    for y in model_data.index_y, z in model_data.index_z
        # Reserve_req_my only includes net load (no green developers' buildouts and exogenous export/import)
        # Green developers' buildouts and exogenous export/import are included on the supply-side of capacity markets
        # Since we model capacity markets as a whole region, no endogenous export/import is considered
        ipp.Reserve_req_my(y, z, :) .= (1.0 + regulator.r(z, y)) * ipp.Max_Net_Load_my(y, z)
    end

    for y in model_data.index_y
        ipp.Capacity_slope_my(y, :) .=
            -ipp.NetCONE(y) / (ipp.DC_length(y) * sum(ipp.Reserve_req_my(y, z) for z in model_data.index_z))
        ipp.Capacity_intercept_my(y, :) .=
            -ipp.Capacity_slope_my(y) * sum(ipp.Reserve_req_my(y, z) for z in model_data.index_z) + ipp.NetCONE(y)
    end

    #####################################################################

    UCAP_p_star =
        y -> begin
            sum(
                ipp.capacity_credit_E_my(y, z, k) * (
                    ipp.x_E_my(p_star, z, k) - sum(
                        x_R[Symbol(Int(y_symbol)), k, z] for y_symbol in
                        model_data.year(first(model_data.index_y)):model_data.year(y)
                    ) - ipp.x_R_cumu(p_star, k, z)
                ) for k in ipp.index_k_existing, z in model_data.index_z
            ) + sum(
                ipp.capacity_credit_C_my(y, z, k) * (
                    sum(
                        x_C[Symbol(Int(y_symbol)), k, z] for y_symbol in
                        model_data.year(first(model_data.index_y)):model_data.year(y)
                    ) + ipp.x_C_cumu(p_star, k, z)
                ) for k in ipp.index_k_new, z in model_data.index_z
            ) + sum(
                ipp.capacity_credit_stor_E_my(y, z, s) * (
                    ipp.x_stor_E_my(p_star, z, s) - sum(
                        x_stor_R[Symbol(Int(y_symbol)), s, z] for y_symbol in
                        model_data.year(first(model_data.index_y)):model_data.year(y)
                    ) - ipp.x_stor_R_cumu(p_star, s, z)
                ) for s in ipp.index_stor_existing, z in model_data.index_z
            ) + sum(
                ipp.capacity_credit_stor_C_my(y, z, s) * (
                    sum(
                        x_stor_C[Symbol(Int(y_symbol)), s, z] for y_symbol in
                        model_data.year(first(model_data.index_y)):model_data.year(y)
                    ) + ipp.x_stor_C_cumu(p_star, s, z)
                ) for s in ipp.index_stor_new, z in model_data.index_z
            )
        end

    # ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p))) doesn't seem to work, maybe change to filter(x -> x != p_star, ipp.index_p.elements)
    if length(ipp.index_p) >= 2
        UCAP_total =
            y -> begin
                UCAP_p_star(y) +
                sum(
                    ipp.capacity_credit_E_my(y, z, k) * (
                        ipp.x_E_my(p, z, k) - sum(
                            ipp.x_R_my(Symbol(Int(y_symbol)), p, k, z) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_R_cumu(p, k, z)
                    ) for k in ipp.index_k_existing,
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p))), 
                    z in model_data.index_z
                ) +
                sum(
                    ipp.capacity_credit_C_my(y, z, k) * (
                        sum(
                            ipp.x_C_my(Symbol(Int(y_symbol)), p, k, z) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p, k, z)
                    ) for k in ipp.index_k_new,
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p))), 
                    z in model_data.index_z
                ) +
                sum(
                    ipp.capacity_credit_stor_E_my(y, z, s) * (
                        ipp.x_stor_E_my(p, z, s) - sum(
                            ipp.x_stor_R_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_stor_R_cumu(p, s, z)
                    ) for s in ipp.index_stor_existing,
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p))), 
                    z in model_data.index_z
                ) +
                sum(
                    ipp.capacity_credit_stor_C_my(y, z, s) * (
                        sum(
                            ipp.x_stor_C_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_stor_C_cumu(p, s, z)
                    ) for s in ipp.index_stor_new,
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p))), 
                    z in model_data.index_z
                ) +
                # green technology subscription
                sum(
                    ipp.capacity_credit_C_my(y, z, j) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, z, h) for y_symbol in
                    model_data.year(first(model_data.index_y_fix)):model_data.year(y))
                    for j in model_data.index_j, h in model_data.index_h, z in model_data.index_z
                ) - 
                # put exogenous export on the supply-side
                # don't have endogenous export/import because the capacity market clearing here assumes the entire region
                sum(ipp.eximport_my(y, z, Max_Net_Load_my_dict[y, z][1], Max_Net_Load_my_dict[y, z][2]) for z in model_data.index_z)
            end
    else
        UCAP_total = y -> begin
            UCAP_p_star(y) +
            # green technology subscription
            sum(
                ipp.capacity_credit_C_my(y, z, j) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, z, h) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):model_data.year(y))
                for j in model_data.index_j, h in model_data.index_h, z in model_data.index_z
            ) - 
            # put exogenous export on the supply-side
            # don't have endogenous export/import because the capacity market clearing here assumes the entire region
            sum(ipp.eximport_my(y, z, Max_Net_Load_my_dict[y, z][1], Max_Net_Load_my_dict[y, z][2]) for z in model_data.index_z)
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
                ipp.fom_E_my(y, p_star, z, k) * (
                    - sum(
                        x_R[Symbol(Int(y_symbol)), k, z] for y_symbol in
                        model_data.year(first(model_data.index_y)):model_data.year(y)
                    )
                ) for k in ipp.index_k_existing, z in model_data.index_z
            ) -
            # fixed costs
            #   fom * cap new for gen type
            sum(
                ipp.fom_C_my(y, p_star, z, k) *
                x_C[y, k, z] *
                sum(
                    ipp.pvf_onm(Symbol(Int(y_symbol)), p_star) for y_symbol in
                    model_data.year(y):model_data.year(last(model_data.index_y.elements))
                ) for k in ipp.index_k_new, z in model_data.index_z
            ) -
            # fixed costs
            #   capex * cap new for gen type
            ipp.pvf_cap(y, p_star) *
            sum(ipp.CapEx_my(y, p_star, z, k) * x_C[y, k, z] for k in ipp.index_k_new, z in model_data.index_z) -
            # fixed costs
            #   fom * (cap exist - cap retiring) for stor type
            ipp.pvf_onm(y, p_star) * sum(
                ipp.fom_stor_E_my(y, p_star, z, s) * (
                    - sum(
                        x_stor_R[Symbol(Int(y_symbol)), s, z] for y_symbol in
                        model_data.year(first(model_data.index_y)):model_data.year(y)
                    )
                ) for s in ipp.index_stor_existing, z in model_data.index_z
            ) -
            # fixed costs
            #   fom * cap new for stor type
            sum(
                ipp.fom_stor_C_my(y, p_star, z, s) *
                x_stor_C[y, s, z] *
                sum(
                    ipp.pvf_onm(Symbol(Int(y_symbol)), p_star) for y_symbol in
                    model_data.year(y):model_data.year(last(model_data.index_y.elements))
                ) for s in ipp.index_stor_new, z in model_data.index_z
            ) -
            # fixed costs
            #   capex * cap new for stor type
            ipp.pvf_cap(y, p_star) *
            sum(ipp.CapEx_stor_my(y, p_star, z, s) * x_stor_C[y, s, z] for s in ipp.index_stor_new, z in model_data.index_z) +
            # Linearized profit term
            ipp.pvf_onm(y, p_star) * (
                sum(
                    miu[y, z, d, t] * (
                        sum(
                            customers.gamma(z, h) * customers.d_my(y, h, z, d, t) for
                            h in model_data.index_h
                        ) + 
                        ipp.eximport_my(y, z, d, t) -
                        # total DG generation at time t
                        sum(
                            customers.rho_DG(h, m, z, d, t) * customers.total_der_capacity_my(y, z, h, m) for
                            h in model_data.index_h, m in customers.index_m
                        ) +
                        # remove aggregated behind-the-meter pv/storage generation/consumption since they're front-of-the-meter now
                        sum(
                            customers.rho_DG(h, m, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) * 
                            customers.total_pv_stor_capacity_my(y, z, h, m) for h in model_data.index_h, m in (:BTMStorage, :BTMPV)
                        ) -
                        # green technology subscription at time t
                        sum(
                            ipp.rho_C_my(Symbol("ipp1"), j, z, d, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, z, h) for y_symbol in
                            model_data.year(first(model_data.index_y_fix)):model_data.year(y))
                            for j in model_data.index_j, h in model_data.index_h
                        )
                    ) for z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) - ( # the reason I choose to do it this way (remove p_star from all p) is because there may be issue when p_star is the only IPP
                    sum(
                        eta[y, p, k, z, d, t] *
                        ipp.rho_E_my(p, k, z, d, t) * (
                            ipp.x_E_my(p, z, k) - sum(
                                ipp.x_R_my(Symbol(Int(y_symbol)), p, k, z) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) - ipp.x_R_cumu(p, k, z)
                        ) for k in ipp.index_k_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t, 
                        p in ipp.index_p
                    ) - sum(
                        eta[y, p_star, k, z, d, t] *
                        ipp.rho_E_my(p_star, k, z, d, t) * (
                            ipp.x_E_my(p_star, z, k) - sum(
                                ipp.x_R_my(Symbol(Int(y_symbol)), p_star, k, z) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) - ipp.x_R_cumu(p_star, k, z)
                        ) for k in ipp.index_k_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                    )
                ) - (
                    sum(
                        lambda[y, p, k, z, d, t] *
                        ipp.rho_C_my(p, k, z, d, t) * (
                            sum(
                                ipp.x_C_my(Symbol(Int(y_symbol)), p, k, z) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) + ipp.x_C_cumu(p, k, z)
                        ) for k in ipp.index_k_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t, 
                        p in ipp.index_p
                    ) - sum(
                        lambda[y, p_star, k, z, d, t] *
                        ipp.rho_C_my(p_star, k, z, d, t) * (
                            sum(
                                ipp.x_C_my(Symbol(Int(y_symbol)), p_star, k, z) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) + ipp.x_C_cumu(p_star, k, z)
                        ) for k in ipp.index_k_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                    )
                ) + 
                sum(
                    iota_min[y, l, d, t] * ipp.trans_capacity(l, :min)
                    for l in ipp.index_l, d in model_data.index_d, t in model_data.index_t
                ) - 
                sum(
                    iota_max[y, l, d, t] * ipp.trans_capacity(l, :max)
                    for l in ipp.index_l, d in model_data.index_d, t in model_data.index_t
                ) - 
                (
                    sum(
                        theta_E_energy[y, p, s, z, d, t] *
                        ipp.stor_duration_existing(s) * (
                            ipp.x_stor_E_my(p, z, s) - sum(
                                ipp.x_stor_R_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) - ipp.x_stor_R_cumu(p, s, z)
                        ) for s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t, 
                        p in ipp.index_p
                    ) - sum(
                        theta_E_energy[y, p_star, s, z, d, t] *
                        ipp.stor_duration_existing(s) * (
                            ipp.x_stor_E_my(p_star, z, s) - sum(
                                ipp.x_stor_R_my(Symbol(Int(y_symbol)), p_star, s, z) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) - ipp.x_stor_R_cumu(p_star, s, z)
                        ) for s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                    )
                ) - 
                (
                    sum(
                        theta_E_discharge[y, p, s, z, d, t] *
                        ipp.rte_stor_E_my(y, p, z, s) * (
                            ipp.x_stor_E_my(p, z, s) - sum(
                                ipp.x_stor_R_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) - ipp.x_stor_R_cumu(p, s, z)
                        ) for s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t, 
                        p in ipp.index_p
                    ) - sum(
                        theta_E_discharge[y, p_star, s, z, d, t] *
                        ipp.rte_stor_E_my(y, p_star, z, s) * (
                            ipp.x_stor_E_my(p_star, z, s) - sum(
                                ipp.x_stor_R_my(Symbol(Int(y_symbol)), p_star, s, z) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) - ipp.x_stor_R_cumu(p_star, s, z)
                        ) for s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                    )
                ) - 
                (
                    sum(
                        theta_E_charge[y, p, s, z, d, t] * (
                            ipp.x_stor_E_my(p, z, s) - sum(
                                ipp.x_stor_R_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) - ipp.x_stor_R_cumu(p, s, z)
                        ) for s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t, 
                        p in ipp.index_p
                    ) - sum(
                        theta_E_charge[y, p_star, s, z, d, t] * (
                            ipp.x_stor_E_my(p_star, z, s) - sum(
                                ipp.x_stor_R_my(Symbol(Int(y_symbol)), p_star, s, z) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) - ipp.x_stor_R_cumu(p_star, s, z)
                        ) for s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                    )
                ) - 
                (
                    sum(
                        pi_E_charge[y, p, s, z, d, t] *
                        ipp.stor_duration_existing(s) * (
                            ipp.x_stor_E_my(p, z, s) - sum(
                                ipp.x_stor_R_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) - ipp.x_stor_R_cumu(p, s, z)
                        ) for s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t, 
                        p in ipp.index_p
                    ) - 
                    sum(
                        pi_E_charge[y, p_star, s, z, d, t] *
                        ipp.stor_duration_existing(s) * (
                            ipp.x_stor_E_my(p_star, z, s) - sum(
                                ipp.x_stor_R_my(Symbol(Int(y_symbol)), p_star, s, z) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) - ipp.x_stor_R_cumu(p_star, s, z)
                        ) for s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                    )
                ) -
                (
                    sum(
                        kappa_E[y, p, s, z, d, t] * (
                            ipp.x_stor_E_my(p, z, s) - sum(
                                ipp.x_stor_R_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) - ipp.x_stor_R_cumu(p, s, z)
                        ) for s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t, 
                        p in ipp.index_p
                    ) - 
                    sum(
                        kappa_E[y, p_star, s, z, d, t] * (
                            ipp.x_stor_E_my(p_star, z, s) - sum(
                                ipp.x_stor_R_my(Symbol(Int(y_symbol)), p_star, s, z) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) - ipp.x_stor_R_cumu(p_star, s, z)
                        ) for s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                    )
                ) - 
                (
                    sum(
                        theta_C_energy[y, p, s, z, d, t] *
                        ipp.stor_duration_new(s) * (
                            sum(
                                ipp.x_stor_C_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) + ipp.x_stor_C_cumu(p, s, z)
                        ) for s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t, 
                        p in ipp.index_p
                    ) - 
                    sum(
                        theta_C_energy[y, p_star, s, z, d, t] *
                        ipp.stor_duration_new(s) * (
                            sum(
                                ipp.x_stor_C_my(Symbol(Int(y_symbol)), p_star, s, z) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) + ipp.x_stor_C_cumu(p_star, s, z)
                        ) for s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                    )
                ) - 
                (
                    sum(
                        theta_C_discharge[y, p, s, z, d, t] *
                        ipp.rte_stor_C_my(y, p, z, s) * (
                            sum(
                                ipp.x_stor_C_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) + ipp.x_stor_C_cumu(p, s, z)
                        ) for s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t, 
                        p in ipp.index_p
                    ) - 
                    sum(
                        theta_C_discharge[y, p_star, s, z, d, t] *
                        ipp.rte_stor_C_my(y, p_star, z, s) * (
                            sum(
                                ipp.x_stor_C_my(Symbol(Int(y_symbol)), p_star, s, z) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) + ipp.x_stor_C_cumu(p_star, s, z)
                        ) for s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t 
                    )
                ) - 
                (
                    sum(
                        theta_C_charge[y, p, s, z, d, t] * (
                            sum(
                                ipp.x_stor_C_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) + ipp.x_stor_C_cumu(p, s, z)
                        ) for s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t, 
                        p in ipp.index_p
                    ) - 
                    sum(
                        theta_C_charge[y, p_star, s, z, d, t] * (
                            sum(
                                ipp.x_stor_C_my(Symbol(Int(y_symbol)), p_star, s, z) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) + ipp.x_stor_C_cumu(p_star, s, z)
                        ) for s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t 
                    )
                ) - 
                (
                    sum(
                        pi_C_charge[y, p, s, z, d, t] *
                        ipp.stor_duration_new(s) * (
                            sum(
                                ipp.x_stor_C_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) + ipp.x_stor_C_cumu(p, s, z)
                        ) for s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t, 
                        p in ipp.index_p
                    ) - 
                    sum(
                        pi_C_charge[y, p_star, s, z, d, t] *
                        ipp.stor_duration_new(s) * (
                            sum(
                                ipp.x_stor_C_my(Symbol(Int(y_symbol)), p_star, s, z) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) + ipp.x_stor_C_cumu(p_star, s, z)
                        ) for s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                    )
                ) - 
                (
                    sum(
                        kappa_C[y, p, s, z, d, t] * (
                            sum(
                                ipp.x_stor_C_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) + ipp.x_stor_C_cumu(p, s, z)
                        ) for s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t, 
                        p in ipp.index_p
                    ) - 
                    sum(
                        kappa_C[y, p_star, s, z, d, t] * (
                            sum(
                                ipp.x_stor_C_my(Symbol(Int(y_symbol)), p_star, s, z) for y_symbol in
                                model_data.year(first(model_data.index_y)):model_data.year(y)
                            ) + ipp.x_stor_C_cumu(p_star, s, z)
                        ) for s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t 
                    )
                ) -
                # part of the item below is from revenue linearization (the part without p_star), part is from cost of p_star
                (
                    sum(
                        model_data.omega(d) * delta_t *
                        (ipp.v_E_my(y, p, k, z, d, t) * y_E[y, p, k, z, d, t]) for
                        d in model_data.index_d, t in model_data.index_t, z in model_data.index_z, k in ipp.index_k_existing, p in ipp.index_p
                    ) + 
                    sum(
                        model_data.omega(d) * delta_t *
                        (ipp.v_C_my(y, p, k, z, d, t) * y_C[y, p, k, z, d, t]) for
                        d in model_data.index_d, t in model_data.index_t, z in model_data.index_z, k in ipp.index_k_new, p in ipp.index_p
                    )
                )
            )
            for y in model_data.index_y
        )
    end

    @objective(WMDER_IPP, Max, objective_function)

    # @constraint(
    #     WMDER_IPP,
    #     Eq_fix_retirement_p130[y in model_data.index_y],
    #     x_R[Symbol(y), :nuclear, :p130] <= 10.0
    # )
    # @constraint(
    #     WMDER_IPP,
    #     Eq_fix_retirement_p132[y in model_data.index_y],
    #     x_R[Symbol(y), :nuclear, :p132] <= 10.0
    # )
    
    @constraint(
        WMDER_IPP,
        Eq_sigma[y in model_data.index_y, k in ipp.index_k_existing, z in model_data.index_z],
        ipp.x_E_my(p_star, z, k) - sum(
            x_R[Symbol(Int(y_symbol)), k, z] for
            y_symbol in model_data.year(first(model_data.index_y)):model_data.year(y)
        ) - ipp.x_R_cumu(p_star, k, z) >= 0
    )

    @constraint(
        WMDER_IPP,
        Eq_sigma_stor[y in model_data.index_y, s in ipp.index_stor_existing, z in model_data.index_z],
        ipp.x_stor_E_my(p_star, z, s) - sum(
            x_stor_R[Symbol(Int(y_symbol)), s, z] for
            y_symbol in model_data.year(first(model_data.index_y)):model_data.year(y)
        ) - ipp.x_stor_R_cumu(p_star, s, z) >= 0
    )

    # dual feasible constraints
    @constraint(
        WMDER_IPP,
        Eq_dual_feasible_y_E[
            y in model_data.index_y,
            p in ipp.index_p,
            k in ipp.index_k_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        model_data.omega(d) * delta_t * ipp.v_E_my(y, p, k, z, d, t) - miu[y, z, d, t] + eta[y, p, k, z, d, t] >= 0
    )
    @constraint(
        WMDER_IPP,
        Eq_dual_feasible_y_C[
            y in model_data.index_y,
            p in ipp.index_p,
            k in ipp.index_k_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        model_data.omega(d) * delta_t * ipp.v_C_my(y, p, k, z, d, t) - miu[y, z, d, t] + lambda[y, p, k, z, d, t] >= 0
    )
    @constraint(
        WMDER_IPP,
        Eq_dual_feasible_flow[
            y in model_data.index_y,
            l in ipp.index_l,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        sum(miu[y, z, d, t] * ipp.trans_topology(l, z) for z in model_data.index_z) - iota_min[y, l, d, t] + iota_max[y, l, d, t] == 0
    )
    @constraint(
        WMDER_IPP,
        Eq_dual_feasible_discharge_E[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        - miu[y, z, d, t] - psi_E[y, p, s, z, d, t] / ipp.rte_stor_E_my(y, p, z, s) * delta_t + theta_E_discharge[y, p, s, z, d, t] + pi_E_discharge[y, p, s, z, d, t] * delta_t + kappa_E[y, p, s, z, d, t] / ipp.rte_stor_E_my(y, p, z, s) >= 0
    )
    @constraint(
        WMDER_IPP,
        Eq_dual_feasible_charge_E[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        miu[y, z, d, t] + psi_E[y, p, s, z, d, t] * delta_t + theta_E_charge[y, p, s, z, d, t] + pi_E_charge[y, p, s, z, d, t] * delta_t + kappa_E[y, p, s, z, d, t] >= 0
    )
    @constraint(
        WMDER_IPP,
        Eq_dual_feasible_energy_E[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t.elements[1:end-1],
        ],
        psi_E[y, p, s, z, d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)+delta_t), model_data.time.values)][1]] - psi_E[y, p, s, z, d, t] + theta_E_energy[y, p, s, z, d, t] 
        - pi_E_discharge[y, p, s, z, d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)+delta_t), model_data.time.values)][1]] 
        + pi_E_charge[y, p, s, z, d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)+delta_t), model_data.time.values)][1]] >= 0
    )
    @constraint(
        WMDER_IPP,
        Eq_dual_feasible_energy_last_E[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in [model_data.index_t.elements[end]],
        ],
        - psi_E[y, p, s, z, d, t] + theta_E_energy[y, p, s, z, d, t] >= 0
    )
    @constraint(
        WMDER_IPP,
        Eq_dual_feasible_discharge_C[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        - miu[y, z, d, t] - psi_C[y, p, s, z, d, t] / ipp.rte_stor_C_my(y, p, z, s) * delta_t + theta_C_discharge[y, p, s, z, d, t] + pi_C_discharge[y, p, s, z, d, t] * delta_t + kappa_C[y, p, s, z, d, t] / ipp.rte_stor_C_my(y, p, z, s) >= 0
    )
    @constraint(
        WMDER_IPP,
        Eq_dual_feasible_charge_C[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        miu[y, z, d, t] + psi_C[y, p, s, z, d, t] * delta_t + theta_C_charge[y, p, s, z, d, t] + pi_C_charge[y, p, s, z, d, t] * delta_t + kappa_C[y, p, s, z, d, t] >= 0
    )
    @constraint(
        WMDER_IPP,
        Eq_dual_feasible_energy_C[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t.elements[1:end-1],
        ],
        psi_C[y, p, s, z, d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)+delta_t), model_data.time.values)][1]] - psi_C[y, p, s, z, d, t] + theta_C_energy[y, p, s, z, d, t] 
        - pi_C_discharge[y, p, s, z, d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)+delta_t), model_data.time.values)][1]] 
        + pi_C_charge[y, p, s, z, d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)+delta_t), model_data.time.values)][1]] >= 0
    )
    @constraint(
        WMDER_IPP,
        Eq_dual_feasible_energy_last_C[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in [model_data.index_t.elements[end]],
        ],
        - psi_C[y, p, s, z, d, t] + theta_C_energy[y, p, s, z, d, t] >= 0
    )

    # primal feasible constraints
    supply_demand_balance =
        (y, z, d, t) -> begin
            # bulk generation at time t
            sum(y_E[y, p, k, z, d, t] for k in ipp.index_k_existing, p in ipp.index_p) +
            sum(y_C[y, p, k, z, d, t] for k in ipp.index_k_new, p in ipp.index_p) -
            # flow out of zone z
            sum(ipp.trans_topology(l, z) * flow[y, l, d, t] for l in ipp.index_l) +
            # battery discharge
            sum(discharge_E[y, p, s, z, d, t] for s in ipp.index_stor_existing, p in ipp.index_p) +
            sum(discharge_C[y, p, s, z, d, t] for s in ipp.index_stor_new, p in ipp.index_p) -
            # battery charge
            sum(charge_E[y, p, s, z, d, t] for s in ipp.index_stor_existing, p in ipp.index_p) -
            sum(charge_C[y, p, s, z, d, t] for s in ipp.index_stor_new, p in ipp.index_p) -
            # demand at time t
            sum(customers.gamma(z, h) * customers.d_my(y, h, z, d, t) for h in model_data.index_h) - ipp.eximport_my(y, z, d, t) +
            # total DG generation at time t
            sum(
                customers.rho_DG(h, m, z, d, t) * customers.total_der_capacity_my(y, z, h, m) for
                h in model_data.index_h, m in customers.index_m
            ) -
            # remove aggregated behind-the-meter pv/storage generation/consumption since they're front-of-the-meter now
            sum(
                customers.rho_DG(h, m, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) *
                customers.total_pv_stor_capacity_my(y, z, h, m) for h in model_data.index_h, m in (:BTMStorage, :BTMPV)
            ) +
            # green technology subscription at time t
            sum(
                ipp.rho_C_my(Symbol("ipp1"), j, z, d, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, z, h) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):model_data.year(y))
                for j in model_data.index_j, h in model_data.index_h
            )
        end

    @constraint(
        WMDER_IPP,
        Eq_primal_feasible_supplydemandbalance[y in model_data.index_y, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t],
        supply_demand_balance(y, z, d, t) == 0
    )

    @constraint(
        WMDER_IPP,
        Eq_primal_feasible_gen_max_E[
            y in model_data.index_y,
            k in ipp.index_k_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        ipp.rho_E_my(p_star, k, z, d, t) * (
            ipp.x_E_my(p_star, z, k) - sum(
                x_R[Symbol(Int(y_symbol)), k, z] for y_symbol in
                model_data.year(first(model_data.index_y)):model_data.year(y)
            ) - ipp.x_R_cumu(p_star, k, z)
        ) - y_E[y, p_star, k, z, d, t] >= 0
    )
    if length(ipp.index_p) >= 2
        @constraint(
            WMDER_IPP,
            Eq_primal_feasible_other_gen_max_E[
                y in model_data.index_y,
                p in filter(x -> x != p_star, ipp.index_p.elements),
                k in ipp.index_k_existing,
                z in model_data.index_z,
                d in model_data.index_d,
                t in model_data.index_t,
            ],
            ipp.rho_E_my(p, k, z, d, t) * (
                ipp.x_E_my(p, z, k) - sum(
                    ipp.x_R_my(Symbol(Int(y_symbol)), p, k, z) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) - ipp.x_R_cumu(p, k, z)
            ) - y_E[y, p, k, z, d, t] >= 0
        )
    end

    @constraint(
        WMDER_IPP,
        Eq_primal_feasible_gen_max_C[
            y in model_data.index_y,
            k in ipp.index_k_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        ipp.rho_C_my(p_star, k, z, d, t) * (
            sum(
                x_C[Symbol(Int(y_symbol)), k, z] for y_symbol in
                model_data.year(first(model_data.index_y)):model_data.year(y)
            ) + ipp.x_C_cumu(p_star, k, z)
        ) - y_C[y, p_star, k, z, d, t] >= 0
    )
    if length(ipp.index_p) >= 2
        @constraint(
            WMDER_IPP,
            Eq_primal_feasible_other_gen_max_C[
                y in model_data.index_y,
                p in filter(x -> x != p_star, ipp.index_p.elements),
                k in ipp.index_k_new,
                z in model_data.index_z,
                d in model_data.index_d,
                t in model_data.index_t,
            ],
            ipp.rho_C_my(p, k, z, d, t) * (
                sum(
                    ipp.x_C_my(Symbol(Int(y_symbol)), p, k, z) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) + ipp.x_C_cumu(p, k, z)
            ) - y_C[y, p, k, z, d, t] >= 0
        )
    end
    
    @constraint(
        WMDER_IPP,
        Eq_primal_feasible_flow_lower[
            y in model_data.index_y,
            l in ipp.index_l,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        flow[y, l, d, t] - ipp.trans_capacity(l, :min) >= 0
    )

    @constraint(
        WMDER_IPP,
        Eq_primal_feasible_flow_upper[
            y in model_data.index_y,
            l in ipp.index_l,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        flow[y, l, d, t] - ipp.trans_capacity(l, :max) <= 0
    )

    @constraint(
        WMDER_IPP,
        Eq_primal_feasible_energy_E[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t.elements[2:end],
        ],
        energy_E[y, p, s, z, d, t] == energy_E[y, p, s, z, d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)-delta_t), model_data.time.values)][1]] - discharge_E[y, p, s, z, d, t] / ipp.rte_stor_E_my(y, p, z, s) * delta_t +
            charge_E[y, p, s, z, d, t] * delta_t
    )

    @constraint(
        WMDER_IPP,
        Eq_primal_feasible_energy_E_0[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in [model_data.index_t.elements[1]],
        ],
        energy_E[y, p, s, z, d, t] == ipp.initial_energy_existing_my(y, p, s, z, d) - discharge_E[y, p, s, z, d, t] / ipp.rte_stor_E_my(y, p, z, s) * delta_t +
            charge_E[y, p, s, z, d, t] * delta_t
    )

    @constraint(
        WMDER_IPP,
        Eq_primal_feasible_energy_upper_bound_E[
            y in model_data.index_y,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        energy_E[y, p_star, s, z, d, t] <= ipp.stor_duration_existing(s) * (
            ipp.x_stor_E_my(p_star, z, s) - sum(
                x_stor_R[Symbol(Int(y_symbol)), s, z] for y_symbol in
                model_data.year(first(model_data.index_y)):model_data.year(y)
            ) - ipp.x_stor_R_cumu(p_star, s, z)
        )
    )
    if length(ipp.index_p) >= 2
        @constraint(
            WMDER_IPP,
            Eq_primal_feasible_energy_upper_bound_other_E[
                y in model_data.index_y,
                p in filter(x -> x != p_star, ipp.index_p.elements),
                s in ipp.index_stor_existing,
                z in model_data.index_z,
                d in model_data.index_d,
                t in model_data.index_t,
            ],
            energy_E[y, p, s, z, d, t] <= ipp.stor_duration_existing(s) * (
                ipp.x_stor_E_my(p, z, s) - sum(
                    ipp.x_stor_R_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) - ipp.x_stor_R_cumu(p, s, z)
            )
        )
    end

    @constraint(
        WMDER_IPP,
        Eq_primal_feasible_discharge_upper_bound_E[
            y in model_data.index_y,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        discharge_E[y, p_star, s, z, d, t] <= ipp.rte_stor_E_my(y, p_star, z, s) * (
            ipp.x_stor_E_my(p_star, z, s) - sum(
                x_stor_R[Symbol(Int(y_symbol)), s, z] for y_symbol in
                model_data.year(first(model_data.index_y)):model_data.year(y)
            ) - ipp.x_stor_R_cumu(p_star, s, z)
        )
    )
    if length(ipp.index_p) >= 2
        @constraint(
            WMDER_IPP,
            Eq_primal_feasible_discharge_upper_bound_other_E[
                y in model_data.index_y,
                p in filter(x -> x != p_star, ipp.index_p.elements),
                s in ipp.index_stor_existing,
                z in model_data.index_z,
                d in model_data.index_d,
                t in model_data.index_t,
            ],
            discharge_E[y, p, s, z, d, t] <= ipp.rte_stor_E_my(y, p, z, s) * (
                ipp.x_stor_E_my(p, z, s) - sum(
                    ipp.x_stor_R_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) - ipp.x_stor_R_cumu(p, s, z)
            )
        )
    end

    @constraint(
        WMDER_IPP,
        Eq_primal_feasible_charge_upper_bound_E[
            y in model_data.index_y,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        charge_E[y, p_star, s, z, d, t] <= 
            ipp.x_stor_E_my(p_star, z, s) - sum(
                x_stor_R[Symbol(Int(y_symbol)), s, z] for y_symbol in
                model_data.year(first(model_data.index_y)):model_data.year(y)
            ) - ipp.x_stor_R_cumu(p_star, s, z)
    )
    if length(ipp.index_p) >= 2
        @constraint(
            WMDER_IPP,
            Eq_primal_feasible_charge_upper_bound_other_E[
                y in model_data.index_y,
                p in filter(x -> x != p_star, ipp.index_p.elements),
                s in ipp.index_stor_existing,
                z in model_data.index_z,
                d in model_data.index_d,
                t in model_data.index_t,
            ],
            charge_E[y, p, s, z, d, t] <= 
                ipp.x_stor_E_my(p, z, s) - sum(
                    ipp.x_stor_R_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) - ipp.x_stor_R_cumu(p, s, z)
        )
    end

    @constraint(
        WMDER_IPP,
        Eq_primal_feasible_discharge_energy_upper_bound_E[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t.elements[2:end],
        ],
        discharge_E[y, p, s, z, d, t] * delta_t <= 
        energy_E[y, p, s, z, d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)-delta_t), model_data.time.values)][1]]
    )

    @constraint(
        WMDER_IPP,
        Eq_primal_feasible_discharge_energy_upper_bound_E_0[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in [model_data.index_t.elements[1]],
        ],
        discharge_E[y, p, s, z, d, t] * delta_t <= 
        ipp.initial_energy_existing_my(y, p, s, z, d)
    )

    @constraint(
        WMDER_IPP,
        Eq_primal_feasible_charge_energy_upper_bound_E[
            y in model_data.index_y,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t.elements[2:end],
        ],
        charge_E[y, p_star, s, z, d, t] * delta_t <= ipp.stor_duration_existing(s) * (
            ipp.x_stor_E_my(p_star, z, s) - sum(
                x_stor_R[Symbol(Int(y_symbol)), s, z] for y_symbol in
                model_data.year(first(model_data.index_y)):model_data.year(y)
            ) - ipp.x_stor_R_cumu(p_star, s, z)
        ) -
        energy_E[y, p_star, s, z, d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)-delta_t), model_data.time.values)][1]]
    )
    if length(ipp.index_p) >= 2
        @constraint(
            WMDER_IPP,
            Eq_primal_feasible_charge_energy_upper_bound_other_E[
                y in model_data.index_y,
                p in filter(x -> x != p_star, ipp.index_p.elements),
                s in ipp.index_stor_existing,
                z in model_data.index_z,
                d in model_data.index_d,
                t in model_data.index_t.elements[2:end],
            ],
            charge_E[y, p, s, z, d, t] * delta_t <= ipp.stor_duration_existing(s) * (
                ipp.x_stor_E_my(p, z, s) - sum(
                    ipp.x_stor_R_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) - ipp.x_stor_R_cumu(p, s, z)
            ) -
            energy_E[y, p, s, z, d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)-delta_t), model_data.time.values)][1]]
        )
    end

    @constraint(
        WMDER_IPP,
        Eq_primal_feasible_charge_energy_upper_bound_E_0[
            y in model_data.index_y,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in [model_data.index_t.elements[1]],
        ],
        charge_E[y, p_star, s, z, d, t] * delta_t <= ipp.stor_duration_existing(s) * (
            ipp.x_stor_E_my(p_star, z, s) - sum(
                x_stor_R[Symbol(Int(y_symbol)), s, z] for y_symbol in
                model_data.year(first(model_data.index_y)):model_data.year(y)
            ) - ipp.x_stor_R_cumu(p_star, s, z)
        ) -
        ipp.initial_energy_existing_my(y, p_star, s, z, d)
    )
    if length(ipp.index_p) >= 2
        @constraint(
            WMDER_IPP,
            Eq_primal_feasible_charge_energy_upper_bound_other_E_0[
                y in model_data.index_y,
                p in filter(x -> x != p_star, ipp.index_p.elements),
                s in ipp.index_stor_existing,
                z in model_data.index_z,
                d in model_data.index_d,
                t in [model_data.index_t.elements[1]],
            ],
            charge_E[y, p, s, z, d, t] * delta_t <= ipp.stor_duration_existing(s) * (
                ipp.x_stor_E_my(p, z, s) - sum(
                    ipp.x_stor_R_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) - ipp.x_stor_R_cumu(p, s, z)
            ) -
            ipp.initial_energy_existing_my(y, p, s, z, d)
        )
    end

    @constraint(
        WMDER_IPP,
        Eq_primal_feasible_charge_discharge_upper_bound_E[
            y in model_data.index_y,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t.elements,
        ],
        charge_E[y, p_star, s, z, d, t] + discharge_E[y, p_star, s, z, d, t] / ipp.rte_stor_E_my(y, p_star, z, s) <= 
            ipp.x_stor_E_my(p_star, z, s) - sum(
                x_stor_R[Symbol(Int(y_symbol)), s, z] for y_symbol in
                model_data.year(first(model_data.index_y)):model_data.year(y)
            ) - ipp.x_stor_R_cumu(p_star, s, z)
    )
    if length(ipp.index_p) >= 2
        @constraint(
            WMDER_IPP,
            Eq_primal_feasible_charge_discharge_upper_bound_other_E[
                y in model_data.index_y,
                p in filter(x -> x != p_star, ipp.index_p.elements),
                s in ipp.index_stor_existing,
                z in model_data.index_z,
                d in model_data.index_d,
                t in model_data.index_t.elements,
            ],
            charge_E[y, p, s, z, d, t] + discharge_E[y, p, s, z, d, t] / ipp.rte_stor_E_my(y, p, z, s) <= 
                ipp.x_stor_E_my(p, z, s) - sum(
                    ipp.x_stor_R_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) - ipp.x_stor_R_cumu(p, s, z)
        )
    end

    @constraint(
        WMDER_IPP,
        Eq_primal_feasible_energy_C[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t.elements[2:end],
        ],
        energy_C[y, p, s, z, d, t] == energy_C[y, p, s, z, d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)-delta_t), model_data.time.values)][1]] - discharge_C[y, p, s, z, d, t] / ipp.rte_stor_C_my(y, p, z, s) * delta_t +
            charge_C[y, p, s, z, d, t] * delta_t
    )

    @constraint(
        WMDER_IPP,
        Eq_primal_feasible_energy_C_0[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in [model_data.index_t.elements[1]],
        ],
        energy_C[y, p, s, z, d, t] == ipp.initial_energy_new_my(y, p, s, z, d) - discharge_C[y, p, s, z, d, t] / ipp.rte_stor_C_my(y, p, z, s) * delta_t +
            charge_C[y, p, s, z, d, t] * delta_t
    )

    @constraint(
        WMDER_IPP,
        Eq_primal_feasible_energy_upper_bound_C[
            y in model_data.index_y,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        energy_C[y, p_star, s, z, d, t] <= ipp.stor_duration_new(s) * (
            sum(
                x_stor_C[Symbol(Int(y_symbol)), s, z] for y_symbol in
                model_data.year(first(model_data.index_y)):model_data.year(y)
            ) + ipp.x_stor_C_cumu(p_star, s, z)
        )
    )
    if length(ipp.index_p) >= 2
        @constraint(
            WMDER_IPP,
            Eq_primal_feasible_energy_upper_bound_other_C[
                y in model_data.index_y,
                p in filter(x -> x != p_star, ipp.index_p.elements),
                s in ipp.index_stor_new,
                z in model_data.index_z,
                d in model_data.index_d,
                t in model_data.index_t,
            ],
            energy_C[y, p, s, z, d, t] <= ipp.stor_duration_new(s) * (
                sum(
                    ipp.x_stor_C_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) + ipp.x_stor_C_cumu(p, s, z)
            )
        )
    end

    @constraint(
        WMDER_IPP,
        Eq_primal_feasible_discharge_upper_bound_C[
            y in model_data.index_y,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        discharge_C[y, p_star, s, z, d, t] <= ipp.rte_stor_C_my(y, p_star, z, s) * (
            sum(
                x_stor_C[Symbol(Int(y_symbol)), s, z] for y_symbol in
                model_data.year(first(model_data.index_y)):model_data.year(y)
            ) + ipp.x_stor_C_cumu(p_star, s, z)
        )
    )
    if length(ipp.index_p) >= 2
        @constraint(
            WMDER_IPP,
            Eq_primal_feasible_discharge_upper_bound_other_C[
                y in model_data.index_y,
                p in filter(x -> x != p_star, ipp.index_p.elements),
                s in ipp.index_stor_new,
                z in model_data.index_z,
                d in model_data.index_d,
                t in model_data.index_t,
            ],
            discharge_C[y, p, s, z, d, t] <= ipp.rte_stor_C_my(y, p, z, s) * (
                sum(
                    ipp.x_stor_C_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) + ipp.x_stor_C_cumu(p, s, z)
            )
        )
    end

    @constraint(
        WMDER_IPP,
        Eq_primal_feasible_charge_upper_bound_C[
            y in model_data.index_y,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        charge_C[y, p_star, s, z, d, t] <= 
            sum(
                x_stor_C[Symbol(Int(y_symbol)), s, z] for y_symbol in
                model_data.year(first(model_data.index_y)):model_data.year(y)
            ) + ipp.x_stor_C_cumu(p_star, s, z)
    )
    if length(ipp.index_p) >= 2
        @constraint(
            WMDER_IPP,
            Eq_primal_feasible_charge_upper_bound_other_C[
                y in model_data.index_y,
                p in filter(x -> x != p_star, ipp.index_p.elements),
                s in ipp.index_stor_new,
                z in model_data.index_z,
                d in model_data.index_d,
                t in model_data.index_t,
            ],
            charge_C[y, p, s, z, d, t] <= 
                sum(
                    ipp.x_stor_C_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) + ipp.x_stor_C_cumu(p, s, z)
        )
    end

    @constraint(
        WMDER_IPP,
        Eq_primal_feasible_discharge_energy_upper_bound_C[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t.elements[2:end],
        ],
        discharge_C[y, p, s, z, d, t] * delta_t <= 
        energy_C[y, p, s, z, d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)-delta_t), model_data.time.values)][1]]
    )

    @constraint(
        WMDER_IPP,
        Eq_primal_feasible_discharge_energy_upper_bound_C_0[
            y in model_data.index_y,
            p in ipp.index_p,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in [model_data.index_t.elements[1]],
        ],
        discharge_C[y, p, s, z, d, t] * delta_t <= 
        ipp.initial_energy_new_my(y, p, s, z, d)
    )

    @constraint(
        WMDER_IPP,
        Eq_primal_feasible_charge_energy_upper_bound_C[
            y in model_data.index_y,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t.elements[2:end],
        ],
        charge_C[y, p_star, s, z, d, t] * delta_t <= ipp.stor_duration_new(s) * (
            sum(
                x_stor_C[Symbol(Int(y_symbol)), s, z] for y_symbol in
                model_data.year(first(model_data.index_y)):model_data.year(y)
            ) + ipp.x_stor_C_cumu(p_star, s, z)
        ) -
        energy_C[y, p_star, s, z, d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)-delta_t), model_data.time.values)][1]]
    )
    if length(ipp.index_p) >= 2
        @constraint(
            WMDER_IPP,
            Eq_primal_feasible_charge_energy_upper_bound_other_C[
                y in model_data.index_y,
                p in filter(x -> x != p_star, ipp.index_p.elements),
                s in ipp.index_stor_new,
                z in model_data.index_z,
                d in model_data.index_d,
                t in model_data.index_t.elements[2:end],
            ],
            charge_C[y, p, s, z, d, t] * delta_t <= ipp.stor_duration_new(s) * (
                sum(
                    ipp.x_stor_C_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) + ipp.x_stor_C_cumu(p, s, z)
            ) -
            energy_C[y, p, s, z, d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)-delta_t), model_data.time.values)][1]]
        )
    end

    @constraint(
        WMDER_IPP,
        Eq_primal_feasible_charge_energy_upper_bound_C_0[
            y in model_data.index_y,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in [model_data.index_t.elements[1]],
        ],
        charge_C[y, p_star, s, z, d, t] * delta_t <= ipp.stor_duration_new(s) * (
            sum(
                x_stor_C[Symbol(Int(y_symbol)), s, z] for y_symbol in
                model_data.year(first(model_data.index_y)):model_data.year(y)
            ) + ipp.x_stor_C_cumu(p_star, s, z)
        ) -
        ipp.initial_energy_new_my(y, p_star, s, z, d)
    )
    if length(ipp.index_p) >= 2
        @constraint(
            WMDER_IPP,
            Eq_primal_feasible_charge_energy_upper_bound_other_C_0[
                y in model_data.index_y,
                p in filter(x -> x != p_star, ipp.index_p.elements),
                s in ipp.index_stor_new,
                z in model_data.index_z,
                d in model_data.index_d,
                t in [model_data.index_t.elements[1]],
            ],
            charge_C[y, p, s, z, d, t] * delta_t <= ipp.stor_duration_new(s) * (
                sum(
                    ipp.x_stor_C_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) + ipp.x_stor_C_cumu(p, s, z)
            ) -
            ipp.initial_energy_new_my(y, p, s, z, d)
        )
    end

    @constraint(
        WMDER_IPP,
        Eq_primal_feasible_charge_discharge_upper_bound_C[
            y in model_data.index_y,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t.elements,
        ],
        charge_C[y, p_star, s, z, d, t] + discharge_C[y, p_star, s, z, d, t] / ipp.rte_stor_C_my(y, p_star, z, s) <= 
            sum(
                x_stor_C[Symbol(Int(y_symbol)), s, z] for y_symbol in
                model_data.year(first(model_data.index_y)):model_data.year(y)
            ) + ipp.x_stor_C_cumu(p_star, s, z)
    )
    if length(ipp.index_p) >= 2
        @constraint(
            WMDER_IPP,
            Eq_primal_feasible_charge_discharge_upper_bound_other_C[
                y in model_data.index_y,
                p in filter(x -> x != p_star, ipp.index_p.elements),
                s in ipp.index_stor_new,
                z in model_data.index_z,
                d in model_data.index_d,
                t in model_data.index_t.elements,
            ],
            charge_C[y, p, s, z, d, t] + discharge_C[y, p, s, z, d, t] / ipp.rte_stor_C_my(y, p, z, s) <= 
                sum(
                    ipp.x_stor_C_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                    model_data.year(first(model_data.index_y)):model_data.year(y)
                ) + ipp.x_stor_C_cumu(p, s, z)
        )
    end

    # strong-duality constraints
    if length(ipp.index_p) >= 2
        # @constraint(
        #     WMDER_IPP,
        #     Eq_strong_duality[
        #         y in model_data.index_y
        #     ],
        #     sum(
        #         model_data.omega(t) *
        #         (ipp.v_E_my(y, p, k, t) * y_E[y, p, k, t]) for
        #         t in model_data.index_t, k in ipp.index_k_existing, p in ipp.index_p
        #     ) + 
        #     sum(
        #         model_data.omega(t) *
        #         (ipp.v_C_my(y, p, k, t) * y_C[y, p, k, t]) for
        #         t in model_data.index_t, k in ipp.index_k_new, p in ipp.index_p
        #     ) == 
        #     sum(
        #         miu[y, t] * (
        #             sum(
        #                 customers.gamma(h) * customers.d_my(y, h, t) for
        #                 h in model_data.index_h
        #             ) + ipp.eximport_my(y, t) - sum(
        #                 customers.rho_DG(h, m, t) * customers.x_DG_E_my(y, h, m) for
        #                 h in model_data.index_h, m in customers.index_m
        #             ) - sum(
        #                 customers.rho_DG(h, m, t) * sum(
        #                     customers.x_DG_new_my(Symbol(Int(y_symbol)), h, m) for
        #                     y_symbol in
        #                     model_data.year(first(model_data.index_y_fix)):model_data.year[y]
        #                 ) for h in model_data.index_h, m in customers.index_m
        #             ) -
        #             # green technology subscription at time t
        #             sum(
        #                 utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
        #                 model_data.year(first(model_data.index_y_fix)):model_data.year(y))
        #                 for j in model_data.index_j, h in model_data.index_h
        #             )
        #         ) for t in model_data.index_t
        #     ) - 
        #     (
        #         sum(
        #             eta[y, p, k, t] *
        #             ipp.rho_E_my(p, k, t) *
        #             (
        #                 ipp.x_E_my(p, k) - ipp.x_R_cumu(p, k)
        #             ) for t in model_data.index_t, k in ipp.index_k_existing, p in ipp.index_p
        #         ) + sum(
        #             lambda[y, p, k, t] *
        #             ipp.rho_C_my(p, k, t) *
        #             ipp.x_C_cumu(p, k) 
        #             for t in model_data.index_t, k in ipp.index_k_new, p in ipp.index_p
        #         )
        #     ) + sum(
        #         eta[y, p, k, t] *
        #         ipp.rho_E_my(p, k, t) *
        #         sum(
        #             ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
        #             model_data.year(first(model_data.index_y)):model_data.year(y)
        #         ) for t in model_data.index_t, k in ipp.index_k_existing,
        #         p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
        #     ) - sum(
        #         lambda[y, p, k, t] *
        #         ipp.rho_C_my(p, k, t) *
        #         sum(
        #             ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
        #             model_data.year(first(model_data.index_y)):model_data.year(y)
        #         ) for t in model_data.index_t, k in ipp.index_k_new,
        #         p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
        #     ) + sum(
        #         ipp.rho_E_my(p_star, k, t) *
        #         mce_E_p_star[y, k, t] for t in model_data.index_t, k in ipp.index_k_existing
        #     ) - sum(
        #         ipp.rho_C_my(p_star, k, t) *
        #         mce_C_p_star[y, k, t] for t in model_data.index_t, k in ipp.index_k_new
        #     )
        # )
    else
        @constraint(
            WMDER_IPP,
            Eq_strong_duality[
                y in model_data.index_y
            ],
            sum(
                model_data.omega(d) * delta_t *
                (ipp.v_E_my(y, p, k, z, d, t) * y_E[y, p, k, z, d, t]) for
                d in model_data.index_d, t in model_data.index_t, z in model_data.index_z, k in ipp.index_k_existing, p in ipp.index_p
            ) + 
            sum(
                model_data.omega(d) * delta_t *
                (ipp.v_C_my(y, p, k, z, d, t) * y_C[y, p, k, z, d, t]) for
                d in model_data.index_d, t in model_data.index_t, z in model_data.index_z, k in ipp.index_k_new, p in ipp.index_p
            ) == 
            sum(
                miu[y, z, d, t] * (
                    sum(
                        customers.gamma(z, h) * customers.d_my(y, h, z, d, t) for h in model_data.index_h
                    ) + ipp.eximport_my(y, z, d, t) -
                    # total DG generation at time t
                    sum(
                        customers.rho_DG(h, m, z, d, t) * customers.total_der_capacity_my(y, z, h, m) for
                        h in model_data.index_h, m in customers.index_m
                    ) + 
                    # remove aggregated behind-the-meter pv/storage generation/consumption since they're front-of-the-meter now
                    sum(
                        customers.rho_DG(h, m, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) *
                        customers.total_pv_stor_capacity_my(y, z, h, m) for h in model_data.index_h, m in (:BTMStorage, :BTMPV)
                    ) -
                    # green technology subscription at time t
                    sum(
                        ipp.rho_C_my(Symbol("ipp1"), j, z, d, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, z, h) for y_symbol in
                        model_data.year(first(model_data.index_y_fix)):model_data.year(y))
                        for j in model_data.index_j, h in model_data.index_h
                    )
                ) for z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
            ) - 
            (
                sum(
                    eta[y, p, k, z, d, t] *
                    ipp.rho_E_my(p, k, z, d, t) * (
                        ipp.x_E_my(p, z, k) - ipp.x_R_cumu(p, k, z)
                    ) for p in ipp.index_p, k in ipp.index_k_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) -
                sum(
                    ipp.rho_E_my(p_star, k, z, d, t) * mce_eta_x_R_p_star[y, k, z, d, t] 
                    for k in ipp.index_k_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    lambda[y, p, k, z, d, t] *
                    ipp.rho_C_my(p, k, z, d, t) * ipp.x_C_cumu(p, k, z) 
                    for p in ipp.index_p, k in ipp.index_k_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) +
                sum(
                    ipp.rho_C_my(p_star, k, z, d, t) * mce_lambda_x_C_p_star[y, k, z, d, t] 
                    for k in ipp.index_k_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) + 
            (
                sum(iota_min[y, l, d, t] * ipp.trans_capacity(l, :min) for l in ipp.index_l, d in model_data.index_d, t in model_data.index_t) - 
                sum(iota_max[y, l, d, t] * ipp.trans_capacity(l, :max) for l in ipp.index_l, d in model_data.index_d, t in model_data.index_t)
            ) - 
            (
                sum(
                    theta_E_energy[y, p, s, z, d, t] *
                    ipp.stor_duration_existing(s) * (
                        ipp.x_stor_E_my(p, z, s) - ipp.x_stor_R_cumu(p, s, z)
                    ) for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) -
                sum(
                    ipp.stor_duration_existing(s) * mce_theta_E_energy_x_stor_R_p_star[y, s, z, d, t] 
                    for s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    theta_E_discharge[y, p, s, z, d, t] *
                    ipp.rte_stor_E_my(y, p, z, s) * (
                        ipp.x_stor_E_my(p, z, s) - ipp.x_stor_R_cumu(p, s, z)
                    ) for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) -
                sum(
                    ipp.rte_stor_E_my(y, p_star, z, s) * mce_theta_E_discharge_x_stor_R_p_star[y, s, z, d, t] 
                    for s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    theta_E_charge[y, p, s, z, d, t] *
                    (
                        ipp.x_stor_E_my(p, z, s) - ipp.x_stor_R_cumu(p, s, z)
                    ) for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) -
                sum(
                    mce_theta_E_charge_x_stor_R_p_star[y, s, z, d, t] 
                    for s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    pi_E_charge[y, p, s, z, d, t] *
                    ipp.stor_duration_existing(s) * (
                        ipp.x_stor_E_my(p, z, s) - ipp.x_stor_R_cumu(p, s, z)
                    ) for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) -
                sum(
                    ipp.stor_duration_existing(s) * mce_pi_E_charge_x_stor_R_p_star[y, s, z, d, t] 
                    for s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    kappa_E[y, p, s, z, d, t] *
                    (
                        ipp.x_stor_E_my(p, z, s) - ipp.x_stor_R_cumu(p, s, z)
                    ) for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) -
                sum(
                    mce_kappa_E_x_stor_R_p_star[y, s, z, d, t] 
                    for s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    theta_C_energy[y, p, s, z, d, t] *
                    ipp.stor_duration_new(s) * ipp.x_stor_C_cumu(p, s, z) 
                    for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) +
                sum(
                    ipp.stor_duration_new(s) * mce_theta_C_energy_x_stor_C_p_star[y, s, z, d, t] 
                    for s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    theta_C_discharge[y, p, s, z, d, t] *
                    ipp.rte_stor_C_my(y, p, z, s) * ipp.x_stor_C_cumu(p, s, z) 
                    for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) +
                sum(
                    ipp.rte_stor_C_my(y, p_star, z, s) * mce_theta_C_discharge_x_stor_C_p_star[y, s, z, d, t] 
                    for s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    theta_C_charge[y, p, s, z, d, t] * ipp.x_stor_C_cumu(p, s, z) 
                    for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) +
                sum(
                    mce_theta_C_charge_x_stor_C_p_star[y, s, z, d, t] 
                    for s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    pi_C_charge[y, p, s, z, d, t] *
                    ipp.stor_duration_new(s) * ipp.x_stor_C_cumu(p, s, z) 
                    for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) +
                sum(
                    ipp.stor_duration_new(s) * mce_pi_C_charge_x_stor_C_p_star[y, s, z, d, t] 
                    for s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    kappa_C[y, p, s, z, d, t] * ipp.x_stor_C_cumu(p, s, z) 
                    for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) +
                sum(
                    mce_kappa_C_x_stor_C_p_star[y, s, z, d, t] 
                    for s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) + 
            sum(psi_E[y, p, s, z, d, model_data.index_t.elements[1]] * ipp.initial_energy_existing_my(y, p, s, z, d)
                for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d) - 
            sum(pi_E_discharge[y, p, s, z, d, model_data.index_t.elements[1]] * ipp.initial_energy_existing_my(y, p, s, z, d)
                for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d) + 
            sum(pi_E_charge[y, p, s, z, d, model_data.index_t.elements[1]] * ipp.initial_energy_existing_my(y, p, s, z, d)
                for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d) + 
            sum(psi_C[y, p, s, z, d, model_data.index_t.elements[1]] * ipp.initial_energy_new_my(y, p, s, z, d)
                for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d) - 
            sum(pi_C_discharge[y, p, s, z, d, model_data.index_t.elements[1]] * ipp.initial_energy_new_my(y, p, s, z, d)
                for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d) + 
            sum(pi_C_charge[y, p, s, z, d, model_data.index_t.elements[1]] * ipp.initial_energy_new_my(y, p, s, z, d)
                for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d)    
        )
    end

    # linking constraints on total retirement of p_star at year y
    @constraint(
        WMDER_IPP,
        Eq_cumu_R_mce[
            y in model_data.index_y,
            k in ipp.index_k_existing,
            z in model_data.index_z
        ],
        x_R_mce[y, k, z] == sum(
            x_R[Symbol(Int(y_symbol)), k, z] for
            y_symbol in model_data.year(first(model_data.index_y)):model_data.year(y)
        )
    )
    @constraint(
        WMDER_IPP,
        Eq_stor_cumu_R_mce[
            y in model_data.index_y,
            s in ipp.index_stor_existing,
            z in model_data.index_z
        ],
        x_stor_R_mce[y, s, z] == sum(
            x_stor_R[Symbol(Int(y_symbol)), s, z] for
            y_symbol in model_data.year(first(model_data.index_y)):model_data.year(y)
        )
    )

    # linking constraints on total investment of p_star at year y
    @constraint(
        WMDER_IPP,
        Eq_cumu_C_mce[
            y in model_data.index_y,
            k in ipp.index_k_new,
            z in model_data.index_z
        ],
        x_C_mce[y, k, z] == sum(
            x_C[Symbol(Int(y_symbol)), k, z] for
            y_symbol in model_data.year(first(model_data.index_y)):model_data.year(y)
        )
    )
    @constraint(
        WMDER_IPP,
        Eq_stor_cumu_C_mce[
            y in model_data.index_y,
            k in ipp.index_stor_new,
            z in model_data.index_z
        ],
        x_stor_C_mce[y, k, z] == sum(
            x_stor_C[Symbol(Int(y_symbol)), k, z] for
            y_symbol in model_data.year(first(model_data.index_y)):model_data.year(y)
        )
    )

    # McCormick-envelope relaxation
    # apply constraint_scaling here
    @constraint(
        WMDER_IPP,
        Eq_mce_eta_x_R_p_star_LL[
            y in model_data.index_y,
            k in ipp.index_k_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_eta_x_R_p_star[y, k, z, d, t] / constraint_scaling >= (X_R_cumu_L(y, k, z) * eta[y, p_star, k, z, d, t] + 
        eta_L(y, k, z, d, t) * x_R_mce[y, k, z] - eta_L(y, k, z, d, t) * X_R_cumu_L(y, k, z)) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_eta_x_R_p_star_UU[
            y in model_data.index_y,
            k in ipp.index_k_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_eta_x_R_p_star[y, k, z, d, t] / constraint_scaling >= (- eta_U(y, k, z, d, t) * X_R_cumu_U(y, k, z) +
        eta_U(y, k, z, d, t) * x_R_mce[y, k, z] + eta[y, p_star, k, z, d, t] * X_R_cumu_U(y, k, z)) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_eta_x_R_p_star_LU[
            y in model_data.index_y,
            k in ipp.index_k_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_eta_x_R_p_star[y, k, z, d, t] / constraint_scaling <= (X_R_cumu_U(y, k, z) * eta[y, p_star, k, z, d, t] -
        eta_L(y, k, z, d, t) * X_R_cumu_U(y, k, z) + eta_L(y, k, z, d, t) * x_R_mce[y, k, z]) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_eta_x_R_p_star_UL[
            y in model_data.index_y,
            k in ipp.index_k_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_eta_x_R_p_star[y, k, z, d, t] / constraint_scaling <= (eta_U(y, k, z, d, t) * x_R_mce[y, k, z] -
        eta_U(y, k, z, d, t) * X_R_cumu_L(y, k, z) + X_R_cumu_L(y, k, z) * eta[y, p_star, k, z, d, t]) / constraint_scaling
    )

    @constraint(
        WMDER_IPP,
        Eq_mce_lambda_x_C_p_star_LL[
            y in model_data.index_y,
            k in ipp.index_k_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_lambda_x_C_p_star[y, k, z, d, t] / constraint_scaling >= (X_C_cumu_L(y, k, z) * lambda[y, p_star, k, z, d, t] + 
        lambda_L(y, k, z, d, t) * x_C_mce[y, k, z] - lambda_L(y, k, z, d, t) * X_C_cumu_L(y, k, z)) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_lambda_x_C_p_star_UU[
            y in model_data.index_y,
            k in ipp.index_k_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_lambda_x_C_p_star[y, k, z, d, t] / constraint_scaling >= (- lambda_U(y, k, z, d, t) * X_C_cumu_U(y, k, z) +
        lambda_U(y, k, z, d, t) * x_C_mce[y, k, z] + lambda[y, p_star, k, z, d, t] * X_C_cumu_U(y, k, z)) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_lambda_x_C_p_star_LU[
            y in model_data.index_y,
            k in ipp.index_k_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_lambda_x_C_p_star[y, k, z, d, t] / constraint_scaling <= (X_C_cumu_U(y, k, z) * lambda[y, p_star, k, z, d, t] -
        lambda_L(y, k, z, d, t) * X_C_cumu_U(y, k, z) + lambda_L(y, k, z, d, t) * x_C_mce[y, k, z]) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_lambda_x_C_p_star_UL[
            y in model_data.index_y,
            k in ipp.index_k_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_lambda_x_C_p_star[y, k, z, d, t] / constraint_scaling <= (lambda_U(y, k, z, d, t) * x_C_mce[y, k, z] -
        lambda_U(y, k, z, d, t) * X_C_cumu_L(y, k, z) + X_C_cumu_L(y, k, z) * lambda[y, p_star, k, z, d, t]) / constraint_scaling
    )

    @constraint(
        WMDER_IPP,
        Eq_mce_theta_E_energy_x_stor_R_p_star_LL[
            y in model_data.index_y,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_theta_E_energy_x_stor_R_p_star[y, s, z, d, t] / constraint_scaling >= (X_R_stor_cumu_L(y, s, z) * theta_E_energy[y, p_star, s, z, d, t] + 
        theta_E_energy_L(y, s, z, d, t) * x_stor_R_mce[y, s, z] - theta_E_energy_L(y, s, z, d, t) * X_R_stor_cumu_L(y, s, z)) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_theta_E_energy_x_stor_R_p_star_UU[
            y in model_data.index_y,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_theta_E_energy_x_stor_R_p_star[y, s, z, d, t] / constraint_scaling >= (- theta_E_energy_U(y, s, z, d, t) * X_R_stor_cumu_U(y, s, z) +
        theta_E_energy_U(y, s, z, d, t) * x_stor_R_mce[y, s, z] + theta_E_energy[y, p_star, s, z, d, t] * X_R_stor_cumu_U(y, s, z)) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_theta_E_energy_x_stor_R_p_star_LU[
            y in model_data.index_y,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_theta_E_energy_x_stor_R_p_star[y, s, z, d, t] / constraint_scaling <= (X_R_stor_cumu_U(y, s, z) * theta_E_energy[y, p_star, s, z, d, t] -
        theta_E_energy_L(y, s, z, d, t) * X_R_stor_cumu_U(y, s, z) + theta_E_energy_L(y, s, z, d, t) * x_stor_R_mce[y, s, z]) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_theta_E_energy_x_stor_R_p_star_UL[
            y in model_data.index_y,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_theta_E_energy_x_stor_R_p_star[y, s, z, d, t] / constraint_scaling <= (theta_E_energy_U(y, s, z, d, t) * x_stor_R_mce[y, s, z] -
        theta_E_energy_U(y, s, z, d, t) * X_R_stor_cumu_L(y, s, z) + X_R_stor_cumu_L(y, s, z) * theta_E_energy[y, p_star, s, z, d, t]) / constraint_scaling
    )

    @constraint(
        WMDER_IPP,
        Eq_mce_theta_E_discharge_x_stor_R_p_star_LL[
            y in model_data.index_y,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_theta_E_discharge_x_stor_R_p_star[y, s, z, d, t] / constraint_scaling >= (X_R_stor_cumu_L(y, s, z) * theta_E_discharge[y, p_star, s, z, d, t] + 
        theta_E_discharge_L(y, s, z, d, t) * x_stor_R_mce[y, s, z] - theta_E_discharge_L(y, s, z, d, t) * X_R_stor_cumu_L(y, s, z)) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_theta_E_discharge_x_stor_R_p_star_UU[
            y in model_data.index_y,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_theta_E_discharge_x_stor_R_p_star[y, s, z, d, t] / constraint_scaling >= (- theta_E_discharge_U(y, s, z, d, t) * X_R_stor_cumu_U(y, s, z) +
        theta_E_discharge_U(y, s, z, d, t) * x_stor_R_mce[y, s, z] + theta_E_discharge[y, p_star, s, z, d, t] * X_R_stor_cumu_U(y, s, z)) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_theta_E_discharge_x_stor_R_p_star_LU[
            y in model_data.index_y,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_theta_E_discharge_x_stor_R_p_star[y, s, z, d, t] / constraint_scaling <= (X_R_stor_cumu_U(y, s, z) * theta_E_discharge[y, p_star, s, z, d, t] -
        theta_E_discharge_L(y, s, z, d, t) * X_R_stor_cumu_U(y, s, z) + theta_E_discharge_L(y, s, z, d, t) * x_stor_R_mce[y, s, z]) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_theta_E_discharge_x_stor_R_p_star_UL[
            y in model_data.index_y,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_theta_E_discharge_x_stor_R_p_star[y, s, z, d, t] / constraint_scaling <= (theta_E_discharge_U(y, s, z, d, t) * x_stor_R_mce[y, s, z] -
        theta_E_discharge_U(y, s, z, d, t) * X_R_stor_cumu_L(y, s, z) + X_R_stor_cumu_L(y, s, z) * theta_E_discharge[y, p_star, s, z, d, t]) / constraint_scaling
    )

    @constraint(
        WMDER_IPP,
        Eq_mce_theta_E_charge_x_stor_R_p_star_LL[
            y in model_data.index_y,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_theta_E_charge_x_stor_R_p_star[y, s, z, d, t] / constraint_scaling >= (X_R_stor_cumu_L(y, s, z) * theta_E_charge[y, p_star, s, z, d, t] + 
        theta_E_charge_L(y, s, z, d, t) * x_stor_R_mce[y, s, z] - theta_E_charge_L(y, s, z, d, t) * X_R_stor_cumu_L(y, s, z)) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_theta_E_charge_x_stor_R_p_star_UU[
            y in model_data.index_y,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_theta_E_charge_x_stor_R_p_star[y, s, z, d, t] / constraint_scaling >= (- theta_E_charge_U(y, s, z, d, t) * X_R_stor_cumu_U(y, s, z) +
        theta_E_charge_U(y, s, z, d, t) * x_stor_R_mce[y, s, z] + theta_E_charge[y, p_star, s, z, d, t] * X_R_stor_cumu_U(y, s, z)) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_theta_E_charge_x_stor_R_p_star_LU[
            y in model_data.index_y,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_theta_E_charge_x_stor_R_p_star[y, s, z, d, t] / constraint_scaling <= (X_R_stor_cumu_U(y, s, z) * theta_E_charge[y, p_star, s, z, d, t] -
        theta_E_charge_L(y, s, z, d, t) * X_R_stor_cumu_U(y, s, z) + theta_E_charge_L(y, s, z, d, t) * x_stor_R_mce[y, s, z]) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_theta_E_charge_x_stor_R_p_star_UL[
            y in model_data.index_y,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_theta_E_charge_x_stor_R_p_star[y, s, z, d, t] / constraint_scaling <= (theta_E_charge_U(y, s, z, d, t) * x_stor_R_mce[y, s, z] -
        theta_E_charge_U(y, s, z, d, t) * X_R_stor_cumu_L(y, s, z) + X_R_stor_cumu_L(y, s, z) * theta_E_charge[y, p_star, s, z, d, t]) / constraint_scaling
    )

    @constraint(
        WMDER_IPP,
        Eq_mce_pi_E_charge_x_stor_R_p_star_LL[
            y in model_data.index_y,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_pi_E_charge_x_stor_R_p_star[y, s, z, d, t] / constraint_scaling >= (X_R_stor_cumu_L(y, s, z) * pi_E_charge[y, p_star, s, z, d, t] + 
        pi_E_charge_L(y, s, z, d, t) * x_stor_R_mce[y, s, z] - pi_E_charge_L(y, s, z, d, t) * X_R_stor_cumu_L(y, s, z)) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_pi_E_charge_x_stor_R_p_star_UU[
            y in model_data.index_y,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_pi_E_charge_x_stor_R_p_star[y, s, z, d, t] / constraint_scaling >= (- pi_E_charge_U(y, s, z, d, t) * X_R_stor_cumu_U(y, s, z) +
        pi_E_charge_U(y, s, z, d, t) * x_stor_R_mce[y, s, z] + pi_E_charge[y, p_star, s, z, d, t] * X_R_stor_cumu_U(y, s, z)) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_pi_E_charge_x_stor_R_p_star_LU[
            y in model_data.index_y,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_pi_E_charge_x_stor_R_p_star[y, s, z, d, t] / constraint_scaling <= (X_R_stor_cumu_U(y, s, z) * pi_E_charge[y, p_star, s, z, d, t] -
        pi_E_charge_L(y, s, z, d, t) * X_R_stor_cumu_U(y, s, z) + pi_E_charge_L(y, s, z, d, t) * x_stor_R_mce[y, s, z]) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_pi_E_charge_x_stor_R_p_star_UL[
            y in model_data.index_y,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_pi_E_charge_x_stor_R_p_star[y, s, z, d, t] / constraint_scaling <= (pi_E_charge_U(y, s, z, d, t) * x_stor_R_mce[y, s, z] -
        pi_E_charge_U(y, s, z, d, t) * X_R_stor_cumu_L(y, s, z) + X_R_stor_cumu_L(y, s, z) * pi_E_charge[y, p_star, s, z, d, t]) / constraint_scaling
    )
    
    @constraint(
        WMDER_IPP,
        Eq_mce_kappa_E_x_stor_R_p_star_LL[
            y in model_data.index_y,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_kappa_E_x_stor_R_p_star[y, s, z, d, t] / constraint_scaling >= (X_R_stor_cumu_L(y, s, z) * kappa_E[y, p_star, s, z, d, t] + 
        kappa_E_L(y, s, z, d, t) * x_stor_R_mce[y, s, z] - kappa_E_L(y, s, z, d, t) * X_R_stor_cumu_L(y, s, z)) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_kappa_E_x_stor_R_p_star_UU[
            y in model_data.index_y,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_kappa_E_x_stor_R_p_star[y, s, z, d, t] / constraint_scaling >= (- kappa_E_U(y, s, z, d, t) * X_R_stor_cumu_U(y, s, z) +
        kappa_E_U(y, s, z, d, t) * x_stor_R_mce[y, s, z] + kappa_E[y, p_star, s, z, d, t] * X_R_stor_cumu_U(y, s, z)) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_kappa_E_x_stor_R_p_star_LU[
            y in model_data.index_y,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_kappa_E_x_stor_R_p_star[y, s, z, d, t] / constraint_scaling <= (X_R_stor_cumu_U(y, s, z) * kappa_E[y, p_star, s, z, d, t] -
        kappa_E_L(y, s, z, d, t) * X_R_stor_cumu_U(y, s, z) + kappa_E_L(y, s, z, d, t) * x_stor_R_mce[y, s, z]) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_kappa_E_x_stor_R_p_star_UL[
            y in model_data.index_y,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_kappa_E_x_stor_R_p_star[y, s, z, d, t] / constraint_scaling <= (kappa_E_U(y, s, z, d, t) * x_stor_R_mce[y, s, z] -
        kappa_E_U(y, s, z, d, t) * X_R_stor_cumu_L(y, s, z) + X_R_stor_cumu_L(y, s, z) * kappa_E[y, p_star, s, z, d, t]) / constraint_scaling
    )

    @constraint(
        WMDER_IPP,
        Eq_mce_theta_C_energy_x_stor_C_p_star_LL[
            y in model_data.index_y,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_theta_C_energy_x_stor_C_p_star[y, s, z, d, t] / constraint_scaling >= (X_C_stor_cumu_L(y, s, z) * theta_C_energy[y, p_star, s, z, d, t] + 
        theta_C_energy_L(y, s, z, d, t) * x_stor_C_mce[y, s, z] - theta_C_energy_L(y, s, z, d, t) * X_C_stor_cumu_L(y, s, z)) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_theta_C_energy_x_stor_C_p_star_UU[
            y in model_data.index_y,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_theta_C_energy_x_stor_C_p_star[y, s, z, d, t] / constraint_scaling >= (- theta_C_energy_U(y, s, z, d, t) * X_C_stor_cumu_U(y, s, z) +
        theta_C_energy_U(y, s, z, d, t) * x_stor_C_mce[y, s, z] + theta_C_energy[y, p_star, s, z, d, t] * X_C_stor_cumu_U(y, s, z)) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_theta_C_energy_x_stor_C_p_star_LU[
            y in model_data.index_y,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_theta_C_energy_x_stor_C_p_star[y, s, z, d, t] / constraint_scaling <= (X_C_stor_cumu_U(y, s, z) * theta_C_energy[y, p_star, s, z, d, t] -
        theta_C_energy_L(y, s, z, d, t) * X_C_stor_cumu_U(y, s, z) + theta_C_energy_L(y, s, z, d, t) * x_stor_C_mce[y, s, z]) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_theta_C_energy_x_stor_C_p_star_UL[
            y in model_data.index_y,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_theta_C_energy_x_stor_C_p_star[y, s, z, d, t] / constraint_scaling <= (theta_C_energy_U(y, s, z, d, t) * x_stor_C_mce[y, s, z] -
        theta_C_energy_U(y, s, z, d, t) * X_C_stor_cumu_L(y, s, z) + X_C_stor_cumu_L(y, s, z) * theta_C_energy[y, p_star, s, z, d, t]) / constraint_scaling
    )

    @constraint(
        WMDER_IPP,
        Eq_mce_theta_C_discharge_x_stor_C_p_star_LL[
            y in model_data.index_y,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_theta_C_discharge_x_stor_C_p_star[y, s, z, d, t] / constraint_scaling >= (X_C_stor_cumu_L(y, s, z) * theta_C_discharge[y, p_star, s, z, d, t] + 
        theta_C_discharge_L(y, s, z, d, t) * x_stor_C_mce[y, s, z] - theta_C_discharge_L(y, s, z, d, t) * X_C_stor_cumu_L(y, s, z)) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_theta_C_discharge_x_stor_C_p_star_UU[
            y in model_data.index_y,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_theta_C_discharge_x_stor_C_p_star[y, s, z, d, t] / constraint_scaling >= (- theta_C_discharge_U(y, s, z, d, t) * X_C_stor_cumu_U(y, s, z) +
        theta_C_discharge_U(y, s, z, d, t) * x_stor_C_mce[y, s, z] + theta_C_discharge[y, p_star, s, z, d, t] * X_C_stor_cumu_U(y, s, z)) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_theta_C_discharge_x_stor_C_p_star_LU[
            y in model_data.index_y,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_theta_C_discharge_x_stor_C_p_star[y, s, z, d, t] / constraint_scaling <= (X_C_stor_cumu_U(y, s, z) * theta_C_discharge[y, p_star, s, z, d, t] -
        theta_C_discharge_L(y, s, z, d, t) * X_C_stor_cumu_U(y, s, z) + theta_C_discharge_L(y, s, z, d, t) * x_stor_C_mce[y, s, z]) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_theta_C_discharge_x_stor_C_p_star_UL[
            y in model_data.index_y,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_theta_C_discharge_x_stor_C_p_star[y, s, z, d, t] / constraint_scaling <= (theta_C_discharge_U(y, s, z, d, t) * x_stor_C_mce[y, s, z] -
        theta_C_discharge_U(y, s, z, d, t) * X_C_stor_cumu_L(y, s, z) + X_C_stor_cumu_L(y, s, z) * theta_C_discharge[y, p_star, s, z, d, t]) / constraint_scaling
    )

    @constraint(
        WMDER_IPP,
        Eq_mce_theta_C_charge_x_stor_C_p_star_LL[
            y in model_data.index_y,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_theta_C_charge_x_stor_C_p_star[y, s, z, d, t] / constraint_scaling >= (X_C_stor_cumu_L(y, s, z) * theta_C_charge[y, p_star, s, z, d, t] + 
        theta_C_charge_L(y, s, z, d, t) * x_stor_C_mce[y, s, z] - theta_C_charge_L(y, s, z, d, t) * X_C_stor_cumu_L(y, s, z)) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_theta_C_charge_x_stor_C_p_star_UU[
            y in model_data.index_y,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_theta_C_charge_x_stor_C_p_star[y, s, z, d, t] / constraint_scaling >= (- theta_C_charge_U(y, s, z, d, t) * X_C_stor_cumu_U(y, s, z) +
        theta_C_charge_U(y, s, z, d, t) * x_stor_C_mce[y, s, z] + theta_C_charge[y, p_star, s, z, d, t] * X_C_stor_cumu_U(y, s, z)) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_theta_C_charge_x_stor_C_p_star_LU[
            y in model_data.index_y,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_theta_C_charge_x_stor_C_p_star[y, s, z, d, t] / constraint_scaling <= (X_C_stor_cumu_U(y, s, z) * theta_C_charge[y, p_star, s, z, d, t] -
        theta_C_charge_L(y, s, z, d, t) * X_C_stor_cumu_U(y, s, z) + theta_C_charge_L(y, s, z, d, t) * x_stor_C_mce[y, s, z]) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_theta_C_charge_x_stor_C_p_star_UL[
            y in model_data.index_y,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_theta_C_charge_x_stor_C_p_star[y, s, z, d, t] / constraint_scaling <= (theta_C_charge_U(y, s, z, d, t) * x_stor_C_mce[y, s, z] -
        theta_C_charge_U(y, s, z, d, t) * X_C_stor_cumu_L(y, s, z) + X_C_stor_cumu_L(y, s, z) * theta_C_charge[y, p_star, s, z, d, t]) / constraint_scaling
    )

    @constraint(
        WMDER_IPP,
        Eq_mce_pi_C_charge_x_stor_C_p_star_LL[
            y in model_data.index_y,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_pi_C_charge_x_stor_C_p_star[y, s, z, d, t] / constraint_scaling >= (X_C_stor_cumu_L(y, s, z) * pi_C_charge[y, p_star, s, z, d, t] + 
        pi_C_charge_L(y, s, z, d, t) * x_stor_C_mce[y, s, z] - pi_C_charge_L(y, s, z, d, t) * X_C_stor_cumu_L(y, s, z)) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_pi_C_charge_x_stor_C_p_star_UU[
            y in model_data.index_y,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_pi_C_charge_x_stor_C_p_star[y, s, z, d, t] / constraint_scaling >= (- pi_C_charge_U(y, s, z, d, t) * X_C_stor_cumu_U(y, s, z) +
        pi_C_charge_U(y, s, z, d, t) * x_stor_C_mce[y, s, z] + pi_C_charge[y, p_star, s, z, d, t] * X_C_stor_cumu_U(y, s, z)) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_pi_C_charge_x_stor_C_p_star_LU[
            y in model_data.index_y,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_pi_C_charge_x_stor_C_p_star[y, s, z, d, t] / constraint_scaling <= (X_C_stor_cumu_U(y, s, z) * pi_C_charge[y, p_star, s, z, d, t] -
        pi_C_charge_L(y, s, z, d, t) * X_C_stor_cumu_U(y, s, z) + pi_C_charge_L(y, s, z, d, t) * x_stor_C_mce[y, s, z]) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_pi_C_charge_x_stor_C_p_star_UL[
            y in model_data.index_y,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_pi_C_charge_x_stor_C_p_star[y, s, z, d, t] / constraint_scaling <= (pi_C_charge_U(y, s, z, d, t) * x_stor_C_mce[y, s, z] -
        pi_C_charge_U(y, s, z, d, t) * X_C_stor_cumu_L(y, s, z) + X_C_stor_cumu_L(y, s, z) * pi_C_charge[y, p_star, s, z, d, t]) / constraint_scaling
    )

    @constraint(
        WMDER_IPP,
        Eq_mce_kappa_C_x_stor_C_p_star_LL[
            y in model_data.index_y,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_kappa_C_x_stor_C_p_star[y, s, z, d, t] / constraint_scaling >= (X_C_stor_cumu_L(y, s, z) * kappa_C[y, p_star, s, z, d, t] + 
        kappa_C_L(y, s, z, d, t) * x_stor_C_mce[y, s, z] - kappa_C_L(y, s, z, d, t) * X_C_stor_cumu_L(y, s, z)) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_kappa_C_x_stor_C_p_star_UU[
            y in model_data.index_y,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_kappa_C_x_stor_C_p_star[y, s, z, d, t] / constraint_scaling >= (- kappa_C_U(y, s, z, d, t) * X_C_stor_cumu_U(y, s, z) +
        kappa_C_U(y, s, z, d, t) * x_stor_C_mce[y, s, z] + kappa_C[y, p_star, s, z, d, t] * X_C_stor_cumu_U(y, s, z)) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_kappa_C_x_stor_C_p_star_LU[
            y in model_data.index_y,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_kappa_C_x_stor_C_p_star[y, s, z, d, t] / constraint_scaling <= (X_C_stor_cumu_U(y, s, z) * kappa_C[y, p_star, s, z, d, t] -
        kappa_C_L(y, s, z, d, t) * X_C_stor_cumu_U(y, s, z) + kappa_C_L(y, s, z, d, t) * x_stor_C_mce[y, s, z]) / constraint_scaling
    )
    @constraint(
        WMDER_IPP,
        Eq_mce_kappa_C_x_stor_C_p_star_UL[
            y in model_data.index_y,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t
        ],
        mce_kappa_C_x_stor_C_p_star[y, s, z, d, t] / constraint_scaling <= (kappa_C_U(y, s, z, d, t) * x_stor_C_mce[y, s, z] -
        kappa_C_U(y, s, z, d, t) * X_C_stor_cumu_L(y, s, z) + X_C_stor_cumu_L(y, s, z) * kappa_C[y, p_star, s, z, d, t]) / constraint_scaling
    )

    if length(ipp.index_p) >= 2
        # planning_reserves_cap =
        #     y -> begin
        #         # bulk generation available capacity at time t
        #         sum(
        #             ipp.capacity_credit_E_my(y, k) * (
        #                 ipp.x_E_my(p_star, k) - sum(
        #                     x_R[Symbol(Int(y_symbol)), k] for y_symbol in
        #                     model_data.year(first(model_data.index_y)):model_data.year(y)
        #                 ) - ipp.x_R_cumu(p_star, k)
        #             ) for k in ipp.index_k_existing
        #         ) +
        #         sum(
        #             ipp.capacity_credit_C_my(y, k) * (
        #                 sum(
        #                     x_C[Symbol(Int(y_symbol)), k] for y_symbol in
        #                     model_data.year(first(model_data.index_y)):model_data.year(y)
        #                 ) + ipp.x_C_cumu(p_star, k)
        #             ) for k in ipp.index_k_new
        #         ) +
        #         sum(
        #             ipp.capacity_credit_E_my(y, k) * (
        #                 ipp.x_E_my(p, k) - sum(
        #                     ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
        #                     model_data.year(first(model_data.index_y)):model_data.year(y)
        #                 ) - ipp.x_R_cumu(p, k)
        #             ) for k in ipp.index_k_existing,
        #             p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
        #         ) +
        #         sum(
        #             ipp.capacity_credit_C_my(y, k) * (
        #                 sum(
        #                     ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
        #                     model_data.year(first(model_data.index_y)):model_data.year(y)
        #                 ) + ipp.x_C_cumu(p, k)
        #             ) for k in ipp.index_k_new,
        #             p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
        #         ) +
        #         # green technology subscription
        #         sum(
        #             ipp.capacity_credit_C_my(y, j) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
        #             model_data.year(first(model_data.index_y_fix)):model_data.year(y))
        #             for j in model_data.index_j, h in model_data.index_h
        #         ) -
        #         # net_load plus planning reserve
        #         ipp.Reserve_req_my(y)
        #     end
    else
        planning_reserves_cap =
        (y, z) -> begin
            # bulk generation available capacity at time t
            sum(
                ipp.capacity_credit_E_my(y, z, k) * (
                    ipp.x_E_my(p_star, z, k) - sum(
                        x_R[Symbol(Int(y_symbol)), k, z] for y_symbol in
                        model_data.year(first(model_data.index_y)):model_data.year(y)
                    ) - ipp.x_R_cumu(p_star, k, z)
                ) for k in ipp.index_k_existing
            ) + sum(
                ipp.capacity_credit_C_my(y, z, k) * (
                    sum(
                        x_C[Symbol(Int(y_symbol)), k, z] for y_symbol in
                        model_data.year(first(model_data.index_y)):model_data.year(y)
                    ) + ipp.x_C_cumu(p_star, k, z)
                ) for k in ipp.index_k_new
            ) +
            sum(
                ipp.capacity_credit_stor_E_my(y, z, s) * (
                    ipp.x_stor_E_my(p_star, z, s) - sum(
                        x_stor_R[Symbol(Int(y_symbol)), s, z] for y_symbol in
                        model_data.year(first(model_data.index_y)):model_data.year(y)
                    ) - ipp.x_stor_R_cumu(p_star, s, z)
                ) for s in ipp.index_stor_existing
            ) + sum(
                ipp.capacity_credit_stor_C_my(y, z, s) * (
                    sum(
                        x_stor_C[Symbol(Int(y_symbol)), s, z] for y_symbol in
                        model_data.year(first(model_data.index_y)):model_data.year(y)
                    ) + ipp.x_stor_C_cumu(p_star, s, z)
                ) for s in ipp.index_stor_new
            ) +
            # green technology subscription
            sum(
                ipp.capacity_credit_C_my(y, z, j) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, z, h) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):model_data.year(y))
                for j in model_data.index_j, h in model_data.index_h
            ) -
            # flow out of zone z
            sum(ipp.trans_topology(l, z) * flow_cap[y, l] for l in ipp.index_l) -
            ipp.eximport_my(y, z, Max_Net_Load_my_dict[y, z][1], Max_Net_Load_my_dict[y, z][2]) - 
            # net_load plus planning reserve
            ipp.Reserve_req_my(y, z)
        end

    end
    @constraint(
        WMDER_IPP,
        Eq_xi_cap[y in model_data.index_y, z in model_data.index_z],
        planning_reserves_cap(y, z) >= 0
    )

    @constraint(
        WMDER_IPP,
        Eq_cap_flow_lower[
            y in model_data.index_y,
            l in ipp.index_l,
        ],
        flow_cap[y, l] - ipp.trans_capacity(l, :min) >= 0
    )

    @constraint(
        WMDER_IPP,
        Eq_cap_flow_upper[
            y in model_data.index_y,
            l in ipp.index_l,
        ],
        flow_cap[y, l] - ipp.trans_capacity(l, :max) <= 0
    )

    # RPS constraint
    index_rps_existing = deepcopy(ipp.index_rps)
    push!(index_rps_existing.elements, Symbol("dera_pv"))

    @constraint(
        WMDER_IPP,
        Eq_rps[y in model_data.index_y],
        (sum(
            model_data.omega(d) * delta_t * y_E[y, p, rps, z, d, t] for
            p in ipp.index_p, rps in index_rps_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
        ) + 
        sum(
            model_data.omega(d) * delta_t * y_C[y, p, rps, z, d, t] for
            p in ipp.index_p, rps in ipp.index_rps, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
        ) -
        ipp.RPS(y) *
        sum(model_data.omega(d) * delta_t * ipp.Net_Load_my(y, z, d, t) for z in model_data.index_z, d in model_data.index_d, t in model_data.index_t)) / constraint_scaling >=
        0
    )

    jump_model[end] = WMDER_IPP

    iteration_year = first(model_data.index_y)
    @info("iteration year is $(iteration_year)")
    @info("MPPDC upper level")
    TimerOutputs.@timeit HEM_TIMER "optimize! MPPDC with transmission and storage" begin
        optimize!(WMDER_IPP)
    end

    objective_value(WMDER_IPP)

    # compute_conflict!(WMDER_IPP)
    # iis_model, _ = copy_conflict(WMDER_IPP)
    # print(iis_model)

    # f = open("iis_model.txt","w"); print(f, iis_model); close(f)


    for y in model_data.index_y, k in ipp.index_k_existing, z in model_data.index_z
        ipp.x_R_my(y, p_star, k, z, :) .= value.(x_R[y, k, z])
    end
    for y in model_data.index_y, k in ipp.index_k_new, z in model_data.index_z
        ipp.x_C_my(y, p_star, k, z, :) .= value.(x_C[y, k, z])
    end
    for y in model_data.index_y, s in ipp.index_stor_existing, z in model_data.index_z
        ipp.x_stor_R_my(y, p_star, s, z, :) .= value.(x_stor_R[y, s, z])
    end
    for y in model_data.index_y, s in ipp.index_stor_new, z in model_data.index_z
        ipp.x_stor_C_my(y, p_star, s, z, :) .= value.(x_stor_C[y, s, z])
    end
    for y in model_data.index_y,
        p in ipp.index_p,
        k in ipp.index_k_existing,
        z in model_data.index_z,
        d in model_data.index_d,
        t in model_data.index_t

        ipp.y_E_my(y, p, k, z, d, t, :) .= value.(y_E[y, p, k, z, d, t])
    end
    for y in model_data.index_y,
        p in ipp.index_p,
        k in ipp.index_k_new,
        z in model_data.index_z,
        d in model_data.index_d,
        t in model_data.index_t

        ipp.y_C_my(y, p, k, z, d, t, :) .= value.(y_C[y, p, k, z, d, t])
    end
    for y in model_data.index_y, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
        ipp.miu_my(y, z, d, t, :) .= value.(miu[y, z, d, t])
    end
    for y in model_data.index_y, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
        ipp.LMP_my(y, z, d, t, :) .= ipp.miu_my(y, z, d, t) / (model_data.omega(d) * delta_t)
    end

    for y in model_data.index_y,
        p in ipp.index_p,
        s in ipp.index_stor_existing,
        z in model_data.index_z,
        d in model_data.index_d,
        t in model_data.index_t

        ipp.charge_E_my(y, p, s, z, d, t, :) .= value.(charge_E[y, p, s, z, d, t])
    end
    for y in model_data.index_y,
        p in ipp.index_p,
        s in ipp.index_stor_existing,
        z in model_data.index_z,
        d in model_data.index_d,
        t in model_data.index_t

        ipp.discharge_E_my(y, p, s, z, d, t, :) .= value.(discharge_E[y, p, s, z, d, t])
    end
    for y in model_data.index_y,
        p in ipp.index_p,
        s in ipp.index_stor_new,
        z in model_data.index_z,
        d in model_data.index_d,
        t in model_data.index_t

        ipp.charge_C_my(y, p, s, z, d, t, :) .= value.(charge_C[y, p, s, z, d, t])
    end
    for y in model_data.index_y,
        p in ipp.index_p,
        s in ipp.index_stor_new,
        z in model_data.index_z,
        d in model_data.index_d,
        t in model_data.index_t

        ipp.discharge_C_my(y, p, s, z, d, t, :) .= value.(discharge_C[y, p, s, z, d, t])
    end
    for y in model_data.index_y,
        p in ipp.index_p,
        s in ipp.index_stor_existing,
        z in model_data.index_z,
        d in model_data.index_d,
        t in model_data.index_t

        ipp.energy_E_my(y, p, s, z, d, t, :) .= value.(energy_E[y, p, s, z, d, t])
    end
    for y in model_data.index_y,
        p in ipp.index_p,
        s in ipp.index_stor_new,
        z in model_data.index_z,
        d in model_data.index_d,
        t in model_data.index_t

        ipp.energy_C_my(y, p, s, z, d, t, :) .= value.(energy_C[y, p, s, z, d, t])
    end
    for y in model_data.index_y,
        l in ipp.index_l,
        d in model_data.index_d,
        t in model_data.index_t

        ipp.flow_my(y, l, d, t, :) .= value.(flow[y, l, d, t])
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
    UCAP_p_star = make_keyed_array(model_data.index_y)
    for y in model_data.index_y
        UCAP_p_star(y, :) .= 
            sum(
                ipp.capacity_credit_E_my(y, z, k) * (
                    ipp.x_E_my(p_star, z, k) - sum(
                        ipp.x_R_my(Symbol(Int(y_symbol)), p_star, k, z) for y_symbol in
                        model_data.year(first(model_data.index_y)):model_data.year(y)
                    ) - ipp.x_R_cumu(p_star, k, z)
                ) for k in ipp.index_k_existing, z in model_data.index_z
            ) + sum(
                ipp.capacity_credit_C_my(y, z, k) * (
                    sum(
                        ipp.x_C_my(Symbol(Int(y_symbol)), p_star, k, z) for y_symbol in
                        model_data.year(first(model_data.index_y)):model_data.year(y)
                    ) + ipp.x_C_cumu(p_star, k, z)
                ) for k in ipp.index_k_new, z in model_data.index_z
            ) + sum(
                ipp.capacity_credit_stor_E_my(y, z, s) * (
                    ipp.x_stor_E_my(p_star, z, s) - sum(
                        ipp.x_stor_R_my(Symbol(Int(y_symbol)), p_star, s, z) for y_symbol in
                        model_data.year(first(model_data.index_y)):model_data.year(y)
                    ) - ipp.x_stor_R_cumu(p_star, s, z)
                ) for s in ipp.index_stor_existing, z in model_data.index_z
            ) + sum(
                ipp.capacity_credit_stor_C_my(y, z, s) * (
                    sum(
                        ipp.x_stor_C_my(Symbol(Int(y_symbol)), p_star, s, z) for y_symbol in
                        model_data.year(first(model_data.index_y)):model_data.year(y)
                    ) + ipp.x_stor_C_cumu(p_star, s, z)
                ) for s in ipp.index_stor_new, z in model_data.index_z
            )
    end

    UCAP_total = make_keyed_array(model_data.index_y)
    if length(ipp.index_p) >= 2
        for y in model_data.index_y
            UCAP_total(y, :) .= 
                UCAP_p_star(y) +
                sum(
                    ipp.capacity_credit_E_my(y, z, k) * (
                        ipp.x_E_my(p, z, k) - sum(
                            ipp.x_R_my(Symbol(Int(y_symbol)), p, k, z) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_R_cumu(p, k, z)
                    ) for k in ipp.index_k_existing,
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p))), 
                    z in model_data.index_z
                ) +
                sum(
                    ipp.capacity_credit_C_my(y, z, k) * (
                        sum(
                            ipp.x_C_my(Symbol(Int(y_symbol)), p, k, z) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_C_cumu(p, k, z)
                    ) for k in ipp.index_k_new,
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p))), 
                    z in model_data.index_z
                ) +
                sum(
                    ipp.capacity_credit_stor_E_my(y, z, s) * (
                        ipp.x_stor_E_my(p, z, s) - sum(
                            ipp.x_stor_R_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) - ipp.x_stor_R_cumu(p, s, z)
                    ) for s in ipp.index_stor_existing,
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p))), 
                    z in model_data.index_z
                ) +
                sum(
                    ipp.capacity_credit_stor_C_my(y, z, s) * (
                        sum(
                            ipp.x_stor_C_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                            model_data.year(first(model_data.index_y)):model_data.year(y)
                        ) + ipp.x_stor_C_cumu(p, s, z)
                    ) for s in ipp.index_stor_new,
                    p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p))), 
                    z in model_data.index_z
                ) +
                # green technology subscription
                sum(
                    ipp.capacity_credit_C_my(y, z, j) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, z, h) for y_symbol in
                    model_data.year(first(model_data.index_y_fix)):model_data.year(y))
                    for j in model_data.index_j, h in model_data.index_h, z in model_data.index_z
                ) - 
                # put exogenous export on the supply-side
                # don't have endogenous export/import because the capacity market clearing here assumes the entire region
                sum(ipp.eximport_my(y, z, Max_Net_Load_my_dict[y, z][1], Max_Net_Load_my_dict[y, z][2]) for z in model_data.index_z)
        end
    else
        for y in model_data.index_y
            UCAP_total(y, :) .= 
                UCAP_p_star(y) +
                # green technology subscription
                sum(
                    ipp.capacity_credit_C_my(y, z, j) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, z, h) for y_symbol in
                    model_data.year(first(model_data.index_y_fix)):model_data.year(y))
                    for j in model_data.index_j, h in model_data.index_h, z in model_data.index_z
                ) - 
                # put exogenous export on the supply-side
                # don't have endogenous export/import because the capacity market clearing here assumes the entire region
                sum(ipp.eximport_my(y, z, Max_Net_Load_my_dict[y, z][1], Max_Net_Load_my_dict[y, z][2]) for z in model_data.index_z)
        end
    end

    for y in model_data.index_y
        ipp.capacity_price(y, :) .=
            ipp.Capacity_intercept_my(y) + ipp.Capacity_slope_my(y) * UCAP_total(y)
        ipp.ucap(y, p_star, :) .= UCAP_p_star(y)
        ipp.ucap_total(y, :) .= UCAP_total(y)
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
            model_data.omega(d) * delta_t *
            (ipp.v_E_my(y, p, k, z, d, t) * value.(y_E[y, p, k, z, d, t])) for
            d in model_data.index_d, t in model_data.index_t, z in model_data.index_z, k in ipp.index_k_existing, p in ipp.index_p
        ) + 
        sum(
            model_data.omega(d) * delta_t *
            (ipp.v_C_my(y, p, k, z, d, t) * value.(y_C[y, p, k, z, d, t])) for
            d in model_data.index_d, t in model_data.index_t, z in model_data.index_z, k in ipp.index_k_new, p in ipp.index_p
        )

        if length(ipp.index_p) >= 2
            # lower_level_dual_obj[1, y] = 
            # sum(
            #     value.(miu[y, t]) * (
            #         sum(
            #             customers.gamma(h) * customers.d_my(y, h, t) for
            #             h in model_data.index_h
            #         ) + ipp.eximport_my(y, t) - sum(
            #             customers.rho_DG(h, m, t) * customers.x_DG_E_my(y, h, m) for
            #             h in model_data.index_h, m in customers.index_m
            #         ) - sum(
            #             customers.rho_DG(h, m, t) * sum(
            #                 customers.x_DG_new_my(Symbol(Int(y_symbol)), h, m) for
            #                 y_symbol in
            #                 model_data.year(first(model_data.index_y_fix)):model_data.year(y)
            #             ) for h in model_data.index_h, m in customers.index_m
            #         ) -
            #         # green technology subscription at time t
            #         sum(
            #             utility.rho_C_my(j, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, h) for y_symbol in
            #             model_data.year(first(model_data.index_y_fix)):model_data.year(y))
            #             for j in model_data.index_j, h in model_data.index_h
            #         )
            #     ) for t in model_data.index_t
            # ) - 
            # (
            #     sum(
            #         value.(eta[y, p, k, t]) *
            #         ipp.rho_E_my(p, k, t) *
            #         (
            #             ipp.x_E_my(p, k) - ipp.x_R_cumu(p, k)
            #         ) for t in model_data.index_t, k in ipp.index_k_existing, p in ipp.index_p
            #     ) + sum(
            #         value.(lambda[y, p, k, t]) *
            #         ipp.rho_C_my(p, k, t) *
            #         ipp.x_C_cumu(p, k) 
            #         for t in model_data.index_t, k in ipp.index_k_new, p in ipp.index_p
            #     )
            # ) + sum(
            #     value.(eta[y, p, k, t]) *
            #     ipp.rho_E_my(p, k, t) *
            #     sum(
            #         ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
            #         model_data.year(first(model_data.index_y)):model_data.year(y)
            #     ) for t in model_data.index_t, k in ipp.index_k_existing,
            #     p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
            # ) - sum(
            #     value.(lambda[y, p, k, t]) *
            #     ipp.rho_C_my(p, k, t) *
            #     sum(
            #         ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
            #         model_data.year(first(model_data.index_y)):model_data.year(y)
            #     ) for t in model_data.index_t, k in ipp.index_k_new,
            #     p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
            # ) + sum(
            #     ipp.rho_E_my(p_star, k, t) *
            #     value.(eta[y, p_star, k, t]) *
            #     value.(x_R_mce[y, k]) for t in model_data.index_t, k in ipp.index_k_existing
            # ) - sum(
            #     ipp.rho_C_my(p_star, k, t) *
            #     value.(lambda[y, p_star, k, t]) *
            #     value.(x_C_mce[y, k]) for t in model_data.index_t, k in ipp.index_k_new
            # )
        else
            lower_level_dual_obj[1, y] = 
            sum(
                value.(miu[y, z, d, t]) * (
                    sum(
                        customers.gamma(z, h) * customers.d_my(y, h, z, d, t) for h in model_data.index_h
                    ) + ipp.eximport_my(y, z, d, t) -
                    # total DG generation at time t
                    sum(
                        customers.rho_DG(h, m, z, d, t) * customers.total_der_capacity_my(y, z, h, m) for
                        h in model_data.index_h, m in customers.index_m
                    ) + 
                    # remove aggregated behind-the-meter pv/storage generation/consumption since they're front-of-the-meter now
                    sum(
                        customers.rho_DG(h, m, z, d, t) * der_aggregator.aggregation_level(reg_year_index_dera, z) *
                        customers.total_pv_stor_capacity_my(y, z, h, m) for h in model_data.index_h, m in (:BTMStorage, :BTMPV)
                    ) -
                    # green technology subscription at time t
                    sum(
                        ipp.rho_C_my(Symbol("ipp1"), j, z, d, t) * sum(green_developer.green_tech_buildout_my(Symbol(Int(y_symbol)), j, z, h) for y_symbol in
                        model_data.year(first(model_data.index_y_fix)):model_data.year(y))
                        for j in model_data.index_j, h in model_data.index_h
                    )
                ) for z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
            ) - 
            (
                sum(
                    value.(eta[y, p, k, z, d, t]) *
                    ipp.rho_E_my(p, k, z, d, t) * (
                        ipp.x_E_my(p, z, k) - ipp.x_R_cumu(p, k, z)
                    ) for p in ipp.index_p, k in ipp.index_k_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) -
                sum(
                    ipp.rho_E_my(p_star, k, z, d, t) * value.(eta[y, p_star, k, z, d, t]) * value.(x_R_mce[y, k, z]) 
                    for k in ipp.index_k_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    value.(lambda[y, p, k, z, d, t]) *
                    ipp.rho_C_my(p, k, z, d, t) * ipp.x_C_cumu(p, k, z) 
                    for p in ipp.index_p, k in ipp.index_k_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) +
                sum(
                    ipp.rho_C_my(p_star, k, z, d, t) * value.(lambda[y, p_star, k, z, d, t]) * value.(x_C_mce[y, k, z]) 
                    for k in ipp.index_k_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) + 
            (
                sum(value.(iota_min[y, l, d, t]) * ipp.trans_capacity(l, :min) for l in ipp.index_l, d in model_data.index_d, t in model_data.index_t) - 
                sum(value.(iota_max[y, l, d, t]) * ipp.trans_capacity(l, :max) for l in ipp.index_l, d in model_data.index_d, t in model_data.index_t)
            ) - 
            (
                sum(
                    value.(theta_E_energy[y, p, s, z, d, t]) *
                    ipp.stor_duration_existing(s) * (
                        ipp.x_stor_E_my(p, z, s) - ipp.x_stor_R_cumu(p, s, z)
                    ) for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) -
                sum(
                    ipp.stor_duration_existing(s) * value.(theta_E_energy[y, p_star, s, z, d, t]) * value.(x_stor_R_mce[y, s, z]) 
                    for s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    value.(theta_E_discharge[y, p, s, z, d, t]) *
                    ipp.rte_stor_E_my(y, p, z, s) * (
                        ipp.x_stor_E_my(p, z, s) - ipp.x_stor_R_cumu(p, s, z)
                    ) for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) -
                sum(
                    ipp.rte_stor_E_my(y, p_star, z, s) * value.(theta_E_discharge[y, p_star, s, z, d, t]) * value.(x_stor_R_mce[y, s, z]) 
                    for s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    value.(theta_E_charge[y, p, s, z, d, t]) *
                    (
                        ipp.x_stor_E_my(p, z, s) - ipp.x_stor_R_cumu(p, s, z)
                    ) for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) -
                sum(
                    value.(theta_E_charge[y, p_star, s, z, d, t]) * value.(x_stor_R_mce[y, s, z]) 
                    for s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    value.(pi_E_charge[y, p, s, z, d, t]) *
                    ipp.stor_duration_existing(s) * (
                        ipp.x_stor_E_my(p, z, s) - ipp.x_stor_R_cumu(p, s, z)
                    ) for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) -
                sum(
                    ipp.stor_duration_existing(s) * value.(pi_E_charge[y, p_star, s, z, d, t]) * value.(x_stor_R_mce[y, s, z]) 
                    for s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    value.(kappa_E[y, p, s, z, d, t]) *
                    (
                        ipp.x_stor_E_my(p, z, s) - ipp.x_stor_R_cumu(p, s, z)
                    ) for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) -
                sum(
                    value.(kappa_E[y, p_star, s, z, d, t]) * value.(x_stor_R_mce[y, s, z]) 
                    for s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    value.(theta_C_energy[y, p, s, z, d, t]) *
                    ipp.stor_duration_new(s) * ipp.x_stor_C_cumu(p, s, z) 
                    for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) +
                sum(
                    ipp.stor_duration_new(s) * value.(theta_C_energy[y, p_star, s, z, d, t]) * value.(x_stor_C_mce[y, s, z])
                    for s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    value.(theta_C_discharge[y, p, s, z, d, t]) *
                    ipp.rte_stor_C_my(y, p, z, s) * ipp.x_stor_C_cumu(p, s, z) 
                    for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) +
                sum(
                    ipp.rte_stor_C_my(y, p_star, z, s) * value.(theta_C_discharge[y, p_star, s, z, d, t]) * value.(x_stor_C_mce[y, s, z])
                    for s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    value.(theta_C_charge[y, p, s, z, d, t]) * ipp.x_stor_C_cumu(p, s, z) 
                    for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) +
                sum(
                    value.(theta_C_charge[y, p_star, s, z, d, t]) * value.(x_stor_C_mce[y, s, z])
                    for s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    value.(pi_C_charge[y, p, s, z, d, t]) *
                    ipp.stor_duration_new(s) * ipp.x_stor_C_cumu(p, s, z) 
                    for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) +
                sum(
                    ipp.stor_duration_new(s) * value.(pi_C_charge[y, p_star, s, z, d, t]) * value.(x_stor_C_mce[y, s, z])
                    for s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    value.(kappa_C[y, p, s, z, d, t]) * ipp.x_stor_C_cumu(p, s, z) 
                    for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) +
                sum(
                    value.(kappa_C[y, p_star, s, z, d, t]) * value.(x_stor_C_mce[y, s, z])
                    for s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) + 
            sum(value.(psi_E[y, p, s, z, d, model_data.index_t.elements[1]]) * ipp.initial_energy_existing_my(y, p, s, z, d)
                for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d) - 
            sum(value.(pi_E_discharge[y, p, s, z, d, model_data.index_t.elements[1]]) * ipp.initial_energy_existing_my(y, p, s, z, d)
                for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d) + 
            sum(value.(pi_E_charge[y, p, s, z, d, model_data.index_t.elements[1]]) * ipp.initial_energy_existing_my(y, p, s, z, d)
                for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d) + 
            sum(value.(psi_C[y, p, s, z, d, model_data.index_t.elements[1]]) * ipp.initial_energy_new_my(y, p, s, z, d)
                for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d) - 
            sum(value.(pi_C_discharge[y, p, s, z, d, model_data.index_t.elements[1]]) * ipp.initial_energy_new_my(y, p, s, z, d)
                for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d) + 
            sum(value.(pi_C_charge[y, p, s, z, d, model_data.index_t.elements[1]]) * ipp.initial_energy_new_my(y, p, s, z, d)
                for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d)
        end
    end
    lower_level_duality_gap = (lower_level_primal_obj .- lower_level_dual_obj) ./ abs.(lower_level_primal_obj)
    @info "lower level primal obj is $(lower_level_primal_obj)"
    @info "lower level dual obj is $(lower_level_dual_obj)"
    @info "lower level duality gap is $(lower_level_duality_gap)"

    # return compute_difference_percentage_one_norm([
    #     (x_R_before, ipp.x_R_my),
    #     (x_C_before, ipp.x_C_my),
    # ])
    return compute_difference_percentage_maximum_one_norm([
        (x_R_before, ipp.x_R_my),
        (x_C_before, ipp.x_C_my),
    ])
end


function solve_agent_problem!(
    ipp::IPPGroup,
    ipp_opts::IPPOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{WM},
    agent_store::AgentStore,
    w_iter,
    window_length,
    jump_model,
    export_file_path,
    update_results::Bool
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
                window_length,
                jump_model
            )
        end
    end

    # report change in key variables from previous iteration to this one
    return diff
end

function save_results(
    ipps::IPPGroup,
    ipp_opts::AgentOptions,
    hem_opts::HEMOptions{WM},
    export_file_path::AbstractString,
)
    # Primal Variables
    save_param(
        ipps.y_E_my.values,
        [:Year, :IPP, :GenTech, :Zone, :Day, :Time],
        :Generation_MWh,
        joinpath(export_file_path, "y_E.csv"),
    )
    save_param(
        ipps.y_C_my.values,
        [:Year, :IPP, :GenTech, :Zone, :Day, :Time],
        :Generation_MWh,
        joinpath(export_file_path, "y_C.csv"),
    )
    save_param(
        ipps.x_R_my.values,
        [:Year, :IPP, :GenTech, :Zone],
        :Capacity_MW,
        joinpath(export_file_path, "x_R.csv"),
    )
    save_param(
        ipps.x_C_my.values,
        [:Year, :IPP, :GenTech, :Zone],
        :Capacity_MW,
        joinpath(export_file_path, "x_C.csv"),
    )
    save_param(
        ipps.x_stor_R_my.values,
        [:Year, :IPP, :StorTech, :Zone],
        :Capacity_MW,
        joinpath(export_file_path, "x_stor_R.csv"),
    )
    save_param(
        ipps.x_stor_C_my.values,
        [:Year, :IPP, :StorTech, :Zone],
        :Capacity_MW,
        joinpath(export_file_path, "x_stor_C.csv"),
    )
    save_param(
        ipps.LMP_my.values,
        [:Year, :Zone, :Day, :Time],
        :MarginalCost,
        joinpath(export_file_path, "LMP.csv"),
    )
    save_param(
        ipps.charge_E_my.values,
        [:Year, :IPP, :StorTech, :Zone, :Day, :Time],
        :Charge_MWh,
        joinpath(export_file_path, "charge_E.csv"),
    )
    save_param(
        ipps.discharge_E_my.values,
        [:Year, :IPP, :StorTech, :Zone, :Day, :Time],
        :Discharge_MWh,
        joinpath(export_file_path, "discharge_E.csv"),
    )
    save_param(
        ipps.charge_C_my.values,
        [:Year, :IPP, :StorTech, :Zone, :Day, :Time],
        :Charge_MWh,
        joinpath(export_file_path, "charge_C.csv"),
    )
    save_param(
        ipps.discharge_C_my.values,
        [:Year, :IPP, :StorTech, :Zone, :Day, :Time],
        :Discharge_MWh,
        joinpath(export_file_path, "discharge_C.csv"),
    )
    save_param(
        ipps.energy_E_my.values,
        [:Year, :IPP, :StorTech, :Zone, :Day, :Time],
        :Energy_MWh,
        joinpath(export_file_path, "energy_E.csv"),
    )
    save_param(
        ipps.energy_C_my.values,
        [:Year, :IPP, :StorTech, :Zone, :Day, :Time],
        :Energy_MWh,
        joinpath(export_file_path, "energy_C.csv"),
    )
    save_param(
        ipps.flow_my.values,
        [:Year, :Line, :Day, :Time],
        :Flow_MWh,
        joinpath(export_file_path, "flow.csv"),
    )
end

function welfare_calculation!(
    ipp::IPPGroup,
    ipp_opts::AgentOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{WM},
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
    for p in ipp.index_p, k in ipp.index_k_existing, z in model_data.index_z
        ipp.x_R_cumu(p, k, z, :) .= min(ipp.x_R_cumu(p, k, z) + ipp.x_R_my(first(model_data.index_y), p, k, z), ipp.x_E_my(p, z, k))
    end

    for p in ipp.index_p, k in ipp.index_k_new, z in model_data.index_z
        ipp.x_C_cumu(p, k, z, :) .= ipp.x_C_cumu(p, k, z) + ipp.x_C_my(first(model_data.index_y), p, k, z)
    end

    for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z
        ipp.x_stor_R_cumu(p, s, z, :) .= min(ipp.x_stor_R_cumu(p, s, z) + ipp.x_stor_R_my(first(model_data.index_y), p, s, z), ipp.x_stor_E_my(p, z, s))
    end

    for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z
        ipp.x_stor_C_cumu(p, s, z, :) .= ipp.x_stor_C_cumu(p, s, z) + ipp.x_stor_C_my(first(model_data.index_y), p, s, z)
    end
end
