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
    current_year::Symbol

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
    BIGM::ParamScalar{<:Integer}
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
    ITC_new_my::ParamArray # ITC of new capacity (%)
    ITCStor_new_my::ParamArray # ITC of new capacity (%)    
    rte_stor_E_my::ParamArray # round trip efficiency of existing storage ($/MW)
    rte_stor_C_my::ParamArray # round trip efficiency of new storage ($/MW)
    v_E_my::ParamArray # variable cost of existing capacity ($/MWh)
    v_C_my::ParamArray # variable cost of new capacity ($/MWh)
    PTC_existing::ParamArray # PTC of existing capacity ($/MWh)
    PTC_new_my::ParamArray # PTC of new capacity ($/MWh)
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

    # McCormic bounds data
    eta_param_vec::Vector{}
    lambda_param_vec::Vector{}
    theta_E_energy_param_vec::Vector{}
    theta_E_discharge_param_vec::Vector{}
    theta_E_charge_param_vec::Vector{}
    pi_E_charge_param_vec::Vector{}
    kappa_E_param_vec::Vector{}
    theta_C_energy_param_vec::Vector{}
    theta_C_discharge_param_vec::Vector{}
    theta_C_charge_param_vec::Vector{}
    pi_C_charge_param_vec::Vector{}
    kappa_C_param_vec::Vector{}

    # duality gap
    lower_level_duality_gap::ParamArray
end

mutable struct McCormickBounds
    eta_L::ParamArray
    eta_U::ParamArray
    lambda_L::ParamArray
    lambda_U::ParamArray
    theta_E_energy_L::ParamArray
    theta_E_energy_U::ParamArray
    theta_E_discharge_L::ParamArray
    theta_E_discharge_U::ParamArray
    theta_E_charge_L::ParamArray
    theta_E_charge_U::ParamArray
    pi_E_charge_L::ParamArray
    pi_E_charge_U::ParamArray
    kappa_E_L::ParamArray
    kappa_E_U::ParamArray
    theta_C_energy_L::ParamArray
    theta_C_energy_U::ParamArray
    theta_C_discharge_L::ParamArray
    theta_C_discharge_U::ParamArray
    theta_C_charge_L::ParamArray
    theta_C_charge_U::ParamArray
    pi_C_charge_L::ParamArray
    pi_C_charge_U::ParamArray
    kappa_C_L::ParamArray
    kappa_C_U::ParamArray
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
    # TODO: Remove if not being used. Also strange that this is combining FOM and CapEx
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
        first(model_data.index_y),
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
        ParamScalar("BIGM", 1000000),
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
            "ITC_new_my",
            input_filename,
            "ITCNewmy",
            index_k_new,
            [model_data.index_y],
        ),
        read_param(
            "ITCStor_new_my",
            input_filename,
            "ITCStorNewmy",
            index_stor_new,
            [model_data.index_y],
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
        read_param("PTC_existing_my", input_filename, "PTCOld", index_k_existing),
        read_param(
            "PTC_new_my",
            input_filename,
            "PTCNewmy",
            index_k_new,
            [model_data.index_y],
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
        [], [], [], [], [], [], [], [], [], [], [], [], 
        initialize_param("lower_level_duality_gap", model_data.index_y, value = 100.0),
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

function get_index_y(ipp, model_data, window_length)
    current_year, _ = get_current_year(ipp, model_data)
    return [Symbol(Int(yr)) for yr in current_year:(current_year + window_length - 1)]
end

"""
Lower level optimization results are used to set variable bounds for McCormick-envelope Relaxation
"""
function ipp_cap_lower(
    ipp, ipp_opts, model_data, delta_t, window_length, 
    customers, der_aggregator, green_developer, solver
)
    MPPDCMER_lower = get_new_jump_model(solver)

    index_y = get_index_y(ipp, model_data, window_length)

    if index_y != model_data.index_y.elements
        cust_year = customers.previous_year
        dera_year = der_aggregator.previous_year
        gd_year = green_developer.previous_year
    else
        cust_year = customers.current_year
        dera_year = der_aggregator.current_year
        gd_year = green_developer.current_year
    end

    # Variables
    @variable(
        MPPDCMER_lower,
        y_E_bounds[index_y, ipp.index_p, ipp.index_k_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower,
        y_C_bounds[index_y, ipp.index_p, ipp.index_k_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower,
        charge_E_bounds[index_y, ipp.index_p, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower,
        discharge_E_bounds[index_y, ipp.index_p, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower,
        charge_C_bounds[index_y, ipp.index_p, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower,
        discharge_C_bounds[index_y, ipp.index_p, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower,
        energy_E_bounds[index_y, ipp.index_p, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower,
        energy_C_bounds[index_y, ipp.index_p, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower,
        flow_bounds[index_y, ipp.index_l, model_data.index_d, model_data.index_t]
    )

    # Objective Function
    objective_function_lower = begin
        sum(
            sum(
                model_data.omega(d) * delta_t *
                ((ipp.v_E_my(y, p, k, z, d, t) - ipp.PTC_existing(k)) * y_E_bounds[y, p, k, z, d, t]) for
                d in model_data.index_d, t in model_data.index_t, z in model_data.index_z, k in ipp.index_k_existing, p in ipp.index_p
            ) + 
            sum(
                model_data.omega(d) * delta_t *
                ((ipp.v_C_my(y, p, k, z, d, t) - ipp.PTC_new_my(y, k)) * y_C_bounds[y, p, k, z, d, t]) for
                d in model_data.index_d, t in model_data.index_t, z in model_data.index_z, k in ipp.index_k_new, p in ipp.index_p
            )
            for y in index_y
        )
    end
    @objective(MPPDCMER_lower, Min, objective_function_lower)

    # Constraints
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
                customers.rho_DG(h, m, z, d, t) * customers.total_der_capacity_my_delay_update(y, z, h, m) for
                h in model_data.index_h, m in customers.index_m
            ) - 
            # remove aggregated behind-the-meter pv/storage generation/consumption since they're front-of-the-meter now
            sum(
                customers.rho_DG(h, m, z, d, t) * der_aggregator.aggregation_level(dera_year, z) *
                customers.total_pv_stor_capacity_my(cust_year, z, h, m) # this actually needs to be a delayed update as well, but all the if-else in customer_group function makes this difficult. However, as long as the initial capacity does not change year over year for the simulation period this is fine (such as in this case).
                for h in model_data.index_h, m in (:BTMStorage, :BTMPV)
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
        Eq_primal_feasible_supplydemandbalance_lower[y in index_y, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t],
        supply_demand_balance_lower(y, z, d, t) == 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_gen_max_E_lower[
            y in index_y,
            p in ipp.index_p,
            k in ipp.index_k_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        # adding this max to avoid rounding errors
        max(0.0, ipp.rho_E_my(p, k, z, d, t) * (
            ipp.x_E_my(p, z, k) - sum(
                ipp.x_R_my(Symbol(Int(y_symbol)), p, k, z) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):model_data.year(y)
            )
        )) - y_E_bounds[y, p, k, z, d, t] >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_gen_max_C_lower[
            y in index_y,
            p in ipp.index_p,
            k in ipp.index_k_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        ipp.rho_C_my(p, k, z, d, t) * (
            sum(
                ipp.x_C_my(Symbol(Int(y_symbol)), p, k, z) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):model_data.year(y)
            )
        ) - y_C_bounds[y, p, k, z, d, t] >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_flow_lower_lower[
            y in index_y,
            l in ipp.index_l,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        flow_bounds[y, l, d, t] - ipp.trans_capacity(l, :min) >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_flow_upper_lower[
            y in index_y,
            l in ipp.index_l,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        ipp.trans_capacity(l, :max) - flow_bounds[y, l, d, t] >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_energy_E_lower[
            y in index_y,
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
            y in index_y,
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
            y in index_y,
            p in ipp.index_p,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        ipp.stor_duration_existing(s) * (
            ipp.x_stor_E_my(p, z, s) - sum(
                ipp.x_stor_R_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):model_data.year(y)
            )
        ) - energy_E_bounds[y, p, s, z, d, t] >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_discharge_upper_bound_E_lower[
            y in index_y,
            p in ipp.index_p,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        ipp.rte_stor_E_my(y, p, z, s) * (
            ipp.x_stor_E_my(p, z, s) - sum(
                ipp.x_stor_R_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):model_data.year(y)
            )
        ) - discharge_E_bounds[y, p, s, z, d, t] >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_charge_upper_bound_E_lower[
            y in index_y,
            p in ipp.index_p,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        ipp.x_stor_E_my(p, z, s) - sum(
            ipp.x_stor_R_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
            model_data.year(first(model_data.index_y_fix)):model_data.year(y)
        ) - charge_E_bounds[y, p, s, z, d, t] >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_discharge_energy_upper_bound_E_lower[
            y in index_y,
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
            y in index_y,
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
            y in index_y,
            p in ipp.index_p,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t.elements[2:end],
        ],
        ipp.stor_duration_existing(s) * (
            ipp.x_stor_E_my(p, z, s) - sum(
                ipp.x_stor_R_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):model_data.year(y)
            )
        ) -
        energy_E_bounds[y, p, s, z, d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)-delta_t), model_data.time.values)][1]] -
        charge_E_bounds[y, p, s, z, d, t] * delta_t >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_charge_energy_upper_bound_E_0_lower[
            y in index_y,
            p in ipp.index_p,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in [model_data.index_t.elements[1]],
        ],
        ipp.stor_duration_existing(s) * (
            ipp.x_stor_E_my(p, z, s) - sum(
                ipp.x_stor_R_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):model_data.year(y)
            )
        ) -
        ipp.initial_energy_existing_my(y, p, s, z, d) - charge_E_bounds[y, p, s, z, d, t] * delta_t >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_charge_discharge_upper_bound_E_lower[
            y in index_y,
            p in ipp.index_p,
            s in ipp.index_stor_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t.elements,
        ],
        ipp.x_stor_E_my(p, z, s) - sum(
                ipp.x_stor_R_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):model_data.year(y)
        ) - 
        charge_E_bounds[y, p, s, z, d, t] - discharge_E_bounds[y, p, s, z, d, t] / ipp.rte_stor_E_my(y, p, z, s) >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_energy_C_lower[
            y in index_y,
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
            y in index_y,
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
            y in index_y,
            p in ipp.index_p,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        ipp.stor_duration_new(s) * (
            sum(
                ipp.x_stor_C_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):model_data.year(y)
            )
        ) - energy_C_bounds[y, p, s, z, d, t] >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_discharge_upper_bound_C_lower[
            y in index_y,
            p in ipp.index_p,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        ipp.rte_stor_C_my(y, p, z, s) * (
            sum(
                ipp.x_stor_C_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):model_data.year(y)
            )
        ) - discharge_C_bounds[y, p, s, z, d, t] >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_charge_upper_bound_C_lower[
            y in index_y,
            p in ipp.index_p,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        sum(
            ipp.x_stor_C_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
            model_data.year(first(model_data.index_y_fix)):model_data.year(y)
        ) - charge_C_bounds[y, p, s, z, d, t] >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_discharge_energy_upper_bound_C_lower[
            y in index_y,
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
            y in index_y,
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
            y in index_y,
            p in ipp.index_p,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t.elements[2:end],
        ],
        ipp.stor_duration_new(s) * (
            sum(
                ipp.x_stor_C_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):model_data.year(y)
            )
        ) -
        energy_C_bounds[y, p, s, z, d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)-delta_t), model_data.time.values)][1]] -
        charge_C_bounds[y, p, s, z, d, t] * delta_t >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_charge_energy_upper_bound_C_0_lower[
            y in index_y,
            p in ipp.index_p,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in [model_data.index_t.elements[1]],
        ],
        ipp.stor_duration_new(s) * (
            sum(
                ipp.x_stor_C_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                model_data.year(first(model_data.index_y_fix)):model_data.year(y)
            )
        ) -
        ipp.initial_energy_new_my(y, p, s, z, d) - charge_C_bounds[y, p, s, z, d, t] * delta_t >= 0
    )

    @constraint(
        MPPDCMER_lower,
        Eq_primal_feasible_charge_discharge_upper_bound_C_lower[
            y in index_y,
            p in ipp.index_p,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t.elements,
        ],
        sum(
            ipp.x_stor_C_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
            model_data.year(first(model_data.index_y_fix)):model_data.year(y)
        ) - charge_C_bounds[y, p, s, z, d, t] - discharge_C_bounds[y, p, s, z, d, t] / ipp.rte_stor_C_my(y, p, z, s) >= 0
    )

    @info("MPPDC lower level primal")
    TimerOutputs.@timeit HEM_TIMER "optimize! MPPDC McCormic Envelope Relaxation lower level primal" begin
        optimize!(MPPDCMER_lower)
    end

    return MPPDCMER_lower
end

"""
Lower level dual optimization results are used to set variable bounds for McCormick-envelope Relaxation
"""
function ipp_cap_lower_dual(
    ipp, ipp_opts, model_data, delta_t, window_length, 
    customers, der_aggregator, green_developer, solver
)
    MPPDCMER_lower_dual = get_new_jump_model(solver)

    index_y = get_index_y(ipp, model_data, window_length)

    if index_y != model_data.index_y.elements
        cust_year = customers.previous_year
        dera_year = der_aggregator.previous_year
        gd_year = green_developer.previous_year
    else
        cust_year = customers.current_year
        dera_year = der_aggregator.current_year
        gd_year = green_developer.current_year
    end

    # Variables
    @variable(
        MPPDCMER_lower_dual, 
        miu_lower[index_y, model_data.index_z, model_data.index_d, model_data.index_t]
    )
    @variable(
        MPPDCMER_lower_dual,
        eta_lower[index_y, ipp.index_p, ipp.index_k_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower_dual,
        lambda_lower[index_y, ipp.index_p, ipp.index_k_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )

    @variable(
        MPPDCMER_lower_dual,
        iota_min_lower[index_y, ipp.index_l, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower_dual,
        iota_max_lower[index_y, ipp.index_l, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower_dual,
        psi_E_lower[index_y, ipp.index_p, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t]
    )
    @variable(
        MPPDCMER_lower_dual,
        theta_E_energy_lower[index_y, ipp.index_p, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower_dual,
        theta_E_discharge_lower[index_y, ipp.index_p, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower_dual,
        theta_E_charge_lower[index_y, ipp.index_p, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower_dual,
        pi_E_discharge_lower[index_y, ipp.index_p, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower_dual,
        pi_E_charge_lower[index_y, ipp.index_p, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower_dual,
        kappa_E_lower[index_y, ipp.index_p, ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower_dual,
        psi_C_lower[index_y, ipp.index_p, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t]
    )
    @variable(
        MPPDCMER_lower_dual,
        theta_C_energy_lower[index_y, ipp.index_p, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower_dual,
        theta_C_discharge_lower[index_y, ipp.index_p, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower_dual,
        theta_C_charge_lower[index_y, ipp.index_p, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower_dual,
        pi_C_discharge_lower[index_y, ipp.index_p, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower_dual,
        pi_C_charge_lower[index_y, ipp.index_p, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        MPPDCMER_lower_dual,
        kappa_C_lower[index_y, ipp.index_p, ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )

    # Objective
    objective_function_lower_dual = begin
        sum(
            sum(
                miu_lower[y, z, d, t] * (
                    sum(
                        customers.gamma(z, h) * customers.d_my(y, h, z, d, t) for h in model_data.index_h
                    ) + ipp.eximport_my(y, z, d, t) -
                    # total DG generation at time t
                    sum(
                        customers.rho_DG(h, m, z, d, t) * customers.total_der_capacity_my_delay_update(y, z, h, m) for
                        h in model_data.index_h, m in customers.index_m
                    ) +
                    # remove aggregated behind-the-meter pv/storage generation/consumption since they're front-of-the-meter now
                    sum(
                        customers.rho_DG(h, m, z, d, t) * der_aggregator.aggregation_level(dera_year, z) * 
                        customers.total_pv_stor_capacity_my(cust_year, z, h, m) for h in model_data.index_h, m in (:BTMStorage, :BTMPV)
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
                    # adding this max to avoid rounding errors
                    max(0.0, ipp.rho_E_my(p, k, z, d, t) * (
                        ipp.x_E_my(p, z, k) - sum(
                            ipp.x_R_my(Symbol(Int(y_symbol)), p, k, z) for y_symbol in
                            model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                        )
                    )) for p in ipp.index_p, k in ipp.index_k_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    lambda_lower[y, p, k, z, d, t] *
                    ipp.rho_C_my(p, k, z, d, t) * (
                        sum(
                            ipp.x_C_my(Symbol(Int(y_symbol)), p, k, z) for y_symbol in
                            model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                        )
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
                        ipp.x_stor_E_my(p, z, s) - sum(
                            ipp.x_stor_R_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                            model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                        )
                    ) for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    theta_E_discharge_lower[y, p, s, z, d, t] *
                    ipp.rte_stor_E_my(y, p, z, s) * (
                        ipp.x_stor_E_my(p, z, s) - sum(
                            ipp.x_stor_R_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                            model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                        )
                    ) for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    theta_E_charge_lower[y, p, s, z, d, t] *
                    (
                        ipp.x_stor_E_my(p, z, s) - sum(
                            ipp.x_stor_R_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                            model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                        )
                    ) for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    pi_E_charge_lower[y, p, s, z, d, t] *
                    ipp.stor_duration_existing(s) * (
                        ipp.x_stor_E_my(p, z, s) - sum(
                            ipp.x_stor_R_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                            model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                        )
                    ) for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    kappa_E_lower[y, p, s, z, d, t] *
                    (
                        ipp.x_stor_E_my(p, z, s) - sum(
                            ipp.x_stor_R_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                            model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                        )
                    ) for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    theta_C_energy_lower[y, p, s, z, d, t] *
                    ipp.stor_duration_new(s) * (
                        sum(
                            ipp.x_stor_C_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                            model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                        )
                    )
                    for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    theta_C_discharge_lower[y, p, s, z, d, t] *
                    ipp.rte_stor_C_my(y, p, z, s) * (
                        sum(
                            ipp.x_stor_C_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                            model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                        )
                    ) 
                    for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    theta_C_charge_lower[y, p, s, z, d, t] * (
                        sum(
                            ipp.x_stor_C_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                            model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                        )
                    ) 
                    for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    pi_C_charge_lower[y, p, s, z, d, t] *
                    ipp.stor_duration_new(s) * (
                        sum(
                            ipp.x_stor_C_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                            model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                        )
                    ) 
                    for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    kappa_C_lower[y, p, s, z, d, t] * (
                        sum(
                            ipp.x_stor_C_my(Symbol(Int(y_symbol)), p, s, z) for y_symbol in
                            model_data.year(first(model_data.index_y_fix)):model_data.year(y)
                        )
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
            for y in index_y
        )
    end
    @objective(MPPDCMER_lower_dual, Max, objective_function_lower_dual)

    # Constraints
    @constraint(
        MPPDCMER_lower_dual,
        Eq_dual_feasible_y_E_lower[
            y in index_y,
            p in ipp.index_p,
            k in ipp.index_k_existing,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        model_data.omega(d) * delta_t * (ipp.v_E_my(y, p, k, z, d, t) - ipp.PTC_existing(k)) - miu_lower[y, z, d, t] + eta_lower[y, p, k, z, d, t] >= 0
    )
    @constraint(
        MPPDCMER_lower_dual,
        Eq_dual_feasible_y_C_lower[
            y in index_y,
            p in ipp.index_p,
            k in ipp.index_k_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        model_data.omega(d) * delta_t * (ipp.v_C_my(y, p, k, z, d, t) - ipp.PTC_new_my(y, k)) - miu_lower[y, z, d, t] + lambda_lower[y, p, k, z, d, t] >= 0
    )
    @constraint(
        MPPDCMER_lower_dual,
        Eq_dual_feasible_flow_lower[
            y in index_y,
            l in ipp.index_l,
            d in model_data.index_d,
            t in model_data.index_t,
        ],
        sum(miu_lower[y, z, d, t] * ipp.trans_topology(l, z) for z in model_data.index_z) - iota_min_lower[y, l, d, t] + iota_max_lower[y, l, d, t] == 0
    )
    @constraint(
        MPPDCMER_lower_dual,
        Eq_dual_feasible_discharge_E_lower[
            y in index_y,
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
            y in index_y,
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
            y in index_y,
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
            y in index_y,
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
            y in index_y,
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
            y in index_y,
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
            y in index_y,
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
            y in index_y,
            p in ipp.index_p,
            s in ipp.index_stor_new,
            z in model_data.index_z,
            d in model_data.index_d,
            t in [model_data.index_t.elements[end]],
        ],
        - psi_C_lower[y, p, s, z, d, t] + theta_C_energy_lower[y, p, s, z, d, t] >= 0
    )

    @info("MPPDC lower level dual")
    TimerOutputs.@timeit HEM_TIMER "optimize! MPPDC McCormic Envelope Relaxation lower level dual problem" begin
        optimize!(MPPDCMER_lower_dual)
    end

    return MPPDCMER_lower_dual
end

function get_mccormick_bounds(
    MPPDCMER_lower, ipp, p_star, model_data, w_iter; 
    bound_size=0.1, first_update=true, skip_lower=false
)

    if !skip_lower && ((termination_status(MPPDCMER_lower) == OPTIMAL) || (termination_status(MPPDCMER_lower) == LOCALLY_SOLVED))
        @info("Lower-level problem completed successfully. Calculate and save McCormick bound parameters.")

        diff = model_data.year(axes(MPPDCMER_lower[:y_E_bounds])[1][1]) - model_data.year(first(model_data.index_y))

        eta_param = initialize_param("eta_param", deepcopy(model_data.index_y), ipp.index_k_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        fill!(eta_param, NaN)
        for y in model_data.index_y, k in ipp.index_k_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
            y_minus = Symbol(model_data.year(y) + diff)
            eta_param(y, k, z, d, t, :) .= abs(dual.(MPPDCMER_lower[:Eq_primal_feasible_gen_max_E_lower][y_minus, p_star, k, z, d, t]))
        end
        first_update ? push!(ipp.eta_param_vec, eta_param) : ipp.eta_param_vec[end] = eta_param

        lambda_param = initialize_param("lambda_param", deepcopy(model_data.index_y), ipp.index_k_new, model_data.index_z, model_data.index_d, model_data.index_t)
        fill!(lambda_param, NaN)
        for y in model_data.index_y, k in ipp.index_k_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
            y_minus = Symbol(model_data.year(y) + diff)
            lambda_param(y, k, z, d, t, :) .= abs(dual.(MPPDCMER_lower[:Eq_primal_feasible_gen_max_C_lower][y_minus, p_star, k, z, d, t]))
        end
        first_update ? push!(ipp.lambda_param_vec, lambda_param) : ipp.lambda_param_vec[end] = lambda_param

        theta_E_energy_param = initialize_param("theta_E_energy_param", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        fill!(theta_E_energy_param, NaN)
        for y in model_data.index_y, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
            y_minus = Symbol(model_data.year(y) + diff)
            theta_E_energy_param(y, s, z, d, t, :) .= abs(dual.(MPPDCMER_lower[:Eq_primal_feasible_energy_upper_bound_E_lower][y_minus, p_star, s, z, d, t]))
        end
        first_update ? push!(ipp.theta_E_energy_param_vec, theta_E_energy_param) : ipp.theta_E_energy_param_vec[end] = theta_E_energy_param

        theta_E_discharge_param = initialize_param("theta_E_discharge_param", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        fill!(theta_E_discharge_param, NaN)
        for y in model_data.index_y, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
            y_minus = Symbol(model_data.year(y) + diff)
            theta_E_discharge_param(y, s, z, d, t, :) .= abs(dual.(MPPDCMER_lower[:Eq_primal_feasible_discharge_upper_bound_E_lower][y_minus, p_star, s, z, d, t]))
        end
        first_update ? push!(ipp.theta_E_discharge_param_vec, theta_E_discharge_param) : ipp.theta_E_discharge_param_vec[end] = theta_E_discharge_param

        theta_E_charge_param = initialize_param("theta_E_charge_param", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        fill!(theta_E_charge_param, NaN)
        for y in model_data.index_y, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
            y_minus = Symbol(model_data.year(y) + diff)
            theta_E_charge_param(y, s, z, d, t, :) .= abs(dual.(MPPDCMER_lower[:Eq_primal_feasible_charge_upper_bound_E_lower][y_minus, p_star, s, z, d, t]))
        end
        first_update ? push!(ipp.theta_E_charge_param_vec, theta_E_charge_param) : ipp.theta_E_charge_param_vec[end] = theta_E_charge_param

        pi_E_charge_param = initialize_param("pi_E_charge_param", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        fill!(pi_E_charge_param, NaN)
        for y in model_data.index_y, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t.elements[2:end]
            y_minus = Symbol(model_data.year(y) + diff)
            pi_E_charge_param(y, s, z, d, t, :) .= abs(dual.(MPPDCMER_lower[:Eq_primal_feasible_charge_energy_upper_bound_E_lower][y_minus, p_star, s, z, d, t]))
        end
        for y in model_data.index_y, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in [model_data.index_t.elements[1]]
            y_minus = Symbol(model_data.year(y) + diff)
            pi_E_charge_param(y, s, z, d, t, :) .= abs(dual.(MPPDCMER_lower[:Eq_primal_feasible_charge_energy_upper_bound_E_0_lower][y_minus, p_star, s, z, d, t]))
        end
        first_update ? push!(ipp.pi_E_charge_param_vec, pi_E_charge_param) : ipp.pi_E_charge_param_vec[end] = pi_E_charge_param

        kappa_E_param = initialize_param("kappa_E_param", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        fill!(kappa_E_param, NaN)
        for y in model_data.index_y, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
            y_minus = Symbol(model_data.year(y) + diff)
            kappa_E_param(y, s, z, d, t, :) .= abs(dual.(MPPDCMER_lower[:Eq_primal_feasible_charge_discharge_upper_bound_E_lower][y_minus, p_star, s, z, d, t]))
        end
        first_update ? push!(ipp.kappa_E_param_vec, kappa_E_param) : ipp.kappa_E_param_vec[end] = kappa_E_param

        theta_C_energy_param = initialize_param("theta_C_energy_param", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        fill!(theta_C_energy_param, NaN)
        for y in model_data.index_y, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
            y_minus = Symbol(model_data.year(y) + diff)
            theta_C_energy_param(y, s, z, d, t, :) .= abs(dual.(MPPDCMER_lower[:Eq_primal_feasible_energy_upper_bound_C_lower][y_minus, p_star, s, z, d, t]))
        end
        first_update ? push!(ipp.theta_C_energy_param_vec, theta_C_energy_param) : ipp.theta_C_energy_param_vec[end] = theta_C_energy_param

        theta_C_discharge_param = initialize_param("theta_C_discharge_param", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        fill!(theta_C_discharge_param, NaN)
        for y in model_data.index_y, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
            y_minus = Symbol(model_data.year(y) + diff)
            theta_C_discharge_param(y, s, z, d, t, :) .= abs(dual.(MPPDCMER_lower[:Eq_primal_feasible_discharge_upper_bound_C_lower][y_minus, p_star, s, z, d, t]))
        end
        first_update ? push!(ipp.theta_C_discharge_param_vec, theta_C_discharge_param) : ipp.theta_C_discharge_param_vec[end] = theta_C_discharge_param

        theta_C_charge_param = initialize_param("theta_C_charge_param", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        fill!(theta_C_charge_param, NaN)
        for y in model_data.index_y, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
            y_minus = Symbol(model_data.year(y) + diff)
            theta_C_charge_param(y, s, z, d, t, :) .= abs(dual.(MPPDCMER_lower[:Eq_primal_feasible_charge_upper_bound_C_lower][y_minus, p_star, s, z, d, t]))
        end
        first_update ? push!(ipp.theta_C_charge_param_vec, theta_C_charge_param) : ipp.theta_C_charge_param_vec[end] = theta_C_charge_param

        pi_C_charge_param = initialize_param("pi_C_charge_param", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        fill!(pi_C_charge_param, NaN)
        for y in model_data.index_y, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t.elements[2:end]
            y_minus = Symbol(model_data.year(y) + diff)
            pi_C_charge_param(y, s, z, d, t, :) .= abs(dual.(MPPDCMER_lower[:Eq_primal_feasible_charge_energy_upper_bound_C_lower][y_minus, p_star, s, z, d, t]))
        end
        for y in model_data.index_y, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in [model_data.index_t.elements[1]]
            y_minus = Symbol(model_data.year(y) + diff)
            pi_C_charge_param(y, s, z, d, t, :) .= abs(dual.(MPPDCMER_lower[:Eq_primal_feasible_charge_energy_upper_bound_C_0_lower][y_minus, p_star, s, z, d, t]))
        end
        first_update ? push!(ipp.pi_C_charge_param_vec, pi_C_charge_param) : ipp.pi_C_charge_param_vec[end] = pi_C_charge_param

        kappa_C_param = initialize_param("kappa_C_param", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        fill!(kappa_C_param, NaN)
        for y in model_data.index_y, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
            y_minus = Symbol(model_data.year(y) + diff)
            kappa_C_param(y, s, z, d, t, :) .= abs(dual.(MPPDCMER_lower[:Eq_primal_feasible_charge_discharge_upper_bound_C_lower][y_minus, p_star, s, z, d, t]))
        end
        first_update ? push!(ipp.kappa_C_param_vec, kappa_C_param) : ipp.kappa_C_param_vec[end] = kappa_C_param

        first_update = false
    else
        why = skip_lower ? "skipped" : "failed"
        @info("Lower-level problem $(why). Use previous McCormick bound parameters.")

        eta_param = initialize_param(
            "eta_param", deepcopy(model_data.index_y), 
            ipp.index_k_existing, model_data.index_z, model_data.index_d, model_data.index_t
        )

        for (i, y) in enumerate(model_data.index_y)
            y_before = ipp.eta_param_vec[end].dims[1][i]

            for k in ipp.index_k_existing
                for z in model_data.index_z
                    for d in model_data.index_d
                        for t in model_data.index_t
                            eta_param(y, k, z, d, t, :) .= ipp.eta_param_vec[end](y_before, k, z, d, t)
                        end
                    end
                end
            end
        end

        lambda_param = initialize_param("lambda_param", deepcopy(model_data.index_y), ipp.index_k_new, model_data.index_z, model_data.index_d, model_data.index_t)

        for (i, y) in enumerate(model_data.index_y)
            y_before = ipp.eta_param_vec[end].dims[1][i]

            for k in ipp.index_k_new
                for z in model_data.index_z
                    for d in model_data.index_d
                        for t in model_data.index_t
                            lambda_param(y, k, z, d, t, :) .= ipp.lambda_param_vec[end](y_before, k, z, d, t)
                        end
                    end
                end
            end
        end

        theta_E_energy_param = initialize_param("theta_E_energy_param", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        theta_E_discharge_param = initialize_param("theta_E_discharge_param", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        theta_E_charge_param = initialize_param("theta_E_charge_param", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        pi_E_charge_param = initialize_param("pi_E_charge_param", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        kappa_E_param = initialize_param("kappa_E_param", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z, model_data.index_d, model_data.index_t)
        
        for (i, y) in enumerate(model_data.index_y)
            y_before = ipp.eta_param_vec[end].dims[1][i]

            for k in ipp.index_stor_existing
                for z in model_data.index_z
                    for d in model_data.index_d
                        for t in model_data.index_t
                            theta_E_energy_param(y, k, z, d, t, :) .= ipp.theta_E_energy_param_vec[end](y_before, k, z, d, t)
                            theta_E_discharge_param(y, k, z, d, t, :) .= ipp.theta_E_discharge_param_vec[end](y_before, k, z, d, t)
                            theta_E_charge_param(y, k, z, d, t, :) .= ipp.theta_E_charge_param_vec[end](y_before, k, z, d, t)
                            pi_E_charge_param(y, k, z, d, t, :) .= ipp.pi_E_charge_param_vec[end](y_before, k, z, d, t)
                            kappa_E_param(y, k, z, d, t, :) .= ipp.kappa_E_param_vec[end](y_before, k, z, d, t)
                        end
                    end
                end
            end
        end

        theta_C_energy_param = initialize_param("theta_C_energy_param", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        theta_C_discharge_param = initialize_param("theta_C_discharge_param", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        theta_C_charge_param = initialize_param("theta_C_charge_param", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        pi_C_charge_param = initialize_param("pi_C_charge_param", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)
        kappa_C_param = initialize_param("kappa_C_param", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z, model_data.index_d, model_data.index_t)

        for (i, y) in enumerate(model_data.index_y)
            y_before = ipp.eta_param_vec[end].dims[1][i]

            for k in ipp.index_stor_new
                for z in model_data.index_z
                    for d in model_data.index_d
                        for t in model_data.index_t
                            theta_C_energy_param(y, k, z, d, t, :) .= ipp.theta_C_energy_param_vec[end](y_before, k, z, d, t)
                            theta_C_discharge_param(y, k, z, d, t, :) .= ipp.theta_C_discharge_param_vec[end](y_before, k, z, d, t)
                            theta_C_charge_param(y, k, z, d, t, :) .= ipp.theta_C_charge_param_vec[end](y_before, k, z, d, t)
                            pi_C_charge_param(y, k, z, d, t, :) .= ipp.pi_C_charge_param_vec[end](y_before, k, z, d, t)
                            kappa_C_param(y, k, z, d, t, :) .= ipp.kappa_C_param_vec[end](y_before, k, z, d, t)
                        end
                    end
                end
            end
        end
    end

    # a bigger bound_size makes the problem easier to solve, but also hurts the duality gap.
    @info("Calcuate McCormick Bounds using bound_size $(bound_size)")

    U_if_zero = 100.0
    L_if_zero = 0.0

    eta_L = initialize_param("eta_L", deepcopy(model_data.index_y), ipp.index_k_existing, model_data.index_z, model_data.index_d, model_data.index_t)
    eta_U = initialize_param("eta_U", deepcopy(model_data.index_y), ipp.index_k_existing, model_data.index_z, model_data.index_d, model_data.index_t)

    for y in model_data.index_y
        for k in ipp.index_k_existing
            for z in model_data.index_z
                for d in model_data.index_d
                    for t in model_data.index_t
                        eta_U(y, k, z, d, t, :) .= eta_param(y, k, z, d, t) == 0.0 ? U_if_zero : eta_param(y, k, z, d, t) * (1.0 + bound_size)
                        eta_L(y, k, z, d, t, :) .= eta_param(y, k, z, d, t) == 0.0 ? L_if_zero : eta_param(y, k, z, d, t) * (1.0 - bound_size)
                    end
                end
            end
        end
    end

    lambda_L = initialize_param("lambda_L", deepcopy(model_data.index_y), ipp.index_k_new, model_data.index_z, model_data.index_d, model_data.index_t)
    lambda_U = initialize_param("lambda_U", deepcopy(model_data.index_y), ipp.index_k_new, model_data.index_z, model_data.index_d, model_data.index_t)
    
    for y in model_data.index_y
        for k in ipp.index_k_new
            for z in model_data.index_z
                for d in model_data.index_d
                    for t in model_data.index_t
                        lambda_U(y, k, z, d, t, :) .= lambda_param(y, k, z, d, t) == 0.0 ? U_if_zero : lambda_param(y, k, z, d, t) * (1.0 + bound_size)
                        lambda_L(y, k, z, d, t, :) .= lambda_param(y, k, z, d, t) == 0.0 ? L_if_zero : lambda_param(y, k, z, d, t) * (1.0 - bound_size)
                    end
                end
            end
        end
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
        for k in ipp.index_stor_existing
            for z in model_data.index_z
                for d in model_data.index_d
                    for t in model_data.index_t
                        theta_E_energy_U(y, k, z, d, t, :) .= theta_E_energy_param(y, k, z, d, t) == 0.0 ? U_if_zero : theta_E_energy_param(y, k, z, d, t) * (1.0 + bound_size)
                        theta_E_energy_L(y, k, z, d, t, :) .= theta_E_energy_param(y, k, z, d, t) == 0.0 ? L_if_zero : theta_E_energy_param(y, k, z, d, t) * (1.0 - bound_size)
                        theta_E_discharge_U(y, k, z, d, t, :) .= theta_E_discharge_param(y, k, z, d, t) == 0.0 ? U_if_zero : theta_E_discharge_param(y, k, z, d, t) * (1.0 + bound_size)
                        theta_E_discharge_L(y, k, z, d, t, :) .= theta_E_discharge_param(y, k, z, d, t) == 0.0 ? L_if_zero : theta_E_discharge_param(y, k, z, d, t) * (1.0 - bound_size)
                        theta_E_charge_U(y, k, z, d, t, :) .= theta_E_charge_param(y, k, z, d, t) == 0.0 ? U_if_zero : theta_E_charge_param(y, k, z, d, t) * (1.0 + bound_size)
                        theta_E_charge_L(y, k, z, d, t, :) .= theta_E_charge_param(y, k, z, d, t) == 0.0 ? L_if_zero : theta_E_charge_param(y, k, z, d, t) * (1.0 - bound_size)
                        pi_E_charge_U(y, k, z, d, t, :) .= pi_E_charge_param(y, k, z, d, t) == 0.0 ? U_if_zero : pi_E_charge_param(y, k, z, d, t) * (1.0 + bound_size)
                        pi_E_charge_L(y, k, z, d, t, :) .= pi_E_charge_param(y, k, z, d, t) == 0.0 ? L_if_zero : pi_E_charge_param(y, k, z, d, t) * (1.0 - bound_size)
                        kappa_E_U(y, k, z, d, t, :) .= kappa_E_param(y, k, z, d, t) == 0.0 ? U_if_zero : kappa_E_param(y, k, z, d, t) * (1.0 + bound_size)
                        kappa_E_L(y, k, z, d, t, :) .= kappa_E_param(y, k, z, d, t) == 0.0 ? L_if_zero : kappa_E_param(y, k, z, d, t) * (1.0 - bound_size)
                    end
                end
            end
        end
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
        for k in ipp.index_stor_new
            for z in model_data.index_z
                for d in model_data.index_d
                    for t in model_data.index_t
                        theta_C_energy_U(y, k, z, d, t, :) .= theta_C_energy_param(y, k, z, d, t) == 0.0 ? U_if_zero : theta_C_energy_param(y, k, z, d, t) * (1.0 + bound_size)
                        theta_C_energy_L(y, k, z, d, t, :) .= theta_C_energy_param(y, k, z, d, t) == 0.0 ? L_if_zero : theta_C_energy_param(y, k, z, d, t) * (1.0 - bound_size)
                        theta_C_discharge_U(y, k, z, d, t, :) .= theta_C_discharge_param(y, k, z, d, t) == 0.0 ? U_if_zero : theta_C_discharge_param(y, k, z, d, t) * (1.0 + bound_size)
                        theta_C_discharge_L(y, k, z, d, t, :) .= theta_C_discharge_param(y, k, z, d, t) == 0.0 ? L_if_zero : theta_C_discharge_param(y, k, z, d, t) * (1.0 - bound_size)
                        theta_C_charge_U(y, k, z, d, t, :) .= theta_C_charge_param(y, k, z, d, t) == 0.0 ? U_if_zero : theta_C_charge_param(y, k, z, d, t) * (1.0 + bound_size)
                        theta_C_charge_L(y, k, z, d, t, :) .= theta_C_charge_param(y, k, z, d, t) == 0.0 ? L_if_zero : theta_C_charge_param(y, k, z, d, t) * (1.0 - bound_size)
                        pi_C_charge_U(y, k, z, d, t, :) .= pi_C_charge_param(y, k, z, d, t) == 0.0 ? U_if_zero : pi_C_charge_param(y, k, z, d, t) * (1.0 + bound_size)
                        pi_C_charge_L(y, k, z, d, t, :) .= pi_C_charge_param(y, k, z, d, t) == 0.0 ? L_if_zero : pi_C_charge_param(y, k, z, d, t) * (1.0 - bound_size)
                        kappa_C_U(y, k, z, d, t, :) .= kappa_C_param(y, k, z, d, t) == 0.0 ? U_if_zero : kappa_C_param(y, k, z, d, t) * (1.0 + bound_size)
                        kappa_C_L(y, k, z, d, t, :) .= kappa_C_param(y, k, z, d, t) == 0.0 ? L_if_zero : kappa_C_param(y, k, z, d, t) * (1.0 - bound_size)
                    end
                end
            end
        end
    end

    return McCormickBounds(
        eta_L, eta_U,
        lambda_L, lambda_U, 
        theta_E_energy_L, theta_E_energy_U, 
        theta_E_discharge_L, theta_E_discharge_U, 
        theta_E_charge_L, theta_E_charge_U, 
        pi_E_charge_L, pi_E_charge_U, 
        kappa_E_L, kappa_E_U, 
        theta_C_energy_L, theta_C_energy_U, 
        theta_C_discharge_L, theta_C_discharge_U, 
        theta_C_charge_L, theta_C_charge_U, 
        pi_C_charge_L, pi_C_charge_U, 
        kappa_C_L, kappa_C_U
    ), first_update
end

function ipp_cap_upper(
    mcbnds, ipp, ipp_opts, p_star, model_data, delta_t,
    regulator, customers, der_aggregator, green_developer; 
    max_build=500.0, constraint_scaling = 1.0
)
    reg_year = regulator.current_year
    cust_year = customers.current_year
    dera_year = der_aggregator.current_year
    gd_year = green_developer.current_year

    WMDER_IPP = get_new_jump_model(ipp_opts.solvers["solve_agent_problem_ipp_mppdc"])

    # Prepare cumulative parameters
    X_R_cumu_L = initialize_param("x_R_L", deepcopy(model_data.index_y), ipp.index_k_existing, model_data.index_z)
    X_R_cumu_U = initialize_param("x_R_U", deepcopy(model_data.index_y), ipp.index_k_existing, model_data.index_z)
    fill!(X_R_cumu_U, NaN)
    X_R_stor_cumu_L = initialize_param("x_R_stor_L", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z)
    X_R_stor_cumu_U = initialize_param("x_R_stor_U", deepcopy(model_data.index_y), ipp.index_stor_existing, model_data.index_z)
    fill!(X_R_stor_cumu_U, NaN)
    for y in model_data.index_y
        for k in ipp.index_k_existing
            for z in model_data.index_z
                # need to have constraints in the optimization as well
                if k == Symbol("nuclear")
                    # X_R_cumu_U(y, k, z, :) .= min(ipp.x_E_my(p_star, z, k) - ipp.x_R_cumu(p_star, k, z), 1.0)
                    X_R_cumu_U(y, k, z, :) .= max(ipp.x_E_my(p_star, z, k) - ipp.x_R_cumu(p_star, k, z), 0.0)
                else
                    X_R_cumu_U(y, k, z, :) .= max(ipp.x_E_my(p_star, z, k) - ipp.x_R_cumu(p_star, k, z), 0.0)
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

    X_C_cumu_L = initialize_param("x_C_L", deepcopy(model_data.index_y), ipp.index_k_new, model_data.index_z)
    X_C_cumu_U = initialize_param("x_C_U", deepcopy(model_data.index_y), ipp.index_k_new, model_data.index_z)
    X_C_stor_cumu_L = initialize_param("x_C_stor_L", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z)
    X_C_stor_cumu_U = initialize_param("x_C_stor_U", deepcopy(model_data.index_y), ipp.index_stor_new, model_data.index_z)
    for y in model_data.index_y
        for k in ipp.index_k_new
            for z in model_data.index_z
                # need to have constraints in the optimization as well
                # this hard-coded number needs to change
                X_C_cumu_U(y, k, z, :) .= max_build
            end
        end
    end
    for y in model_data.index_y
        for s in ipp.index_stor_new
            for z in model_data.index_z
                # need to have constraints in the optimization as well
                # this hard-coded number needs to change
                X_C_stor_cumu_U(y, s, z, :) .= max_build
            end
        end
    end

    # Positive Variables
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
    @variable(
        WMDER_IPP, 
        miu[model_data.index_y, model_data.index_z, model_data.index_d, model_data.index_t]
    )
    @variable(
        WMDER_IPP,
        eta[model_data.index_y, ipp.index_p, ipp.index_k_existing, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )
    @variable(
        WMDER_IPP,
        lambda[model_data.index_y, ipp.index_p, ipp.index_k_new, model_data.index_z, model_data.index_d, model_data.index_t] >= 0
    )

    # Variables that Represent Bilinear Terms
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

    # Cumulative Retirement of p_star
    @variable(WMDER_IPP, x_R_mce[model_data.index_y, ipp.index_k_existing, model_data.index_z] >= 0)
    @variable(WMDER_IPP, x_stor_R_mce[model_data.index_y, ipp.index_stor_existing, model_data.index_z] >= 0)
    # cumulative investment of p_star
    @variable(WMDER_IPP, x_C_mce[model_data.index_y, ipp.index_k_new, model_data.index_z] >= 0)
    @variable(WMDER_IPP, x_stor_C_mce[model_data.index_y, ipp.index_stor_new, model_data.index_z] >= 0)

    # Bounds on Retirements and Builds
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

    # More Variables
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

    # Calculate Parameters

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
                customers.rho_DG(h, m, z, d, t) * customers.total_der_capacity_my(cust_year, z, h, m) for
                h in model_data.index_h, m in customers.index_m
            ) +
            # remove aggregated behind-the-meter pv/storage generation/consumption since they're front-of-the-meter now
            sum(
                customers.rho_DG(h, m, z, d, t) * der_aggregator.aggregation_level(dera_year, z) *
                customers.total_pv_stor_capacity_my(cust_year, z, h, m) for h in model_data.index_h, m in (:BTMStorage, :BTMPV)
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

    # Objective Function
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
            ipp.pvf_cap(y, p_star) * sum(
                (1.0 - ipp.ITC_new_my(y, k)) * ipp.CapEx_my(y, p_star, z, k) * x_C[y, k, z] 
                for k in ipp.index_k_new, z in model_data.index_z) -
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
            ipp.pvf_cap(y, p_star) * sum(
                (1.0 - ipp.ITCStor_new_my(y, s)) * ipp.CapEx_stor_my(y, p_star, z, s) * x_stor_C[y, s, z] 
                for s in ipp.index_stor_new, z in model_data.index_z) +
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
                            customers.rho_DG(h, m, z, d, t) * customers.total_der_capacity_my(cust_year, z, h, m) for
                            h in model_data.index_h, m in customers.index_m
                        ) +
                        # remove aggregated behind-the-meter pv/storage generation/consumption since they're front-of-the-meter now
                        sum(
                            customers.rho_DG(h, m, z, d, t) * der_aggregator.aggregation_level(dera_year, z) * 
                            customers.total_pv_stor_capacity_my(cust_year, z, h, m) for h in model_data.index_h, m in (:BTMStorage, :BTMPV)
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
                        ((ipp.v_E_my(y, p, k, z, d, t) - ipp.PTC_existing(k)) * y_E[y, p, k, z, d, t]) for
                        d in model_data.index_d, t in model_data.index_t, z in model_data.index_z, k in ipp.index_k_existing, p in ipp.index_p
                    ) + 
                    sum(
                        model_data.omega(d) * delta_t *
                        ((ipp.v_C_my(y, p, k, z, d, t) - ipp.PTC_new_my(y, k)) * y_C[y, p, k, z, d, t]) for
                        d in model_data.index_d, t in model_data.index_t, z in model_data.index_z, k in ipp.index_k_new, p in ipp.index_p
                    )
                )
            )
            for y in model_data.index_y
        )
    end

    @objective(WMDER_IPP, Max, objective_function)

    # Constraints
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
        model_data.omega(d) * delta_t * (ipp.v_E_my(y, p, k, z, d, t) - ipp.PTC_existing(k)) - miu[y, z, d, t] + eta[y, p, k, z, d, t] >= 0
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
        model_data.omega(d) * delta_t * (ipp.v_C_my(y, p, k, z, d, t) - ipp.PTC_new_my(y, k)) - miu[y, z, d, t] + lambda[y, p, k, z, d, t] >= 0
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
                customers.rho_DG(h, m, z, d, t) * customers.total_der_capacity_my(cust_year, z, h, m) for
                h in model_data.index_h, m in customers.index_m
            ) -
            # remove aggregated behind-the-meter pv/storage generation/consumption since they're front-of-the-meter now
            sum(
                customers.rho_DG(h, m, z, d, t) * der_aggregator.aggregation_level(dera_year, z) *
                customers.total_pv_stor_capacity_my(cust_year, z, h, m) for h in model_data.index_h, m in (:BTMStorage, :BTMPV)
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
        #         ((ipp.v_E_my(y, p, k, t) - ipp.PTC_existing(k)) * y_E[y, p, k, t]) for
        #         t in model_data.index_t, k in ipp.index_k_existing, p in ipp.index_p
        #     ) + 
        #     sum(
        #         model_data.omega(t) *
        #         ((ipp.v_C_my(y, p, k, t) - ipp.PTC_new_my(y, k)) * y_C[y, p, k, t]) for
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
                ((ipp.v_E_my(y, p, k, z, d, t) - ipp.PTC_existing(k)) * y_E[y, p, k, z, d, t]) for
                d in model_data.index_d, t in model_data.index_t, z in model_data.index_z, k in ipp.index_k_existing, p in ipp.index_p
            ) + 
            sum(
                model_data.omega(d) * delta_t *
                ((ipp.v_C_my(y, p, k, z, d, t) - ipp.PTC_new_my(y, k)) * y_C[y, p, k, z, d, t]) for
                d in model_data.index_d, t in model_data.index_t, z in model_data.index_z, k in ipp.index_k_new, p in ipp.index_p
            ) == 
            sum(
                miu[y, z, d, t] * (
                    sum(
                        customers.gamma(z, h) * customers.d_my(y, h, z, d, t) for h in model_data.index_h
                    ) + ipp.eximport_my(y, z, d, t) -
                    # total DG generation at time t
                    sum(
                        customers.rho_DG(h, m, z, d, t) * customers.total_der_capacity_my(cust_year, z, h, m) for
                        h in model_data.index_h, m in customers.index_m
                    ) + 
                    # remove aggregated behind-the-meter pv/storage generation/consumption since they're front-of-the-meter now
                    sum(
                        customers.rho_DG(h, m, z, d, t) * der_aggregator.aggregation_level(dera_year, z) *
                        customers.total_pv_stor_capacity_my(cust_year, z, h, m) for h in model_data.index_h, m in (:BTMStorage, :BTMPV)
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
        mcbnds.eta_L(y, k, z, d, t) * x_R_mce[y, k, z] - mcbnds.eta_L(y, k, z, d, t) * X_R_cumu_L(y, k, z)) / constraint_scaling
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
        mce_eta_x_R_p_star[y, k, z, d, t] / constraint_scaling >= (- mcbnds.eta_U(y, k, z, d, t) * X_R_cumu_U(y, k, z) +
        mcbnds.eta_U(y, k, z, d, t) * x_R_mce[y, k, z] + eta[y, p_star, k, z, d, t] * X_R_cumu_U(y, k, z)) / constraint_scaling
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
        mcbnds.eta_L(y, k, z, d, t) * X_R_cumu_U(y, k, z) + mcbnds.eta_L(y, k, z, d, t) * x_R_mce[y, k, z]) / constraint_scaling
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
        mce_eta_x_R_p_star[y, k, z, d, t] / constraint_scaling <= (mcbnds.eta_U(y, k, z, d, t) * x_R_mce[y, k, z] -
        mcbnds.eta_U(y, k, z, d, t) * X_R_cumu_L(y, k, z) + X_R_cumu_L(y, k, z) * eta[y, p_star, k, z, d, t]) / constraint_scaling
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
        mcbnds.lambda_L(y, k, z, d, t) * x_C_mce[y, k, z] - mcbnds.lambda_L(y, k, z, d, t) * X_C_cumu_L(y, k, z)) / constraint_scaling
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
        mce_lambda_x_C_p_star[y, k, z, d, t] / constraint_scaling >= (- mcbnds.lambda_U(y, k, z, d, t) * X_C_cumu_U(y, k, z) +
        mcbnds.lambda_U(y, k, z, d, t) * x_C_mce[y, k, z] + lambda[y, p_star, k, z, d, t] * X_C_cumu_U(y, k, z)) / constraint_scaling
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
        mcbnds.lambda_L(y, k, z, d, t) * X_C_cumu_U(y, k, z) + mcbnds.lambda_L(y, k, z, d, t) * x_C_mce[y, k, z]) / constraint_scaling
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
        mce_lambda_x_C_p_star[y, k, z, d, t] / constraint_scaling <= (mcbnds.lambda_U(y, k, z, d, t) * x_C_mce[y, k, z] -
        mcbnds.lambda_U(y, k, z, d, t) * X_C_cumu_L(y, k, z) + X_C_cumu_L(y, k, z) * lambda[y, p_star, k, z, d, t]) / constraint_scaling
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
        mcbnds.theta_E_energy_L(y, s, z, d, t) * x_stor_R_mce[y, s, z] - mcbnds.theta_E_energy_L(y, s, z, d, t) * X_R_stor_cumu_L(y, s, z)) / constraint_scaling
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
        mce_theta_E_energy_x_stor_R_p_star[y, s, z, d, t] / constraint_scaling >= (- mcbnds.theta_E_energy_U(y, s, z, d, t) * X_R_stor_cumu_U(y, s, z) +
        mcbnds.theta_E_energy_U(y, s, z, d, t) * x_stor_R_mce[y, s, z] + theta_E_energy[y, p_star, s, z, d, t] * X_R_stor_cumu_U(y, s, z)) / constraint_scaling
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
        mcbnds.theta_E_energy_L(y, s, z, d, t) * X_R_stor_cumu_U(y, s, z) + mcbnds.theta_E_energy_L(y, s, z, d, t) * x_stor_R_mce[y, s, z]) / constraint_scaling
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
        mce_theta_E_energy_x_stor_R_p_star[y, s, z, d, t] / constraint_scaling <= (mcbnds.theta_E_energy_U(y, s, z, d, t) * x_stor_R_mce[y, s, z] -
        mcbnds.theta_E_energy_U(y, s, z, d, t) * X_R_stor_cumu_L(y, s, z) + X_R_stor_cumu_L(y, s, z) * theta_E_energy[y, p_star, s, z, d, t]) / constraint_scaling
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
        mcbnds.theta_E_discharge_L(y, s, z, d, t) * x_stor_R_mce[y, s, z] - mcbnds.theta_E_discharge_L(y, s, z, d, t) * X_R_stor_cumu_L(y, s, z)) / constraint_scaling
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
        mce_theta_E_discharge_x_stor_R_p_star[y, s, z, d, t] / constraint_scaling >= (- mcbnds.theta_E_discharge_U(y, s, z, d, t) * X_R_stor_cumu_U(y, s, z) +
        mcbnds.theta_E_discharge_U(y, s, z, d, t) * x_stor_R_mce[y, s, z] + theta_E_discharge[y, p_star, s, z, d, t] * X_R_stor_cumu_U(y, s, z)) / constraint_scaling
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
        mcbnds.theta_E_discharge_L(y, s, z, d, t) * X_R_stor_cumu_U(y, s, z) + mcbnds.theta_E_discharge_L(y, s, z, d, t) * x_stor_R_mce[y, s, z]) / constraint_scaling
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
        mce_theta_E_discharge_x_stor_R_p_star[y, s, z, d, t] / constraint_scaling <= (mcbnds.theta_E_discharge_U(y, s, z, d, t) * x_stor_R_mce[y, s, z] -
        mcbnds.theta_E_discharge_U(y, s, z, d, t) * X_R_stor_cumu_L(y, s, z) + X_R_stor_cumu_L(y, s, z) * theta_E_discharge[y, p_star, s, z, d, t]) / constraint_scaling
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
        mcbnds.theta_E_charge_L(y, s, z, d, t) * x_stor_R_mce[y, s, z] - mcbnds.theta_E_charge_L(y, s, z, d, t) * X_R_stor_cumu_L(y, s, z)) / constraint_scaling
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
        mce_theta_E_charge_x_stor_R_p_star[y, s, z, d, t] / constraint_scaling >= (- mcbnds.theta_E_charge_U(y, s, z, d, t) * X_R_stor_cumu_U(y, s, z) +
        mcbnds.theta_E_charge_U(y, s, z, d, t) * x_stor_R_mce[y, s, z] + theta_E_charge[y, p_star, s, z, d, t] * X_R_stor_cumu_U(y, s, z)) / constraint_scaling
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
        mcbnds.theta_E_charge_L(y, s, z, d, t) * X_R_stor_cumu_U(y, s, z) + mcbnds.theta_E_charge_L(y, s, z, d, t) * x_stor_R_mce[y, s, z]) / constraint_scaling
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
        mce_theta_E_charge_x_stor_R_p_star[y, s, z, d, t] / constraint_scaling <= (mcbnds.theta_E_charge_U(y, s, z, d, t) * x_stor_R_mce[y, s, z] -
        mcbnds.theta_E_charge_U(y, s, z, d, t) * X_R_stor_cumu_L(y, s, z) + X_R_stor_cumu_L(y, s, z) * theta_E_charge[y, p_star, s, z, d, t]) / constraint_scaling
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
        mcbnds.pi_E_charge_L(y, s, z, d, t) * x_stor_R_mce[y, s, z] - mcbnds.pi_E_charge_L(y, s, z, d, t) * X_R_stor_cumu_L(y, s, z)) / constraint_scaling
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
        mce_pi_E_charge_x_stor_R_p_star[y, s, z, d, t] / constraint_scaling >= (- mcbnds.pi_E_charge_U(y, s, z, d, t) * X_R_stor_cumu_U(y, s, z) +
        mcbnds.pi_E_charge_U(y, s, z, d, t) * x_stor_R_mce[y, s, z] + pi_E_charge[y, p_star, s, z, d, t] * X_R_stor_cumu_U(y, s, z)) / constraint_scaling
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
        mcbnds.pi_E_charge_L(y, s, z, d, t) * X_R_stor_cumu_U(y, s, z) + mcbnds.pi_E_charge_L(y, s, z, d, t) * x_stor_R_mce[y, s, z]) / constraint_scaling
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
        mce_pi_E_charge_x_stor_R_p_star[y, s, z, d, t] / constraint_scaling <= (mcbnds.pi_E_charge_U(y, s, z, d, t) * x_stor_R_mce[y, s, z] -
        mcbnds.pi_E_charge_U(y, s, z, d, t) * X_R_stor_cumu_L(y, s, z) + X_R_stor_cumu_L(y, s, z) * pi_E_charge[y, p_star, s, z, d, t]) / constraint_scaling
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
        mcbnds.kappa_E_L(y, s, z, d, t) * x_stor_R_mce[y, s, z] - mcbnds.kappa_E_L(y, s, z, d, t) * X_R_stor_cumu_L(y, s, z)) / constraint_scaling
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
        mce_kappa_E_x_stor_R_p_star[y, s, z, d, t] / constraint_scaling >= (- mcbnds.kappa_E_U(y, s, z, d, t) * X_R_stor_cumu_U(y, s, z) +
        mcbnds.kappa_E_U(y, s, z, d, t) * x_stor_R_mce[y, s, z] + kappa_E[y, p_star, s, z, d, t] * X_R_stor_cumu_U(y, s, z)) / constraint_scaling
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
        mcbnds.kappa_E_L(y, s, z, d, t) * X_R_stor_cumu_U(y, s, z) + mcbnds.kappa_E_L(y, s, z, d, t) * x_stor_R_mce[y, s, z]) / constraint_scaling
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
        mce_kappa_E_x_stor_R_p_star[y, s, z, d, t] / constraint_scaling <= (mcbnds.kappa_E_U(y, s, z, d, t) * x_stor_R_mce[y, s, z] -
        mcbnds.kappa_E_U(y, s, z, d, t) * X_R_stor_cumu_L(y, s, z) + X_R_stor_cumu_L(y, s, z) * kappa_E[y, p_star, s, z, d, t]) / constraint_scaling
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
        mcbnds.theta_C_energy_L(y, s, z, d, t) * x_stor_C_mce[y, s, z] - mcbnds.theta_C_energy_L(y, s, z, d, t) * X_C_stor_cumu_L(y, s, z)) / constraint_scaling
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
        mce_theta_C_energy_x_stor_C_p_star[y, s, z, d, t] / constraint_scaling >= (- mcbnds.theta_C_energy_U(y, s, z, d, t) * X_C_stor_cumu_U(y, s, z) +
        mcbnds.theta_C_energy_U(y, s, z, d, t) * x_stor_C_mce[y, s, z] + theta_C_energy[y, p_star, s, z, d, t] * X_C_stor_cumu_U(y, s, z)) / constraint_scaling
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
        mcbnds.theta_C_energy_L(y, s, z, d, t) * X_C_stor_cumu_U(y, s, z) + mcbnds.theta_C_energy_L(y, s, z, d, t) * x_stor_C_mce[y, s, z]) / constraint_scaling
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
        mce_theta_C_energy_x_stor_C_p_star[y, s, z, d, t] / constraint_scaling <= (mcbnds.theta_C_energy_U(y, s, z, d, t) * x_stor_C_mce[y, s, z] -
        mcbnds.theta_C_energy_U(y, s, z, d, t) * X_C_stor_cumu_L(y, s, z) + X_C_stor_cumu_L(y, s, z) * theta_C_energy[y, p_star, s, z, d, t]) / constraint_scaling
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
        mcbnds.theta_C_discharge_L(y, s, z, d, t) * x_stor_C_mce[y, s, z] - mcbnds.theta_C_discharge_L(y, s, z, d, t) * X_C_stor_cumu_L(y, s, z)) / constraint_scaling
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
        mce_theta_C_discharge_x_stor_C_p_star[y, s, z, d, t] / constraint_scaling >= (- mcbnds.theta_C_discharge_U(y, s, z, d, t) * X_C_stor_cumu_U(y, s, z) +
        mcbnds.theta_C_discharge_U(y, s, z, d, t) * x_stor_C_mce[y, s, z] + theta_C_discharge[y, p_star, s, z, d, t] * X_C_stor_cumu_U(y, s, z)) / constraint_scaling
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
        mcbnds.theta_C_discharge_L(y, s, z, d, t) * X_C_stor_cumu_U(y, s, z) + mcbnds.theta_C_discharge_L(y, s, z, d, t) * x_stor_C_mce[y, s, z]) / constraint_scaling
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
        mce_theta_C_discharge_x_stor_C_p_star[y, s, z, d, t] / constraint_scaling <= (mcbnds.theta_C_discharge_U(y, s, z, d, t) * x_stor_C_mce[y, s, z] -
        mcbnds.theta_C_discharge_U(y, s, z, d, t) * X_C_stor_cumu_L(y, s, z) + X_C_stor_cumu_L(y, s, z) * theta_C_discharge[y, p_star, s, z, d, t]) / constraint_scaling
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
        mcbnds.theta_C_charge_L(y, s, z, d, t) * x_stor_C_mce[y, s, z] - mcbnds.theta_C_charge_L(y, s, z, d, t) * X_C_stor_cumu_L(y, s, z)) / constraint_scaling
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
        mce_theta_C_charge_x_stor_C_p_star[y, s, z, d, t] / constraint_scaling >= (- mcbnds.theta_C_charge_U(y, s, z, d, t) * X_C_stor_cumu_U(y, s, z) +
        mcbnds.theta_C_charge_U(y, s, z, d, t) * x_stor_C_mce[y, s, z] + theta_C_charge[y, p_star, s, z, d, t] * X_C_stor_cumu_U(y, s, z)) / constraint_scaling
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
        mcbnds.theta_C_charge_L(y, s, z, d, t) * X_C_stor_cumu_U(y, s, z) + mcbnds.theta_C_charge_L(y, s, z, d, t) * x_stor_C_mce[y, s, z]) / constraint_scaling
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
        mce_theta_C_charge_x_stor_C_p_star[y, s, z, d, t] / constraint_scaling <= (mcbnds.theta_C_charge_U(y, s, z, d, t) * x_stor_C_mce[y, s, z] -
        mcbnds.theta_C_charge_U(y, s, z, d, t) * X_C_stor_cumu_L(y, s, z) + X_C_stor_cumu_L(y, s, z) * theta_C_charge[y, p_star, s, z, d, t]) / constraint_scaling
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
        mcbnds.pi_C_charge_L(y, s, z, d, t) * x_stor_C_mce[y, s, z] - mcbnds.pi_C_charge_L(y, s, z, d, t) * X_C_stor_cumu_L(y, s, z)) / constraint_scaling
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
        mce_pi_C_charge_x_stor_C_p_star[y, s, z, d, t] / constraint_scaling >= (- mcbnds.pi_C_charge_U(y, s, z, d, t) * X_C_stor_cumu_U(y, s, z) +
        mcbnds.pi_C_charge_U(y, s, z, d, t) * x_stor_C_mce[y, s, z] + pi_C_charge[y, p_star, s, z, d, t] * X_C_stor_cumu_U(y, s, z)) / constraint_scaling
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
        mcbnds.pi_C_charge_L(y, s, z, d, t) * X_C_stor_cumu_U(y, s, z) + mcbnds.pi_C_charge_L(y, s, z, d, t) * x_stor_C_mce[y, s, z]) / constraint_scaling
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
        mce_pi_C_charge_x_stor_C_p_star[y, s, z, d, t] / constraint_scaling <= (mcbnds.pi_C_charge_U(y, s, z, d, t) * x_stor_C_mce[y, s, z] -
        mcbnds.pi_C_charge_U(y, s, z, d, t) * X_C_stor_cumu_L(y, s, z) + X_C_stor_cumu_L(y, s, z) * pi_C_charge[y, p_star, s, z, d, t]) / constraint_scaling
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
        mcbnds.kappa_C_L(y, s, z, d, t) * x_stor_C_mce[y, s, z] - mcbnds.kappa_C_L(y, s, z, d, t) * X_C_stor_cumu_L(y, s, z)) / constraint_scaling
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
        mce_kappa_C_x_stor_C_p_star[y, s, z, d, t] / constraint_scaling >= (- mcbnds.kappa_C_U(y, s, z, d, t) * X_C_stor_cumu_U(y, s, z) +
        mcbnds.kappa_C_U(y, s, z, d, t) * x_stor_C_mce[y, s, z] + kappa_C[y, p_star, s, z, d, t] * X_C_stor_cumu_U(y, s, z)) / constraint_scaling
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
        mcbnds.kappa_C_L(y, s, z, d, t) * X_C_stor_cumu_U(y, s, z) + mcbnds.kappa_C_L(y, s, z, d, t) * x_stor_C_mce[y, s, z]) / constraint_scaling
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
        mce_kappa_C_x_stor_C_p_star[y, s, z, d, t] / constraint_scaling <= (mcbnds.kappa_C_U(y, s, z, d, t) * x_stor_C_mce[y, s, z] -
        mcbnds.kappa_C_U(y, s, z, d, t) * X_C_stor_cumu_L(y, s, z) + X_C_stor_cumu_L(y, s, z) * kappa_C[y, p_star, s, z, d, t]) / constraint_scaling
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

    @info("MPPDC upper level")
    TimerOutputs.@timeit HEM_TIMER "optimize! MPPDC with transmission and storage" begin
        optimize!(WMDER_IPP)
    end

    return WMDER_IPP
end

function ipp_calc_duality_gap(
    WMDER_IPP, ipp, p_star, model_data, delta_t,
    customers, der_aggregator, green_developer
)

    cust_year = customers.current_year
    dera_year = der_aggregator.current_year
    gd_year = green_developer.current_year

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
            ((ipp.v_E_my(y, p, k, z, d, t) - ipp.PTC_existing(k)) * value.(WMDER_IPP[:y_E][y, p, k, z, d, t])) for
            d in model_data.index_d, t in model_data.index_t, z in model_data.index_z, k in ipp.index_k_existing, p in ipp.index_p
        ) + 
        sum(
            model_data.omega(d) * delta_t *
            ((ipp.v_C_my(y, p, k, z, d, t) - ipp.PTC_new_my(y, k)) * value.(WMDER_IPP[:y_C][y, p, k, z, d, t])) for
            d in model_data.index_d, t in model_data.index_t, z in model_data.index_z, k in ipp.index_k_new, p in ipp.index_p
        )

        if length(ipp.index_p) >= 2
            # lower_level_dual_obj[1, y] = 
            # sum(
            #     value.(WMDER_IPP[:miu][y, t]) * (
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
            #         value.(WMDER_IPP[:eta][y, p, k, t]) *
            #         ipp.rho_E_my(p, k, t) *
            #         (
            #             ipp.x_E_my(p, k) - ipp.x_R_cumu(p, k)
            #         ) for t in model_data.index_t, k in ipp.index_k_existing, p in ipp.index_p
            #     ) + sum(
            #         value.(WMDER_IPP[:lambda][y, p, k, t]) *
            #         ipp.rho_C_my(p, k, t) *
            #         ipp.x_C_cumu(p, k) 
            #         for t in model_data.index_t, k in ipp.index_k_new, p in ipp.index_p
            #     )
            # ) + sum(
            #     value.(WMDER_IPP[:eta][y, p, k, t]) *
            #     ipp.rho_E_my(p, k, t) *
            #     sum(
            #         ipp.x_R_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
            #         model_data.year(first(model_data.index_y)):model_data.year(y)
            #     ) for t in model_data.index_t, k in ipp.index_k_existing,
            #     p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
            # ) - sum(
            #     value.(WMDER_IPP[:lambda][y, p, k, t]) *
            #     ipp.rho_C_my(p, k, t) *
            #     sum(
            #         ipp.x_C_my(Symbol(Int(y_symbol)), p, k) for y_symbol in
            #         model_data.year(first(model_data.index_y)):model_data.year(y)
            #     ) for t in model_data.index_t, k in ipp.index_k_new,
            #     p in ipp.index_p(Not(findall(x -> x == p_star, ipp.index_p)))
            # ) + sum(
            #     ipp.rho_E_my(p_star, k, t) *
            #     value.(WMDER_IPP[:eta][y, p_star, k, t]) *
            #     value.(WMDER_IPP[:x_R_mce][y, k]) for t in model_data.index_t, k in ipp.index_k_existing
            # ) - sum(
            #     ipp.rho_C_my(p_star, k, t) *
            #     value.(WMDER_IPP[:lambda][y, p_star, k, t]) *
            #     value.(WMDER_IPP[:x_C_mce][y, k]) for t in model_data.index_t, k in ipp.index_k_new
            # )
        else
            lower_level_dual_obj[1, y] = 
            sum(
                value.(WMDER_IPP[:miu][y, z, d, t]) * (
                    sum(
                        customers.gamma(z, h) * customers.d_my(y, h, z, d, t) for h in model_data.index_h
                    ) + ipp.eximport_my(y, z, d, t) -
                    # total DG generation at time t
                    sum(
                        customers.rho_DG(h, m, z, d, t) * customers.total_der_capacity_my(cust_year, z, h, m) for
                        h in model_data.index_h, m in customers.index_m
                    ) + 
                    # remove aggregated behind-the-meter pv/storage generation/consumption since they're front-of-the-meter now
                    sum(
                        customers.rho_DG(h, m, z, d, t) * der_aggregator.aggregation_level(dera_year, z) *
                        customers.total_pv_stor_capacity_my(cust_year, z, h, m) for h in model_data.index_h, m in (:BTMStorage, :BTMPV)
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
                    value.(WMDER_IPP[:eta][y, p, k, z, d, t]) *
                    ipp.rho_E_my(p, k, z, d, t) * (
                        ipp.x_E_my(p, z, k) - ipp.x_R_cumu(p, k, z)
                    ) for p in ipp.index_p, k in ipp.index_k_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) -
                sum(
                    ipp.rho_E_my(p_star, k, z, d, t) * value.(WMDER_IPP[:eta][y, p_star, k, z, d, t]) * value.(WMDER_IPP[:x_R_mce][y, k, z]) 
                    for k in ipp.index_k_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    value.(WMDER_IPP[:lambda][y, p, k, z, d, t]) *
                    ipp.rho_C_my(p, k, z, d, t) * ipp.x_C_cumu(p, k, z) 
                    for p in ipp.index_p, k in ipp.index_k_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) +
                sum(
                    ipp.rho_C_my(p_star, k, z, d, t) * value.(WMDER_IPP[:lambda][y, p_star, k, z, d, t]) * value.(WMDER_IPP[:x_C_mce][y, k, z]) 
                    for k in ipp.index_k_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) + 
            (
                sum(value.(WMDER_IPP[:iota_min][y, l, d, t]) * ipp.trans_capacity(l, :min) for l in ipp.index_l, d in model_data.index_d, t in model_data.index_t) - 
                sum(value.(WMDER_IPP[:iota_max][y, l, d, t]) * ipp.trans_capacity(l, :max) for l in ipp.index_l, d in model_data.index_d, t in model_data.index_t)
            ) - 
            (
                sum(
                    value.(WMDER_IPP[:theta_E_energy][y, p, s, z, d, t]) *
                    ipp.stor_duration_existing(s) * (
                        ipp.x_stor_E_my(p, z, s) - ipp.x_stor_R_cumu(p, s, z)
                    ) for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) -
                sum(
                    ipp.stor_duration_existing(s) * value.(WMDER_IPP[:theta_E_energy][y, p_star, s, z, d, t]) * value.(WMDER_IPP[:x_stor_R_mce][y, s, z]) 
                    for s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    value.(WMDER_IPP[:theta_E_discharge][y, p, s, z, d, t]) *
                    ipp.rte_stor_E_my(y, p, z, s) * (
                        ipp.x_stor_E_my(p, z, s) - ipp.x_stor_R_cumu(p, s, z)
                    ) for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) -
                sum(
                    ipp.rte_stor_E_my(y, p_star, z, s) * value.(WMDER_IPP[:theta_E_discharge][y, p_star, s, z, d, t]) * value.(WMDER_IPP[:x_stor_R_mce][y, s, z]) 
                    for s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    value.(WMDER_IPP[:theta_E_charge][y, p, s, z, d, t]) *
                    (
                        ipp.x_stor_E_my(p, z, s) - ipp.x_stor_R_cumu(p, s, z)
                    ) for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) -
                sum(
                    value.(WMDER_IPP[:theta_E_charge][y, p_star, s, z, d, t]) * value.(WMDER_IPP[:x_stor_R_mce][y, s, z]) 
                    for s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    value.(WMDER_IPP[:pi_E_charge][y, p, s, z, d, t]) *
                    ipp.stor_duration_existing(s) * (
                        ipp.x_stor_E_my(p, z, s) - ipp.x_stor_R_cumu(p, s, z)
                    ) for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) -
                sum(
                    ipp.stor_duration_existing(s) * value.(WMDER_IPP[:pi_E_charge][y, p_star, s, z, d, t]) * value.(WMDER_IPP[:x_stor_R_mce][y, s, z]) 
                    for s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    value.(WMDER_IPP[:kappa_E][y, p, s, z, d, t]) *
                    (
                        ipp.x_stor_E_my(p, z, s) - ipp.x_stor_R_cumu(p, s, z)
                    ) for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) -
                sum(
                    value.(WMDER_IPP[:kappa_E][y, p_star, s, z, d, t]) * value.(WMDER_IPP[:x_stor_R_mce][y, s, z]) 
                    for s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    value.(WMDER_IPP[:theta_C_energy][y, p, s, z, d, t]) *
                    ipp.stor_duration_new(s) * ipp.x_stor_C_cumu(p, s, z) 
                    for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) +
                sum(
                    ipp.stor_duration_new(s) * value.(WMDER_IPP[:theta_C_energy][y, p_star, s, z, d, t]) * value.(WMDER_IPP[:x_stor_C_mce][y, s, z])
                    for s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    value.(WMDER_IPP[:theta_C_discharge][y, p, s, z, d, t]) *
                    ipp.rte_stor_C_my(y, p, z, s) * ipp.x_stor_C_cumu(p, s, z) 
                    for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) +
                sum(
                    ipp.rte_stor_C_my(y, p_star, z, s) * value.(WMDER_IPP[:theta_C_discharge][y, p_star, s, z, d, t]) * value.(WMDER_IPP[:x_stor_C_mce][y, s, z])
                    for s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    value.(WMDER_IPP[:theta_C_charge][y, p, s, z, d, t]) * ipp.x_stor_C_cumu(p, s, z) 
                    for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) +
                sum(
                    value.(WMDER_IPP[:theta_C_charge][y, p_star, s, z, d, t]) * value.(WMDER_IPP[:x_stor_C_mce][y, s, z])
                    for s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    value.(WMDER_IPP[:pi_C_charge][y, p, s, z, d, t]) *
                    ipp.stor_duration_new(s) * ipp.x_stor_C_cumu(p, s, z) 
                    for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) +
                sum(
                    ipp.stor_duration_new(s) * value.(WMDER_IPP[:pi_C_charge][y, p_star, s, z, d, t]) * value.(WMDER_IPP[:x_stor_C_mce][y, s, z])
                    for s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) - 
            (
                sum(
                    value.(WMDER_IPP[:kappa_C][y, p, s, z, d, t]) * ipp.x_stor_C_cumu(p, s, z) 
                    for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                ) +
                sum(
                    value.(WMDER_IPP[:kappa_C][y, p_star, s, z, d, t]) * value.(WMDER_IPP[:x_stor_C_mce][y, s, z])
                    for s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
                )
            ) + 
            sum(value.(WMDER_IPP[:psi_E][y, p, s, z, d, model_data.index_t.elements[1]]) * ipp.initial_energy_existing_my(y, p, s, z, d)
                for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d) - 
            sum(value.(WMDER_IPP[:pi_E_discharge][y, p, s, z, d, model_data.index_t.elements[1]]) * ipp.initial_energy_existing_my(y, p, s, z, d)
                for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d) + 
            sum(value.(WMDER_IPP[:pi_E_charge][y, p, s, z, d, model_data.index_t.elements[1]]) * ipp.initial_energy_existing_my(y, p, s, z, d)
                for p in ipp.index_p, s in ipp.index_stor_existing, z in model_data.index_z, d in model_data.index_d) + 
            sum(value.(WMDER_IPP[:psi_C][y, p, s, z, d, model_data.index_t.elements[1]]) * ipp.initial_energy_new_my(y, p, s, z, d)
                for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d) - 
            sum(value.(WMDER_IPP[:pi_C_discharge][y, p, s, z, d, model_data.index_t.elements[1]]) * ipp.initial_energy_new_my(y, p, s, z, d)
                for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d) + 
            sum(value.(WMDER_IPP[:pi_C_charge][y, p, s, z, d, model_data.index_t.elements[1]]) * ipp.initial_energy_new_my(y, p, s, z, d)
                for p in ipp.index_p, s in ipp.index_stor_new, z in model_data.index_z, d in model_data.index_d)
        end
    end
    lower_level_duality_gap = (lower_level_primal_obj .- lower_level_dual_obj) ./ abs.(lower_level_primal_obj)
    @info "lower level primal obj is $(lower_level_primal_obj)"
    @info "lower level dual obj is $(lower_level_dual_obj)"
    @info "lower level duality gap is $(lower_level_duality_gap)"

    ipp.current_year = first(model_data.index_y)

    return maximum.(vec(Matrix(lower_level_duality_gap)))[1]
end

function ipp_cap_save_results(
    WMDER_IPP, ipp, p_star, model_data, delta_t,
    customers, der_aggregator, green_developer
)
    cust_year = customers.current_year
    dera_year = der_aggregator.current_year
    gd_year = green_developer.current_year

    # Save Results
    for y in model_data.index_y, k in ipp.index_k_existing, z in model_data.index_z
        ipp.x_R_my(y, p_star, k, z, :) .= value.(WMDER_IPP[:x_R][y, k, z])
    end
    for y in model_data.index_y, k in ipp.index_k_new, z in model_data.index_z
        ipp.x_C_my(y, p_star, k, z, :) .= value.(WMDER_IPP[:x_C][y, k, z])
    end
    for y in model_data.index_y, s in ipp.index_stor_existing, z in model_data.index_z
        ipp.x_stor_R_my(y, p_star, s, z, :) .= value.(WMDER_IPP[:x_stor_R][y, s, z])
    end
    for y in model_data.index_y, s in ipp.index_stor_new, z in model_data.index_z
        ipp.x_stor_C_my(y, p_star, s, z, :) .= value.(WMDER_IPP[:x_stor_C][y, s, z])
    end
    for y in model_data.index_y,
        p in ipp.index_p,
        k in ipp.index_k_existing,
        z in model_data.index_z,
        d in model_data.index_d,
        t in model_data.index_t

        ipp.y_E_my(y, p, k, z, d, t, :) .= value.(WMDER_IPP[:y_E][y, p, k, z, d, t])
    end
    for y in model_data.index_y,
        p in ipp.index_p,
        k in ipp.index_k_new,
        z in model_data.index_z,
        d in model_data.index_d,
        t in model_data.index_t

        ipp.y_C_my(y, p, k, z, d, t, :) .= value.(WMDER_IPP[:y_C][y, p, k, z, d, t])
    end
    for y in model_data.index_y, z in model_data.index_z, d in model_data.index_d, t in model_data.index_t
        ipp.miu_my(y, z, d, t, :) .= value.(WMDER_IPP[:miu][y, z, d, t])
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

        ipp.charge_E_my(y, p, s, z, d, t, :) .= value.(WMDER_IPP[:charge_E][y, p, s, z, d, t])
    end
    for y in model_data.index_y,
        p in ipp.index_p,
        s in ipp.index_stor_existing,
        z in model_data.index_z,
        d in model_data.index_d,
        t in model_data.index_t

        ipp.discharge_E_my(y, p, s, z, d, t, :) .= value.(WMDER_IPP[:discharge_E][y, p, s, z, d, t])
    end
    for y in model_data.index_y,
        p in ipp.index_p,
        s in ipp.index_stor_new,
        z in model_data.index_z,
        d in model_data.index_d,
        t in model_data.index_t

        ipp.charge_C_my(y, p, s, z, d, t, :) .= value.(WMDER_IPP[:charge_C][y, p, s, z, d, t])
    end
    for y in model_data.index_y,
        p in ipp.index_p,
        s in ipp.index_stor_new,
        z in model_data.index_z,
        d in model_data.index_d,
        t in model_data.index_t

        ipp.discharge_C_my(y, p, s, z, d, t, :) .= value.(WMDER_IPP[:discharge_C][y, p, s, z, d, t])
    end
    for y in model_data.index_y,
        p in ipp.index_p,
        s in ipp.index_stor_existing,
        z in model_data.index_z,
        d in model_data.index_d,
        t in model_data.index_t

        ipp.energy_E_my(y, p, s, z, d, t, :) .= value.(WMDER_IPP[:energy_E][y, p, s, z, d, t])
    end
    for y in model_data.index_y,
        p in ipp.index_p,
        s in ipp.index_stor_new,
        z in model_data.index_z,
        d in model_data.index_d,
        t in model_data.index_t

        ipp.energy_C_my(y, p, s, z, d, t, :) .= value.(WMDER_IPP[:energy_C][y, p, s, z, d, t])
    end
    for y in model_data.index_y,
        l in ipp.index_l,
        d in model_data.index_d,
        t in model_data.index_t

        ipp.flow_my(y, l, d, t, :) .= value.(WMDER_IPP[:flow][y, l, d, t])
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

    # Calculated prior to WMDER_IPP
    Max_Net_Load_my_dict = ipp.Max_Net_Load_my_dict

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

end

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
    # utility = get_agent(Utility, agent_store)
    regulator = get_agent(Regulator, agent_store)
    customers = get_agent(CustomerGroup, agent_store)
    der_aggregator = get_agent(DERAggregator, agent_store)
    green_developer = get_agent(GreenDeveloper, agent_store)

    reg_year, reg_year_index = get_reg_year(model_data)
    #### !!!! we also need to do delayed update for customers.total_pv_stor_capacity_my, but need to figure out all the if-else in customer_group !!!! ####
    # as long as the existing capacity for customers.total_pv_stor_capacity_my does not change from year to year for the simulation period this is fine (such as in this case).
    for z in model_data.index_z, h in model_data.index_h, m in (:BTMStorage, :BTMPV)
        update_total_capacity!(customers.total_der_capacity_my_delay_update, customers.x_DG_new, model_data, reg_year, z, h, m)
    end

    x_R_before = ParamArray(ipp.x_R_my)
    x_C_before = ParamArray(ipp.x_C_my)
    delta_t = get_delta_t(model_data)

    iteration_year = model_data.index_y_fix.elements[w_iter]

    max_iter = 10
    preferred_duality_gap = 0.05
    acceptable_duality_gap = 0.2
    bound_adjust = 0.5

    first_update = true
    bound_size = 0.3
    lower_level_duality_gap = 1.0
    skip_lower = false
    WMDER_IPP = nothing

    cust_year = customers.current_year
    dera_year = der_aggregator.current_year

    cust_pre_year = customers.previous_year
    dera_pre_year = der_aggregator.previous_year

    # TODO: Only update here on first year? Or find a better place to initialize?
    for z in model_data.index_z
        # simply assign DERAggregator to a random ipp (ipp1)
        ipp.x_stor_E_my(:ipp1, z, Symbol("der_aggregator"), :) .= sum(
            der_aggregator.dera_stor_my(dera_pre_year, z, h) for h in model_data.index_h
        )
        ipp.x_E_my(:ipp1, z, Symbol("dera_pv"), :) .= sum(
            der_aggregator.dera_pv_my(dera_pre_year, z, h) for h in model_data.index_h
        )
    end

    # TODO: Only update in customers
    for z in model_data.index_z, h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
        customers.rho_DG(h, :BTMStorage, z, d, t, :) .= customers.rho_DG_my(cust_pre_year, h, :BTMStorage, z, d, t)
    end

    for iter in 1:max_iter
        @info "IPP Problem Iteration $(iteration_year) - $(iter)"
        @info "" eta_param_vec_length=length(ipp.eta_param_vec) first_update bound_size lower_level_duality_gap skip_lower

        MPPDCMER_lower = nothing
        if !skip_lower
            lower_level_solver = []
            push!(lower_level_solver, ipp_opts.solvers["solve_agent_problem_ipp_mppdc_mccormic_lower"])
            solver = lower_level_solver[1]
            MPPDCMER_lower = ipp_cap_lower(
                ipp, ipp_opts, model_data, delta_t, window_length,
                customers, der_aggregator, green_developer, solver
            )
            push!(jump_model, MPPDCMER_lower)
            if termination_status(MPPDCMER_lower) != OPTIMAL
                lower_level_solver[1] = ipp_opts.solvers["solve_agent_problem_ipp_mppdc_mccormic_lower_presolve"]
                solver = lower_level_solver[1]
                MPPDCMER_lower = ipp_cap_lower(
                    ipp, ipp_opts, model_data, delta_t, window_length,
                    customers, der_aggregator, green_developer, solver
                )
                jump_model[end] = MPPDCMER_lower
                # if termination_status(MPPDCMER_lower) != OPTIMAL
                #     error("lower-level LP failed")
                # end
            end

            # TEMPORARY CODE FOR TESTING
            # objective_value(MPPDCMER_lower) # calling this errors the program if the optimization failed
            # dual_model = dualize(MPPDCMER_lower; dual_names = DualNames("dual", ""))
            # f = open("lower_level_dual.txt","w"); print(f, dual_model); close(f)

            ############ lower-level dual problem ############
            solver = lower_level_solver[1]
            MPPDCMER_lower_dual = ipp_cap_lower_dual(
                ipp, ipp_opts, model_data, delta_t, window_length,
                customers, der_aggregator, green_developer, solver
            )
            # jump_model[end] = MPPDCMER_lower_dual
        end

        # TEMPORARY CODE FOR TESTING
        # objective_value(MPPDCMER_lower_dual) # calling this errors the program if the optimization failed
        # f = open("lower_level_dual_my_version.txt","w"); print(f, MPPDCMER_lower_dual); close(f)
        # abs.(value.(miu_lower).data) .- abs.(dual.(Eq_primal_feasible_supplydemandbalance_lower).data)

        for z in model_data.index_z
            # simply assign DERAggregator to a random ipp (ipp1)
            ipp.x_stor_E_my(:ipp1, z, Symbol("der_aggregator"), :) .= sum(der_aggregator.dera_stor_my(dera_year, z, h) for h in model_data.index_h)
            ipp.x_E_my(:ipp1, z, Symbol("dera_pv"), :) .= sum(der_aggregator.dera_pv_my(dera_year, z, h) for h in model_data.index_h)
        end
    
        for z in model_data.index_z, h in model_data.index_h, d in model_data.index_d, t in model_data.index_t
            customers.rho_DG(h, :BTMStorage, z, d, t, :) .= customers.rho_DG_my(cust_year, h, :BTMStorage, z, d, t)
        end

        mcbnds, first_update = get_mccormick_bounds(
            MPPDCMER_lower, ipp, p_star, model_data, w_iter; 
            bound_size=bound_size, first_update=first_update, skip_lower=skip_lower
        )

        WMDER_IPP = ipp_cap_upper(
            mcbnds, ipp, ipp_opts, p_star, model_data, delta_t,
            regulator, customers, der_aggregator, green_developer
        )
        push!(jump_model, WMDER_IPP)

        if (termination_status(WMDER_IPP) != OPTIMAL) && (termination_status(WMDER_IPP) != LOCALLY_SOLVED)
            bound_size = min(bound_size * (1.0 + bound_adjust), 0.9)
            @info("Upper-level problem terminated as $(termination_status(WMDER_IPP)). Loosening McCormick bounds to $(bound_size).")
            skip_lower = true
            continue
        end

        objective_value(WMDER_IPP)
        lower_level_duality_gap = ipp_calc_duality_gap(
            WMDER_IPP, ipp, p_star, model_data, delta_t,
            customers, der_aggregator, green_developer
        )

        if lower_level_duality_gap < ipp.lower_level_duality_gap(reg_year_index)
            ipp.lower_level_duality_gap(reg_year_index, :) .= lower_level_duality_gap

            ipp_cap_save_results(
                WMDER_IPP, ipp, p_star, model_data, delta_t,
                customers, der_aggregator, green_developer
            )
        end

        if ipp.lower_level_duality_gap(reg_year_index) > preferred_duality_gap
            # tighten the bound
            bound_size = bound_size * bound_adjust
            @info("Duality gap $(ipp.lower_level_duality_gap(reg_year_index)) is greater than $(preferred_duality_gap). "*
                  "Re-running lower level problem and applying McCormick bounds with $(bound_size)")
            skip_lower = false
        else
            break
        end
    end

    if (ipp.lower_level_duality_gap(reg_year_index) > acceptable_duality_gap)
        msg = "MPPDC upper level never solved to optimality with an acceptable optimality gap, i.e., $(ipp.lower_level_duality_gap(reg_year_index)) > $(acceptable_duality_gap))."
        if (termination_status(WMDER_IPP) != OPTIMAL)
            @error "$(msg) In addition, MPPDC upper level did not solve to optimality and terminated as $(termination_status(WMDER_IPP)). Calling compute_conflict! and outputting information to iis_model_$(iteration_year).txt."
            compute_conflict!(WMDER_IPP)
            iis_model, _ = copy_conflict(WMDER_IPP)
            print(iis_model)
            f = open("iis_model_$(iteration_year).txt","w"); print(f, iis_model); close(f)
            objective_value(WMDER_IPP)
        else
            error(msg)
        end
    end

    @info "IPP Problem solved with duality gap $(ipp.lower_level_duality_gap(reg_year_index)) < $(preferred_duality_gap)."

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
                model_data.omega(t) * ((ipp.v_E_my(y, p, k, t) - ipp.PTC_existing(k)) * ipp.y_E_my(y, p, k, t)) for
                t in model_data.index_t, k in ipp.index_k_existing
            ) + sum(
                model_data.omega(t) * ((ipp.v_C_my(y, p, k, t) - ipp.PTC_new_my(y, k)) * ipp.y_C_my(y, p, k, t)) for
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
