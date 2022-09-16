# This file defines the data and functions associated with the utility.

abstract type AbstractUtility <: Agent end

mutable struct Utility <: AbstractUtility
    id::String
    # Sets
    "existing bulk generation technologies"
    index_k_existing::Dimension
    "potential bulk generation technologies"
    index_k_new::Dimension
    "RPS-qualified technologies"
    index_rps::Dimension

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
    "net export (MWh)"
    eximport::ParamAxisArray
    "peak export"
    Peak_eximport::ParamScalar

    # Primal Variables
    y_E::ParamAxisArray
    y_C::ParamAxisArray
    x_R::ParamAxisArray
    x_C::ParamAxisArray
    # Dual Variables
    miu::ParamAxisArray

    # Finance Related Parameters
    CapEx_existing::ParamAxisArray # capital cost of existing capacity ($/MW)
    CumuTaxDepre_existing::ParamAxisArray # cumulative tax depreciation of existing capacity (%)
    CumuAccoutDepre_existing::ParamAxisArray # cumulative accounting depreciation of existing capacity (%)
    AnnualAccoutDepre_existing::ParamAxisArray # annual accounting depreciation of existing capacity (%)
    ITC_existing::ParamAxisArray # ITC of existing capacity (%)
    CumuITCAmort_existing::ParamAxisArray # ITC amortization of existing capacity (%)
    ADIT_existing::ParamAxisArray # accumulated deferred income taxes ($/MW)
    RateBaseNoWC_existing::ParamAxisArray # rate base (excluding working capital) ($/MW)

    CapEx_new::ParamAxisArray # capital cost of new capacity ($/MW)
    CumuTaxDepre_new::ParamAxisArray # cumulative tax depreciation of new capacity (%)
    CumuAccoutDepre_new::ParamAxisArray # cumulative accounting depreciation of new capacity (%)
    ITC_new::ParamAxisArray # ITC of new capacity (%)
    CumuITCAmort_new::ParamAxisArray # ITC amortization of new capacity (%)
    ADIT_new::ParamAxisArray # accumulated deferred income taxes ($/MW)
    FOM_new::ParamAxisArray # fixed O&M of new capacity ($/MW-yr)
    Lifetime_new::ParamAxisArray # lifetime of new capacity (yrs)
    RateBaseNoWC_new::ParamAxisArray # rate base (excluding working capital) ($/MW)

    DebtRatio::ParamScalar
    COD::ParamScalar
    COE::ParamScalar
    Tax::ParamScalar
    DaysofWC::ParamScalar

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
    CRF_default::ParamScalar
    Peak_eximport_my::ParamAxisArray

    # Primal Variables (multi-year)
    y_E_my::ParamAxisArray
    y_C_my::ParamAxisArray
    x_R_my::ParamAxisArray
    x_C_my::ParamAxisArray
    # Dual Variables (multi-year)
    miu_my::ParamAxisArray

    # Finance Related Parameters (multi-year)
    CapEx_existing_my::ParamAxisArray # capital cost of existing capacity ($/MW)
    CumuTaxDepre_existing_my::ParamAxisArray # cumulative tax depreciation of existing capacity (%)
    CumuAccoutDepre_existing_my::ParamAxisArray # cumulative accounting depreciation of existing capacity (%)
    AnnualAccoutDepre_existing_my::ParamAxisArray # annual accounting depreciation of existing capacity (%)
    AnnualTaxDepre_existing_my::ParamAxisArray # annual tax depreciation of existing capacity (%)
    ITC_existing_my::ParamAxisArray # ITC of existing capacity (%)
    CumuITCAmort_existing_my::ParamAxisArray # ITC amortization of existing capacity (%)
    AnnualITCAmort_existing_my::ParamAxisArray # ITC amortization of existing capacity (%)
    ADIT_existing_my::ParamAxisArray # accumulated deferred income taxes ($/MW)
    RateBaseNoWC_existing_my::ParamAxisArray # rate base (excluding working capital) ($/MW)

    CumuTaxDepre_new_my::ParamAxisArray # cumulative tax depreciation of new capacity (%)
    CumuAccoutDepre_new_my::ParamAxisArray # cumulative accounting depreciation of new capacity (%)
    ITC_new_my::ParamAxisArray # ITC of new capacity (%)
    CumuITCAmort_new_my::ParamAxisArray # ITC amortization of new capacity (%)
    AnnualITCAmort_new_my::ParamAxisArray # ITC amortization of new capacity (%)
    AnnualAccoutDepre_new_my::ParamAxisArray # annual accounting depreciation of new capacity (%)
    AnnualTaxDepre_new_my::ParamAxisArray # annual tax depreciation of new capacity (%) 

    x_R_cumu::ParamAxisArray
    x_C_cumu::ParamAxisArray

    # capacity
    capacity_credit_E_my::ParamAxisArray # capacity credit of existing resources
    capacity_credit_C_my::ParamAxisArray # capacity credit of new resources
    Net_Load_my::ParamAxisArray
    Max_Net_Load_my::ParamAxisArray
    Reserve_req_my::ParamAxisArray

    # RPS
    RPS::ParamAxisArray
    rec_my::ParamAxisArray

    # distribution loss
    loss_dist::ParamScalar

    # emission rate
    emission_rate_E_my::ParamAxisArray
    emission_rate_C_my::ParamAxisArray

    # Lagrange decomposition
    x_R_feasible::ParamAxisArray
    x_C_feasible::ParamAxisArray
    obj_feasible::ParamScalar
    obj_upper_bound::ParamScalar
    obj_lower_bound::ParamScalar
    L_R_my::ParamAxisArray
    L_C_my::ParamAxisArray
    x_R_my_decomp::ParamAxisArray
    x_C_my_decomp::ParamAxisArray
    obj_my::ParamAxisArray
    obj_my_feasible::ParamAxisArray
end

function Utility(
    input_filename::String,
    model_data::HEMData,
    regulator::Regulator;
    id = DEFAULT_ID,
)
    index_k_existing = read_set(
        input_filename,
        "index_k_existing",
        "index_k_existing",
        prose_name = "existing bulk generation technologies",
    )

    index_k_new = read_set(
        input_filename,
        "index_k_new",
        "index_k_new",
        prose_name = "potential bulk generation technologies",
    )

    index_rps = read_set(
        input_filename,
        "index_rps",
        "index_rps",
        prose_name = "RPS-qualified technologies",
    )

    eximport = read_param(
        "eximport",
        input_filename,
        "Export",
        model_data.index_t,
        description = "net export (MWh)",
    )

    peak_eximport =
        ParamScalar("Peak_eximport", findmax(eximport)[1], description = "peak export")

    FOMNew = read_param("FOM_new", input_filename, "FOMNew", index_k_new)
    CapExNew = read_param("CapEx_new", input_filename, "CapExNew", index_k_new)
    LifetimeNew = read_param("Lifetime_new", input_filename, "LifetimeNew", index_k_new)
    debt_ratio = ParamScalar("DebtRatio", 0.6, description = "debt ratio")
    cost_of_debt = ParamScalar("COD", 0.06, description = "cost of debt")
    cost_of_equity = regulator.z
    tax_rate = ParamScalar("Tax", 0.26, description = "tax rate")
    atwacc = debt_ratio * cost_of_debt * (1 - tax_rate) + (1 - debt_ratio) * cost_of_equity
    CRF = Dict(
        k => atwacc * (1 + atwacc)^LifetimeNew[k] / ((1 + atwacc)^LifetimeNew[k] - 1)
        for k in index_k_new
    )
    FixedCostNew = AxisArray(
        [FOMNew[k] + CapExNew[k] * CRF[k] for k in index_k_new],
        index_k_new.elements,
    )

    CapExOld = read_param("CapEx_existing", input_filename, "CapExOld", index_k_existing)
    CumuTaxDepreOld = read_param(
        "CumuTaxDepre_existing",
        input_filename,
        "CumuTaxDepreOld",
        index_k_existing,
    )
    CumuAccoutDepreOld = read_param(
        "CumuAccoutDepre_existing",
        input_filename,
        "CumuAccoutDepreOld",
        index_k_existing,
    )
    AnnualAccoutDepreOld = read_param(
        "AnnualAccoutDepre_existing",
        input_filename,
        "AnnualAccoutDepreOld",
        index_k_existing,
    )
    ITCOld = read_param("ITC_existing", input_filename, "ITCOld", index_k_existing)
    CumuITCAmortOld = read_param(
        "CumuITCAmort_existing",
        input_filename,
        "CumuITCAmortOld",
        index_k_existing,
    )
    ADITOld = AxisArray(
        [
            CapExOld[k] * (CumuTaxDepreOld[k] - CumuAccoutDepreOld[k]) * tax_rate +
            ITCOld[k] * CapExOld[k] * (1 - CumuITCAmortOld[k]) for k in index_k_existing
        ],
        index_k_existing.elements,
    )
    RateBaseNoWCOld = AxisArray(
        [CapExOld[k] * (1 - CumuAccoutDepreOld[k]) - ADITOld[k] for k in index_k_existing],
        index_k_existing.elements,
    )

    CumuTaxDepreNew =
        read_param("CumuTaxDepre_new", input_filename, "CumuTaxDepreNew", index_k_new)
    CumuAccoutDepreNew =
        read_param("CumuAccoutDepre_new", input_filename, "CumuAccoutDepreNew", index_k_new)
    ITCNew = read_param("ITC_new", input_filename, "ITCNew", index_k_new)
    CumuITCAmortNew =
        read_param("CumuITCAmort_new", input_filename, "CumuITCAmortNew", index_k_new)
    ADITNew = AxisArray(
        [
            CapExNew[k] * (CumuTaxDepreNew[k] - CumuAccoutDepreNew[k]) * tax_rate +
            ITCNew[k] * CapExNew[k] * (1 - CumuITCAmortNew[k]) for k in index_k_new
        ],
        index_k_new.elements,
    )
    RateBaseNoWCNew = AxisArray(
        [CapExNew[k] * (1 - CumuAccoutDepreNew[k]) - ADITNew[k] for k in index_k_new],
        index_k_new.elements,
    )

    eximport_my = read_param(
        "eximport_my",
        input_filename,
        "Exportmy",
        model_data.index_t,
        [model_data.index_y],
    )
    array_eximport_my = zeros(length(model_data.index_y), length(model_data.index_t))
    for y in model_data.index_y, t in model_data.index_t
        array_eximport_my[
            Int(model_data.year[y] - model_data.year[first(model_data.index_y)] + 1),
            Int(model_data.hour[t]),
        ] = eximport_my[y, t]
    end
    peak_eximport_my = AxisArray(
        [
            findmax(array_eximport_my, dims = 2)[1][Int(
                model_data.year[y] - model_data.year[first(model_data.index_y)] + 1,
            )] for y in model_data.index_y
        ],
        model_data.index_y.elements,
    )

    CRF_default = atwacc * (1 + atwacc)^20 / ((1 + atwacc)^20 - 1)
    pvf_cap = AxisArray(
        [
            1 / (1 + atwacc)^(model_data.year[y] - model_data.year_start) for
            y in model_data.index_y
        ],
        model_data.index_y.elements,
    )
    pvf_onm = AxisArray(
        [
            1 / (1 + atwacc)^(model_data.year[y] - model_data.year_start) for
            y in model_data.index_y
        ],
        model_data.index_y.elements,
    )

    # capital expense of existing units ($/MW)
    CapExOld_my =
        read_param("CapEx_existing_my", input_filename, "CapExOld", index_k_existing)
    # cumulative tax depreciation of existing units (%)
    CumuTaxDepreOld_my = read_param(
        "CumuTaxDepre_existing_my",
        input_filename,
        "CumuTaxDepreOldmy",
        index_k_existing,
        [model_data.index_y],
    )
    # cumulative accounting depreciation of existing units (%)
    CumuAccoutDepreOld_my = read_param(
        "CumuAccoutDepre_existing_my",
        input_filename,
        "CumuAccoutDepreOldmy",
        index_k_existing,
        [model_data.index_y],
    )
    # annual accounting depreciation of existing units (%)
    AnnualAccoutDepreOld_my = read_param(
        "AnnualAccoutDepre_existing_my",
        input_filename,
        "AnnualAccoutDepreOldmy",
        index_k_existing,
        [model_data.index_y],
    )
    # annual tax depreciation of existing units (%)
    AnnualTaxDepreOld_my = read_param(
        "AnnualTaxDepre_existing_my",
        input_filename,
        "AnnualTaxDepreOldmy",
        index_k_existing,
        [model_data.index_y],
    )
    # ITC of existing units (%)
    ITCOld_my = read_param("ITC_existing_my", input_filename, "ITCOld", index_k_existing)
    # cumulative ITC ammortization of existing units (%)
    CumuITCAmortOld_my = read_param(
        "CumuITCAmort_existing_my",
        input_filename,
        "CumuITCAmortOldmy",
        index_k_existing,
        [model_data.index_y],
    )
    AnnualITCAmortOld_my = read_param(
        "AnnualITCAmort_existing_my",
        input_filename,
        "AnnualITCAmortOldmy",
        index_k_existing,
        [model_data.index_y],
    )
    # accumulative deferred income tax of existing units ($/MW)
    ADITOld_my = make_axis_array(model_data.index_y, index_k_existing)
    for y in model_data.index_y, k in index_k_existing
        ADITOld_my[y, k] =
            CapExOld_my[k] *
            (CumuTaxDepreOld_my[y, k] - CumuAccoutDepreOld_my[y, k]) *
            tax_rate + ITCOld_my[k] * CapExOld_my[k] * (1 - CumuITCAmortOld_my[y, k])
    end
    # rate base (without working capital) ($/MW)
    RateBaseNoWCOld_my = make_axis_array(model_data.index_y, index_k_existing)
    for y in model_data.index_y, k in index_k_existing
        RateBaseNoWCOld_my[y, k] =
            CapExOld_my[k] * (1 - CumuAccoutDepreOld_my[y, k]) - ADITOld_my[y, k]
    end

    # cumulative tax depreciation of new units (for each schedule year) (%)
    CumuTaxDepreNew_my = read_param(
        "CumuTaxDepre_new_my",
        input_filename,
        "CumuTaxDepreNewmy",
        index_k_new,
        [model_data.index_s],
    )
    # cumulative accounting depreciation of new units (for each schedule year) (%)
    CumuAccoutDepreNew_my = read_param(
        "CumuAccoutDepre_new_my",
        input_filename,
        "CumuAccoutDepreNewmy",
        index_k_new,
        [model_data.index_s],
    )
    # ITC of new units (%)
    ITCNew_my = read_param(
        "ITC_new_my",
        input_filename,
        "ITCNewmy",
        index_k_new,
        [model_data.index_y],
    )
    # cumulative ITC ammortization of new units (for each schedule year) (%)
    CumuITCAmortNew_my = read_param(
        "CumuITCAmort_new_my",
        input_filename,
        "CumuITCAmortNewmy",
        index_k_new,
        [model_data.index_s],
    )
    AnnualITCAmortNew_my = read_param(
        "AnnualITCAmort_new_my",
        input_filename,
        "AnnualITCAmortNewmy",
        index_k_new,
        [model_data.index_s],
    )
    # annual accounting depreciation of new units (%)
    AnnualAccoutDepreNew_my = read_param(
        "AnnualAccoutDepre_new_my",
        input_filename,
        "AnnualAccoutDepreNewmy",
        index_k_new,
        [model_data.index_s],
    )
    # annual tax depreciation of new units (%)
    AnnualTaxDepreNew_my = read_param(
        "AnnualTaxDepre_new_my",
        input_filename,
        "AnnualTaxDepreNewmy",
        index_k_new,
        [model_data.index_s],
    )

    return Utility(
        id,
        index_k_existing,
        index_k_new,
        index_rps,
        read_param(
            "x_E",
            input_filename,
            "ExistingCapacity",
            index_k_existing,
            description = "existing capacity (MW)",
        ),
        read_param(
            "f_E",
            input_filename,
            "FixedCostOld",
            index_k_existing,
            description = "fixed cost of existing capacity (\$/MW-yr)",
        ),
        ParamAxisArray("f_C", (index_k_new,), FixedCostNew),
        read_param(
            "v_E",
            input_filename,
            "VariableCostOld",
            model_data.index_t,
            [index_k_existing],
            description = "variable cost of existing capacity (\$/MWh)",
        ),
        read_param(
            "v_C",
            input_filename,
            "VariableCostNew",
            model_data.index_t,
            [index_k_new],
            description = "variable cost of new capacity (\$/Mwh)",
        ),
        read_param(
            "rho_E",
            input_filename,
            "AvailabilityOld",
            model_data.index_t,
            [index_k_existing],
            description = "availability of existing capacity (fraction)",
        ),
        read_param(
            "rho_C",
            input_filename,
            "AvailabilityNew",
            model_data.index_t,
            [index_k_new],
            description = "availability of new capacity (fraction)",
        ),
        eximport,
        peak_eximport,
        initialize_param(
            "y_E",
            index_k_existing,
            model_data.index_t,
            description = "existing capacity generation in each time period",
        ),
        initialize_param(
            "y_C",
            index_k_new,
            model_data.index_t,
            description = "new capacity generation in each time period",
        ),
        initialize_param("x_R", index_k_existing),
        initialize_param("x_C", index_k_new),
        initialize_param("miu", model_data.index_t),
        CapExOld,
        CumuTaxDepreOld,
        CumuAccoutDepreOld,
        AnnualAccoutDepreOld,
        ITCOld,
        CumuITCAmortOld,
        ParamAxisArray("ADIT_existing", (index_k_existing,), ADITOld),
        ParamAxisArray("RateBaseNoWC_existing", (index_k_existing,), RateBaseNoWCOld),
        CapExNew,
        CumuTaxDepreNew,
        CumuAccoutDepreNew,
        ITCNew,
        CumuITCAmortNew,
        ParamAxisArray("ADIT_new", (index_k_existing,), ADITNew),
        FOMNew,
        LifetimeNew,
        ParamAxisArray("RateBaseNoWC_new", (index_k_existing,), RateBaseNoWCNew),
        debt_ratio,
        cost_of_debt,
        cost_of_equity,
        tax_rate,
        ParamScalar("DaysofWC", 45.0, description = "number of days of working capital"),
        read_param("x_E_my", input_filename, "ExistingCapacity", index_k_existing),
        read_param(
            "fom_E_my",
            input_filename,
            "FixedCostOldmy",
            index_k_existing,
            [model_data.index_y],
        ),
        read_param(
            "fom_C_my",
            input_filename,
            "FOMNewmy",
            index_k_new,
            [model_data.index_y],
        ),
        read_param(
            "CapEx_my",
            input_filename,
            "CapExNewmy",
            index_k_new,
            [model_data.index_y],
        ),
        read_param(
            "v_E_my",
            input_filename,
            "VariableCostOldmy",
            model_data.index_t,
            [model_data.index_y, index_k_existing],
        ),
        read_param(
            "v_C_my",
            input_filename,
            "VariableCostNewmy",
            model_data.index_t,
            [model_data.index_y, index_k_new],
        ),
        read_param(
            "rho_E_my",
            input_filename,
            "AvailabilityOld",
            model_data.index_t,
            [index_k_existing],
        ),
        read_param(
            "rho_C_my",
            input_filename,
            "AvailabilityNew",
            model_data.index_t,
            [index_k_new],
        ),
        eximport_my,
        ParamAxisArray("pvf_cap", (model_data.index_y,), pvf_cap),
        ParamAxisArray("pvf_onm", (model_data.index_y,), pvf_onm),
        ParamScalar("CRF_default", CRF_default, description = "capital recovery factor"),
        ParamAxisArray("Peak_eximport_my", (model_data.index_y,), peak_eximport_my),
        initialize_param(
            "y_E_my",
            model_data.index_y,
            index_k_existing,
            model_data.index_t,
        ),
        initialize_param("y_C_my", model_data.index_y, index_k_new, model_data.index_t),
        initialize_param("x_R_my", model_data.index_y, index_k_existing),
        initialize_param("x_C_my", model_data.index_y, index_k_new),
        initialize_param("miu_my", model_data.index_y, model_data.index_t),
        CapExOld_my,
        CumuTaxDepreOld_my,
        CumuAccoutDepreOld_my,
        AnnualAccoutDepreOld_my,
        AnnualTaxDepreOld_my,
        ITCOld_my,
        CumuITCAmortOld_my,
        AnnualITCAmortOld_my,
        ParamAxisArray(
            "ADIT_existing_my",
            Tuple(push!(copy([model_data.index_y]), index_k_existing)),
            ADITOld_my,
        ),
        ParamAxisArray(
            "RateBaseNoWC_existing_my",
            Tuple(push!(copy([model_data.index_y]), index_k_existing)),
            RateBaseNoWCOld_my,
        ),
        CumuTaxDepreNew_my,
        CumuAccoutDepreNew_my,
        ITCNew_my,
        CumuITCAmortNew_my,
        AnnualITCAmortNew_my,
        AnnualAccoutDepreNew_my,
        AnnualTaxDepreNew_my,
        initialize_param("x_R_cumu", index_k_existing),
        initialize_param("x_C_cumu", index_k_new),
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
        initialize_param("Net_Load_my", model_data.index_y, model_data.index_t),
        initialize_param("Max_Net_Load_my", model_data.index_y),
        initialize_param("Reserve_req_my", model_data.index_y),
        read_param("RPS", input_filename, "RPS", model_data.index_y),
        initialize_param("rec_my", model_data.index_y),
        ParamScalar("loss_dist", 0.053, description = "distribution system loss factor"),
        read_param(
            "emission_rate_E_my",
            input_filename,
            "EmissionRateOldmy",
            index_k_existing,
            [model_data.index_y],
        ),
        read_param(
            "emission_rate_C_my",
            input_filename,
            "EmissionRateNewmy",
            index_k_new,
            [model_data.index_y],
        ),
        initialize_param("x_R_feasible", model_data.index_y, index_k_existing),
        initialize_param("x_C_feasible", model_data.index_y, index_k_new),
        ParamScalar("obj_feasible", 1.0, description = "feasible objective value"),
        ParamScalar("obj_upper_bound", Inf, description = "upper bound of objective value"),
        ParamScalar(
            "obj_lower_bound",
            -Inf,
            description = "lower bound of objective value",
        ),
        initialize_param(
            "L_R_my",
            model_data.index_y,
            model_data.index_y,
            index_k_existing,
        ),
        initialize_param("L_C_my", model_data.index_y, model_data.index_y, index_k_new),
        initialize_param(
            "x_R_my_decomp",
            model_data.index_y,
            model_data.index_y,
            index_k_existing,
        ),
        initialize_param(
            "x_C_my_decomp",
            model_data.index_y,
            model_data.index_y,
            index_k_new,
        ),
        initialize_param("obj_my", model_data.index_y),
        initialize_param("obj_my_feasible", model_data.index_y),
    )
end

get_id(x::Utility) = x.id

function solve_agent_problem!(
    utility::Utility,
    utility_opts::AgentOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{WholesaleMarket},
    agent_store::AgentStore,
    w_iter,
)
    return 0.0
end

function solve_agent_problem!(
    utility::Utility,
    utility_opts::AgentOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{VerticallyIntegratedUtility},
    agent_store::AgentStore,
    w_iter,
)
    regulator = get_agent(Regulator, agent_store)
    customers = get_agent(CustomerGroup, agent_store)
    green_developer = get_agent(GreenDeveloper, agent_store)

    VIUDER_Utility = get_new_jump_model(hem_opts.MIP_solver)

    # Define positive variables
    @variable(VIUDER_Utility, x_C[model_data.index_y, utility.index_k_new] >= 0)
    @variable(VIUDER_Utility, x_R[model_data.index_y, utility.index_k_existing] >= 0)
    @variable(
        VIUDER_Utility,
        y_E[model_data.index_y, utility.index_k_existing, model_data.index_t] >= 0
    )
    @variable(
        VIUDER_Utility,
        y_C[model_data.index_y, utility.index_k_new, model_data.index_t] >= 0
    )

    for y in model_data.index_y
        if y == last(model_data.index_y.elements)
            utility.pvf_onm[y] = utility.pvf_cap[y] / utility.CRF_default
        else
            utility.pvf_onm[y] = utility.pvf_cap[y]
        end
    end

    fill!(utility.Net_Load_my, NaN)
    for y in model_data.index_y, t in model_data.index_t
        utility.Net_Load_my[y, t] =
            sum(customers.gamma[h] * customers.d_my[y, h, t] for h in model_data.index_h) +
            utility.eximport_my[y, t] - sum(
                customers.rho_DG[h, m, t] * customers.x_DG_E_my[y, h, m] for
                h in model_data.index_h, m in customers.index_m
            ) - sum(
                customers.rho_DG[h, m, t] * sum(
                    customers.x_DG_new_my[Symbol(Int(y_symbol)), h, m] for y_symbol in
                    model_data.year[first(model_data.index_y_fix)]:model_data.year[y]
                ) for h in model_data.index_h, m in customers.index_m
            )
    end

    fill!(utility.Max_Net_Load_my, NaN)
    for y in model_data.index_y
        utility.Max_Net_Load_my[y] =
            findmax(Dict(t => utility.Net_Load_my[y, t] for t in model_data.index_t))[1]
    end

    Max_Net_Load_my_index = Dict(
        y =>
            findmax(Dict(t => utility.Net_Load_my[y, t] for t in model_data.index_t))[2]
        for y in model_data.index_y
    )

    fill!(utility.capacity_credit_E_my, NaN)
    for y in model_data.index_y, k in utility.index_k_existing
        utility.capacity_credit_E_my[y, k] = utility.rho_E_my[k, Max_Net_Load_my_index[y]]
    end
    fill!(utility.capacity_credit_C_my, NaN)
    for y in model_data.index_y, k in utility.index_k_new
        utility.capacity_credit_C_my[y, k] = utility.rho_C_my[k, Max_Net_Load_my_index[y]]
    end

    fill!(utility.Reserve_req_my, NaN)
    for y in model_data.index_y
        utility.Reserve_req_my[y] = (1 + regulator.r) * utility.Max_Net_Load_my[y]
    end

    objective_function = begin
        sum(
            # generation costs
            #   num hrs * ((fuel + vom) * gen existing + (fuel + vom) * gen new) for t and gen type
            utility.pvf_onm[y] * (
                sum(
                    model_data.omega[t] * (utility.v_E_my[y, k, t] * y_E[y, k, t]) for
                    t in model_data.index_t, k in utility.index_k_existing
                ) + sum(
                    model_data.omega[t] * (utility.v_C_my[y, k, t] * y_C[y, k, t]) for
                    t in model_data.index_t, k in utility.index_k_new
                )
            ) +
            # fixed o&m costs
            #   fom * (cap exist - cap retiring) + fom * cap new for gen type
            # the discount factor is applied to fom of existing capacity remaining at year y
            utility.pvf_onm[y] * sum(
                utility.fom_E_my[y, k] * (
                    utility.x_E_my[k] - sum(
                        x_R[Symbol(Int(y_symbol)), k] for y_symbol in
                        model_data.year[first(model_data.index_y)]:model_data.year[y]
                    )
                ) for k in utility.index_k_existing
            ) +
            # the discount factor is applied to fom of new capacity for every year since year y (when it is built)
            sum(
                utility.fom_C_my[y, k] *
                x_C[y, k] *
                sum(
                    utility.pvf_onm[Symbol(Int(y_symbol))] for y_symbol in
                    model_data.year[y]:model_data.year[last(model_data.index_y.elements)]
                ) for k in utility.index_k_new
            ) +
            # capital costs
            # the discout factor is applied to new capacity for the year it is built
            utility.pvf_cap[y] *
            sum(utility.CapEx_my[y, k] * x_C[y, k] for k in utility.index_k_new)

            for y in model_data.index_y
        )
    end

    @objective(VIUDER_Utility, Min, objective_function)

    supply_demand_balance =
        (y, t) -> begin
            # bulk generation at time t
            sum(y_E[y, k, t] for k in utility.index_k_existing) +
            sum(y_C[y, k, t] for k in utility.index_k_new) -
            # demand at time t
            sum(customers.gamma[h] * customers.d_my[y, h, t] for h in model_data.index_h) - utility.eximport_my[y, t] +
            # existing DG generation at time t
            sum(
                customers.rho_DG[h, m, t] * customers.x_DG_E_my[y, h, m] for
                h in model_data.index_h, m in customers.index_m
            ) +
            # new DG generation at time t
            sum(
                customers.rho_DG[h, m, t] * sum(
                    customers.x_DG_new_my[Symbol(Int(y_symbol)), h, m] for y_symbol in
                    model_data.year[first(model_data.index_y_fix)]:model_data.year[y]
                ) for h in model_data.index_h, m in customers.index_m
            ) +
            # green technology subscription at time t
            sum(
                utility.rho_C_my[j, t] * sum(green_developer.green_tech_buildout_my[Symbol(Int(y_symbol)), j, h] for y_symbol in
                model_data.year[first(model_data.index_y_fix)]:model_data.year[y])
                for j in model_data.index_j, h in model_data.index_h
            )
        end

    @constraint(
        VIUDER_Utility,
        Eq_miu[y in model_data.index_y, t in model_data.index_t],
        supply_demand_balance(y, t) >= 0
    )

    # HERE -- once running try defining function over two indices
    # y_E must be less than available capacity
    @constraint(
        VIUDER_Utility,
        Eq_eta[
            y in model_data.index_y,
            k in utility.index_k_existing,
            t in model_data.index_t,
        ],
        utility.rho_E_my[k, t] * (
            utility.x_E_my[k] - sum(
                x_R[Symbol(Int(y_symbol)), k] for
                y_symbol in model_data.year[first(model_data.index_y)]:model_data.year[y]
            ) - utility.x_R_cumu[k]
        ) - y_E[y, k, t] >= 0
    )
    # y_C must be less than available capacity
    @constraint(
        VIUDER_Utility,
        Eq_lambda[
            y in model_data.index_y,
            k in utility.index_k_new,
            t in model_data.index_t,
        ],
        utility.rho_C_my[k, t] * (
            sum(
                x_C[Symbol(Int(y_symbol)), k] for
                y_symbol in model_data.year[first(model_data.index_y)]:model_data.year[y]
            ) + utility.x_C_cumu[k]
        ) - y_C[y, k, t] >= 0
    )
    # retiring capacity is bounded by existing capacity
    @constraint(
        VIUDER_Utility,
        Eq_sigma[y in model_data.index_y, k in utility.index_k_existing],
        utility.x_E[k] - sum(
            x_R[Symbol(Int(y_symbol)), k] for
            y_symbol in model_data.year[first(model_data.index_y)]:model_data.year[y]
        ) - utility.x_R_cumu[k] >= 0
    )

    planning_reserves =
        (y, t) -> begin
            # bulk generation available capacity at time t
            sum(
                utility.rho_E_my[k, t] * (
                    utility.x_E_my[k] - sum(
                        x_R[Symbol(Int(y_symbol)), k] for y_symbol in
                        model_data.year[first(model_data.index_y)]:model_data.year[y]
                    ) - utility.x_R_cumu[k]
                ) for k in utility.index_k_existing
            ) + sum(
                utility.rho_C_my[k, t] * (
                    sum(
                        x_C[Symbol(Int(y_symbol)), k] for y_symbol in
                        model_data.year[first(model_data.index_y)]:model_data.year[y]
                    ) + utility.x_C_cumu[k]
                ) for k in utility.index_k_new
            ) +
            # green technology subscription
            sum(
                utility.rho_C_my[j, t] * sum(green_developer.green_tech_buildout_my[Symbol(Int(y_symbol)), j, h] for y_symbol in
                model_data.year[first(model_data.index_y_fix)]:model_data.year[y])
                for j in model_data.index_j, h in model_data.index_h
            ) -
            # net_load plus planning reserve
            (1 + regulator.r) * (
                sum(
                    customers.gamma[h] * customers.d_my[y, h, t] for
                    h in model_data.index_h
                ) + utility.eximport_my[y, t] - sum(
                    customers.rho_DG[h, m, t] * customers.x_DG_E_my[y, h, m] for
                    h in model_data.index_h, m in customers.index_m
                ) - sum(
                    customers.rho_DG[h, m, t] * sum(
                        customers.x_DG_new_my[Symbol(Int(y_symbol)), h, m] for
                        y_symbol in
                        model_data.year[first(model_data.index_y_fix)]:model_data.year[y]
                    ) for h in model_data.index_h, m in customers.index_m
                )
            )
        end
    @constraint(
        VIUDER_Utility,
        Eq_xi[y in model_data.index_y, t in model_data.index_t],
        planning_reserves(y, t) >= 0
    )

    planning_reserves_cap =
        y -> begin
            # bulk generation available capacity at time t
            sum(
                utility.capacity_credit_E_my[y, k] * (
                    utility.x_E_my[k] - sum(
                        x_R[Symbol(Int(y_symbol)), k] for y_symbol in
                        model_data.year[first(model_data.index_y)]:model_data.year[y]
                    ) - utility.x_R_cumu[k]
                ) for k in utility.index_k_existing
            ) + sum(
                utility.capacity_credit_C_my[y, k] * (
                    sum(
                        x_C[Symbol(Int(y_symbol)), k] for y_symbol in
                        model_data.year[first(model_data.index_y)]:model_data.year[y]
                    ) + utility.x_C_cumu[k]
                ) for k in utility.index_k_new
            ) +
            # green technology subscription
            sum(
                utility.capacity_credit_C_my[y, j] * sum(green_developer.green_tech_buildout_my[Symbol(Int(y_symbol)), j, h] for y_symbol in
                model_data.year[first(model_data.index_y_fix)]:model_data.year[y])
                for j in model_data.index_j, h in model_data.index_h
            ) -
            # net_load plus planning reserve
            utility.Reserve_req_my[y]
        end
    @constraint(
        VIUDER_Utility,
        Eq_xi_cap[y in model_data.index_y],
        planning_reserves_cap(y) >= 0
    )

    # RPS constraint
    @constraint(
        VIUDER_Utility,
        Eq_rps[y in model_data.index_y],
        sum(
            model_data.omega[t] * (y_E[y, rps, t] + y_C[y, rps, t]) for
            rps in utility.index_rps, t in model_data.index_t
        ) -
        utility.RPS[y] *
        sum(model_data.omega[t] * utility.Net_Load_my[y, t] for t in model_data.index_t) >=
        0
    )

    TimerOutputs.@timeit HEM_TIMER "optimize! VIUDER_Utility 1" begin
        optimize!(VIUDER_Utility)
    end

    # record current primary variable values
    for y in model_data.index_y, k in utility.index_k_existing, t in model_data.index_t
        utility.y_E_my[y, k, t] = value.(y_E[y, k, t])
    end

    for y in model_data.index_y, k in utility.index_k_new, t in model_data.index_t
        utility.y_C_my[y, k, t] = value.(y_C[y, k, t])
    end

    x_R_before = ParamAxisArray(utility.x_R_my)
    x_C_before = ParamAxisArray(utility.x_C_my)
    for y in model_data.index_y, k in utility.index_k_existing
        utility.x_R_my[y, k] = value.(x_R[y, k])
    end

    for y in model_data.index_y, k in utility.index_k_new
        utility.x_C_my[y, k] = value.(x_C[y, k])
    end

    for y in model_data.index_y, t in model_data.index_t
        utility.miu_my[y, t] = dual.(Eq_miu[y, t])
    end

    for y in model_data.index_y
        utility.rec_my[y] = dual.(Eq_rps[y]) ./ utility.pvf_onm[y]  # exact REC needs to be carefully evaluated
    end

    # @info "Original built capacity" x_C_before
    # @info "New built capacity" utility.x_C_my

    # report change in key variables from previous iteration to this one
    return compute_difference_percentage_one_norm([
        (x_R_before, utility.x_R_my),
        (x_C_before, utility.x_C_my),
    ])
end

function save_results(
    utility::Utility,
    utility_opts::AgentOptions,
    hem_opts::HEMOptions{VerticallyIntegratedUtility},
    export_file_path::AbstractString,
    fileprefix::AbstractString,
)
    # Primal Variables
    save_param(
        utility.y_E_my.values,
        [:Year, :GenTech, :Time],
        :Generation_MWh,
        joinpath(export_file_path, "$(fileprefix)_y_E.csv"),
    )
    save_param(
        utility.y_C_my.values,
        [:Year, :GenTech, :Time],
        :Generation_MWh,
        joinpath(export_file_path, "$(fileprefix)_y_C.csv"),
    )
    save_param(
        utility.x_R_my.values,
        [:Year, :GenTech],
        :Capacity_MW,
        joinpath(export_file_path, "$(fileprefix)_x_R.csv"),
    )
    save_param(
        utility.x_C_my.values,
        [:Year, :GenTech],
        :Capacity_MW,
        joinpath(export_file_path, "$(fileprefix)_x_C.csv"),
    )
end

function welfare_calculation!(
    utility::Utility,
    utility_opts::AgentOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{VerticallyIntegratedUtility},
    agent_store::AgentStore,
)
    regulator = get_agent(Regulator, agent_store)

    Utility_Revenue_my = make_axis_array(model_data.index_y_fix)
    Utility_Cost_my = make_axis_array(model_data.index_y_fix)
    Utility_debt_interest_my = make_axis_array(model_data.index_y_fix)
    Utility_income_tax_my = make_axis_array(model_data.index_y_fix)
    Utility_operational_cost_my = make_axis_array(model_data.index_y_fix)
    Utility_depreciation_my = make_axis_array(model_data.index_y_fix)
    Utility_depreciation_tax_my = make_axis_array(model_data.index_y_fix)
    Utility_total_emission_my = make_axis_array(model_data.index_y_fix)

    for y in model_data.index_y_fix
        Utility_Revenue_my[y] = regulator.revenue_req_my[y] + regulator.othercost[y]
        Utility_Cost_my[y] = regulator.cost_my[y] + regulator.othercost[y]
        Utility_debt_interest_my[y] = regulator.debt_interest_my[y]
        Utility_income_tax_my[y] = regulator.income_tax_my[y]
        Utility_operational_cost_my[y] = regulator.operational_cost_my[y]
        Utility_depreciation_my[y] = regulator.depreciation_my[y]
        Utility_depreciation_tax_my[y] = regulator.depreciation_tax_my[y]
        Utility_total_emission_my[y] =
            sum(
                model_data.omega[t] * (
                    sum(
                        utility.y_E_my[y, k, t] * utility.emission_rate_E_my[y, k] for
                        k in utility.index_k_existing
                    ) + sum(
                        utility.y_C_my[y, k, t] * utility.emission_rate_C_my[y, k] for
                        k in utility.index_k_new
                    )
                ) for t in model_data.index_t
            ) * 0.000453592
    end

    return Utility_Revenue_my,
    Utility_Cost_my,
    Utility_debt_interest_my,
    Utility_income_tax_my,
    Utility_operational_cost_my,
    Utility_depreciation_my,
    Utility_depreciation_tax_my,
    Utility_total_emission_my
end

function solve_agent_problem_decomposition_by_year(
    utility::Utility,
    utility_opts::AgentOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{VerticallyIntegratedUtility},
    agent_store::AgentStore,
    w_iter,
    year,
)
    regulator = get_agent(Regulator, agent_store)
    customers = get_agent(CustomerGroup, agent_store)

    VIUDER_Utility = get_new_jump_model(hem_opts.MIP_solver)
    set_optimizer_attribute(VIUDER_Utility, "OUTPUTLOG", 0)

    Number_of_years =
        Int(model_data.year[year] - model_data.year[first(model_data.index_y)] + 1)
    # Define positive variables
    @variable(
        VIUDER_Utility,
        0 <=
        x_C[model_data.index_y.elements[1:Number_of_years], utility.index_k_new] <=
        5000
    )
    @variable(
        VIUDER_Utility,
        x_R[model_data.index_y.elements[1:Number_of_years], utility.index_k_existing] >= 0
    )
    @variable(VIUDER_Utility, y_E[utility.index_k_existing, model_data.index_t] >= 0)
    @variable(VIUDER_Utility, y_C[utility.index_k_new, model_data.index_t] >= 0)

    for y in model_data.index_y
        if y == last(model_data.index_y.elements)
            utility.pvf_onm[y] = utility.pvf_cap[y] / utility.CRF_default
        else
            utility.pvf_onm[y] = utility.pvf_cap[y]
        end
    end

    fill!(utility.Net_Load_my, NaN)
    for y in model_data.index_y, t in model_data.index_t
        utility.Net_Load_my[y, t] =
            sum(customers.gamma[h] * customers.d_my[y, h, t] for h in model_data.index_h) +
            utility.eximport_my[y, t] - sum(
                customers.rho_DG[h, m, t] * customers.x_DG_E_my[y, h, m] for
                h in model_data.index_h, m in customers.index_m
            ) - sum(
                customers.rho_DG[h, m, t] * sum(
                    customers.x_DG_new_my[Symbol(Int(y_symbol)), h, m] for y_symbol in
                    model_data.year[first(model_data.index_y_fix)]:model_data.year[y]
                ) for h in model_data.index_h, m in customers.index_m
            )
    end
    fill!(utility.Max_Net_Load_my, NaN)
    for y in model_data.index_y
        utility.Max_Net_Load_my[y] =
            findmax((t => utility.Net_Load_my[y, t] for t in model_data.index_t))[1]
    end

    Max_Net_Load_my_index = Dict(
        y =>
            findmax(Dict(t => utility.Net_Load_my[y, t] for t in model_data.index_t))[2]
        for y in model_data.index_y
    )
    fill!(utility.capacity_credit_E_my, NaN)
    for y in model_data.index_y, k in utility.index_k_existing
        utility.capacity_credit_E_my[y, k] = utility.rho_E_my[k, Max_Net_Load_my_index[y]]
    end
    fill!(utility.capacity_credit_C_my, NaN)
    for y in model_data.index_y, k in utility.index_k_new
        utility.capacity_credit_C_my[y, k] = utility.rho_C_my[k, Max_Net_Load_my_index[y]]
    end
    fill!(utility.Reserve_req_my, NaN)
    for y in model_data.index_y
        utility.Reserve_req_my[y] = (1 + regulator.r) * utility.Max_Net_Load_my[y]
    end

    y = year

    # define Lagrange terms
    if (model_data.year[year] == model_data.year[first(model_data.index_y)]) &&
       (model_data.year[year] == model_data.year[last(model_data.index_y.elements)])
        Lagrange_terms = begin
            0
        end
    elseif (model_data.year[year] == model_data.year[first(model_data.index_y)]) &&
           (model_data.year[year] < model_data.year[last(model_data.index_y.elements)])
        Lagrange_terms = begin
            sum(
                utility.L_R_my[year, Symbol(Int(y_decision)), k] * x_R[year, k] for
                y_decision in
                (model_data.year[first(model_data.index_y)] + 1):model_data.year[last(
                    model_data.index_y.elements,
                )], k in utility.index_k_existing
            ) + sum(
                utility.L_C_my[year, Symbol(Int(y_decision)), k] * x_C[year, k] for
                y_decision in
                (model_data.year[first(model_data.index_y)] + 1):model_data.year[last(
                    model_data.index_y.elements,
                )], k in utility.index_k_new
            )
        end
    elseif (model_data.year[year] == model_data.year[last(model_data.index_y.elements)]) &&
           (model_data.year[year] > model_data.year[first(model_data.index_y)])
        Lagrange_terms = begin
            -sum(
                utility.L_R_my[Symbol(Int(y_inv_ret)), year, k] *
                x_R[Symbol(Int(y_inv_ret)), k] for y_inv_ret in
                model_data.year[first(model_data.index_y)]:(model_data.year[year] - 1),
                k in utility.index_k_existing
            ) - sum(
                utility.L_C_my[Symbol(Int(y_inv_ret)), year, k] *
                x_C[Symbol(Int(y_inv_ret)), k] for y_inv_ret in
                model_data.year[first(model_data.index_y)]:(model_data.year[year] - 1),
                k in utility.index_k_new
            )
        end
    else
        Lagrange_terms = begin
            sum(
                utility.L_R_my[year, Symbol(Int(y_decision)), k] * x_R[year, k] for
                y_decision in
                (model_data.year[year] + 1):model_data.year[last(
                    model_data.index_y.elements,
                )], k in utility.index_k_existing
            ) + sum(
                utility.L_C_my[year, Symbol(Int(y_decision)), k] * x_C[year, k] for
                y_decision in
                (model_data.year[year] + 1):model_data.year[last(
                    model_data.index_y.elements,
                )], k in utility.index_k_new
            ) - sum(
                utility.L_R_my[Symbol(Int(y_inv_ret)), year, k] *
                x_R[Symbol(Int(y_inv_ret)), k] for y_inv_ret in
                model_data.year[first(model_data.index_y)]:(model_data.year[year] - 1),
                k in utility.index_k_existing
            ) - sum(
                utility.L_C_my[Symbol(Int(y_inv_ret)), year, k] *
                x_C[Symbol(Int(y_inv_ret)), k] for y_inv_ret in
                model_data.year[first(model_data.index_y)]:(model_data.year[year] - 1),
                k in utility.index_k_new
            )
        end
    end

    objective_function = begin
        # generation costs
        #   num hrs * ((fuel + vom) * gen existing + (fuel + vom) * gen new) for t and gen type
        utility.pvf_onm[y] * (
            sum(
                model_data.omega[t] * (utility.v_E_my[y, k, t] * y_E[k, t]) for
                t in model_data.index_t, k in utility.index_k_existing
            ) + sum(
                model_data.omega[t] * (utility.v_C_my[y, k, t] * y_C[k, t]) for
                t in model_data.index_t, k in utility.index_k_new
            )
        ) +
        # fixed o&m costs
        #   fom * (cap exist - cap retiring) + fom * cap new for gen type
        # the discount factor is applied to fom of existing capacity remaining at year y
        utility.pvf_onm[y] * sum(
            utility.fom_E_my[y, k] * (
                utility.x_E_my[k] - sum(
                    x_R[Symbol(Int(y_symbol)), k] for y_symbol in
                    model_data.year[first(model_data.index_y)]:model_data.year[y]
                )
            ) for k in utility.index_k_existing
        ) +
        # the discount factor is applied to fom of new capacity for every year since year y (when it is built)
        sum(
            utility.fom_C_my[y, k] *
            x_C[y, k] *
            sum(
                utility.pvf_onm[Symbol(Int(y_symbol))] for y_symbol in
                model_data.year[y]:model_data.year[last(model_data.index_y.elements)]
            ) for k in utility.index_k_new
        ) +
        # capital costs
        # the discout factor is applied to new capacity for the year it is built
        utility.pvf_cap[y] *
        sum(utility.CapEx_my[y, k] * x_C[y, k] for k in utility.index_k_new) +
        # lagrange multiplier
        Lagrange_terms
    end

    @objective(VIUDER_Utility, Min, objective_function)

    supply_demand_balance =
        t -> begin
            # bulk generation at time t
            sum(y_E[k, t] for k in utility.index_k_existing) +
            sum(y_C[k, t] for k in utility.index_k_new) -
            # demand at time t
            sum(customers.gamma[h] * customers.d_my[y, h, t] for h in model_data.index_h) - utility.eximport_my[y, t] +
            # existing DG generation at time t
            sum(
                customers.rho_DG[h, m, t] * customers.x_DG_E_my[y, h, m] for
                h in model_data.index_h, m in customers.index_m
            ) +
            # new DG generation at time t
            sum(
                customers.rho_DG[h, m, t] * sum(
                    customers.x_DG_new_my[Symbol(Int(y_symbol)), h, m] for y_symbol in
                    model_data.year[first(model_data.index_y_fix)]:model_data.year[y]
                ) for h in model_data.index_h, m in customers.index_m
            )
        end

    @constraint(
        VIUDER_Utility,
        Eq_miu[t in model_data.index_t],
        supply_demand_balance(t) >= 0
    )

    # HERE -- once running try defining function over two indices
    # y_E must be less than available capacity
    @constraint(
        VIUDER_Utility,
        Eq_eta[k in utility.index_k_existing, t in model_data.index_t],
        utility.rho_E_my[k, t] * (
            utility.x_E_my[k] - sum(
                x_R[Symbol(Int(y_symbol)), k] for
                y_symbol in model_data.year[first(model_data.index_y)]:model_data.year[y]
            ) - utility.x_R_cumu[k]
        ) - y_E[k, t] >= 0
    )
    # y_C must be less than available capacity
    @constraint(
        VIUDER_Utility,
        Eq_lambda[k in utility.index_k_new, t in model_data.index_t],
        utility.rho_C_my[k, t] * (
            sum(
                x_C[Symbol(Int(y_symbol)), k] for
                y_symbol in model_data.year[first(model_data.index_y)]:model_data.year[y]
            ) + utility.x_C_cumu[k]
        ) - y_C[k, t] >= 0
    )
    # retiring capacity is bounded by existing capacity
    @constraint(
        VIUDER_Utility,
        Eq_sigma[k in utility.index_k_existing],
        utility.x_E[k] - sum(
            x_R[Symbol(Int(y_symbol)), k] for
            y_symbol in model_data.year[first(model_data.index_y)]:model_data.year[y]
        ) - utility.x_R_cumu[k] >= 0
    )

    planning_reserves =
        t -> begin
            # bulk generation available capacity at time t
            sum(
                utility.rho_E_my[k, t] * (
                    utility.x_E_my[k] - sum(
                        x_R[Symbol(Int(y_symbol)), k] for y_symbol in
                        model_data.year[first(model_data.index_y)]:model_data.year[y]
                    ) - utility.x_R_cumu[k]
                ) for k in utility.index_k_existing
            ) + sum(
                utility.rho_C_my[k, t] * (
                    sum(
                        x_C[Symbol(Int(y_symbol)), k] for y_symbol in
                        model_data.year[first(model_data.index_y)]:model_data.year[y]
                    ) + utility.x_C_cumu[k]
                ) for k in utility.index_k_new
            ) -
            # net_load plus planning reserve
            (1 + regulator.r) * (
                sum(
                    customers.gamma[h] * customers.d_my[y, h, t] for
                    h in model_data.index_h
                ) + utility.eximport_my[y, t] - sum(
                    customers.rho_DG[h, m, t] * customers.x_DG_E_my[y, h, m] for
                    h in model_data.index_h, m in customers.index_m
                ) - sum(
                    customers.rho_DG[h, m, t] * sum(
                        customers.x_DG_new_my[Symbol(Int(y_symbol)), h, m] for
                        y_symbol in
                        model_data.year[first(model_data.index_y_fix)]:model_data.year[y]
                    ) for h in model_data.index_h, m in customers.index_m
                )
            )
        end
    @constraint(VIUDER_Utility, Eq_xi[t in model_data.index_t], planning_reserves(t) >= 0)

    planning_reserves_cap = begin
        # bulk generation available capacity at time t
        sum(
            utility.capacity_credit_E_my[y, k] * (
                utility.x_E_my[k] - sum(
                    x_R[Symbol(Int(y_symbol)), k] for y_symbol in
                    model_data.year[first(model_data.index_y)]:model_data.year[y]
                ) - utility.x_R_cumu[k]
            ) for k in utility.index_k_existing
        ) + sum(
            utility.capacity_credit_C_my[y, k] * (
                sum(
                    x_C[Symbol(Int(y_symbol)), k] for y_symbol in
                    model_data.year[first(model_data.index_y)]:model_data.year[y]
                ) + utility.x_C_cumu[k]
            ) for k in utility.index_k_new
        ) -
        # net_load plus planning reserve
        utility.Reserve_req_my[y]
    end
    @constraint(VIUDER_Utility, Eq_xi_cap, planning_reserves_cap >= 0)

    # RPS constraint
    @constraint(
        VIUDER_Utility,
        Eq_rps,
        sum(
            model_data.omega[t] * (y_E[rps, t] + y_C[rps, t]) for rps in utility.index_rps,
            t in model_data.index_t
        ) -
        utility.RPS[y] *
        sum(model_data.omega[t] * utility.Net_Load_my[y, t] for t in model_data.index_t) >=
        0
    )

    TimerOutputs.@timeit HEM_TIMER "optimize! VIUDER_Utility 2" begin
        optimize!(VIUDER_Utility)
    end

    for y_inv_ret in model_data.year[first(model_data.index_y)]:model_data.year[year],
        k in utility.index_k_existing

        utility.x_R_my_decomp[Symbol(Int(y_inv_ret)), year, k] =
            value.(x_R[Symbol(Int(y_inv_ret)), k])
    end

    for y_inv_ret in model_data.year[first(model_data.index_y)]:model_data.year[year],
        k in utility.index_k_new

        utility.x_C_my_decomp[Symbol(Int(y_inv_ret)), year, k] =
            value.(x_C[Symbol(Int(y_inv_ret)), k])
    end

    utility.obj_my[year] = objective_value(VIUDER_Utility)

    # return VIUDER_Utility
end

function solve_agent_problem_decomposition_by_year_feasible(
    utility::Utility,
    utility_opts::AgentOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{VerticallyIntegratedUtility},
    agent_store::AgentStore,
    w_iter,
    year,
)
    regulator = get_agent(Regulator, agent_store)
    customers = get_agent(CustomerGroup, agent_store)

    VIUDER_Utility = get_new_jump_model(hem_opts.MIP_solver)
    set_optimizer_attribute(VIUDER_Utility, "OUTPUTLOG", 0)

    Number_of_years =
        Int(model_data.year[year] - model_data.year[first(model_data.index_y)] + 1)
    # Define positive variables
    @variable(
        VIUDER_Utility,
        0 <=
        x_C[model_data.index_y.elements[1:Number_of_years], utility.index_k_new] <=
        5000
    )
    @variable(
        VIUDER_Utility,
        x_R[model_data.index_y.elements[1:Number_of_years], utility.index_k_existing] >= 0
    )
    @variable(VIUDER_Utility, y_E[utility.index_k_existing, model_data.index_t] >= 0)
    @variable(VIUDER_Utility, y_C[utility.index_k_new, model_data.index_t] >= 0)
    # fix previous year solution
    if model_data.year[year] > model_data.year[first(model_data.index_y)]
        for y in model_data.year[first(model_data.index_y)]:(model_data.year[year] - 1),
            k in utility.index_k_existing

            fix(
                x_R[Symbol(Int(y)), k],
                utility.x_R_feasible[Symbol(Int(y)), k];
                force = true,
            )
        end
        for y in model_data.year[first(model_data.index_y)]:(model_data.year[year] - 1),
            k in utility.index_k_new

            fix(
                x_C[Symbol(Int(y)), k],
                utility.x_C_feasible[Symbol(Int(y)), k];
                force = true,
            )
        end
    end

    for y in model_data.index_y
        if y == last(model_data.index_y.elements)
            utility.pvf_onm[y] = utility.pvf_cap[y] / utility.CRF_default
        else
            utility.pvf_onm[y] = utility.pvf_cap[y]
        end
    end

    fill!(utility.Net_Load_my, NaN)
    for y in model_data.index_y, t in model_data.index_t
        utility.Net_Load_my[y, t] =
            sum(customers.gamma[h] * customers.d_my[y, h, t] for h in model_data.index_h) +
            utility.eximport_my[y, t] - sum(
                customers.rho_DG[h, m, t] * customers.x_DG_E_my[y, h, m] for
                h in model_data.index_h, m in customers.index_m
            ) - sum(
                customers.rho_DG[h, m, t] * sum(
                    customers.x_DG_new_my[Symbol(Int(y_symbol)), h, m] for y_symbol in
                    model_data.year[first(model_data.index_y_fix)]:model_data.year[y]
                ) for h in model_data.index_h, m in customers.index_m
            )
    end
    fill!(utility.Max_Net_Load_my, NaN)
    for y in model_data.index_y
        utility.Max_Net_Load_my[y] =
            findmax((t => utility.Net_Load_my[y, t] for t in model_data.index_t))[1]
    end

    Max_Net_Load_my_index = Dict(
        y =>
            findmax(Dict(t => utility.Net_Load_my[y, t] for t in model_data.index_t))[2]
        for y in model_data.index_y
    )
    fill!(utility.capacity_credit_E_my, NaN)
    for y in model_data.index_y, k in utility.index_k_existing
        utility.capacity_credit_E_my[y, k] = utility.rho_E_my[k, Max_Net_Load_my_index[y]]
    end
    fill!(utility.capacity_credit_C_my, NaN)
    for y in model_data.index_y, k in utility.index_k_new
        utility.capacity_credit_C_my[y, k] = utility.rho_C_my[k, Max_Net_Load_my_index[y]]
    end
    fill!(utility.Reserve_req_my, NaN)
    for y in model_data.index_y
        utility.Reserve_req_my[y] = (1 + regulator.r) * utility.Max_Net_Load_my[y]
    end

    y = year

    # define Lagrange terms
    if (model_data.year[year] == model_data.year[first(model_data.index_y)]) &&
       (model_data.year[year] == model_data.year[last(model_data.index_y.elements)])
        Lagrange_terms = begin
            0
        end
    elseif (model_data.year[year] == model_data.year[first(model_data.index_y)]) &&
           (model_data.year[year] < model_data.year[last(model_data.index_y.elements)])
        Lagrange_terms = begin
            sum(
                utility.L_R_my[year, Symbol(Int(y_decision)), k] * x_R[year, k] for
                y_decision in
                (model_data.year[first(model_data.index_y)] + 1):model_data.year[last(
                    model_data.index_y.elements,
                )], k in utility.index_k_existing
            ) + sum(
                utility.L_C_my[year, Symbol(Int(y_decision)), k] * x_C[year, k] for
                y_decision in
                (model_data.year[first(model_data.index_y)] + 1):model_data.year[last(
                    model_data.index_y.elements,
                )], k in utility.index_k_new
            )
        end
    elseif (model_data.year[year] == model_data.year[last(model_data.index_y.elements)]) &&
           (model_data.year[year] > model_data.year[first(model_data.index_y)])
        Lagrange_terms = begin
            -sum(
                utility.L_R_my[Symbol(Int(y_inv_ret)), year, k] *
                x_R[Symbol(Int(y_inv_ret)), k] for y_inv_ret in
                model_data.year[first(model_data.index_y)]:(model_data.year[year] - 1),
                k in utility.index_k_existing
            ) - sum(
                utility.L_C_my[Symbol(Int(y_inv_ret)), year, k] *
                x_C[Symbol(Int(y_inv_ret)), k] for y_inv_ret in
                model_data.year[first(model_data.index_y)]:(model_data.year[year] - 1),
                k in utility.index_k_new
            )
        end
    else
        Lagrange_terms = begin
            sum(
                utility.L_R_my[year, Symbol(Int(y_decision)), k] * x_R[year, k] for
                y_decision in
                (model_data.year[year] + 1):model_data.year[last(
                    model_data.index_y.elements,
                )], k in utility.index_k_existing
            ) + sum(
                utility.L_C_my[year, Symbol(Int(y_decision)), k] * x_C[year, k] for
                y_decision in
                (model_data.year[year] + 1):model_data.year[last(
                    model_data.index_y.elements,
                )], k in utility.index_k_new
            ) - sum(
                utility.L_R_my[Symbol(Int(y_inv_ret)), year, k] *
                x_R[Symbol(Int(y_inv_ret)), k] for y_inv_ret in
                model_data.year[first(model_data.index_y)]:(model_data.year[year] - 1),
                k in utility.index_k_existing
            ) - sum(
                utility.L_C_my[Symbol(Int(y_inv_ret)), year, k] *
                x_C[Symbol(Int(y_inv_ret)), k] for y_inv_ret in
                model_data.year[first(model_data.index_y)]:(model_data.year[year] - 1),
                k in utility.index_k_new
            )
        end
    end

    objective_function = begin
        # generation costs
        #   num hrs * ((fuel + vom) * gen existing + (fuel + vom) * gen new) for t and gen type
        utility.pvf_onm[y] * (
            sum(
                model_data.omega[t] * (utility.v_E_my[y, k, t] * y_E[k, t]) for
                t in model_data.index_t, k in utility.index_k_existing
            ) + sum(
                model_data.omega[t] * (utility.v_C_my[y, k, t] * y_C[k, t]) for
                t in model_data.index_t, k in utility.index_k_new
            )
        ) +
        # fixed o&m costs
        #   fom * (cap exist - cap retiring) + fom * cap new for gen type
        # the discount factor is applied to fom of existing capacity remaining at year y
        utility.pvf_onm[y] * sum(
            utility.fom_E_my[y, k] * (
                utility.x_E_my[k] - sum(
                    x_R[Symbol(Int(y_symbol)), k] for y_symbol in
                    model_data.year[first(model_data.index_y)]:model_data.year[y]
                )
            ) for k in utility.index_k_existing
        ) +
        # the discount factor is applied to fom of new capacity for every year since year y (when it is built)
        sum(
            utility.fom_C_my[y, k] *
            x_C[y, k] *
            sum(
                utility.pvf_onm[Symbol(Int(y_symbol))] for y_symbol in
                model_data.year[y]:model_data.year[last(model_data.index_y.elements)]
            ) for k in utility.index_k_new
        ) +
        # capital costs
        # the discout factor is applied to new capacity for the year it is built
        utility.pvf_cap[y] *
        sum(utility.CapEx_my[y, k] * x_C[y, k] for k in utility.index_k_new) +
        # lagrange multiplier
        Lagrange_terms
    end

    @objective(VIUDER_Utility, Min, objective_function)

    supply_demand_balance =
        t -> begin
            # bulk generation at time t
            sum(y_E[k, t] for k in utility.index_k_existing) +
            sum(y_C[k, t] for k in utility.index_k_new) -
            # demand at time t
            sum(customers.gamma[h] * customers.d_my[y, h, t] for h in model_data.index_h) - utility.eximport_my[y, t] +
            # existing DG generation at time t
            sum(
                customers.rho_DG[h, m, t] * customers.x_DG_E_my[y, h, m] for
                h in model_data.index_h, m in customers.index_m
            ) +
            # new DG generation at time t
            sum(
                customers.rho_DG[h, m, t] * sum(
                    customers.x_DG_new_my[Symbol(Int(y_symbol)), h, m] for y_symbol in
                    model_data.year[first(model_data.index_y_fix)]:model_data.year[y]
                ) for h in model_data.index_h, m in customers.index_m
            )
        end

    @constraint(
        VIUDER_Utility,
        Eq_miu[t in model_data.index_t],
        supply_demand_balance(t) >= 0
    )

    # HERE -- once running try defining function over two indices
    # y_E must be less than available capacity
    @constraint(
        VIUDER_Utility,
        Eq_eta[k in utility.index_k_existing, t in model_data.index_t],
        utility.rho_E_my[k, t] * (
            utility.x_E_my[k] - sum(
                x_R[Symbol(Int(y_symbol)), k] for
                y_symbol in model_data.year[first(model_data.index_y)]:model_data.year[y]
            ) - utility.x_R_cumu[k]
        ) - y_E[k, t] >= 0
    )
    # y_C must be less than available capacity
    @constraint(
        VIUDER_Utility,
        Eq_lambda[k in utility.index_k_new, t in model_data.index_t],
        utility.rho_C_my[k, t] * (
            sum(
                x_C[Symbol(Int(y_symbol)), k] for
                y_symbol in model_data.year[first(model_data.index_y)]:model_data.year[y]
            ) + utility.x_C_cumu[k]
        ) - y_C[k, t] >= 0
    )
    # retiring capacity is bounded by existing capacity
    @constraint(
        VIUDER_Utility,
        Eq_sigma[k in utility.index_k_existing],
        utility.x_E[k] - sum(
            x_R[Symbol(Int(y_symbol)), k] for
            y_symbol in model_data.year[first(model_data.index_y)]:model_data.year[y]
        ) - utility.x_R_cumu[k] >= 0
    )

    planning_reserves =
        t -> begin
            # bulk generation available capacity at time t
            sum(
                utility.rho_E_my[k, t] * (
                    utility.x_E_my[k] - sum(
                        x_R[Symbol(Int(y_symbol)), k] for y_symbol in
                        model_data.year[first(model_data.index_y)]:model_data.year[y]
                    ) - utility.x_R_cumu[k]
                ) for k in utility.index_k_existing
            ) + sum(
                utility.rho_C_my[k, t] * (
                    sum(
                        x_C[Symbol(Int(y_symbol)), k] for y_symbol in
                        model_data.year[first(model_data.index_y)]:model_data.year[y]
                    ) + utility.x_C_cumu[k]
                ) for k in utility.index_k_new
            ) -
            # net_load plus planning reserve
            (1 + regulator.r) * (
                sum(
                    customers.gamma[h] * customers.d_my[y, h, t] for
                    h in model_data.index_h
                ) + utility.eximport_my[y, t] - sum(
                    customers.rho_DG[h, m, t] * customers.x_DG_E_my[y, h, m] for
                    h in model_data.index_h, m in customers.index_m
                ) - sum(
                    customers.rho_DG[h, m, t] * sum(
                        customers.x_DG_new_my[Symbol(Int(y_symbol)), h, m] for
                        y_symbol in
                        model_data.year[first(model_data.index_y_fix)]:model_data.year[y]
                    ) for h in model_data.index_h, m in customers.index_m
                )
            )
        end
    @constraint(VIUDER_Utility, Eq_xi[t in model_data.index_t], planning_reserves(t) >= 0)

    planning_reserves_cap = begin
        # bulk generation available capacity at time t
        sum(
            utility.capacity_credit_E_my[y, k] * (
                utility.x_E_my[k] - sum(
                    x_R[Symbol(Int(y_symbol)), k] for y_symbol in
                    model_data.year[first(model_data.index_y)]:model_data.year[y]
                ) - utility.x_R_cumu[k]
            ) for k in utility.index_k_existing
        ) + sum(
            utility.capacity_credit_C_my[y, k] * (
                sum(
                    x_C[Symbol(Int(y_symbol)), k] for y_symbol in
                    model_data.year[first(model_data.index_y)]:model_data.year[y]
                ) + utility.x_C_cumu[k]
            ) for k in utility.index_k_new
        ) -
        # net_load plus planning reserve
        utility.Reserve_req_my[y]
    end
    @constraint(VIUDER_Utility, Eq_xi_cap, planning_reserves_cap >= 0)

    # RPS constraint
    @constraint(
        VIUDER_Utility,
        Eq_rps,
        sum(
            model_data.omega[t] * (y_E[rps, t] + y_C[rps, t]) for rps in utility.index_rps,
            t in model_data.index_t
        ) -
        utility.RPS[y] *
        sum(model_data.omega[t] * utility.Net_Load_my[y, t] for t in model_data.index_t) >=
        0
    )

    TimerOutputs.@timeit HEM_TIMER "optimize! VIUDER_Utility 3" begin
        optimize!(VIUDER_Utility)
    end

    for k in utility.index_k_existing
        utility.x_R_feasible[year, k] = value.(x_R[year, k])
    end

    for k in utility.index_k_new
        utility.x_C_feasible[year, k] = value.(x_C[year, k])
    end

    utility.obj_my_feasible[year] = objective_value(VIUDER_Utility)

    # return VIUDER_Utility
end

function solve_agent_problem_decomposition_by_year_feasible_obj(
    utility::Utility,
    utility_opts::AgentOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{VerticallyIntegratedUtility},
    agent_store::AgentStore,
    w_iter,
)

    # solve the original problem with feasible investment and retirement fixed
    customers = get_agent(CustomerGroup, agent_store)

    VIUDER_Utility_feasible = get_new_jump_model(hem_opts.MIP_solver)
    set_optimizer_attribute(VIUDER_Utility_feasible, "OUTPUTLOG", 0)

    # Define positive variables
    @variable(
        VIUDER_Utility_feasible,
        y_E[model_data.index_y, utility.index_k_existing, model_data.index_t] >= 0
    )
    @variable(
        VIUDER_Utility_feasible,
        y_C[model_data.index_y, utility.index_k_new, model_data.index_t] >= 0
    )

    objective_function = begin
        sum(
            # generation costs
            #   num hrs * ((fuel + vom) * gen existing + (fuel + vom) * gen new) for t and gen type
            utility.pvf_onm[y] * (
                sum(
                    model_data.omega[t] * (utility.v_E_my[y, k, t] * y_E[y, k, t]) for
                    t in model_data.index_t, k in utility.index_k_existing
                ) + sum(
                    model_data.omega[t] * (utility.v_C_my[y, k, t] * y_C[y, k, t]) for
                    t in model_data.index_t, k in utility.index_k_new
                )
            ) +
            # fixed o&m costs
            #   fom * (cap exist - cap retiring) + fom * cap new for gen type
            # the discount factor is applied to fom of existing capacity remaining at year y
            utility.pvf_onm[y] * sum(
                utility.fom_E_my[y, k] * (
                    utility.x_E_my[k] - sum(
                        utility.x_R_feasible[Symbol(Int(y_symbol)), k] for y_symbol in
                        model_data.year[first(model_data.index_y)]:model_data.year[y]
                    )
                ) for k in utility.index_k_existing
            ) +
            # the discount factor is applied to fom of new capacity for every year since year y (when it is built)
            sum(
                utility.fom_C_my[y, k] *
                utility.x_C_feasible[y, k] *
                sum(
                    utility.pvf_onm[Symbol(Int(y_symbol))] for y_symbol in
                    model_data.year[y]:model_data.year[last(model_data.index_y.elements)]
                ) for k in utility.index_k_new
            ) +
            # capital costs
            # the discout factor is applied to new capacity for the year it is built
            utility.pvf_cap[y] * sum(
                utility.CapEx_my[y, k] * utility.x_C_feasible[y, k] for
                k in utility.index_k_new
            )

            for y in model_data.index_y
        )
    end

    @objective(VIUDER_Utility_feasible, Min, objective_function)

    supply_demand_balance =
        (y, t) -> begin
            # bulk generation at time t
            sum(y_E[y, k, t] for k in utility.index_k_existing) +
            sum(y_C[y, k, t] for k in utility.index_k_new) -
            # demand at time t
            sum(customers.gamma[h] * customers.d_my[y, h, t] for h in model_data.index_h) - utility.eximport_my[y, t] +
            # existing DG generation at time t
            sum(
                customers.rho_DG[h, m, t] * customers.x_DG_E_my[y, h, m] for
                h in model_data.index_h, m in customers.index_m
            ) +
            # new DG generation at time t
            sum(
                customers.rho_DG[h, m, t] * sum(
                    customers.x_DG_new_my[Symbol(Int(y_symbol)), h, m] for y_symbol in
                    model_data.year[first(model_data.index_y_fix)]:model_data.year[y]
                ) for h in model_data.index_h, m in customers.index_m
            )
        end

    @constraint(
        VIUDER_Utility_feasible,
        Eq_miu[y in model_data.index_y, t in model_data.index_t],
        supply_demand_balance(y, t) >= 0
    )

    # HERE -- once running try defining function over two indices
    # y_E must be less than available capacity
    @constraint(
        VIUDER_Utility_feasible,
        Eq_eta[
            y in model_data.index_y,
            k in utility.index_k_existing,
            t in model_data.index_t,
        ],
        utility.rho_E_my[k, t] * (
            utility.x_E_my[k] - sum(
                utility.x_R_feasible[Symbol(Int(y_symbol)), k] for
                y_symbol in model_data.year[first(model_data.index_y)]:model_data.year[y]
            ) - utility.x_R_cumu[k]
        ) - y_E[y, k, t] >= 0
    )
    # y_C must be less than available capacity
    @constraint(
        VIUDER_Utility_feasible,
        Eq_lambda[
            y in model_data.index_y,
            k in utility.index_k_new,
            t in model_data.index_t,
        ],
        utility.rho_C_my[k, t] * (
            sum(
                utility.x_C_feasible[Symbol(Int(y_symbol)), k] for
                y_symbol in model_data.year[first(model_data.index_y)]:model_data.year[y]
            ) + utility.x_C_cumu[k]
        ) - y_C[y, k, t] >= 0
    )

    # RPS constraint
    @constraint(
        VIUDER_Utility_feasible,
        Eq_rps[y in model_data.index_y],
        sum(
            model_data.omega[t] * (y_E[y, rps, t] + y_C[y, rps, t]) for
            rps in utility.index_rps, t in model_data.index_t
        ) -
        utility.RPS[y] *
        sum(model_data.omega[t] * utility.Net_Load_my[y, t] for t in model_data.index_t) >=
        0
    )

    TimerOutputs.@timeit HEM_TIMER "optimize! VIUDER_Utility_feasible" begin
        optimize!(VIUDER_Utility_feasible)
    end

    update!(utility.obj_feasible, objective_value(VIUDER_Utility_feasible))
end

function solve_agent_problem_decomposition_by_year_master(
    utility::Utility,
    utility_opts::AgentOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{VerticallyIntegratedUtility},
    agent_store::AgentStore,
    w_iter,
)
    global optimality_gap = 100.0
    global tol = 0.000001  #0.001
    global rho = 1.0
    global alpha = 1.0
    global lagrange_iter = 1
    global max_iter_utility = 50   #100

    # i = 1
    # global dual_optimum = zeros(max_iter_ipp)

    update!(utility.obj_upper_bound, Inf)
    update!(utility.obj_lower_bound, -Inf)

    # x_R_before = copy(ipp.x_R_my)
    # x_C_before = copy(ipp.x_C_my)

    while optimality_gap >= tol && lagrange_iter <= max_iter_utility
        for y in model_data.index_y
            solve_agent_problem_decomposition_by_year(
                utility,
                utility_opts,
                model_data,
                hem_opts,
                agent_store,
                w_iter,
                y,
            )
            solve_agent_problem_decomposition_by_year_feasible(
                utility,
                utility_opts,
                model_data,
                hem_opts,
                agent_store,
                w_iter,
                y,
            )
        end
        solve_agent_problem_decomposition_by_year_feasible_obj(
            utility,
            utility_opts,
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

        # if the supposedly lower bound is greater than the upper bound, modify lagrange multiplier to move away from this solution
        while sum(utility.obj_my[y] for y in model_data.index_y) > utility.obj_feasible
            for k in utility.index_k_existing,
                y_inv_ret in model_data.index_y,
                y_decision in model_data.index_y

                utility.L_R_my[y_inv_ret, y_decision, k] =
                    utility.L_R_my[y_inv_ret, y_decision, k] + 1e-9
            end
            for k in utility.index_k_new,
                y_inv_ret in model_data.index_y,
                y_decision in model_data.index_y

                utility.L_C_my[y_inv_ret, y_decision, k] =
                    utility.L_C_my[y_inv_ret, y_decision, k] + 1e-9
            end
            for y in model_data.index_y
                solve_agent_problem_decomposition_by_year(
                    utility,
                    utility_opts,
                    model_data,
                    hem_opts,
                    agent_store,
                    w_iter,
                    y,
                )
                solve_agent_problem_decomposition_by_year_feasible(
                    utility,
                    utility_opts,
                    model_data,
                    hem_opts,
                    agent_store,
                    w_iter,
                    y,
                )
            end
            solve_agent_problem_decomposition_by_year_feasible_obj(
                utility,
                utility_opts,
                model_data,
                hem_opts,
                agent_store,
                w_iter,
            )
        end

        # update lower bound
        if sum(utility.obj_my[y] for y in model_data.index_y) < utility.obj_upper_bound &&
           sum(utility.obj_my[y] for y in model_data.index_y) > utility.obj_lower_bound
            update!(
                utility.obj_lower_bound,
                sum(utility.obj_my[y] for y in model_data.index_y),
            )
            # utility.x_R_year_1_dual = utility.x_R_year_1
            # utility.x_C_year_1_dual = utility.x_C_year_1
            # utility.x_R_year_2_dual = utility.x_R_year_2
            # utility.x_C_year_2_dual = utility.x_C_year_2
        end

        feasible = utility.obj_feasible

        # update upper bound
        if utility.obj_feasible < utility.obj_upper_bound
            @info "Iteration $lagrange_iter feasibility $feasible"
            update!(utility.obj_upper_bound, utility.obj_feasible.value)
            for y in model_data.index_y, k in utility.index_k_existing
                utility.x_R_my[y, k] = utility.x_R_feasible[y, k]
            end
            for y in model_data.index_y, k in utility.index_k_new
                utility.x_C_my[y, k] = utility.x_C_feasible[y, k]
            end
            # for y in model_data.index_y, p in ipp.index_p, k in ipp.index_k_existing, t in model_data.index_t
            #     ipp.y_E_my[y,p,k,t] = ipp.y_E_my_temp[y,p,k,t]
            # end
            # for y in model_data.index_y, p in ipp.index_p, k in ipp.index_k_new, t in model_data.index_t
            #     ipp.y_C_my[y,p,k,t] = ipp.y_C_my_temp[y,p,k,t]
            # end
            # for y in model_data.index_y, t in model_data.index_t
            #     ipp.miu_my[y,t] = ipp.miu_my_temp[y,t]
            # end
            # for y in model_data.index_y, t in model_data.index_t
            #     ipp.LMP_my[y,t] = ipp.LMP_my_temp[y,t]
            # end
            # for y in model_data.index_y
            #     ipp.capacity_price[y] = ipp.capacity_price_my_temp[y]
            #     ipp.ucap[y,p_star] = ipp.ucap_temp[y,p_star]
            # end
        end

        optimality_gap =
            (utility.obj_upper_bound - utility.obj_lower_bound) / utility.obj_upper_bound
        UB = utility.obj_upper_bound
        LB = utility.obj_lower_bound

        alpha_denominator = 0
        for y_inv_ret in
            model_data.year[first(model_data.index_y)]:(model_data.year[last(
            model_data.index_y.elements,
        )] - 1)
            for y_decision in
                (y_inv_ret + 1):model_data.year[last(model_data.index_y.elements)]
                alpha_denominator =
                    alpha_denominator +
                    sum(
                        (
                            utility.x_C_my_decomp[
                                Symbol(Int(y_inv_ret)),
                                Symbol(Int(y_inv_ret)),
                                k,
                            ] - utility.x_C_my_decomp[
                                Symbol(Int(y_inv_ret)),
                                Symbol(Int(y_decision)),
                                k,
                            ]
                        )^2 for k in utility.index_k_existing
                    ) +
                    sum(
                        (
                            utility.x_R_my_decomp[
                                Symbol(Int(y_inv_ret)),
                                Symbol(Int(y_inv_ret)),
                                k,
                            ] - utility.x_R_my_decomp[
                                Symbol(Int(y_inv_ret)),
                                Symbol(Int(y_decision)),
                                k,
                            ]
                        )^2 for k in utility.index_k_new
                    )
            end
        end

        alpha =
            rho * (utility.obj_upper_bound - utility.obj_lower_bound) / alpha_denominator

        for y_inv_ret in
            model_data.year[first(model_data.index_y)]:(model_data.year[last(
            model_data.index_y.elements,
        )] - 1)
            for y_decision in
                (y_inv_ret + 1):model_data.year[last(model_data.index_y.elements)]
                for k in utility.index_k_existing
                    utility.L_R_my[Symbol(Int(y_inv_ret)), Symbol(Int(y_decision)), k] =
                        utility.L_R_my[Symbol(Int(y_inv_ret)), Symbol(Int(y_decision)), k] +
                        alpha * (
                            utility.x_R_my_decomp[
                                Symbol(Int(y_inv_ret)),
                                Symbol(Int(y_inv_ret)),
                                k,
                            ] - utility.x_R_my_decomp[
                                Symbol(Int(y_inv_ret)),
                                Symbol(Int(y_decision)),
                                k,
                            ]
                        )
                end
                for k in utility.index_k_new
                    utility.L_C_my[Symbol(Int(y_inv_ret)), Symbol(Int(y_decision)), k] =
                        utility.L_C_my[Symbol(Int(y_inv_ret)), Symbol(Int(y_decision)), k] +
                        alpha * (
                            utility.x_C_my_decomp[
                                Symbol(Int(y_inv_ret)),
                                Symbol(Int(y_inv_ret)),
                                k,
                            ] - utility.x_C_my_decomp[
                                Symbol(Int(y_inv_ret)),
                                Symbol(Int(y_decision)),
                                k,
                            ]
                        )
                end
            end
        end

        lagrange_iter += 1
        @info "Lagrange iteration $lagrange_iter optimality gap: $optimality_gap"
        @info "Iteration $lagrange_iter upper bound: $UB"
        @info "Iteration $lagrange_iter lower bound: $LB"

        #=
        @info "Iteration $lagrange_iter stage 1 objective: $objective_st1"
        @info "Iteration $lagrange_iter stage 2 objective: $objective_st2"
        =#

        # dual_optimum[i] = utility.obj_year_1 + utility.obj_year_2
        # i = i + 1

    end

    # return compute_difference_one_norm([
    #     (x_R_before, ipp.x_R_my),
    #     (x_C_before, ipp.x_C_my)
    # ])

end

"""
Update Utility cumulative parameters
"""
function update_cumulative!(model_data::HEMData, utility::Utility)
    for k in utility.index_k_existing
        utility.x_R_cumu[k] = utility.x_R_cumu[k] + utility.x_R_my[first(model_data.index_y),k]
    end

    for k in utility.index_k_new
        utility.x_C_cumu[k] = utility.x_C_cumu[k] + utility.x_C_my[first(model_data.index_y),k]
    end
end
