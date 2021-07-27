# This module defines the data and functions associated with the Independent Power Producer

abstract type AbstractIPP <: Agent end

struct IPPCommon
    # TODO: populate this struct with fields common to all IPPs.
    # Everything in this struct should probably be immutable.
end

mutable struct IPP <: AbstractIPP
    id::String
    # Sets
    common::IPPCommon
    index_k::Dimension # bulk generation technologies

    # Parameters
    "existing capacity (MW)"
    x_E::ParamVector
    "fixed cost of existing capacity (\$/MW-yr)"
    f_E::ParamVector
    "fixed cost of new capacity (\$/MW-yr)"
    f_C::ParamVector
    "variable cost of existing capacity (\$/MWh)"
    v_E::ParamArray
    "variable cost of new capacity (\$/MWh)"
    v_C::ParamArray
    "availability of existing capacity (fraction)"
    rho_E::ParamArray
    "availability of new capacity (fraction)"
    rho_C::ParamArray
    "Big M Parameter"
    B1GM::ParamScalar{<:Integer}

    # Primal Variables
    y_E::ParamArray
    y_C::ParamArray
    x_R::ParamVector
    x_C::ParamVector
    miu::ParamVector
end

function IPP(input_filename::String, model_data::HEMData, id = DEFAULT_ID)
    index_k = read_set(
        input_filename,
        "index_k",
        "index_k",
        prose_name = "bulk generation technologies k",
    )

    return IPP(
        id,
        IPPCommon(),
        index_k,
        read_param(
            "x_E",
            input_filename,
            "ExistingCapacity",
            index_k,
            description = "existing capacity (MW)",
        ),
        read_param(
            "f_E",
            input_filename,
            "FixedCostOld",
            index_k,
            description = "fixed cost of existing capacity (\$/MW-yr)",
        ),
        read_param(
            "f_C",
            input_filename,
            "FixedCostNew",
            index_k,
            description = "fixed cost of new capacity (\$/MW-yr)",
        ),
        read_param(
            "v_E",
            input_filename,
            "VariableCostOld",
            model_data.index_t,
            [index_k],
            description = "variable cost of existing capacity (\$/MWh)",
        ),
        read_param(
            "v_C",
            input_filename,
            "VariableCostNew",
            model_data.index_t,
            [index_k],
            description = "variable cost of new capacity (\$/Mwh)",
        ),
        read_param(
            "rho_E",
            input_filename,
            "AvailabilityOld",
            model_data.index_t,
            [index_k],
            description = "availability of existing capacity (fraction)",
        ),
        read_param(
            "rho_C",
            input_filename,
            "AvailabilityNew",
            model_data.index_t,
            [index_k],
            description = "availability of new capacity (fraction)",
        ),
        ParamScalar("B1GM", 1000000),
        initialize_param(
            "y_E",
            index_k,
            model_data.index_t,
            description = "existing capacity generation in each time period",
        ),
        initialize_param(
            "y_C",
            index_k,
            model_data.index_t,
            description = "new capacity generation in each time period",
        ),
        initialize_param("x_R", index_k),
        initialize_param("x_C", index_k),
        initialize_param("miu", model_data.index_t),
    )
end

get_id(x::IPP) = x.id

function solve_agent_problem!(
    ipp::IPP,
    ipp_opts::AgentOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{VerticallyIntegratedUtility},
    agent_store::AgentStore,
)
    return 0.0
end

function solve_agent_problem!(
    ipp::IPP,
    ipp_opts::AgentOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{WholesaleMarket},
    agent_store::AgentStore,
)
    regulator = get_agent(Regulator, agent_store)
    customers = get_agent(Customers, agent_store)

    WMDER_IPP = get_new_jump_model(hem_opts.solver)

    # Define positive variables
    @variable(WMDER_IPP, x_C[ipp.index_k] >= 0)
    @variable(WMDER_IPP, x_R[ipp.index_k] >= 0)
    @variable(WMDER_IPP, y_E[ipp.index_k, model_data.index_t] >= 0)
    @variable(WMDER_IPP, y_C[ipp.index_k, model_data.index_t] >= 0)
    @variable(WMDER_IPP, miu[model_data.index_t] >= 0)
    @variable(WMDER_IPP, eta[ipp.index_k, model_data.index_t] >= 0)
    @variable(WMDER_IPP, lambda[ipp.index_k, model_data.index_t] >= 0)

    @variable(WMDER_IPP, u_y_E[ipp.index_k, model_data.index_t], Bin)
    @variable(WMDER_IPP, u_y_C[ipp.index_k, model_data.index_t], Bin)
    @variable(WMDER_IPP, u_miu[model_data.index_t], Bin)
    @variable(WMDER_IPP, u_eta[ipp.index_k, model_data.index_t], Bin)
    @variable(WMDER_IPP, u_lambda[ipp.index_k, model_data.index_t], Bin)

    objective_function = begin
        # Linearized revenue term 
        sum(
            miu[t] * (
                sum(customers.gamma[h] * customers.d[h, t] for h in model_data.index_h) - sum(
                    customers.rho_DG[h, m, t] * customers.x_DG_E[h, m] for
                    h in model_data.index_h, m in customers.index_m
                ) - sum(
                    customers.rho_DG[h, m, t] * customers.x_DG_new[h, m] for
                    h in model_data.index_h, m in customers.index_m
                )
            ) for t in model_data.index_t
        ) -
        # generation costs
        #   num hrs * ((fuel + vom) * gen existing + (fuel + vom) * gen new) for t and gen type
        sum(
            model_data.omega[t] * (ipp.v_E[k, t] * y_E[k, t] + ipp.v_C[k, t] * y_C[k, t]) for
            t in model_data.index_t, k in ipp.index_k
        ) -
        # fixed costs
        #   fom * (cap exist - cap retiring) + fom * cap new for gen type
        sum(ipp.f_E[k] * (ipp.x_E[k] - x_R[k]) + ipp.f_C[k] * x_C[k] for k in ipp.index_k)
    end

    @objective(WMDER_IPP, Max, objective_function)

    @constraint(
        WMDER_IPP,
        Eq_y_E_1[k in ipp.index_k, t in model_data.index_t],
        y_E[k, t] <= ipp.B1GM * u_y_E[k, t]
    )
    @constraint(
        WMDER_IPP,
        Eq_y_E_2[k in ipp.index_k, t in model_data.index_t],
        model_data.omega[t] * ipp.v_E[k, t] - miu[t] + eta[k, t] <=
        ipp.B1GM * (1 - u_y_E[k, t])
    )
    @constraint(
        WMDER_IPP,
        Eq_y_E_3[k in ipp.index_k, t in model_data.index_t],
        model_data.omega[t] * ipp.v_E[k, t] - miu[t] + eta[k, t] >= 0
    )

    @constraint(
        WMDER_IPP,
        Eq_y_C_1[k in ipp.index_k, t in model_data.index_t],
        y_C[k, t] <= ipp.B1GM * u_y_C[k, t]
    )
    @constraint(
        WMDER_IPP,
        Eq_y_C_2[k in ipp.index_k, t in model_data.index_t],
        model_data.omega[t] * ipp.v_C[k, t] - miu[t] + lambda[k, t] <=
        ipp.B1GM * (1 - u_y_C[k, t])
    )
    @constraint(
        WMDER_IPP,
        Eq_y_C_3[k in ipp.index_k, t in model_data.index_t],
        model_data.omega[t] * ipp.v_C[k, t] - miu[t] + lambda[k, t] >= 0
    )

    supply_demand_balance =
        t -> begin
            # bulk generation at time t
            sum(y_E[k, t] + y_C[k, t] for k in ipp.index_k) -
            # demand at time t
            sum(customers.gamma[h] * customers.d[h, t] for h in model_data.index_h) +
            # existing DG generation at time t
            sum(
                customers.rho_DG[h, m, t] * customers.x_DG_E[h, m] for
                h in model_data.index_h, m in customers.index_m
            ) +
            # new DG generation at time t
            sum(
                customers.rho_DG[h, m, t] * customers.x_DG_new[h, m] for
                h in model_data.index_h, m in customers.index_m
            )
        end

    @constraint(WMDER_IPP, Eq_miu_1[t in model_data.index_t], miu[t] <= ipp.B1GM * u_miu[t])
    @constraint(
        WMDER_IPP,
        Eq_miu_2[t in model_data.index_t],
        supply_demand_balance(t) <= ipp.B1GM * (1 - u_miu[t])
    )
    @constraint(WMDER_IPP, Eq_miu_3[t in model_data.index_t], supply_demand_balance(t) >= 0)

    @constraint(
        WMDER_IPP,
        Eq_eta_1[k in ipp.index_k, t in model_data.index_t],
        eta[k, t] <= ipp.B1GM * u_eta[k, t]
    )
    @constraint(
        WMDER_IPP,
        Eq_eta_2[k in ipp.index_k, t in model_data.index_t],
        ipp.rho_E[k, t] * (ipp.x_E[k] - x_R[k]) - y_E[k, t] <= ipp.B1GM * (1 - u_eta[k, t])
    )
    @constraint(
        WMDER_IPP,
        Eq_eta_3[k in ipp.index_k, t in model_data.index_t],
        ipp.rho_E[k, t] * (ipp.x_E[k] - x_R[k]) - y_E[k, t] >= 0
    )

    @constraint(
        WMDER_IPP,
        Eq_lambda_1[k in ipp.index_k, t in model_data.index_t],
        lambda[k, t] <= ipp.B1GM * u_lambda[k, t]
    )
    @constraint(
        WMDER_IPP,
        Eq_lambda_2[k in ipp.index_k, t in model_data.index_t],
        ipp.rho_C[k, t] * x_C[k] - y_C[k, t] <= ipp.B1GM * (1 - u_lambda[k, t])
    )
    @constraint(
        WMDER_IPP,
        Eq_lambda_3[k in ipp.index_k, t in model_data.index_t],
        ipp.rho_C[k, t] * x_C[k] - y_C[k, t] >= 0
    )

    @constraint(WMDER_IPP, Eq_sigma[k in ipp.index_k], ipp.x_E[k] - x_R[k] >= 0)

    planning_reserves =
        t -> begin
            # bulk generation available capacity at time t
            sum(
                ipp.rho_E[k, t] * (ipp.x_E[k] - x_R[k]) + ipp.rho_C[k, t] * x_C[k] for
                k in ipp.index_k
            ) -
            # net_load plus planning reserve
            (1 + regulator.r) * (
                sum(customers.gamma[h] * customers.d[h, t] for h in model_data.index_h) - sum(
                    customers.rho_DG[h, m, t] * customers.x_DG_E[h, m] for
                    h in model_data.index_h, m in customers.index_m
                ) - sum(
                    customers.rho_DG[h, m, t] * customers.x_DG_new[h, m] for
                    h in model_data.index_h, m in customers.index_m
                )
            )
        end
    @constraint(WMDER_IPP, Eq_xi[t in model_data.index_t], planning_reserves(t) >= 0)

    optimize!(WMDER_IPP)

    # record current primary variable values
    for k in ipp.index_k, t in model_data.index_t
        ipp.y_E[k, t] = value.(y_E[k, t])
        ipp.y_C[k, t] = value.(y_C[k, t])
    end

    x_R_before = ParamVector(ipp.x_R)
    x_C_before = ParamVector(ipp.x_C)
    for k in ipp.index_k
        ipp.x_R[k] = value.(x_R[k])
        ipp.x_C[k] = value.(x_C[k])
    end

    for t in model_data.index_t
        ipp.miu[t] = value.(miu[t])
    end

    @info "Original built capacity" x_C_before
    @info "New built capacity" ipp.x_C

    # report change in key variables from previous iteration to this one
    return compute_difference_one_norm([(x_R_before, ipp.x_R), (x_C_before, ipp.x_C)])
end

function save_results(
    ipp::IPP,
    ipp_opts::AgentOptions,
    hem_opts::HEMOptions{WholesaleMarket},
    export_file_path::AbstractString,
    fileprefix::AbstractString,
)
    # Primal Variables
    save_param(
        ipp.y_E.values,
        [:GenTech, :Time],
        :Generation_MWh,
        joinpath(export_file_path, "$(fileprefix)_y_E.csv"),
    )
    save_param(
        ipp.y_C.values,
        [:GenTech, :Time],
        :Generation_MWh,
        joinpath(export_file_path, "$(fileprefix)_y_C.csv"),
    )
    save_param(
        ipp.x_R.values,
        [:GenTech],
        :Capacity_MW,
        joinpath(export_file_path, "$(fileprefix)_x_R.csv"),
    )
    save_param(
        ipp.x_C.values,
        [:GenTech],
        :Capacity_MW,
        joinpath(export_file_path, "$(fileprefix)_x_C.csv"),
    )
end

function welfare_calculation(
    ipp::IPP,
    model_data::HEMData,
    regulator::Agent,
    customers::Agent,
)
    h = :Residential
    IPP_Revenue =
    # Linearized revenue term 
        sum(
            ipp.miu[t] * (
                sum(customers.gamma[h] * customers.d[h, t] for h in model_data.index_h) -
                sum(
                    customers.rho_DG[h, m, t] * customers.x_DG_E[h, m] for
                    h in model_data.index_h, m in customers.index_m
                ) - sum(
                    customers.rho_DG[h, m, t] * customers.x_DG_new[h, m] for
                    h in model_data.index_h, m in customers.index_m
                )
            ) for t in model_data.index_t
        )
    IPP_Cost = # generation costs
    #   num hrs * ((fuel + vom) * gen existing + (fuel + vom) * gen new) for t and gen type
        sum(
            model_data.omega[t] *
            (ipp.v_E[k, t] * ipp.y_E[k, t] + ipp.v_C[k, t] * ipp.y_C[k, t]) for
            t in model_data.index_t, k in ipp.index_k
        ) +
        # fixed costs
        #   fom * (cap exist - cap retiring) + fom * cap new for gen type
        sum(
            ipp.f_E[k] * (ipp.x_E[k] - ipp.x_R[k]) + ipp.f_C[k] * ipp.x_C[k] for
            k in ipp.index_k
        )

    return IPP_Revenue, IPP_Cost
end
