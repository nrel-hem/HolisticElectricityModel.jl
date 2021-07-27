# This file defines the data and functions associated with the utility.

mutable struct Utility <: Agent
    # Sets
    # bulk generation technologies
    index_k::Set1D

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

    # Primal Variables
    y_E::ParamArray
    y_C::ParamArray
    x_R::ParamVector
    x_C::ParamVector
    # Dual Variables
    miu::ParamVector
end

function Utility(input_filename::String, model_data::HEMData)
    index_k = read_set(input_filename, "index_k", "index_k",
                       prose_name = "bulk generation technologies k")

    return Utility(
        index_k,
        read_param("x_E", input_filename, "ExistingCapacity", index_k,
            description = "existing capacity (MW)"),
        read_param("f_E", input_filename, "FixedCostOld", index_k,
            description = "fixed cost of existing capacity (\$/MW-yr)"),
        read_param("f_C", input_filename, "FixedCostNew", index_k,
            description = "fixed cost of new capacity (\$/MW-yr)"),
        read_param("v_E", input_filename, "VariableCostOld", model_data.index_t, 
            row_indices = [index_k],
            description = "variable cost of existing capacity (\$/MWh)"),
        read_param("v_C", input_filename, "VariableCostNew", model_data.index_t, 
            row_indices = [index_k],
            description = "variable cost of new capacity (\$/Mwh)"),
        read_param("rho_E", input_filename, "AvailabilityOld", model_data.index_t,
            row_indices = [index_k],
            description = "availability of existing capacity (fraction)"),
        read_param("rho_C", input_filename, "AvailabilityNew", model_data.index_t, 
            row_indices = [index_k],
            description = "availability of new capacity (fraction)"),
        initialize_param("y_E", index_k, model_data.index_t,
            description = "existing capacity generation in each time period"),
        initialize_param("y_C", index_k, model_data.index_t,
            description = "new capacity generation in each time period"),
        initialize_param("x_R", index_k),
        initialize_param("x_C", index_k),
        initialize_param("miu", model_data.index_t)
    )
end

function solve_agent_problem(
        utility::Utility,
        utility_opts::AgentOptions,
        model_data::HEMData,
        hem_opts::HEMOptions{WholesaleMarket},
        other_agents::Vector{Agent}
    )

    return 0.0
end

function solve_agent_problem(
        utility::Utility, 
        utility_opts::AgentOptions,
        model_data::HEMData, 
        hem_opts::HEMOptions{VerticallyIntegratedUtility},
        other_agents::Vector{Agent}
    )

    regulator = get_agent(other_agents, Regulator)
    customers = get_agent(other_agents, Customers)

    VIUDER_Utility = get_new_jump_model(hem_opts.solver)

    # Define positive variables
    @variable(VIUDER_Utility, x_C[utility.index_k] >= 0)
    @variable(VIUDER_Utility, x_R[utility.index_k] >= 0)
    @variable(VIUDER_Utility, y_E[utility.index_k, model_data.index_t] >= 0)
    @variable(VIUDER_Utility, y_C[utility.index_k, model_data.index_t] >= 0)

    objective_function = begin
        # generation costs
        #   num hrs * ((fuel + vom) * gen existing + (fuel + vom) * gen new) for t and gen type
        sum(model_data.omega[t]*(utility.v_E[k,t]*y_E[k,t] + utility.v_C[k,t]*y_C[k,t]) 
            for t in model_data.index_t, k in utility.index_k) + 
        # fixed costs
        #   fom * (cap exist - cap retiring) + fom * cap new for gen type
        sum(utility.f_E[k]*(utility.x_E[k]-x_R[k]) + utility.f_C[k]*x_C[k] 
            for k in utility.index_k)
    end

    @objective(VIUDER_Utility, Min, objective_function)

    supply_demand_balance = t->begin
        # bulk generation at time t
        sum(y_E[k,t] + y_C[k,t] for k in utility.index_k) - 
        # demand at time t
        sum(customers.gamma[h]*customers.d[h,t] for h in model_data.index_h) + 
        # existing DG generation at time t
        sum(customers.rho_DG[h,m,t]*customers.x_DG_E[h,m] for h in model_data.index_h, m in customers.index_m) + 
        # new DG generation at time t
        sum(customers.rho_DG[h,m,t]*customers.x_DG_new[h,m] for h in model_data.index_h, m in customers.index_m)
    end
    
    @constraint(VIUDER_Utility, Eq_miu[t in model_data.index_t], supply_demand_balance(t) >= 0)

    # HERE -- once running try defining function over two indices
    # y_E must be less than available capacity
    @constraint(VIUDER_Utility, 
                Eq_eta[k in utility.index_k, t in model_data.index_t], 
                utility.rho_E[k,t]*(utility.x_E[k]-x_R[k])-y_E[k,t] >= 0)
    # y_C must be less than available capacity
    @constraint(VIUDER_Utility, 
                Eq_lambda[k in utility.index_k, t in model_data.index_t], 
                utility.rho_C[k,t]*x_C[k]-y_C[k,t] >= 0)
    # retiring capacity is bounded by existing capacity
    @constraint(VIUDER_Utility, 
                Eq_sigma[k in utility.index_k], 
                utility.x_E[k] - x_R[k] >= 0)

    planning_reserves = t->begin
        # bulk generation available capacity at time t
        sum(utility.rho_E[k,t]*(utility.x_E[k] - x_R[k]) + utility.rho_C[k,t]*x_C[k] 
            for k in utility.index_k) - 
        # net_load plus planning reserve
        (1+regulator.r)*(
            sum(customers.gamma[h]*customers.d[h,t] for h in model_data.index_h) - 
            sum(customers.rho_DG[h,m,t]*customers.x_DG_E[h,m] for h in model_data.index_h, m in customers.index_m) - 
            sum(customers.rho_DG[h,m,t]*customers.x_DG_new[h,m] for h in model_data.index_h, m in customers.index_m)
        )
    end
    @constraint(VIUDER_Utility, Eq_xi[t in model_data.index_t], planning_reserves(t) >= 0)
    
    optimize!(VIUDER_Utility)
        
    # record current primary variable values
    for k in utility.index_k, t in model_data.index_t
        utility.y_E[k,t] = value.(y_E[k,t])
        utility.y_C[k,t] = value.(y_C[k,t])
    end
    
    x_R_before = copy(utility.x_R)
    x_C_before = copy(utility.x_C)
    for k in utility.index_k
        utility.x_R[k] = value.(x_R[k])
        utility.x_C[k] = value.(x_C[k])
    end

    for t in model_data.index_t
        utility.miu[t] = dual.(Eq_miu[t])
    end

    @info "Original built capacity" x_C_before
    @info "New built capacity" utility.x_C

    # report change in key variables from previous iteration to this one
    return compute_difference_one_norm([
        (x_R_before, utility.x_R),
        (x_C_before, utility.x_C)
    ])
end

function save_results(
        utility::Utility, 
        utility_opts::AgentOptions,
        hem_opts::HEMOptions{VerticallyIntegratedUtility},
        exportfilepath::AbstractString, 
        fileprefix::AbstractString)
    # Primal Variables
    save_param(utility.y_E.values, [:GenTech, :Time], :Generation_MWh, 
               joinpath(exportfilepath, "$(fileprefix)_y_E.csv"))
    save_param(utility.y_C.values, [:GenTech, :Time], :Generation_MWh, 
               joinpath(exportfilepath, "$(fileprefix)_y_C.csv"))
    save_param(utility.x_R.values, [:GenTech], :Capacity_MW, 
               joinpath(exportfilepath, "$(fileprefix)_x_R.csv"))
    save_param(utility.x_C.values, [:GenTech], :Capacity_MW, 
               joinpath(exportfilepath, "$(fileprefix)_x_C.csv"))
end

function welfare_calculation(
    utility::Utility, 
    model_data::HEMData, 
    regulator::Agent,    
    customers::Agent)
    
    h = :Residential
    Utility_Revenue = sum(model_data.omega[t] * regulator.p[h,t] * sum(utility.y_E[k,t]+utility.y_C[k,t] for k in utility.index_k) 
        for t in model_data.index_t)
    Utility_Cost = # generation costs
        #   num hrs * ((fuel + vom) * gen existing + (fuel + vom) * gen new) for t and gen type
        sum(model_data.omega[t]*(utility.v_E[k,t]*utility.y_E[k,t] + utility.v_C[k,t]*utility.y_C[k,t]) 
        for t in model_data.index_t, k in utility.index_k) + 
        # fixed costs
        #   fom * (cap exist - cap retiring) + fom * cap new for gen type
        sum(utility.f_E[k]*(utility.x_E[k]-utility.x_R[k]) + utility.f_C[k]*utility.x_C[k] 
        for k in utility.index_k)

    return Utility_Revenue, Utility_Cost

end

