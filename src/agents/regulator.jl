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


struct RegulatorOptions{T <: RateDesign, U <: NetMeteringPolicy} <: AgentOptions
    rate_design::T
    net_metering_policy::U
end


mutable struct Regulator <: Agent
    # Parameters
    r::ParamScalar # planning reserve (fraction)
    z::ParamScalar # allowed return on investment (fraction)

    # Primal Variables
    p::ParamArray # retail price
    p_ex::ParamArray # DER excess generation rate
end

function Regulator(input_filename::String, model_data::HEMData)
    return Regulator(
        ParamScalar("r", 0.14, description = "planning reserve (fraction)"),
        ParamScalar("z", 0.09, description = "allowed return on investment (fraction)"),
        initialize_param("p", model_data.index_h, model_data.index_t, value=10.0,
            description = "retail price"),
        initialize_param("p_ex", model_data.index_h, model_data.index_t, value=10.0,
            description = "DER excess generation rate")
    )
end

# although Customer is subtype of Agent, 
# Vector{Customer} is not subtype of Vector{Agent}
# But if a vector of customers c1, c2, c3 is defined 
# using the syntax Agent[c1, c2, c3], calling 
# this function will work. Can also:
# Vector{Agent}([c1, c2, c3])

function solve_agent_problem(
        regulator::Regulator,    
        regulator_opts::RegulatorOptions,
        model_data::HEMData, 
        hem_opts::HEMOptions{VerticallyIntegratedUtility},
        other_agents::Vector{Agent})

    utility = get_agent(other_agents, Utility)
    customers = get_agent(other_agents, Customers)
    
    # pure volumetric rate
    energy_cost = sum(
        model_data.omega[t]*(utility.v_E[k,t]*utility.y_E[k,t]+utility.v_C[k,t]*utility.y_C[k,t]) 
        for k in utility.index_k, t in model_data.index_t)
    capital_cost = sum(
        utility.f_E[k]*(utility.x_E[k]-utility.x_R[k]) + utility.f_C[k]*utility.x_C[k] 
        for k in utility.index_k)
    der_excess_cost = sum(model_data.omega[t]*regulator.p_ex[h,t] * max(0, sum(customers.rho_DG[h,m,t]*customers.Opti_DG[h,m] for m in customers.index_m)-customers.d[h,t]) *
        (customers.rho_DG[h,m,t]*customers.Opti_DG[h,m])/customers.DERGen[h,t] * (customers.x_DG_E[h,m]+customers.x_DG_new[h,m])/customers.Opti_DG[h,m]
        for t in model_data.index_t, h in model_data.index_h, m in customers.index_m)
    net_demand = (
        # demand
        sum(customers.gamma[h]*model_data.omega[t]*customers.d[h,t] 
            for h in model_data.index_h, t in model_data.index_t) - 
        # DG
        sum(min(sum(customers.rho_DG[h,m,t]*customers.Opti_DG[h,m] for m in customers.index_m), customers.d[h,t]) *
        (customers.rho_DG[h,m,t]*customers.Opti_DG[h,m])/customers.DERGen[h,t] * (customers.x_DG_E[h,m]+customers.x_DG_new[h,m])/customers.Opti_DG[h,m]
            for h in model_data.index_h, t in model_data.index_t, m in customers.index_m)
    )

    energy_cost_t = Dict(t => sum(utility.v_E[k,t]*utility.y_E[k,t]+utility.v_C[k,t]*utility.y_C[k,t] for k in utility.index_k) for t in model_data.index_t)
    der_excess_cost_t = Dict(t => sum(regulator.p_ex[h,t] * max(0, sum(customers.rho_DG[h,m,t]*customers.Opti_DG[h,m] for m in customers.index_m)-customers.d[h,t]) *
        (customers.rho_DG[h,m,t]*customers.Opti_DG[h,m])/customers.DERGen[h,t] * (customers.x_DG_E[h,m]+customers.x_DG_new[h,m])/customers.Opti_DG[h,m]
        for h in model_data.index_h, m in customers.index_m) for t in model_data.index_t)
    net_demand_t = Dict(t => 
    # demand
        sum(customers.gamma[h]*customers.d[h,t] 
            for h in model_data.index_h) - 
    # DG
        sum(min(sum(customers.rho_DG[h,m,t]*customers.Opti_DG[h,m] for m in customers.index_m), customers.d[h,t]) *
        (customers.rho_DG[h,m,t]*customers.Opti_DG[h,m])/customers.DERGen[h,t] * (customers.x_DG_E[h,m]+customers.x_DG_new[h,m])/customers.Opti_DG[h,m]
            for h in model_data.index_h, m in customers.index_m)
            for t in model_data.index_t
    )

    # compute the retail price
    p_before = copy(regulator.p, "p_before")
    p_ex_before = copy(regulator.p_ex, "p_ex_before")

    # TODO: Call a function instead of using if-then
    if regulator_opts.rate_design isa FlatRate
        regulator.p = update(regulator.p, Dict(
            (h,t) => (energy_cost + der_excess_cost + (1 + regulator.z) * capital_cost) / net_demand
            for h in model_data.index_h, t in model_data.index_t))
    elseif regulator_opts.rate_design isa TOU
        regulator.p = update(regulator.p, Dict(
            (h,t) => (energy_cost_t[t] + der_excess_cost_t[t]) / net_demand_t[t] + (1.0 + regulator.z) * capital_cost / net_demand
            for h in model_data.index_h, t in model_data.index_t))
    end

    # TODO: Call a function instead of using if-then
    if regulator_opts.net_metering_policy isa ExcessRetailRate
        regulator.p_ex = update(regulator.p, regulator.p.values)
    elseif regulator_opts.net_metering_policy isa ExcessMarginalCost
        regulator.p_ex = update(regulator.p_ex, 
            Dict((h,t) => utility.miu[t]/model_data.omega[t] for h in model_data.index_h, t in model_data.index_t))
    elseif regulator_opts.net_metering_policy isa ExcessZero
        regulator.p_ex = update(regulator.p_ex,
            Dict((h,t) => 0.0 for h in model_data.index_h, t in model_data.index_t))
    end

    @info "Original retail price" p_before
    @info "Original DER excess rate" p_ex_before
    @info "New retail price" regulator.p
    @info "New DER excess rate" regulator.p_ex
                   
    return compute_difference_one_norm([
        (p_before.values, regulator.p.values),
        (p_ex_before.values, regulator.p_ex.values)
    ])    
end




function solve_agent_problem(
        regulator::Regulator,    
        regulator_opts::RegulatorOptions,
        model_data::HEMData, 
        hem_opts::HEMOptions{WholesaleMarket},
        other_agents::Vector{Agent})

    customers = get_agent(other_agents, Customers)
    ipp = get_agent(other_agents, IPP)

    # pure volumetric rate
    der_excess_cost = sum(model_data.omega[t]*regulator.p_ex[h,t] * max(0, sum(customers.rho_DG[h,m,t]*customers.Opti_DG[h,m] for m in customers.index_m)-customers.d[h,t]) *
        (customers.rho_DG[h,m,t]*customers.Opti_DG[h,m])/customers.DERGen[h,t] * (customers.x_DG_E[h,m]+customers.x_DG_new[h,m])/customers.Opti_DG[h,m]
        for t in model_data.index_t, h in model_data.index_h, m in customers.index_m)
    net_demand = (
        # demand
        sum(customers.gamma[h]*model_data.omega[t]*customers.d[h,t] 
            for h in model_data.index_h, t in model_data.index_t) - 
        # DG
        sum(min(sum(customers.rho_DG[h,m,t]*customers.Opti_DG[h,m] for m in customers.index_m), customers.d[h,t]) *
            (customers.rho_DG[h,m,t]*customers.Opti_DG[h,m])/customers.DERGen[h,t] * (customers.x_DG_E[h,m]+customers.x_DG_new[h,m])/customers.Opti_DG[h,m]
            for h in model_data.index_h, t in model_data.index_t, m in customers.index_m)
    )

    der_excess_cost_t = Dict(t => sum(regulator.p_ex[h,t] * max(0, sum(customers.rho_DG[h,m,t]*customers.Opti_DG[h,m] for m in customers.index_m)-customers.d[h,t]) *
        (customers.rho_DG[h,m,t]*customers.Opti_DG[h,m])/customers.DERGen[h,t] * (customers.x_DG_E[h,m]+customers.x_DG_new[h,m])/customers.Opti_DG[h,m]
        for h in model_data.index_h, m in customers.index_m) for t in model_data.index_t)
    net_demand_t = Dict(t => 
    # demand
        sum(customers.gamma[h]*customers.d[h,t] 
            for h in model_data.index_h) - 
    # DG
        sum(min(sum(customers.rho_DG[h,m,t]*customers.Opti_DG[h,m] for m in customers.index_m), customers.d[h,t]) *
            (customers.rho_DG[h,m,t]*customers.Opti_DG[h,m])/customers.DERGen[h,t] * (customers.x_DG_E[h,m]+customers.x_DG_new[h,m])/customers.Opti_DG[h,m]
            for h in model_data.index_h, m in customers.index_m)
            for t in model_data.index_t
    )

    # compute the retail price
    p_before = copy(regulator.p, "p_before")
    p_ex_before = copy(regulator.p_ex, "p_ex_before")

    # TODO: Call a function instead of using if-then
    if regulator_opts.rate_design isa FlatRate
        regulator.p = Dict(
            (h,t) => (sum(ipp.miu[t] * sum(ipp.y_E[k,t]+ipp.y_C[k,t] for k in ipp.index_k) for t in model_data.index_t) + der_excess_cost) / net_demand
            for h in model_data.index_h, t in model_data.index_t)
    elseif regulator_opts.rate_design isa TOU
        regulator.p = Dict(
            (h,t) => (ipp.miu[t]/model_data.omega[t] * sum(ipp.y_E[k,t]+ipp.y_C[k,t] for k in ipp.index_k) + der_excess_cost_t[t]) / net_demand_t[t]
            for h in model_data.index_h, t in model_data.index_t)
    end

    # TODO: Call a function instead of using if-then
    if regulator_opts.net_metering_policy isa ExcessRetailRate
        regulator.p_ex = regulator.p
    elseif regulator_opts.net_metering_policy isa ExcessMarginalCost
        regulator.p_ex = Dict((h,t) => ipp.miu[t]/model_data.omega[t] for h in model_data.index_h, t in model_data.index_t)
    elseif regulator_opts.net_metering_policy isa ExcessZero
        regulator.p_ex = Dict((h,t) => 0.0 for h in model_data.index_h, t in model_data.index_t)
    end

    @info "Original retail price" p_before
    @info "Original DER excess rate" p_ex_before
    @info "New retail price" regulator.p
    @info "New DER excess rate" regulator.p_ex
                   
    return compute_difference_one_norm([
        (p_before.values, regulator.p.values),
        (p_ex_before.values, regulator.p_ex.values)
    ])    
end

function save_results(
        regulator::Regulator, 
        regulator_opts::RegulatorOptions,
        hem_opts::HEMOptions,
        exportfilepath::AbstractString, 
        fileprefix::AbstractString)
    # Primal Variables
    save_param(regulator.p.values, [:CustomerType, :Time], :Price, 
               joinpath(exportfilepath, "$(fileprefix)_p.csv"))
    save_param(regulator.p_ex.values, [:CustomerType, :Time], :Price, 
               joinpath(exportfilepath, "$(fileprefix)_p_ex.csv"))
end
