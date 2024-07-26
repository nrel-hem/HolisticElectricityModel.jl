abstract type AbstractDERA <: Agent end

struct DERAOptions <: AgentOptions
    solvers::HEMSolver
    # solvers::Union{HEMSolver, Dict{String, <:HEMSolver}}
end

"""
Construct DERAOptions with an MOI.OptimizerWithAttributes instance.
"""
function DERAOptions(attributes::MOI.OptimizerWithAttributes)
    return DERAOptions(AnySolver(attributes))
end

mutable struct DERA <: AbstractDERA
    id::String
    # incentive function for DER aggregation (piece-wise linear function: incentive (x) vs participation (y))
    dera_stor_incentive_function::DataFrame
    # parameters applied to aggregated capacities of DER to represent potential friction of aggregation
    aggregation_friction::ParamArray
    # incentive level by year
    incentive_level::ParamArray
    # aggregation level by year
    aggregation_level::ParamArray
end

function DERA(input_filename::AbstractString, model_data::HEMData; id = DEFAULT_ID)

    # need to have the incentive function for each customer type
    dera_stor_incentive_function = CSV.read(joinpath(input_filename, "dera_stor_incentive_function.csv"), DataFrame)

    return DERA(
        id,
        dera_stor_incentive_function,
        initialize_param("aggregation_friction", model_data.index_y, model_data.index_z, model_data.index_h, value = 0.0),
        initialize_param("incentive_level", model_data.index_y, model_data.index_z, value = 0.0),
        initialize_param("aggregation_level", model_data.index_y, model_data.index_z, value = 0.0),
    )
end

get_id(x::DERA) = x.id

function solve_agent_problem!(
    der_aggregator::DERA,
    dera_opts::AgentOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{WholesaleMarket},
    agent_store::AgentStore,
    w_iter,
    jump_model
)

    reg_year = model_data.year(first(model_data.index_y))
    reg_year_index = Symbol(Int(reg_year))
    delta_t = parse(Int64, chop(string(model_data.index_t.elements[2]), head = 1, tail = 0)) - parse(Int64, chop(string(model_data.index_t.elements[1]), head = 1, tail = 0))

    ipp = get_agent(IPPGroup, agent_store)
    customers = get_agent(CustomerGroup, agent_store)

    # the x-axis of incentive function is assumed to be distributed storage capacity (even though the technology is PV+storage)
    incentive_function_dimension = size(der_aggregator.dera_stor_incentive_function)[1]
    incentive_function_slope = zeros(incentive_function_dimension - 1)
    incentive_function_intercept = zeros(incentive_function_dimension - 1)
    for i in 1:incentive_function_dimension - 1
        incentive_function_slope[i] = (der_aggregator.dera_stor_incentive_function[i+1, "participation"] - der_aggregator.dera_stor_incentive_function[i, "participation"]) /
        (der_aggregator.dera_stor_incentive_function[i+1, "incentive"] - der_aggregator.dera_stor_incentive_function[i, "incentive"])

        incentive_function_intercept[i] = der_aggregator.dera_stor_incentive_function[i+1, "participation"] - incentive_function_slope[i] * der_aggregator.dera_stor_incentive_function[i+1, "incentive"]
    end
    incentive_level_by_segment = Dict(z => zeros(incentive_function_dimension - 1) for z in model_data.index_z)
    participation_by_segment = Dict(z => zeros(incentive_function_dimension - 1) for z in model_data.index_z)
    obj_by_segment = Dict(z => zeros(incentive_function_dimension - 1) for z in model_data.index_z)

    total_der_stor_capacity = make_keyed_array(model_data.index_z, model_data.index_h)
    for z in model_data.index_z, h in model_data.index_h
        if w_iter >= 2
            total_der_stor_capacity(z, h, :) .=
                customers.x_DG_E_my(reg_year_index, h, z, :BTMStorage) + sum(
                    customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMStorage) for
                    y in model_data.year(first(model_data.index_y_fix)):(reg_year - 1)
                )
        else
            total_der_stor_capacity(z, h, :) .= customers.x_DG_E_my(reg_year_index, h, z, :BTMStorage)
        end
    end
    
    for z in model_data.index_z
        for i in 1:incentive_function_dimension - 1
            # x (incentive) should be $/MW (per year)?
            DERA_WM = get_new_jump_model(dera_opts.solvers)
            @variable(DERA_WM, der_aggregator.dera_stor_incentive_function[i, "incentive"] <= incentive <= der_aggregator.dera_stor_incentive_function[i+1, "incentive"])
            @variable(DERA_WM, dera_charge[model_data.index_d, model_data.index_t] >= 0)
            @variable(DERA_WM, dera_discharge[model_data.index_d, model_data.index_t] >= 0)
            @variable(DERA_WM, dera_energy[model_data.index_d, model_data.index_t] >= 0)
            @variable(DERA_WM, dera_pv_gen[model_data.index_d, model_data.index_t] >= 0)

            dera_stor_capacity_h =
                h -> begin

                (incentive_function_intercept[i] + incentive_function_slope[i] * incentive) * total_der_stor_capacity(z, h)

            end

            dera_stor_capacity =
                begin

                sum(dera_stor_capacity_h(h) for h in model_data.index_h)

            end

            # when it comes to the aggregated pv capacity, use aggregated storage capacity, divide by optimal storage size to get the number of households that's participating the dera aggregation,
            # then multiply by the optimal pv size,
            # TODO: I'm using optimal storage and pv capacities for existing resources, need to be more flexible to account for changing optimal tech sizes overtime.
            # right now, assume the optimal sizes are the same for both existing and new techs.
            dera_pv_capacity_h = 
                h -> begin

                dera_stor_capacity_h(h) / customers.Opti_DG_E(z, h, :BTMStorage) * customers.Opti_DG_E(z, h, :BTMPV)

            end

            dera_pv_capacity =
                begin

                sum(dera_pv_capacity_h(h) for h in model_data.index_h)

            end

            objective_revenue = begin
                sum(
                    model_data.omega(d) * delta_t *
                    ipp.LMP_my(reg_year_index, z, d, t) *
                    (dera_discharge[d, t] - dera_charge[d, t] + dera_pv_gen[d, t]) 
                    for d in model_data.index_d, t in model_data.index_t
                ) + ipp.capacity_price(reg_year_index) * (dera_stor_capacity .* ipp.capacity_credit_stor_E_my(reg_year_index, z, Symbol("der_aggregator")) + 
                dera_pv_capacity .* ipp.capacity_credit_E_my(reg_year_index, z, Symbol("dera_pv")))
            end

            objective_cost = begin
                incentive * dera_stor_capacity
            end

            @objective(DERA_WM, Max, objective_revenue - objective_cost)

            # DERA (storage) constraints -- We assume all distributed storage resources have the same initial energy, duration and round-trip-efficiency

            @constraint(
                DERA_WM,
                Eq_stor_energy_balance[
                    d in model_data.index_d,
                    t in model_data.index_t.elements[2:end],
                ],
                dera_energy[d, t] == dera_energy[d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)-delta_t), model_data.time.values)][1]] - dera_discharge[d, t] / customers.rte_dist_stor(z, :Commercial) * delta_t +
                    dera_charge[d, t] * delta_t
            )

            @constraint(
                DERA_WM,
                Eq_stor_energy_balance_0[
                    d in model_data.index_d,
                    t in [model_data.index_t.elements[1]],
                ],
                dera_energy[d, t] == customers.initial_energy_dist_stor(z, :Commercial, d) - dera_discharge[d, t] / customers.rte_dist_stor(z, :Commercial) * delta_t + dera_charge[d, t] * delta_t
            )

            @constraint(
                DERA_WM,
                Eq_stor_energy_upper_bound[
                    d in model_data.index_d,
                    t in model_data.index_t,
                ],
                dera_energy[d, t] <= customers.duration_dist_stor(z, :Commercial) * dera_stor_capacity
            )

            @constraint(
                DERA_WM,
                Eq_stor_discharge_upper_bound[
                    d in model_data.index_d,
                    t in model_data.index_t,
                ],
                dera_discharge[d, t] <= customers.rte_dist_stor(z, :Commercial) * dera_stor_capacity
            )

            @constraint(
                DERA_WM,
                Eq_stor_charge_upper_bound[
                    d in model_data.index_d,
                    t in model_data.index_t,
                ],
                dera_charge[d, t] <= dera_stor_capacity
            )

            @constraint(
                DERA_WM,
                Eq_stor_discharge_energy_upper_bound[
                    d in model_data.index_d,
                    t in model_data.index_t.elements[2:end],
                ],
                dera_discharge[d, t] * delta_t <= 
                dera_energy[d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)-delta_t), model_data.time.values)][1]]
            )

            @constraint(
                DERA_WM,
                Eq_stor_discharge_energy_upper_bound_0[
                    d in model_data.index_d,
                    t in [model_data.index_t.elements[1]],
                ],
                dera_discharge[d, t] * delta_t <= customers.initial_energy_dist_stor(z, :Commercial, d)
            )

            @constraint(
                DERA_WM,
                Eq_stor_charge_energy_upper_bound[
                    d in model_data.index_d,
                    t in model_data.index_t.elements[2:end],
                ],
                dera_charge[d, t] * delta_t <= customers.duration_dist_stor(z, :Commercial) * dera_stor_capacity -
                dera_energy[d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)-delta_t), model_data.time.values)][1]]
            )

            @constraint(
                DERA_WM,
                Eq_stor_charge_energy_upper_bound_0[
                    d in model_data.index_d,
                    t in [model_data.index_t.elements[1]],
                ],
                dera_charge[d, t] * delta_t <= customers.duration_dist_stor(z, :Commercial) * dera_stor_capacity -
                customers.initial_energy_dist_stor(z, :Commercial, d)
            )

            @constraint(
                DERA_WM,
                Eq_stor_charge_discharge_upper_bound[
                    d in model_data.index_d,
                    t in model_data.index_t.elements,
                ],
                dera_charge[d, t] + dera_discharge[d, t] / customers.rte_dist_stor(z, :Commercial) <= 
                dera_stor_capacity
            )

            @constraint(
                DERA_WM,
                Eq_pv_upper_bound[
                    d in model_data.index_d,
                    t in model_data.index_t,
                ],
                dera_pv_gen[d, t] <= ipp.rho_E_my(:ipp1, :dera_pv, z, d, t) * dera_pv_capacity
            )

            # TimerOutputs.@timeit HEM_TIMER "optimize! DER Aggregator dispatch" begin
                optimize!(DERA_WM)
            # end

            incentive_level_by_segment[z][i] = value(incentive)
            participation_by_segment[z][i] = incentive_function_intercept[i] + incentive_function_slope[i] * value(incentive)
            obj_by_segment[z][i] = objective_value(DERA_WM)

        end
    end

    max_seg_index = Dict(z => findfirst(x -> x == maximum(obj_by_segment[z]), obj_by_segment[z]) for z in model_data.index_z)
    for z in model_data.index_z
        der_aggregator.incentive_level(reg_year_index, z, :) .= incentive_level_by_segment[z][max_seg_index[z]]
        der_aggregator.aggregation_level(reg_year_index, z, :) .= participation_by_segment[z][max_seg_index[z]]
    end

    dera_agg_stor_capacity_h = make_keyed_array(model_data.index_z, model_data.index_h)
    dera_agg_pv_capacity_h = make_keyed_array(model_data.index_z, model_data.index_h)
    for z in model_data.index_z, h in model_data.index_h
        dera_agg_stor_capacity_h(z, h, :) .= der_aggregator.aggregation_level(reg_year_index, z) * total_der_stor_capacity(z, h)
        dera_agg_pv_capacity_h(z, h, :) .= dera_agg_stor_capacity_h(z, h) / customers.Opti_DG_E(z, h, :BTMStorage) * customers.Opti_DG_E(z, h, :BTMPV)
    end

    for z in model_data.index_z
        # simply assign DERA to a random ipp (ipp1)
        ipp.x_stor_E_my(:ipp1, z, Symbol("der_aggregator"), :) .= sum(dera_agg_stor_capacity_h(z, h) for h in model_data.index_h)
        ipp.x_E_my(:ipp1, z, Symbol("dera_pv"), :) .= sum(dera_agg_pv_capacity_h(z, h) for h in model_data.index_h)
    end

    # since we moved some BTM storage to transmission level, need to reduce the BTM net load accordingly (in bulk power system, regulator, customers (maybe?)).

    return 0.0
    
end

function solve_agent_problem!(
    der_aggregator::DERA,
    dera_opts::AgentOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{VerticallyIntegratedUtility},
    agent_store::AgentStore,
    w_iter,
    jump_model
)

    reg_year = model_data.year(first(model_data.index_y))
    reg_year_index = Symbol(Int(reg_year))
    delta_t = parse(Int64, chop(string(model_data.index_t.elements[2]), head = 1, tail = 0)) - parse(Int64, chop(string(model_data.index_t.elements[1]), head = 1, tail = 0))

    utility = get_agent(Utility, agent_store)
    customers = get_agent(CustomerGroup, agent_store)
    utility_opts = get_option(Utility, agent_store)

    # the x-axis of incentive function is assumed to be distributed storage capacity (even though the technology is PV+storage)
    incentive_function_dimension = size(der_aggregator.dera_stor_incentive_function)[1]
    incentive_function_slope = zeros(incentive_function_dimension - 1)
    incentive_function_intercept = zeros(incentive_function_dimension - 1)
    for i in 1:incentive_function_dimension - 1
        incentive_function_slope[i] = (der_aggregator.dera_stor_incentive_function[i+1, "participation"] - der_aggregator.dera_stor_incentive_function[i, "participation"]) /
        (der_aggregator.dera_stor_incentive_function[i+1, "incentive"] - der_aggregator.dera_stor_incentive_function[i, "incentive"])

        incentive_function_intercept[i] = der_aggregator.dera_stor_incentive_function[i+1, "participation"] - incentive_function_slope[i] * der_aggregator.dera_stor_incentive_function[i+1, "incentive"]
    end
    incentive_level_by_segment = Dict(z => zeros(incentive_function_dimension - 1) for z in model_data.index_z)
    participation_by_segment = Dict(z => zeros(incentive_function_dimension - 1) for z in model_data.index_z)
    obj_by_segment = Dict(z => zeros(incentive_function_dimension - 1) for z in model_data.index_z)

    total_der_stor_capacity = make_keyed_array(model_data.index_z, model_data.index_h)
    for z in model_data.index_z, h in model_data.index_h
        if w_iter >= 2
            total_der_stor_capacity(z, h, :) .=
                customers.x_DG_E_my(reg_year_index, h, z, :BTMStorage) + sum(
                    customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMStorage) for
                    y in model_data.year(first(model_data.index_y_fix)):(reg_year - 1)
                )
        else
            total_der_stor_capacity(z, h, :) .= customers.x_DG_E_my(reg_year_index, h, z, :BTMStorage)
        end
    end

    # construct revenue curves for VPP based on cost savings of DER aggregation

    diff_one, viu_obj_value = solve_agent_problem!(
                            utility,
                            utility_opts,
                            model_data,
                            hem_opts,
                            agent_store,
                            w_iter,
                            jump_model
    )



















    
    for z in model_data.index_z
        for i in 1:incentive_function_dimension - 1
            # x (incentive) should be $/MW (per year)?
            DERA_WM = get_new_jump_model(dera_opts.solvers)
            @variable(DERA_WM, der_aggregator.dera_stor_incentive_function[i, "incentive"] <= incentive <= der_aggregator.dera_stor_incentive_function[i+1, "incentive"])
            @variable(DERA_WM, dera_charge[model_data.index_d, model_data.index_t] >= 0)
            @variable(DERA_WM, dera_discharge[model_data.index_d, model_data.index_t] >= 0)
            @variable(DERA_WM, dera_energy[model_data.index_d, model_data.index_t] >= 0)
            @variable(DERA_WM, dera_pv_gen[model_data.index_d, model_data.index_t] >= 0)

            dera_stor_capacity_h =
                h -> begin

                (incentive_function_intercept[i] + incentive_function_slope[i] * incentive) * total_der_stor_capacity(z, h)

            end

            dera_stor_capacity =
                begin

                sum(dera_stor_capacity_h(h) for h in model_data.index_h)

            end

            # when it comes to the aggregated pv capacity, use aggregated storage capacity, divide by optimal storage size to get the number of households that's participating the dera aggregation,
            # then multiply by the optimal pv size,
            # TODO: I'm using optimal storage and pv capacities for existing resources, need to be more flexible to account for changing optimal tech sizes overtime.
            # right now, assume the optimal sizes are the same for both existing and new techs.
            dera_pv_capacity_h = 
                h -> begin

                dera_stor_capacity_h(h) / customers.Opti_DG_E(z, h, :BTMStorage) * customers.Opti_DG_E(z, h, :BTMPV)

            end

            dera_pv_capacity =
                begin

                sum(dera_pv_capacity_h(h) for h in model_data.index_h)

            end

            objective_revenue = begin
                sum(
                    model_data.omega(d) * delta_t *
                    ipp.LMP_my(reg_year_index, z, d, t) *
                    (dera_discharge[d, t] - dera_charge[d, t] + dera_pv_gen[d, t]) 
                    for d in model_data.index_d, t in model_data.index_t
                ) + ipp.capacity_price(reg_year_index) * (dera_stor_capacity .* ipp.capacity_credit_stor_E_my(reg_year_index, z, Symbol("der_aggregator")) + 
                dera_pv_capacity .* ipp.capacity_credit_E_my(reg_year_index, z, Symbol("dera_pv")))
            end

            objective_cost = begin
                incentive * dera_stor_capacity
            end

            @objective(DERA_WM, Max, objective_revenue - objective_cost)

            # DERA (storage) constraints -- We assume all distributed storage resources have the same initial energy, duration and round-trip-efficiency

            @constraint(
                DERA_WM,
                Eq_stor_energy_balance[
                    d in model_data.index_d,
                    t in model_data.index_t.elements[2:end],
                ],
                dera_energy[d, t] == dera_energy[d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)-delta_t), model_data.time.values)][1]] - dera_discharge[d, t] / customers.rte_dist_stor(z, :Commercial) * delta_t +
                    dera_charge[d, t] * delta_t
            )

            @constraint(
                DERA_WM,
                Eq_stor_energy_balance_0[
                    d in model_data.index_d,
                    t in [model_data.index_t.elements[1]],
                ],
                dera_energy[d, t] == customers.initial_energy_dist_stor(z, :Commercial, d) - dera_discharge[d, t] / customers.rte_dist_stor(z, :Commercial) * delta_t + dera_charge[d, t] * delta_t
            )

            @constraint(
                DERA_WM,
                Eq_stor_energy_upper_bound[
                    d in model_data.index_d,
                    t in model_data.index_t,
                ],
                dera_energy[d, t] <= customers.duration_dist_stor(z, :Commercial) * dera_stor_capacity
            )

            @constraint(
                DERA_WM,
                Eq_stor_discharge_upper_bound[
                    d in model_data.index_d,
                    t in model_data.index_t,
                ],
                dera_discharge[d, t] <= customers.rte_dist_stor(z, :Commercial) * dera_stor_capacity
            )

            @constraint(
                DERA_WM,
                Eq_stor_charge_upper_bound[
                    d in model_data.index_d,
                    t in model_data.index_t,
                ],
                dera_charge[d, t] <= dera_stor_capacity
            )

            @constraint(
                DERA_WM,
                Eq_stor_discharge_energy_upper_bound[
                    d in model_data.index_d,
                    t in model_data.index_t.elements[2:end],
                ],
                dera_discharge[d, t] * delta_t <= 
                dera_energy[d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)-delta_t), model_data.time.values)][1]]
            )

            @constraint(
                DERA_WM,
                Eq_stor_discharge_energy_upper_bound_0[
                    d in model_data.index_d,
                    t in [model_data.index_t.elements[1]],
                ],
                dera_discharge[d, t] * delta_t <= customers.initial_energy_dist_stor(z, :Commercial, d)
            )

            @constraint(
                DERA_WM,
                Eq_stor_charge_energy_upper_bound[
                    d in model_data.index_d,
                    t in model_data.index_t.elements[2:end],
                ],
                dera_charge[d, t] * delta_t <= customers.duration_dist_stor(z, :Commercial) * dera_stor_capacity -
                dera_energy[d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)-delta_t), model_data.time.values)][1]]
            )

            @constraint(
                DERA_WM,
                Eq_stor_charge_energy_upper_bound_0[
                    d in model_data.index_d,
                    t in [model_data.index_t.elements[1]],
                ],
                dera_charge[d, t] * delta_t <= customers.duration_dist_stor(z, :Commercial) * dera_stor_capacity -
                customers.initial_energy_dist_stor(z, :Commercial, d)
            )

            @constraint(
                DERA_WM,
                Eq_stor_charge_discharge_upper_bound[
                    d in model_data.index_d,
                    t in model_data.index_t.elements,
                ],
                dera_charge[d, t] + dera_discharge[d, t] / customers.rte_dist_stor(z, :Commercial) <= 
                dera_stor_capacity
            )

            @constraint(
                DERA_WM,
                Eq_pv_upper_bound[
                    d in model_data.index_d,
                    t in model_data.index_t,
                ],
                dera_pv_gen[d, t] <= ipp.rho_E_my(:ipp1, :dera_pv, z, d, t) * dera_pv_capacity
            )

            # TimerOutputs.@timeit HEM_TIMER "optimize! DER Aggregator dispatch" begin
                optimize!(DERA_WM)
            # end

            incentive_level_by_segment[z][i] = value(incentive)
            participation_by_segment[z][i] = incentive_function_intercept[i] + incentive_function_slope[i] * value(incentive)
            obj_by_segment[z][i] = objective_value(DERA_WM)

        end
    end

    max_seg_index = Dict(z => findfirst(x -> x == maximum(obj_by_segment[z]), obj_by_segment[z]) for z in model_data.index_z)
    for z in model_data.index_z
        der_aggregator.incentive_level(reg_year_index, z, :) .= incentive_level_by_segment[z][max_seg_index[z]]
        der_aggregator.aggregation_level(reg_year_index, z, :) .= participation_by_segment[z][max_seg_index[z]]
    end

    dera_agg_stor_capacity_h = make_keyed_array(model_data.index_z, model_data.index_h)
    dera_agg_pv_capacity_h = make_keyed_array(model_data.index_z, model_data.index_h)
    for z in model_data.index_z, h in model_data.index_h
        dera_agg_stor_capacity_h(z, h, :) .= der_aggregator.aggregation_level(reg_year_index, z) * total_der_stor_capacity(z, h)
        dera_agg_pv_capacity_h(z, h, :) .= dera_agg_stor_capacity_h(z, h) / customers.Opti_DG_E(z, h, :BTMStorage) * customers.Opti_DG_E(z, h, :BTMPV)
    end

    for z in model_data.index_z
        # simply assign DERA to a random ipp (ipp1)
        ipp.x_stor_E_my(:ipp1, z, Symbol("der_aggregator"), :) .= sum(dera_agg_stor_capacity_h(z, h) for h in model_data.index_h)
        ipp.x_E_my(:ipp1, z, Symbol("dera_pv"), :) .= sum(dera_agg_pv_capacity_h(z, h) for h in model_data.index_h)
    end

    # since we moved some BTM storage to transmission level, need to reduce the BTM net load accordingly (in bulk power system, regulator, customers (maybe?)).

    return 0.0
    
end