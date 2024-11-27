abstract type AbstractDERAggregator <: Agent end

struct DERAggregatorOptions <: AgentOptions
    solvers::HEMSolver
    
    incentive_curve::Int
    frac_viu_cost_savings_as_revenue::AbstractFloat
end

"""
Construct DERAggregatorOptions with an MOI.OptimizerWithAttributes instance.
"""
function DERAggregatorOptions(
    attributes::MOI.OptimizerWithAttributes; 
    incentive_curve::Int = 1, 
    frac_viu_cost_savings_as_revenue::AbstractFloat = 0.1
)
    @assert (incentive_curve >= 1) && (incentive_curve <= 5) "dera_stor_incenctive_function_{x} is available for 1 <= x <= 5, was passed $(incentive_curve)"
    @assert (frac_viu_cost_savings_as_revenue > 0.0) && (frac_viu_cost_savings_as_revenue <= 1.0) "Expected fraction > 0 and <= 1, got $(frac_viu_cost_savings_as_revenue)"
    return DERAggregatorOptions(AnySolver(attributes), incentive_curve, frac_viu_cost_savings_as_revenue)
end

function get_file_prefix(options::DERAggregatorOptions)
    return join(["IncentiveCurve$(options.incentive_curve)",
                 "VIURevenueFrac$(options.frac_viu_cost_savings_as_revenue)"], "_")
end

mutable struct DERAggregator <: AbstractDERAggregator
    id::String
    current_year::Symbol

    # incentive function for DER aggregation (piece-wise linear function: incentive (x) vs participation (y))
    dera_stor_incentive_function::DataFrame
    # parameters applied to aggregated capacities of DER to represent potential friction of aggregation
    aggregation_friction::ParamArray
    # incentive level by year
    incentive_level::ParamArray
    # DER aggregator revenue by year (under VIU)
    revenue::ParamArray
    # aggregation level by year
    aggregation_level::ParamArray
    # aggregation level by year (used for output)
    aggregation_level_output::ParamArray
    # percentage of cost savings as revenues for DERAggregator under VIU
    rev_perc_cost_saving_viu::ParamArray
    # aggregated distributed storage capacity (MW)
    dera_stor_my::ParamArray
    # aggregated distributed pv capacity associated with pv+storage (MW)
    dera_pv_my::ParamArray
end

function DERAggregator(input_filename::AbstractString, model_data::HEMData, opts::DERAggregatorOptions; id = DEFAULT_ID)

    # need to have the incentive function for each customer type
    dera_stor_incentive_function = CSV.read(joinpath(input_filename, "dera_stor_incentive_function_$(opts.incentive_curve).csv"), DataFrame)

    return DERAggregator(
        id,
        first(model_data.index_y),
        dera_stor_incentive_function,
        initialize_param("aggregation_friction", model_data.index_y, model_data.index_z, model_data.index_h, value = 0.0),
        initialize_param("incentive_level", model_data.index_y, model_data.index_z, value = 0.0),
        initialize_param("revenue", model_data.index_y, model_data.index_z, value = 0.0),
        initialize_param("aggregation_level", model_data.index_y, model_data.index_z, value = 0.0),
        initialize_param("aggregation_level_output", model_data.index_y, model_data.index_z, value = 0.0),
        initialize_param("rev_perc_cost_saving_viu", model_data.index_y, model_data.index_z, value = opts.frac_viu_cost_savings_as_revenue),
        initialize_param("dera_stor_my", model_data.index_y, model_data.index_z, model_data.index_h, value = 0.0),
        initialize_param("dera_pv_my", model_data.index_y, model_data.index_z, model_data.index_h, value = 0.0),
    )
end

get_id(x::DERAggregator) = x.id

function solve_agent_problem!(
    der_aggregator::DERAggregator,
    dera_opts::DERAggregatorOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{WM, <:UseCase, <:UseCase, NullUseCase},
    agent_store::AgentStore,
    w_iter,
    window_length,
    jump_model,
    export_file_path,
    update_results::Bool
)

    ipp = get_agent(IPPGroup, agent_store)
    reg_year, reg_year_index = get_reg_year(model_data)

    for z in model_data.index_z
        der_aggregator.incentive_level(reg_year_index, z, :) .= 0.0
        der_aggregator.aggregation_level(reg_year_index, z, :) .= 0.0
    end

    for z in model_data.index_z, h in model_data.index_h
        der_aggregator.dera_stor_my(reg_year_index, z, h, :) .= 0.0
        der_aggregator.dera_pv_my(reg_year_index, z, h, :) .= 0.0
    end

    for z in model_data.index_z
        # simply assign DERA to a random ipp (ipp1)
        ipp.x_stor_E_my(:ipp1, z, Symbol("der_aggregator"), :) .= 0.0
        ipp.x_E_my(:ipp1, z, Symbol("dera_pv"), :) .= 0.0
    end

    der_aggregator.current_year = reg_year_index

    # since we moved some BTM storage to transmission level, need to reduce the BTM net load accordingly (in bulk power system, regulator, customers (maybe?)).

    return 0.0

end

function solve_agent_problem!(
    der_aggregator::DERAggregator,
    dera_opts::AgentOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{VIU, <:UseCase, <:UseCase, NullUseCase},
    agent_store::AgentStore,
    w_iter,
    window_length,
    jump_model,
    export_file_path,
    update_results::Bool
)

    utility = get_agent(Utility, agent_store)
    reg_year, reg_year_index = get_reg_year(model_data)
    
    for z in model_data.index_z
        der_aggregator.incentive_level(reg_year_index, z, :) .= 0.0
        der_aggregator.aggregation_level(reg_year_index, z, :) .= 0.0
    end

    for z in model_data.index_z, h in model_data.index_h
        der_aggregator.dera_stor_my(reg_year_index, z, h, :) .= 0.0
        der_aggregator.dera_pv_my(reg_year_index, z, h, :) .= 0.0
    end

    for z in model_data.index_z
        utility.x_stor_E_my(z, Symbol("der_aggregator"), :) .= 0.0
        utility.x_E_my(z, Symbol("dera_pv"), :) .= 0.0
    end

    der_aggregator.current_year = reg_year_index

    # since we moved some BTM storage to transmission level, need to reduce the BTM net load accordingly (in bulk power system, regulator, customers (maybe?)).

    return 0.0

end

function solve_agent_problem!(
    der_aggregator::DERAggregator,
    dera_opts::AgentOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{WM, <:UseCase, <:UseCase, DERAggregation},
    agent_store::AgentStore,
    w_iter,
    window_length,
    jump_model,
    export_file_path,
    update_results::Bool
)

    reg_year, reg_year_index = get_reg_year(model_data)
    delta_t = get_delta_t(model_data)

    ipp = get_agent(IPPGroup, agent_store)
    customers = get_agent(CustomerGroup, agent_store)

    # the x-axis of incentive function is assumed to be dollars per distributed storage capacity (even though the technology is PV+storage)
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
                    y in model_data.year(first(model_data.index_y_fix)):(reg_year)
                )
        else
            total_der_stor_capacity(z, h, :) .= customers.x_DG_E_my(reg_year_index, h, z, :BTMStorage)
        end
    end
    
    for z in model_data.index_z
        for i in 1:incentive_function_dimension - 1
            # x (incentive) should be $/MW (per year)?
            DERAggregator_WM = get_new_jump_model(dera_opts.solvers)
            @variable(DERAggregator_WM, der_aggregator.dera_stor_incentive_function[i, "incentive"] <= incentive <= der_aggregator.dera_stor_incentive_function[i+1, "incentive"])
            @variable(DERAggregator_WM, dera_charge[model_data.index_d, model_data.index_t] >= 0)
            @variable(DERAggregator_WM, dera_discharge[model_data.index_d, model_data.index_t] >= 0)
            @variable(DERAggregator_WM, dera_energy[model_data.index_d, model_data.index_t] >= 0)
            @variable(DERAggregator_WM, dera_pv_gen[model_data.index_d, model_data.index_t] >= 0)

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

            @objective(DERAggregator_WM, Max, objective_revenue - objective_cost)

            # DERAggregator (storage) constraints -- We assume all distributed storage resources have the same initial energy, duration and round-trip-efficiency

            @constraint(
                DERAggregator_WM,
                Eq_stor_energy_balance[
                    d in model_data.index_d,
                    t in model_data.index_t.elements[2:end],
                ],
                dera_energy[d, t] == dera_energy[d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)-delta_t), model_data.time.values)][1]] - dera_discharge[d, t] / customers.rte_dist_stor(z, :Commercial) * delta_t +
                    dera_charge[d, t] * delta_t
            )

            @constraint(
                DERAggregator_WM,
                Eq_stor_energy_balance_0[
                    d in model_data.index_d,
                    t in [model_data.index_t.elements[1]],
                ],
                dera_energy[d, t] == customers.initial_energy_dist_stor(z, :Commercial, d) - dera_discharge[d, t] / customers.rte_dist_stor(z, :Commercial) * delta_t + dera_charge[d, t] * delta_t
            )

            @constraint(
                DERAggregator_WM,
                Eq_stor_energy_upper_bound[
                    d in model_data.index_d,
                    t in model_data.index_t,
                ],
                dera_energy[d, t] <= customers.duration_dist_stor(z, :Commercial) * dera_stor_capacity
            )

            @constraint(
                DERAggregator_WM,
                Eq_stor_discharge_upper_bound[
                    d in model_data.index_d,
                    t in model_data.index_t,
                ],
                dera_discharge[d, t] <= customers.rte_dist_stor(z, :Commercial) * dera_stor_capacity
            )

            @constraint(
                DERAggregator_WM,
                Eq_stor_charge_upper_bound[
                    d in model_data.index_d,
                    t in model_data.index_t,
                ],
                dera_charge[d, t] <= dera_stor_capacity
            )

            @constraint(
                DERAggregator_WM,
                Eq_stor_discharge_energy_upper_bound[
                    d in model_data.index_d,
                    t in model_data.index_t.elements[2:end],
                ],
                dera_discharge[d, t] * delta_t <= 
                dera_energy[d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)-delta_t), model_data.time.values)][1]]
            )

            @constraint(
                DERAggregator_WM,
                Eq_stor_discharge_energy_upper_bound_0[
                    d in model_data.index_d,
                    t in [model_data.index_t.elements[1]],
                ],
                dera_discharge[d, t] * delta_t <= customers.initial_energy_dist_stor(z, :Commercial, d)
            )

            @constraint(
                DERAggregator_WM,
                Eq_stor_charge_energy_upper_bound[
                    d in model_data.index_d,
                    t in model_data.index_t.elements[2:end],
                ],
                dera_charge[d, t] * delta_t <= customers.duration_dist_stor(z, :Commercial) * dera_stor_capacity -
                dera_energy[d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)-delta_t), model_data.time.values)][1]]
            )

            @constraint(
                DERAggregator_WM,
                Eq_stor_charge_energy_upper_bound_0[
                    d in model_data.index_d,
                    t in [model_data.index_t.elements[1]],
                ],
                dera_charge[d, t] * delta_t <= customers.duration_dist_stor(z, :Commercial) * dera_stor_capacity -
                customers.initial_energy_dist_stor(z, :Commercial, d)
            )

            @constraint(
                DERAggregator_WM,
                Eq_stor_charge_discharge_upper_bound[
                    d in model_data.index_d,
                    t in model_data.index_t.elements,
                ],
                dera_charge[d, t] + dera_discharge[d, t] / customers.rte_dist_stor(z, :Commercial) <= 
                dera_stor_capacity
            )

            @constraint(
                DERAggregator_WM,
                Eq_pv_upper_bound[
                    d in model_data.index_d,
                    t in model_data.index_t,
                ],
                dera_pv_gen[d, t] <= ipp.rho_E_my(:ipp1, :dera_pv, z, d, t) * dera_pv_capacity
            )

            # TimerOutputs.@timeit HEM_TIMER "optimize! DER Aggregator dispatch" begin
                optimize!(DERAggregator_WM)
            # end

            incentive_level_by_segment[z][i] = value(incentive)
            participation_by_segment[z][i] = incentive_function_intercept[i] + incentive_function_slope[i] * value(incentive)
            obj_by_segment[z][i] = objective_value(DERAggregator_WM)

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
        der_aggregator.dera_stor_my(reg_year_index, z, h, :) .= dera_agg_stor_capacity_h(z, h)
        der_aggregator.dera_pv_my(reg_year_index, z, h, :) .= dera_agg_pv_capacity_h(z, h)
    end

    # TODO: I think this might also be happening in the IPPGroup right now. Choose one location?
    for z in model_data.index_z
        # simply assign DERAggregator to a random ipp (ipp1)
        ipp.x_stor_E_my(:ipp1, z, Symbol("der_aggregator"), :) .= sum(dera_agg_stor_capacity_h(z, h) for h in model_data.index_h)
        ipp.x_E_my(:ipp1, z, Symbol("dera_pv"), :) .= sum(dera_agg_pv_capacity_h(z, h) for h in model_data.index_h)
    end

    der_aggregator.current_year = reg_year_index

    # since we moved some BTM storage to transmission level, need to reduce the BTM net load accordingly (in bulk power system, regulator, customers (maybe?)).

    return 0.0    
end

function solve_agent_problem!(
    der_aggregator::DERAggregator,
    dera_opts::DERAggregatorOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{VIU, <:UseCase, <:UseCase, DERAggregation},
    agent_store::AgentStore,
    w_iter,
    window_length,
    jump_model,
    export_file_path,
    update_results::Bool
)

    reg_year, reg_year_index = get_reg_year(model_data)
    reg_year_dera, reg_year_index_dera = get_prev_reg_year(model_data, w_iter)
    delta_t = get_delta_t(model_data)

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
    obj_revenue_by_segment = Dict(z => zeros(incentive_function_dimension - 1) for z in model_data.index_z)

    cem_cost_saving_function_one = copy(der_aggregator.dera_stor_incentive_function)
    rename!(cem_cost_saving_function_one, :incentive => Symbol("cost_savings"))
    cem_cost_saving_function_one[!, "cost_savings"] = convert.(Float64, cem_cost_saving_function_one[!, "cost_savings"])
    cem_cost_saving_function = Dict(z => copy(cem_cost_saving_function_one) for z in model_data.index_z)

    total_der_stor_capacity = make_keyed_array(model_data.index_z, model_data.index_h)
    for z in model_data.index_z, h in model_data.index_h
        if w_iter >= 2
            total_der_stor_capacity(z, h, :) .=
                customers.x_DG_E_my(reg_year_index, h, z, :BTMStorage) + sum(
                    customers.x_DG_new_my(Symbol(Int(y)), h, z, :BTMStorage) for
                    y in model_data.year(first(model_data.index_y_fix)):(reg_year)
                )
        else
            total_der_stor_capacity(z, h, :) .= customers.x_DG_E_my(reg_year_index, h, z, :BTMStorage)
        end
    end

    # construct revenue curves for VPP based on cost savings of DER aggregation
    for z in model_data.index_z
        # simply assign DERAggregator to a random ipp (ipp1)
        utility.x_stor_E_my(z, Symbol("der_aggregator"), :) .= 0.0
        utility.x_E_my(z, Symbol("dera_pv"), :) .= 0.0
        der_aggregator.aggregation_level(reg_year_index_dera, z, :) .= 0.0
    end

    diff_one_base = solve_agent_problem!(
        utility,
        utility_opts,
        model_data,
        hem_opts,
        agent_store,
        w_iter,
        window_length,
        jump_model,
        export_file_path,
        false
    )
    viu_obj_value_base = deepcopy(utility._obj_value)

    for z in model_data.index_z
        if sum(total_der_stor_capacity(z, h) for h in model_data.index_h) != 0.0
            for i in 1:incentive_function_dimension - 1
                # set der aggregation level to points on the curve before running CEM
                der_aggregator.aggregation_level(reg_year_index_dera, z, :) .= der_aggregator.dera_stor_incentive_function[i+1, "participation"]
                utility.x_stor_E_my(z, Symbol("der_aggregator"), :) .= der_aggregator.dera_stor_incentive_function[i+1, "participation"] * sum(total_der_stor_capacity(z, h) for h in model_data.index_h)
                utility.x_E_my(z, Symbol("dera_pv"), :) .= der_aggregator.dera_stor_incentive_function[i+1, "participation"] * 
                sum(
                    customers.Opti_DG_E(z, h, :BTMStorage) == 0 || customers.Opti_DG_E(z, h, :BTMPV) == 0 ? 0.0 : 
                    total_der_stor_capacity(z, h) / customers.Opti_DG_E(z, h, :BTMStorage) * customers.Opti_DG_E(z, h, :BTMPV)
                    for h in model_data.index_h
                )
            

                diff_one = solve_agent_problem!(
                    utility,
                    utility_opts,
                    model_data,
                    hem_opts,
                    agent_store,
                    w_iter,
                    window_length,
                    jump_model,
                    export_file_path,
                    false
                )
                viu_obj_value = deepcopy(utility._obj_value)

                cem_cost_saving_function[z][i+1, "cost_savings"] = max(0.0, viu_obj_value_base - viu_obj_value)
            end
        else
            cem_cost_saving_function[z][!, "cost_savings"] .= 0.0
        end
        CSV.write(joinpath(export_file_path, "cem_cost_saving_function_$(reg_year)_$(z).csv"), cem_cost_saving_function[z])
    end

    cem_cost_saving_function_slope = Dict(z => zeros(incentive_function_dimension - 1) for z in model_data.index_z)
    cem_cost_saving_function_intercept = Dict(z => zeros(incentive_function_dimension - 1) for z in model_data.index_z)

    for z in model_data.index_z
        for i in 1:incentive_function_dimension - 1
            if (cem_cost_saving_function[z][i+1, "cost_savings"] - cem_cost_saving_function[z][i, "cost_savings"]) < 0.0
                @info "cost saving function is decreasing at some points."
            else
                if (cem_cost_saving_function[z][i+1, "participation"] - cem_cost_saving_function[z][i, "participation"]) != 0.0
                    cem_cost_saving_function_slope[z][i] = (cem_cost_saving_function[z][i+1, "cost_savings"] - cem_cost_saving_function[z][i, "cost_savings"]) / 
                    (cem_cost_saving_function[z][i+1, "participation"] - cem_cost_saving_function[z][i, "participation"])
                else
                    cem_cost_saving_function_slope[z][i] = 0.0
                end
                cem_cost_saving_function_intercept[z][i] = cem_cost_saving_function[z][i+1, "cost_savings"] - cem_cost_saving_function_slope[z][i] * cem_cost_saving_function[z][i+1, "participation"]
            end
        end
    end

    for z in model_data.index_z
        for i in 1:incentive_function_dimension - 1
            # x (incentive) should be $/MW (per year)?
            DERAggregator_VIU = get_new_jump_model(dera_opts.solvers)
            @variable(DERAggregator_VIU, der_aggregator.dera_stor_incentive_function[i, "incentive"] <= incentive <= der_aggregator.dera_stor_incentive_function[i+1, "incentive"])
            @variable(DERAggregator_VIU, dera_charge[model_data.index_d, model_data.index_t] >= 0)
            @variable(DERAggregator_VIU, dera_discharge[model_data.index_d, model_data.index_t] >= 0)
            @variable(DERAggregator_VIU, dera_energy[model_data.index_d, model_data.index_t] >= 0)
            @variable(DERAggregator_VIU, dera_pv_gen[model_data.index_d, model_data.index_t] >= 0)

            participation_perc = begin
                incentive_function_intercept[i] + incentive_function_slope[i] * incentive
            end

            dera_stor_capacity_h =
                h -> begin
                participation_perc * total_der_stor_capacity(z, h)
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
                der_aggregator.rev_perc_cost_saving_viu(reg_year_index, z) * (cem_cost_saving_function_intercept[z][i] + cem_cost_saving_function_slope[z][i] * participation_perc)
            end

            objective_cost = begin
                incentive * dera_stor_capacity
            end

            @objective(DERAggregator_VIU, Max, objective_revenue - objective_cost)

            # DERAggregator (storage) constraints -- We assume all distributed storage resources have the same initial energy, duration and round-trip-efficiency

            @constraint(
                DERAggregator_VIU,
                Eq_stor_energy_balance[
                    d in model_data.index_d,
                    t in model_data.index_t.elements[2:end],
                ],
                dera_energy[d, t] == dera_energy[d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)-delta_t), model_data.time.values)][1]] - dera_discharge[d, t] / customers.rte_dist_stor(z, :Commercial) * delta_t +
                    dera_charge[d, t] * delta_t
            )

            @constraint(
                DERAggregator_VIU,
                Eq_stor_energy_balance_0[
                    d in model_data.index_d,
                    t in [model_data.index_t.elements[1]],
                ],
                dera_energy[d, t] == customers.initial_energy_dist_stor(z, :Commercial, d) - dera_discharge[d, t] / customers.rte_dist_stor(z, :Commercial) * delta_t + dera_charge[d, t] * delta_t
            )

            @constraint(
                DERAggregator_VIU,
                Eq_stor_energy_upper_bound[
                    d in model_data.index_d,
                    t in model_data.index_t,
                ],
                dera_energy[d, t] <= customers.duration_dist_stor(z, :Commercial) * dera_stor_capacity
            )

            @constraint(
                DERAggregator_VIU,
                Eq_stor_discharge_upper_bound[
                    d in model_data.index_d,
                    t in model_data.index_t,
                ],
                dera_discharge[d, t] <= customers.rte_dist_stor(z, :Commercial) * dera_stor_capacity
            )

            @constraint(
                DERAggregator_VIU,
                Eq_stor_charge_upper_bound[
                    d in model_data.index_d,
                    t in model_data.index_t,
                ],
                dera_charge[d, t] <= dera_stor_capacity
            )

            @constraint(
                DERAggregator_VIU,
                Eq_stor_discharge_energy_upper_bound[
                    d in model_data.index_d,
                    t in model_data.index_t.elements[2:end],
                ],
                dera_discharge[d, t] * delta_t <= 
                dera_energy[d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)-delta_t), model_data.time.values)][1]]
            )

            @constraint(
                DERAggregator_VIU,
                Eq_stor_discharge_energy_upper_bound_0[
                    d in model_data.index_d,
                    t in [model_data.index_t.elements[1]],
                ],
                dera_discharge[d, t] * delta_t <= customers.initial_energy_dist_stor(z, :Commercial, d)
            )

            @constraint(
                DERAggregator_VIU,
                Eq_stor_charge_energy_upper_bound[
                    d in model_data.index_d,
                    t in model_data.index_t.elements[2:end],
                ],
                dera_charge[d, t] * delta_t <= customers.duration_dist_stor(z, :Commercial) * dera_stor_capacity -
                dera_energy[d, model_data.index_t.elements[findall(x -> x == (model_data.time(t)-delta_t), model_data.time.values)][1]]
            )

            @constraint(
                DERAggregator_VIU,
                Eq_stor_charge_energy_upper_bound_0[
                    d in model_data.index_d,
                    t in [model_data.index_t.elements[1]],
                ],
                dera_charge[d, t] * delta_t <= customers.duration_dist_stor(z, :Commercial) * dera_stor_capacity -
                customers.initial_energy_dist_stor(z, :Commercial, d)
            )

            @constraint(
                DERAggregator_VIU,
                Eq_stor_charge_discharge_upper_bound[
                    d in model_data.index_d,
                    t in model_data.index_t.elements,
                ],
                dera_charge[d, t] + dera_discharge[d, t] / customers.rte_dist_stor(z, :Commercial) <= 
                dera_stor_capacity
            )

            @constraint(
                DERAggregator_VIU,
                Eq_pv_upper_bound[
                    d in model_data.index_d,
                    t in model_data.index_t,
                ],
                dera_pv_gen[d, t] <= utility.rho_E_my(:dera_pv, z, d, t) * dera_pv_capacity
            )

            # TimerOutputs.@timeit HEM_TIMER "optimize! DER Aggregator dispatch" begin
                optimize!(DERAggregator_VIU)
            # end

            incentive_level_by_segment[z][i] = value(incentive)
            participation_by_segment[z][i] = incentive_function_intercept[i] + incentive_function_slope[i] * value(incentive)
            obj_by_segment[z][i] = objective_value(DERAggregator_VIU)
            obj_revenue_by_segment[z][i] = der_aggregator.rev_perc_cost_saving_viu(reg_year_index, z) * (cem_cost_saving_function_intercept[z][i] + cem_cost_saving_function_slope[z][i] * (incentive_function_intercept[i] + incentive_function_slope[i] * value(incentive)))

        end
    end

    max_seg_index = Dict(z => findfirst(x -> x == maximum(obj_by_segment[z]), obj_by_segment[z]) for z in model_data.index_z)
    for z in model_data.index_z
        der_aggregator.incentive_level(reg_year_index, z, :) .= incentive_level_by_segment[z][max_seg_index[z]]
        der_aggregator.aggregation_level(reg_year_index, z, :) .= participation_by_segment[z][max_seg_index[z]]
        der_aggregator.aggregation_level_output(reg_year_index, z, :) .= participation_by_segment[z][max_seg_index[z]]
        der_aggregator.revenue(reg_year_index, z, :) .= obj_revenue_by_segment[z][max_seg_index[z]]
    end

    dera_agg_stor_capacity_h = make_keyed_array(model_data.index_z, model_data.index_h)
    dera_agg_pv_capacity_h = make_keyed_array(model_data.index_z, model_data.index_h)
    for z in model_data.index_z, h in model_data.index_h
        dera_agg_stor_capacity_h(z, h, :) .= der_aggregator.aggregation_level(reg_year_index, z) * total_der_stor_capacity(z, h)
        dera_agg_pv_capacity_h(z, h, :) .= customers.Opti_DG_E(z, h, :BTMStorage) == 0 || customers.Opti_DG_E(z, h, :BTMPV) == 0 ? 
        0.0 : dera_agg_stor_capacity_h(z, h) / customers.Opti_DG_E(z, h, :BTMStorage) * customers.Opti_DG_E(z, h, :BTMPV)   
        der_aggregator.dera_stor_my(reg_year_index, z, h, :) .= dera_agg_stor_capacity_h(z, h)
        der_aggregator.dera_pv_my(reg_year_index, z, h, :) .= dera_agg_pv_capacity_h(z, h)
    end

    for z in model_data.index_z
        utility.x_stor_E_my(z, Symbol("der_aggregator"), :) .= sum(dera_agg_stor_capacity_h(z, h) for h in model_data.index_h)
        utility.x_E_my(z, Symbol("dera_pv"), :) .= sum(dera_agg_pv_capacity_h(z, h) for h in model_data.index_h)
    end

    der_aggregator.current_year = reg_year_index

    # since we moved some BTM storage to transmission level, need to reduce the BTM net load accordingly (in bulk power system, regulator, customers (maybe?)).

    return 0.0
    
end

function save_results(
    der_aggregator::DERAggregator,
    dera_opts::AgentOptions,
    hem_opts::HEMOptions{<:MarketStructure},
    export_file_path::AbstractString,
)
    # Primal Variables
    save_param(
        der_aggregator.incentive_level.values,
        [:Year, :Zone],
        :incentive_dollar_per_MW,
        joinpath(export_file_path, "dera_incentive.csv"),
    )
    save_param(
        der_aggregator.aggregation_level_output.values,
        [:Year, :Zone],
        :aggregation_perc,
        joinpath(export_file_path, "dera_aggregation_perc.csv"),
    )
    save_param(
        der_aggregator.dera_stor_my.values,
        [:Year, :Zone, :CustomerType],
        :aggregation_stor_mw,
        joinpath(export_file_path, "dera_aggregation_stor_mw.csv"),
    )
    save_param(
        der_aggregator.dera_pv_my.values,
        [:Year, :Zone, :CustomerType],
        :aggregation_pv_mw,
        joinpath(export_file_path, "dera_aggregation_pv_mw.csv"),
    )
end
