function get_regulator_options(config::Dict{Any,Any})
    rate_design, net_metering_policy, tou_suffix, planning_reserve_margin,
    allowed_return_on_investment = parse(config, "regulator_options", validators)

    return RegulatorOptions(
        rate_design,
        net_metering_policy;
        tou_suffix=tou_suffix,
        planning_reserve_margin=planning_reserve_margin,
        allowed_return_on_investment=allowed_return_on_investment,
    )
end

function get_ipp_options(config::Dict{Any,Any})
    # Get the optimizer depending on the solver defined the config
    ipp_algorithm, = parse(config, "ipp_options", validators)

    ipp_solvers = Dict()
    addsolvers_ipp!(ipp_solvers, :Ipopt)
    addsolvers_ipp!(ipp_solvers, solver)
    return IPPOptions(ipp_algorithm, ipp_solvers)
end

function get_utility_options()
    return UtilityOptions(JuMP.optimizer_with_attributes(
        () -> get_optimizer_for_solver(solver),
        # "OUTPUTLOG" => 0,
    ))
end

function get_customer_options(config::Dict{Any,Any})
    pv_adoption_type, = parse(config, "customer_options", validators)

    return CustomerOptions(
        pv_adoption_type,
        JuMP.optimizer_with_attributes(
            () -> get_optimizer_for_solver(solver),
            # "OUTPUTLOG" => 0,
        ),
    )
end

function get_green_developer_options()
    return GreenDeveloperOptions(
        JuMP.optimizer_with_attributes(
            () -> get_optimizer_for_solver(solver)
            # "OUTPUTLOG" => 0,
        ),
    )
end

function get_der_aggregator_options(config::Dict{Any,Any})
    incentive_curve, frac_viu_cost_savings_as_revenue =
        parse(config, "der_aggregator_options", validators)

    return DERAggregatorOptions(
        JuMP.optimizer_with_attributes(
            () -> get_optimizer_for_solver(solver),
            # "OUTPUTLOG" => 0,
        );
        incentive_curve=incentive_curve,
        frac_viu_cost_savings_as_revenue=frac_viu_cost_savings_as_revenue
    )
end

function get_agent_options(config::Dict{Any,Any}, ::HEMOptions{VIU})
    return Dict(
        "regulator_options" => get_regulator_options(config),
        "utility_options" => get_utility_options(),
        "customer_options" => get_customer_options(config),
        "green_developer_options" => get_green_developer_options(),
        "der_aggregator_options" => get_der_aggregator_options(config)
    )
end

function get_agent_options(config::Dict{Any,Any}, ::HEMOptions{WM})
    return Dict(
        "regulator_options" => get_regulator_options(config),
        "ipp_options" => get_ipp_options(config),
        "customer_options" => get_customer_options(config),
        "green_developer_options" => get_green_developer_options(),
        "der_aggregator_options" => get_der_aggregator_options(config)
    )
end