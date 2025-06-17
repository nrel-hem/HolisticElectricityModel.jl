function get_regulator_options(config::Dict{Any,Any})
    rate_design, net_metering_policy, tou_suffix, planning_reserve_margin,
    allowed_return_on_investment = parse(config, "Regulator", validators)

    return RegulatorOptions(
        rate_design,
        net_metering_policy;
        tou_suffix=tou_suffix,
        planning_reserve_margin=planning_reserve_margin,
        allowed_return_on_investment=allowed_return_on_investment,
    )
end

function get_ipp_options(config::Dict{Any,Any}, solver::Symbol)
    # Get the optimizer depending on the solver defined the config
    ipp_algorithm, = parse(config, "IPPGroup", validators)

    ipp_solvers = Dict()
    addsolvers_ipp!(ipp_solvers, :Ipopt)
    addsolvers_ipp!(ipp_solvers, solver)
    return IPPOptions(ipp_algorithm, ipp_solvers)
end

function get_utility_options(solver::Symbol)
    return UtilityOptions(JuMP.optimizer_with_attributes(
        () -> get_optimizer_for_solver(solver),
        # "OUTPUTLOG" => 0,
    ))
end

function get_customer_options(config::Dict{Any,Any}, solver::Symbol)
    pv_adoption_type, = parse(config, "CustomerGroup", validators)

    return CustomerOptions(
        pv_adoption_type,
        JuMP.optimizer_with_attributes(
            () -> get_optimizer_for_solver(solver),
            # "OUTPUTLOG" => 0,
        ),
    )
end

function get_green_developer_options(solver::Symbol)
    return GreenDeveloperOptions(
        JuMP.optimizer_with_attributes(
            () -> get_optimizer_for_solver(solver)
            # "OUTPUTLOG" => 0,
        ),
    )
end

function get_der_aggregator_options(config::Dict{Any,Any}, solver::Symbol)
    incentive_curve, frac_viu_cost_savings_as_revenue =
        parse(config, "DERAggregator", validators)

    return DERAggregatorOptions(
        JuMP.optimizer_with_attributes(
            () -> get_optimizer_for_solver(solver),
            # "OUTPUTLOG" => 0,
        );
        incentive_curve=incentive_curve,
        frac_viu_cost_savings_as_revenue=frac_viu_cost_savings_as_revenue
    )
end

function get_agent_options(config::Dict{Any,Any}, ::HEMOptions{VIU}, solver::Symbol)
    return AgentOptionsStore(
        Dict(
            Regulator => get_regulator_options(config),
            Utility => get_utility_options(solver),
            CustomerGroup => get_customer_options(config, solver),
            GreenDeveloper => get_green_developer_options(solver),
            DERAggregator => get_der_aggregator_options(config, solver),
            # DistributionUtility => NullAgentOptions()
        )
    )
end

function get_agent_options(config::Dict{Any,Any}, ::HEMOptions{WM}, solver::Symbol)
    return AgentOptionsStore(
        Dict(
            Regulator => get_regulator_options(config),
            IPPGroup => get_ipp_options(config, solver),
            CustomerGroup => get_customer_options(config, solver),
            GreenDeveloper => get_green_developer_options(solver),
            DERAggregator => get_der_aggregator_options(config, solver),
            # DistributionUtility => NullAgentOptions()
        )
    )
end