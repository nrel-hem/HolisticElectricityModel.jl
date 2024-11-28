# include config.jl before this file

reeds_bas = collect("p$n" for n = 1:134)

null_use_case_identifier = "null_use_case"

market_structure_map = Dict(
    "wholesale_market" => WM(),
    "vertically_integrated_utility" => VIU()
)

der_use_case_map = Dict(
    "der_use_case" => DERAdoption(),
    null_use_case_identifier => NullUseCase()
)

supply_choice_use_case_map = Dict(
    "supply_choice_use_case" => SupplyChoice(),
    null_use_case_identifier => NullUseCase()
)

der_aggregation_use_case_map = Dict(
    "der_aggregation_use_case" => DERAggregation(),
    null_use_case_identifier => NullUseCase()
)

rate_design_map = Dict(
    "flat_rate" => FlatRate(),
    "time_of_use" => TOU()
)

net_metering_policy_map = Dict(
    "excess_retail_rate" => ExcessRetailRate(),
    "excess_marginal_cost" => ExcessMarginalCost(),
    "excess_zero" => ExcessZero()
)

ipp_algorithm_map = Dict(
    "lagrange_decomposition" => LagrangeDecomposition(),
    "mppdcmer_transportation_storage" => MPPDCMERTransStorage(),
    "mppdcmer" => MPPDCMER(),
    "miqp" => MIQP()
)

pv_adoption_type_map = Dict(
    "standalone_pv" => StandalonePVOnly(),
    "solar_plus_storage" => SolarPlusStorageOnly(),
    "compete_der_configs" => CompeteDERConfigs()
)

validators = Dict(
    "data_selection" => [
        FieldValidatorBasic(
            "input_path",
            value -> check_path(value)
        ),
    ],
    "simulation_parameters" => [
        FieldValidatorBasic(
            "solver",
            value -> check_chain(value, [
                val -> check_in_collection(val, ("Gurobi", "Xpress")),
                check_symbol
            ])
        ),
    ],
    "hem_options" => [
        FieldValidatorBasic(
            "market_structure",
            value -> check_and_return_from_map(value, market_structure_map)
        ),
        FieldValidatorBasic(
            "der_use_case",
            value -> check_and_return_from_map(value, der_use_case_map)
        ),
        FieldValidatorBasic(
            "supply_choice_use_case",
            value -> check_and_return_from_map(value, supply_choice_use_case_map)
        ),
        FieldValidatorBasic(
            "der_aggregation_use_case",
            value -> check_and_return_from_map(value, der_aggregation_use_case_map)
        ),
    ],
    "regulator_options" => [
        FieldValidatorBasic(
            "rate_design",
            value -> check_and_return_from_map(value, rate_design_map)
        ),
        FieldValidatorBasic(
            "net_metering_policy",
            value -> check_and_return_from_map(value, net_metering_policy_map)
        ),
        FieldValidatorHasDefault(
            "tou_suffix",
            value -> check_in_collection(value, ("NE2025", "NE2035")),
            "NE2025"
        ),
        FieldValidatorHasDefault(
            "planning_reserve_margin",
            value -> check_float(value; min=0.0, max=0.5),
            0.129
        ),
        FieldValidatorHasDefault(
            "allowed_return_on_investment",
            value -> check_float(value; min=0.0, max=0.5),
            0.112
        )
    ],
    "ipp_options" => [
        FieldValidatorBasic(
            "ipp_algorithm",
            value -> check_and_return_from_map(value, ipp_algorithm_map)
        )
    ],
    "customer_options" => [
        FieldValidatorBasic(
            "pv_adoption_type",
            value -> check_and_return_from_map(value, pv_adoption_type_map)
        )
    ],
    "der_aggregator_options" => [
        FieldValidatorHasDefault(
            "incentive_curve",
            value -> check_integer(value; min=1, max=5),
            1
        ),
        FieldValidatorHasDefault(
            "frac_viu_cost_savings_as_revenue",
            value -> check_float(value; min=0.0, max=1.0),
            0.5
        ),
    ],
)
