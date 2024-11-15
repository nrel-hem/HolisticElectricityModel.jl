# include config.jl before this file

reeds_bas = collect("p$n" for n = 1:134)

null_use_case_identifier = "null_use_case"

market_structure_map = Dict(
    "wholesale_market" => WholesaleMarket(),
    "vertically_integrated_utility" => VerticallyIntegratedUtility()
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

ipp_algorithm_map = Dict(
    "lagrange_decomposition" => LagrangeDecomposition(),
    "mppdcmer_transportation_storage" => MPPDCMERTransStorage(),
    "mppdcmer" => MPPDCMER(),
    "miqp" => MIQP()
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

pv_adoption_type_map = Dict(
    "standalone_pv" => StandalonePVOnly(),
    "solar_plus_storage" => SolarPlusStorageOnly(),
    "compete_standalone_pv_solar_storage" => Compete_StandalonePV_SolarPlusStorage()
)

validators = Dict(
    "data_selection" => [
        FieldValidatorBasic(
            "balancing_areas",
            value -> check_iterable(value, 1, x -> check_in_collection(x, reeds_bas)),
        ),
        FieldValidatorBasic(
            "base_year",
            value -> check_integer(value, 2000, 2050)
        ),
        FieldValidatorBasic(
            "num_future_years",
            value -> check_integer(value, 0, 100)
        ),
        FieldValidatorBasic(
            "num_ipps",
            value -> check_integer(value, 1)
        ),
        FieldValidatorHasDefault(
            "folder_suffix",
            value -> check_string(value),
            String("")
        ),
        FieldValidatorBasic(
            "der_aggregator",
            value -> check_bool(value)
        ),
        FieldValidatorBasic(
            "set_nuclear_varcost_negative",
            value -> check_bool(value)
        ),
        FieldValidatorHasDefault(
            "sectors_with_county_level_load",
            value -> check_iterable(value, 1, x -> check_string(x)),
            Vector{String}()
        ),
        FieldValidatorHasDefault(
            "load_profiles_subdir",
            value -> check_string(value),
            nothing
        ),
        FieldValidatorBasic(
            "skip_parse",
            value -> check_bool(value)
        )
    ],
    "simulation_parameters" => [
        FieldValidatorBasic(
            "solver",
            value -> check_in_collection(value, ("Gurobi", "Xpress"))
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
    "ipp_options" => [
        FieldValidatorBasic(
            "ipp_algorithm",
            value -> check_and_return_from_map(value, ipp_algorithm_map)
        )
    ],
    "regulator_options" => [
        FieldValidatorBasic(
            "rate_design",
            value -> check_and_return_from_map(value, rate_design_map)
        ),
        FieldValidatorBasic(
            "net_metering_policy",
            value -> check_and_return_from_map(value, net_metering_policy_map)
        )
    ],
    "customer_options" => [
        FieldValidatorBasic(
            "pv_adoption_type",
            value -> check_and_return_from_map(value, pv_adoption_type_map)
        )
    ]
)
