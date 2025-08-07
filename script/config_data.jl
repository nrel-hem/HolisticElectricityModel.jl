# include config.jl before this file

reeds_bas = collect("p$n" for n = 1:134)

null_use_case_identifier = "NullUseCase"

market_structure_map = Dict(
    "WM" => WM(),
    "VIU" => VIU()
)

der_use_case_map = Dict(
    "DERAdoption" => DERAdoption(),
    null_use_case_identifier => NullUseCase()
)

supply_choice_use_case_map = Dict(
    "supply_choice_use_case" => SupplyChoice(),
    null_use_case_identifier => NullUseCase()
)

der_aggregation_use_case_map = Dict(
    "DERAggregation" => DERAggregation(),
    null_use_case_identifier => NullUseCase()
)

rate_design_map = Dict(
    "FlatRate" => FlatRate(),
    "TOU" => TOU()
)

net_metering_policy_map = Dict(
    "ExcessRetailRate" => ExcessRetailRate(),
    "ExcessMarginalCost" => ExcessMarginalCost(),
    "ExcessZero" => ExcessZero()
)

ipp_algorithm_map = Dict(
    "LagrangeDecomposition" => LagrangeDecomposition(),
    "MPPDCMERTransStorage" => MPPDCMERTransStorage(),
    "MPPDCMER" => MPPDCMER(),
    "MIQP" => MIQP()
)

pv_adoption_type_map = Dict(
    "StandalonePVOnly" => StandalonePVOnly(),
    "SolarPlusStorageOnly" => SolarPlusStorageOnly(),
    "CompeteDERConfigs" => CompeteDERConfigs()
)

validators = Dict(
    "DataSelection" => [
        FieldValidatorBasic(
            "input_path",
            value -> check_path(value)
        ),
        FieldValidatorHasDefault(
            "year_start",
            value -> check_integer(value),
            2020
        ),
    ],
    "RunOptions" => [
        FieldValidatorHasDefault(
            "output_dir",
            value -> check_string(value),
            nothing
        ),
    ],      
    "SimulationParameters" => [
        FieldValidatorBasic(
            "solver",
            value -> check_chain(value, [
                val -> check_in_collection(val, ("Gurobi", "Xpress")),
                check_symbol
            ])
        ),
        FieldValidatorHasDefault(
            "delta_t",
            value -> check_integer(value; min=1, max=24),
            4
        ),
    ],
    "HEMOptions" => [
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
    "Regulator" => [
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
    "IPPGroup" => [
        FieldValidatorBasic(
            "ipp_algorithm",
            value -> check_and_return_from_map(value, ipp_algorithm_map)
        )
    ],
    "CustomerGroup" => [
        FieldValidatorBasic(
            "pv_adoption_type",
            value -> check_and_return_from_map(value, pv_adoption_type_map)
        )
    ],
    "DERAggregator" => [
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
