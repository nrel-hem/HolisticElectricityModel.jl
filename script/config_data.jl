# include config.jl before this file

reeds_bas = collect("p$n" for n = 1:134)

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
            "skip_parse",
            value -> check_bool(value)
        )
    ],
    "simulation_parameters" => [
        FieldValidatorBasic(
            "solver",
            value -> check_in_collection(value, ("Gurobi", "Xpress"))
        )
    ],
)
