#driver_name = "driver_gurobi.jl"
driver_name = "driver_xpress.jl"

@testset "Test driver" begin
    script_filename = joinpath(BASE_DIR, "script", driver_name)

    include(script_filename)
    scenario_dir = basename(dirname(output_dir))
    results_dir = basename(output_dir)

    expected_results_dir =
        joinpath(BASE_DIR, "test", "driver_outputs", scenario_dir, results_dir)
    for exp_file in readdir(expected_results_dir, join = true)
        !(splitext(exp_file)[2] == ".csv") && continue
        actual = sort!(read_dataframe(joinpath(output_dir, basename(exp_file))))
        expected = sort!(read_dataframe(exp_file))
        for col in names(expected)
            actual_data = actual[!, col]
            expected_data = expected[!, col]
            if eltype(expected_data) === Float64
                @test isapprox(actual_data, expected_data, nans = true)
            else
                @test actual_data == expected_data
            end
        end
    end
end
