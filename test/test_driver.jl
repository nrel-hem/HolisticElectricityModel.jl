@testset "Test driver" begin
    script_filename = joinpath(BASE_DIR, "script", "driver.jl")
    output_dir = joinpath(DATA_DIR, "outputs")
    for file in [x for x in readdir(output_dir, join = true) if splitext(x) == ".csv"]
        isfile(file) && rm(file)
    end

    include(script_filename)

    expected_results_dir = joinpath(BASE_DIR, "test", "driver_outputs")
    for exp_file in readdir(expected_results_dir, join = true)
        actual = sort!(read_dataframe(joinpath(DATA_DIR, "outputs", basename(exp_file))))
        expected = sort!(read_dataframe(exp_file))
        for col in names(expected)
            actual_data = actual[!, col]
            expected_data = expected[!, col]
            if eltype(expected_data) === Float64
                @test isapprox(actual_data, expected_data, nans=true)
            else
                @test actual_data == expected_data
            end
        end
    end
end
