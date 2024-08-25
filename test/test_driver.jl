using CSV
using DataFrames
using Test

using HolisticElectricityModel

driver_name = "driver_gurobi.jl"
#driver_name = "driver_xpress.jl"

BASE_DIR = abspath(joinpath(dirname(Base.find_package("HolisticElectricityModel")), ".."))

@testset "Test driver" begin
    script_filename = joinpath(BASE_DIR, "script", driver_name)

    include(script_filename)
    scenario_dir = basename(dirname(output_dir))
    results_dir = basename(output_dir)

    # expected_results_dir = joinpath(BASE_DIR, "test", "driver_outputs", scenario_dir, results_dir)
    expected_results_dir = joinpath("/scratch/ehale/hem/ba_6_base_2020_future_15_ipps_1_enhanced_test_full_dera_pv_w_EV_2035_EE", results_dir)
    data = DataFrame(filename=String[], column=String[], n=Integer[], N=Integer[], orig_filepath=String[], new_filepath=String[])
    for exp_file in readdir(expected_results_dir, join = true)
        !(splitext(exp_file)[2] == ".csv") && continue
        actual_filename = joinpath(output_dir, basename(exp_file))
        @test isfile(actual_filename) || @warn "Expected file $actual_filename does not exist"
        actual = sort!(read_dataframe(actual_filename))
        expected = sort!(read_dataframe(exp_file))
        msg = "Expected dataframe $exp_file:\n$(first(expected,5))\nActual dataframe $actual_filename:\n$(first(actual,5))"
        for col in names(expected)
            fail = false
            actual_data = actual[!, col]
            expected_data = expected[!, col]
            @test length(actual_data) == length(expected_data) || @warn "When testing $(basename(exp_file)), column $col expected to find $(length(expected_data)) elements, but found $(length(actual_data)).\n$msg"
            if eltype(expected_data) === Float64
                if !isapprox(actual_data, expected_data, nans = true)
                    fail = true
                    n = length(actual_data[.~isapprox.(actual_data, expected_data, nans = true)])
                end
            elseif !(actual_data == expected_data)
                fail = true
                n = length(actual_data[.~(actual_data .== expected_data)])
            end
            if fail
                push!(data, (basename(exp_file), col, n, length(actual_data), exp_file, actual_filename))
            end
        end    
    end
    @test isempty(data) || @warn "Some result files do not match:\n$(data)"
    if !isempty(data)
        CSV.write(joinpath(output_dir, "..", "Mismatched_Data_for_$(results_dir).csv"), data)
    end
end
