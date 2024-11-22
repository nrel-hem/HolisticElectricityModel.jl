using Revise
using CSV
using DataFrames
using Test

using HolisticElectricityModel

driver_name = "hem.jl"

BASE_DIR = abspath(joinpath(dirname(Base.find_package("HolisticElectricityModel")), ".."))

@testset "Test driver" begin
    script_filename = joinpath(BASE_DIR, "script", driver_name)

    include(script_filename)
    scenario_dir = basename(dirname(output_dir))
    results_dir = basename(output_dir)

    expected_results_dir = joinpath(BASE_DIR, "test", "driver_outputs", scenario_dir, results_dir)
    data = DataFrame(filename=String[], column=String[], max_abserr = Float64[], relerr_of_max_abserr = Float64[], max_relerr = Float64[], abserr_of_max_relerr = Float64[], n=Integer[], N=Integer[], orig_filepath=String[], new_filepath=String[])
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
                if !isapprox(actual_data, expected_data, nans = true, rtol = 1e-5, atol=1e-8,)
                    fail = true
                    n = length(actual_data[.~isapprox.(actual_data, expected_data, rtol = 1e-5, atol=1e-8, nans = true)])
                    biggest_diff = findmax(actual_data[.~(actual_data .== expected_data)] .- expected_data[.~(actual_data .== expected_data)])
                    max_abserr = biggest_diff[1]
                    relerr_of_max_abserr = max_abserr/expected_data[.~(actual_data .== expected_data)][biggest_diff[2]]  * 100
                    perc_err = findmax(((actual_data[.~(actual_data .== expected_data)] .- expected_data[.~(actual_data .== expected_data)]) ./ expected_data[.~(actual_data .== expected_data)])*100)
                    max_relerr = perc_err[1]
                    abserr_of_max_relerr = (actual_data[.~(actual_data .== expected_data)] .- expected_data[.~(actual_data .== expected_data)])[perc_err[2]]
                end
            elseif !(actual_data == expected_data)
                fail = true
                n = length(actual_data[.~(actual_data .== expected_data, rtol = 1e-5, atol=1e-8,)])
                biggest_diff = findmax(actual_data[.~(actual_data .== expected_data)] .- expected_data[.~(actual_data .== expected_data)])
                max_abserr = biggest_diff[1]
                relerr_of_max_abserr = max_abserr/expected_data[.~(actual_data .== expected_data)][biggest_diff[2]]  * 100
                perc_err = findmax(((actual_data[.~(actual_data .== expected_data)] .- expected_data[.~(actual_data .== expected_data)]) ./ expected_data[.~(actual_data .== expected_data)])*100)
                max_relerr = perc_err[1]
                abserr_of_max_relerr = (actual_data[.~(actual_data .== expected_data)] .- expected_data[.~(actual_data .== expected_data)])[perc_err[2]]
            end
            if fail
                push!(data, (basename(exp_file), col, max_abserr, relerr_of_max_abserr, max_relerr, abserr_of_max_relerr, n, length(actual_data), exp_file, actual_filename))
            end
        end    
    end
    @test isempty(data) || @warn "Some result files do not match:\n$(data)"
    if !isempty(data)
        CSV.write(joinpath(output_dir, "..", "Mismatched_Data_for_$(results_dir).csv"), data)
    end
end
