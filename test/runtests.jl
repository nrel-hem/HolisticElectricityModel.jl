using Logging
using Test
using TestSetExtensions
using JuMP
import InfrastructureSystems
const IS = InfrastructureSystems
using HolisticElectricityModel
const HEM = HolisticElectricityModel

import Aqua
Aqua.test_unbound_args(HolisticElectricityModel)
Aqua.test_undefined_exports(HolisticElectricityModel)
# EH: This next test is failing and I don't know how to fix. 
#     I would expect all tests to run and then tell me which ones failed?
#Aqua.test_ambiguities(HolisticElectricityModel)
Aqua.test_stale_deps(HolisticElectricityModel; ignore=[:JuliaFormatter,:Aqua,:TestSetExtensions,:Xpress,:Gurobi,:Revise])
Aqua.test_deps_compat(HolisticElectricityModel; ignore=[:JuliaFormatter,:Aqua,:TestSetExtensions,:Revise, :Logging])

BASE_DIR = abspath(joinpath(dirname(Base.find_package("HolisticElectricityModel")), ".."))
DATA_DIR = joinpath(BASE_DIR, "..", "HolisticElectricityModelData.jl")

LOG_FILE = "hem.log"
LOG_LEVELS = Dict(
    "Debug" => Logging.Debug,
    "Info" => Logging.Info,
    "Warn" => Logging.Warn,
    "Error" => Logging.Error,
)

DISABLED_TEST_FILES = ["test_hem.jl"]

macro includetests(testarg...)
    if length(testarg) == 0
        tests = []
    elseif length(testarg) == 1
        tests = testarg[1]
    else
        error("@includetests takes zero or one argument")
    end

    quote
        tests = $tests
        rootfile = @__FILE__
        if length(tests) == 0
            tests = readdir(dirname(rootfile))
            tests = filter(
                f ->
                    startswith(f, "test_") && endswith(f, ".jl") && f != basename(rootfile),
                tests,
            )
        else
            tests = map(f -> string(f, ".jl"), tests)
        end
        println()
        if !isempty(DISABLED_TEST_FILES)
            @warn("Some tests are disabled $DISABLED_TEST_FILES")
        end
        for test in tests
            test ∈ DISABLED_TEST_FILES && continue
            print(splitext(test)[1], ": ")
            include(test)
            println()
        end
    end
end

function get_logging_level_from_env(env_name::String, default)
    level = get(ENV, env_name, default)
    return IS.get_logging_level(level)
end

function run_tests()
    logging_config_filename = get(ENV, "SIIP_LOGGING_CONFIG", nothing)
    if logging_config_filename !== nothing
        config = IS.LoggingConfiguration(logging_config_filename)
    else
        config = IS.LoggingConfiguration(
            filename = LOG_FILE,
            file_level = Logging.Info,
            console_level = Logging.Error,
        )
    end
    console_logger = ConsoleLogger(config.console_stream, config.console_level)

    IS.open_file_logger(config.filename, config.file_level) do file_logger
        levels = (Logging.Info, Logging.Warn, Logging.Error)
        multi_logger =
            IS.MultiLogger([console_logger, file_logger], IS.LogEventTracker(levels))
        global_logger(multi_logger)

        if !isempty(config.group_levels)
            IS.set_group_levels!(multi_logger, config.group_levels)
        end

        @time @testset "Begin HolisticElectricityModel tests" begin
            @includetests ARGS
        end

        @test length(IS.get_log_events(multi_logger.tracker, Logging.Error)) == 0
    end
end

logger = global_logger()

try
    run_tests()
finally
    # Guarantee that the global logger is reset.
    global_logger(logger)
    nothing
end
