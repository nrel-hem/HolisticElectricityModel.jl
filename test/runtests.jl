using Logging
using Test
using TestSetExtensions
import InfrastructureSystems
const IS = InfrastructureSystems
using HolisticElectricityModel
const HEM = HolisticElectricityModel

import Aqua
Aqua.test_unbound_args(HolisticElectricityModel)
Aqua.test_undefined_exports(HolisticElectricityModel)
Aqua.test_ambiguities(HolisticElectricityModel)
Aqua.test_stale_deps(HolisticElectricityModel)
Aqua.test_deps_compat(HolisticElectricityModel)

BASE_DIR = abspath(joinpath(dirname(Base.find_package("HolisticElectricityModel")), ".."))
DATA_DIR = joinpath(BASE_DIR, "..", "HolisticElectricityModel-Data")

LOG_FILE = "hem.log"
LOG_LEVELS = Dict(
    "Debug" => Logging.Debug,
    "Info" => Logging.Info,
    "Warn" => Logging.Warn,
    "Error" => Logging.Error,
)

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
        @info IS.report_log_summary(multi_logger)
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
