"""
Convert configuration files from the old format to the new format.
To run this script, run:
```julia`
julia convert_config.jl <input_file_or_directory> [overwrite]
```
where <input_file_or_directory> is the path to a single YAML file or a directory containing YAML files.
If `overwrite` is true, the converted files will be saved with the same name as the input file.
Otherwise, the converted files will be saved with "_converted" appended to the filename.
"""
using YAML
using HolisticElectricityModel

CONVERSION_MAP = Dict{String, String}(
    "wholesale_market" => "WM",
    "vertically_integrated_utility" => "VIU",
    "der_use_case" => "DERAdoption",
    "supply_choice_use_case" => "SupplyChoice",
    "der_aggregation_use_case" => "DERAggregation",
    "null_use_case" => "NullUseCase",
    "flat_rate" => "FlatRate",
    "time_of_use" => "TOU",
    "excess_retail_rate" => "ExcessRetailRate",
    "excess_marginal_cost" => "ExcessMarginalCost",
    "excess_zero" => "ExcessZero",
    "lagrange_decomposition" => "LagrangeDecomposition",
    "mppdcmer_transportation_storage" => "MPPDCMERTransStorage",
    "mppdcmer" => "MPPDCMER",
    "miqp" => "MIQP",
    "standalone_pv" => "StandalonePVOnly",
    "solar_plus_storage" => "SolarPlusStorageOnly",
    "compete_der_configs" => "CompeteDERConfigs"
)

function convert_config(input_file_name::AbstractString, output_file_name::AbstractString)
    config = YAML.load_file(input_file_name)

    # Convert the configuration to a new format
    new_config = Dict{String, Any}()

    # DataSelection
    if haskey(config, "data_selection")
        new_config["DataSelection"] = Dict{String, Any}()
        for (key, value) in config["data_selection"]
            if key == "input_path"
                new_config["DataSelection"]["input_path"] = value
            end
        end
    end
    if haskey(config, "run_options")
        new_config["RunOptions"] = Dict{String, Any}()
        for (key, value) in config["run_options"]
            if key == "output_dir"
                new_config["RunOptions"]["output_dir"] = value
            end
        end
    end
    # SimulationParameters
    if haskey(config, "simulation_parameters")
        new_config["SimulationParameters"] = Dict{String, Any}()
        for (key, value) in config["simulation_parameters"]
            if key == "solver"
                new_config["SimulationParameters"]["solver"] = value
            end
        end
    end
    # HEMOptions
    if haskey(config, "hem_options")
        new_config["HEMOptions"] = Dict{String, Any}()
        for (key, value) in config["hem_options"]
            new_config["HEMOptions"][key] = get(CONVERSION_MAP, value, value) 
        end
    end
    # Regulator
    if haskey(config, "regulator_options")
        new_config["Regulator"] = Dict{String, Any}()
        for (key, value) in config["regulator_options"]
            new_config["Regulator"][key] = get(CONVERSION_MAP, value, value)
        end
    end
    # IPPGroup
    if haskey(config, "ipp_options")
        new_config["IPPGroup"] = Dict{String, Any}()
        for (key, value) in config["ipp_options"]
            new_config["IPPGroup"][key] = get(CONVERSION_MAP, value, value)
        end
    end
    # CustomerGroup
    if haskey(config, "customer_options")
        new_config["CustomerGroup"] = Dict{String, Any}()
        for (key, value) in config["customer_options"]
            new_config["CustomerGroup"][key] = get(CONVERSION_MAP, value, value)
        end
    end

    # DERAggregator
    if haskey(config, "der_aggregator_options")
        new_config["DERAggregator"] = Dict{String, Any}()
        for (key, value) in config["der_aggregator_options"]
            new_config["DERAggregator"][key] = get(CONVERSION_MAP, value, value)
        end
    end
    # Save the new configuration to a YAML file
    @info "Converting configuration from $(input_file_name) to $(output_file_name)"
    YAML.write_file(output_file_name, new_config)
end

if length(ARGS) > 1
    overwrite = ARGS[2] == "true" || Bool(ARGS[2])
elseif length(ARGS) == 1
    overwrite = false
else
    error("Usage: julia convert_config.jl <input_file_or_directory> [overwrite]")
end

input_file_name = Vector{String}()

if isdir(ARGS[1])
    for file in readdir(ARGS[1]; join=true)
        file_name = basename(file)
        if endswith(file_name, ".yaml") || endswith(file_name, ".yml")
            push!(input_file_name, file)
        end
    end
else
    push!(input_file_name, ARGS[1])
end

for file in input_file_name
    output_file_name = overwrite ? file : replace(file, r"\.yaml$" => "_converted.yaml")
    convert_config(file, output_file_name)
end