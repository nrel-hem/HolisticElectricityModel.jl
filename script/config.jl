
abstract type FieldValidator end

struct FieldValidatorBasic <: FieldValidator
    name::String
    func::Function
end

struct FieldValidatorHasDefault <: FieldValidator
    name::String
    func::Function
    default::Any
end

function validate(section::String, validator::FieldValidator, value::Any)
    ret, msg_or_value = validator.func(value)
    if !ret
        error("Error in $section, $(validator.name): $msg_or_value")
    end
    return msg_or_value
end

function check_in_collection(value, collection)
    if !(value in collection)
        return false, "$value not in $collection"
    end
    return true, value
end

function check_iterable(value, min_length, func)
    if length(value) < min_length
        return false, "only contains $length(value) elements, when at least $min_length are required"
    end
    result = []
    for x in value
        ret, msg_or_value = func(x)
        if !ret
            return ret, msg_or_value
        end
        push!(result, msg_or_value)
    end
    return true, typeof(value)(result)
end

function check_integer(value, min=nothing, max=nothing)
    result = Integer(value)

    if !isnothing(min) && (result < min)
        return false, "$result < minimum value $min"
    end
    if !isnothing(max) && (result > max)
        return false, "$result > maximum value $max"
    end

    return true, result
end

function check_string(value)
    return true, String(value)
end

function check_bool(value)
    return true, Bool(value)
end

function parse(config::Dict{Any,Any}, section::String, validators::Dict{String,Vector{T}}) where T<: FieldValidator
    if !(section in keys(config))
        error("Invalid config section: $section.")
    end

    fields = config[section]
    field_validators = validators[section]

    result = []
    for validator in field_validators
        if !(validator.name in keys(fields))
            if validator isa FieldValidatorBasic
                error("Did not find field $validator.name in $section")
            end
            push!(result, validator.default)
        else
            push!(result, validate(section, validator, fields[validator.name]))
        end
    end
    return result
end

function unpack_config_struct(config, config_field, option_dict)
    # if the config_field isn't supported
    if !(config_field in keys(config))
        error("Invalid config field: $config_field.")
        # elseif the option specified in the yaml for the config_field isn't an acceptable option
    elseif !(config[config_field] in option_dict[config_field])
        error("Invalid option $config[config_field] for config_field $config_field.")
        # all good -- unpack and return
    else
        config[config_field]
    end
end

function check_and_return_from_map(value::AbstractString, mapping::Dict)
    check_in_collection(value, keys(mapping))
    return true, mapping[value]
end

function check_path(value)
    if !(isdir(value))
        return false, "$(value) not a valid path"
    end
    return true, value
end