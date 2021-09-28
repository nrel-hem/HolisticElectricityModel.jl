# Question: Also, when should I specify types?
# Answer: When necessary for dispatch, and often for clarity. 
#         Always use AbstractString, not String.
#         Iterables can be tricky--go with duck typing when you can in that case

"""
    read_set(filename, sheetname)

Reads the column names of Excel filename, sheetname in as symbols.
"""
function read_set(
    filename::AbstractString,
    sheetname::AbstractString,
    name::AbstractString;
    prose_name = "",
    description = "",
)
    result = Dimension(
        name,
        Symbol.(names(DataFrame(XLSX.readtable(filename, sheetname)...))),
        prose_name = prose_name,
        description = description,
    )
    @info "Loaded $sheetname = $result" result
    return result
end

"""
Reads parameter data from the Excel workbook at filename. Assumes the data is
in sheetname, with a single header row, and that the column_index dimension is
expressed in the column names.

Returns the data loaded into a Dict with the Symbol values in column_index as its keys.
"""
function read_param(
    name::AbstractString,
    filename::AbstractString,
    sheetname::AbstractString,
    column_index::Dimension;
    prose_name::AbstractString = "",
    description::AbstractString = "",
)
    table = DataFrame(XLSX.readtable(filename, sheetname)...)
    vals = Dict{Symbol, Float64}()

    for row in eachrow(table)
        # register all of the values in this row
        for i in column_index
            push!(vals, i => row[i])
        end
    end
    result = ParamVector(
        name,
        column_index,
        vals,
        prose_name = prose_name,
        description = description,
    )
    @debug "Loaded $sheetname" result
    return result
end

"""
Reads parameter data from the Excel workbook at filename. Assumes the data is
in sheetname, with a single header row, and that the column_index dimension is
expressed in the column names. All other data dimension values are listed in the 
first length(row_indices) columns of sheetname, corresonding one-to-one and in 
the same order as the row_indices.

Returns the data loaded into a Dict. The Dict keys are tuples of
length(row_indices) + 1, and with values taken from each of the row_indices in
turn, plus a column_index value.
"""
function read_param(
    name::AbstractString,
    filename::AbstractString,
    sheetname::AbstractString,
    column_index::Dimension,
    row_indices::Vector{Dimension};
    prose_name::AbstractString = "",
    description::AbstractString = "",
)
    table = DataFrame(XLSX.readtable(filename, sheetname)...)

    n = length(row_indices)
    dims = Tuple(push!(copy(row_indices), column_index))
    vals = Dict{DimensionKey{n + 1}, Float64}()

    for row in eachrow(table)
        # the Dict key starts with the first n values in row
        preamble = Symbol.(Array(row[1:n]))
        # check that the key values are expected
        for (i, rindex) in enumerate(row_indices)
            @assert preamble[i] in rindex "Entry $i $(preamble[i]) not in $rindex"
        end
        # now register all of the values in this row
        for i in column_index
            key = Tuple(push!(copy(preamble), i))
            push!(vals, key => row[i])
        end
    end
    result =
        ParamArray(name, dims, vals, prose_name = prose_name, description = description)
    @debug "Loaded $sheetname" result
    return result
end

function save_param(
    param::Dict,
    set_names::Array,
    value_name::Symbol,
    filepath::AbstractString,
)
    # TODO: Replace this with more efficient code

    # make sure we always have the same order
    the_keys = sort([k for k in keys(param)])

    # get categorical data
    result = DataFrame()
    for (i, set_name) in enumerate(set_names)
        data = []
        for key in the_keys
            push!(data, isa(key, Tuple) ? key[i] : key)
        end
        result[!, set_name] = data
    end

    # get values
    data = []
    for key in the_keys
        push!(data, param[key])
    end
    result[!, value_name] = data

    # save out
    CSV.write(filepath, result)
    return result
end

function compute_difference_one_norm(before_after_pairs)
    result = 0.0
    for (before, after) in before_after_pairs
        result += sum(abs(after[k] - before[k]) for k in keys(before))
    end
    return result
end

"""
Sets the log level to, e.g., Logging.Debug, Logging.Warn. Currently only 
implemented for ConsoleLogger.
"""
function set_log_level(level::Base.CoreLogging.LogLevel)
    current_logger = global_logger()
    new_logger = nothing
    if typeof(current_logger) == ConsoleLogger
        new_logger = ConsoleLogger(
            current_logger.stream,
            level,
            current_logger.meta_formatter,
            current_logger.show_limited,
            current_logger.right_justify,
            current_logger.message_limits,
        )
    end
    if isnothing(new_logger)
        @warn "Unhandled logger type $(typeof(new_logger)) with fields $(fieldnames(typeof(new_logger)))"
        return
    end
    global_logger(new_logger)
end

"""
Creates console and file loggers and sets the global_logger.

**Note:** Log messages may not be written to the file until flush() or close() is called on
the returned logger.

# Arguments
- `console_level = Logging.Error`: level for console messages
- `file_level = Logging.Info`: level for file messages
- `filename::Union{Nothing, AbstractString} = "hem.log"`: log file; pass nothing
  to disable file logging

# Example
```julia
logger = configure_logging(console_level = Logging.Info)
@info "log message"
close(logger)
```
"""
function configure_logging(;
    console_level = Logging.Error,
    file_level = Logging.Info,
    filename::Union{Nothing, AbstractString} = "hem.log",
)
    return IS.configure_logging(
        console = true,
        console_stream = stderr,
        console_level = console_level,
        file = filename !== nothing,
        filename = filename,
        file_level = file_level,
        file_mode = "w+",
        tracker = nothing,
        set_global = true,
    )
end

"""
Returns a Dict with all values set to value, and keys formed from 
indices. Each key is a Symbol.
"""
function initialize_param(
    name::AbstractString,
    indices::Dimension;
    value = 0.0,
    prose_name = "",
    description = "",
)
    return ParamVector(
        name,
        indices,
        Dict(t => value for t in indices),
        prose_name = prose_name,
        description = description,
    )
end

"""
Returns a Dict with all values set to value, and keys formed from
Iterators.product(indices...).
"""
function initialize_param(
    name::AbstractString,
    indices...;
    value = 0.0,
    prose_name = "",
    description = "",
)
    return ParamArray(
        name,
        indices,
        Dict(t => value for t in Iterators.product(indices...)),
        prose_name = prose_name,
        description = description,
    )
end

function read_dataframe(filename::AbstractString)
    open(filename) do io
        CSV.read(io, DataFrame)
    end
end

# Code from https://stackoverflow.com/a/60907180 that provides a stop() function
# that can be helpful for testing.

function stop()
    throw(StopException("Stop."))
end

struct StopException{T}
    S::T
end

function Base.showerror(io::IO, ex::StopException, bt; backtrace = true)
    Base.with_output_color(get(io, :color, false) ? :green : :nothing, io) do io
        showerror(io, ex.S)
    end
end
