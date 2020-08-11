using DataFrames
using Logging
using XLSX

# Question: Best practice for documenting data types? Why doesn't my type 
# mark-up in the code show up when I look at ?read_set or ?read_param
# Answer: Documentation below is correct--types are not listed explicitly in the docstrings. 
#     using DocStringExtensions is your friend for pulling that information 
#     together. For example, see 
#     https://github.com/NREL-SIIP/PowerSystems.jl/blob/master/src/PowerSystems.jl#L292

# Question: Also, when should I specify types?
# Answer: When necessary for dispatch, and often for clarity. 
#         Always use AbstractString, not String.
#         Iterables can be tricky--go with duck typing when you can in that case

const SetType1D = Array{Symbol,1}
const ParamType1D = Dict{Symbol, AbstractFloat}
const ParamTypeND = Dict{Tuple, AbstractFloat}

"""
    read_set(filename, sheetname)

Reads the column names of Excel filename, sheetname in as symbols.
"""
function read_set(filename::AbstractString, sheetname::AbstractString)
    result = Symbol.(names(DataFrame(XLSX.readtable(filename, sheetname)...)))
    @info "Loaded $sheetname = $result" result
    return result
end


"""
    read_param(filename, sheetname, column_index, row_indices=[])

Reads parameter data from the Excel workbook at filename. Assumes the data is
in sheetname, with a single header row, and that the column_index dimension is
expressed in the column names. All other data dimension values are listed in the 
first length(row_indices) columns of sheetname, corresonding one-to-one and in 
the same order as the row_indices.

Returns the data loaded into a Dict. If isempty(row_indices), the Dict has the 
Symbol values in column_index as its keys. If ~isempty(row_indices), then the 
Dict keys are tuples of length(row_indices) + 1, and with values taken from each
of the row_indices in turn, plus a column_index value.
"""
function read_param(filename::AbstractString, 
                    sheetname::AbstractString, 
                    column_index::Array{Symbol,1}, 
                    row_indices::Array{Array{Symbol,1},1}=Array{Array{Symbol,1},1}())
    # load the sheet as a DataFrame
    table = DataFrame(XLSX.readtable(filename, sheetname)...)

    # determine the return type
    n = 0
    result = Dict{Symbol, AbstractFloat}()
    if ~isempty(row_indices)
        n = length(row_indices)
        result = Dict{Tuple{Vararg{Symbol,length(row_indices)+1}}, AbstractFloat}()
    end

    # process each row
    for row in eachrow(table)
        # the Dict key starts with the first n values in row
        preamble = isempty(row_indices) ? [] : Symbol.(Array(row[1:n]))
        # check that the key values are expected
        for (i, rindex) in enumerate(row_indices)
            @assert preamble[i] in rindex "Entry $i $(preamble[i]) not in $rindex"
        end
        # now register all of the values in this row
        for j in column_index
            key = j
            if ~isempty(row_indices)
                key = Tuple(push!(copy(preamble), j))
            end
            push!(result, key => row[j])
        end
    end
    @debug "Loaded $sheetname" result
    return result
end


function save_param(param::Dict, set_names::Array, value_name::Symbol, filepath::AbstractString)
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
    set_log_level(level)

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
            current_logger.message_limits
        )
    end
    if isnothing(new_logger)
        @warn "Unhandled logger type $(typeof(new_logger)) with fields $(fieldnames(typeof(new_logger)))"
        return
    end
    global_logger(new_logger)
end

"""
    initialize_param(indices...; value=0.0)

Returns a Dict with all values set to value, and keys formed from 
Iterators.product(indices..). If only one list of indices is passed, the key is 
not a tuple, but is instead a bare Symbol.
"""
function initialize_param(indices...; value=0.0)
    if length(indices) == 1
        return Dict(t[1] => value for t in Iterators.product(indices...))
    end
    return Dict(t => value for t in Iterators.product(indices...))
end

# Code from https://stackoverflow.com/a/60907180 that provides a stop() function
# that can be helpful for testing.

function stop()
    throw(StopException("Stop."))
end

struct StopException{T}
    S::T
end

function Base.showerror(io::IO, ex::StopException, bt; backtrace=true)
    Base.with_output_color(get(io, :color, false) ? :green : :nothing, io) do io
        showerror(io, ex.S)
    end
end