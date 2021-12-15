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
    vals = read_record_file(DataFrame, filename, sheetname)
    result = Dimension(
        name,
        Symbol.(names(DataFrame(vals))),
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
    index::Dimension;
    prose_name::AbstractString = "",
    description::AbstractString = "",
)
    vals = read_record_file(AxisArray, filename, sheetname, 1)
    result = ParamAxisArray(
        name,
        (index,),
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
    dims = Tuple(push!(copy(row_indices), column_index))
    vals = read_record_file(AxisArray, filename, sheetname, length(dims))
    ar_axes = AxisArrays.axes(vals)
    if length(ar_axes) != length(dims)
        throw(
            ArgumentError(
                "dimension length of AxisArray ($(length(ar_axes))) does not match passed dimensions ($(length(dims)))",
            ),
        )
    end

    for (i, ax) in enumerate(ar_axes)
        elements = dims[i].elements
        if ax.val != elements
            throw(
                ArgumentError(
                    "dimension elements of AxisArray ($ax) does not match passed dimension ($elements)",
                ),
            )
        end
    end
    result =
        ParamAxisArray(name, dims, vals, prose_name = prose_name, description = description)
    @debug "Loaded $sheetname" result
    return result
end

function read_axis_array(filename::AbstractString, num_dims::Int)
    @debug "read_axis_array" filename num_dims
    file = open(filename) do io
        CSV.File(io)
    end

    isempty(file) && error("$filename is empty")
    return _read_nd_axis_array(file, num_dims)
end

function _read_1d_axis_array(file::CSV.File)
    data = [getproperty(file, x)[1] for x in file.names]
    return AxisArray(data, file.names)
end

function _read_2d_axis_array(file::CSV.File)
    num_rows = length(file)
    num_columns = length(file.columns) - 1
    row_names = Symbol.(file.columns[1])
    column_names = file.names[2:end]

    # TODO DT: There is no place to store this. Do we need it?
    row_variable_name = file.names[1]

    data = Array{Float64,2}(undef, num_rows, num_columns)
    for (i, column) in enumerate(file.columns[2:end])
        data[1:end, i] = column
    end

    return AxisArray(data, row_names, column_names)
end

function _read_3d_axis_array(file::CSV.File)
    i_names = Symbol.(unique(file.columns[1]))
    j_names = Symbol.(unique(file.columns[2]))
    k_names = Symbol.(file.names[3:end])

    num_i = length(i_names)
    num_j = length(j_names)
    num_k = length(file.columns) - 2

    data =
        AxisArray(Array{Float64,3}(undef, num_i, num_j, num_k), i_names, j_names, k_names)
    for i = 1:(file.rows)
        i_val = Symbol(file.columns[1][i])
        j_val = Symbol(file.columns[2][i])
        k_vals = [file.columns[x+2][i] for x = 1:num_k]
        data[i_val, j_val] = k_vals
    end

    # TODO DT: There is no place to store this. Do we need it?
    #row_variable_name = file.names[1]

    return data
end

function _read_4d_axis_array(file::CSV.File)
    i_names = Symbol.(unique(file.columns[1]))
    j_names = Symbol.(unique(file.columns[2]))
    k_names = Symbol.(unique(file.columns[3]))
    m_names = Symbol.(file.names[4:end])

    num_i = length(i_names)
    num_j = length(j_names)
    num_k = length(k_names)
    num_m = length(file.columns) - 3

    data = AxisArray(
        Array{Float64,4}(undef, num_i, num_j, num_k, num_m),
        i_names,
        j_names,
        k_names,
        m_names,
    )
    for i = 1:(file.rows)
        i_val = Symbol(file.columns[1][i])
        j_val = Symbol(file.columns[2][i])
        k_val = Symbol(file.columns[3][i])
        m_vals = [file.columns[x+3][i] for x = 1:num_m]
        data[i_val, j_val, k_val] = m_vals
    end

    # TODO DT: There is no place to store this. Do we need it?
    #row_variable_name = file.names[1]

    return data
end

function _read_nd_axis_array(file::CSV.File, num_dims::Int)
    if num_dims == 1
        return _read_1d_axis_array(file)
    elseif num_dims == 2
        return _read_2d_axis_array(file)
    elseif num_dims == 3
        return _read_3d_axis_array(file)
    elseif num_dims == 4
        return _read_4d_axis_array(file)
    else
        error("read_axis_array does not support $num_dims dimensions")
    end
end

function read_record_file(::Type{DataFrame}, filename, sheetname)
    # TODO DT
    base_dir = joinpath("..", "HolisticElectricityModel-Data", "inputs", "workbooks")
    workbook_dir = joinpath(base_dir, splitext(basename(filename))[1])
    record_file = joinpath(workbook_dir, sheetname * ".csv")
    if !isfile(record_file)
        df = DataFrame(XLSX.readtable(filename, sheetname)...)
        mkpath(workbook_dir)
        to_csv(df, record_file)
        @info "Converted $filename $sheetname to $record_file"
    end

    return read_dataframe(record_file)
end

function read_record_file(::Type{AxisArray}, filename, sheetname, num_dims)
    base_dir = joinpath("..", "HolisticElectricityModel-Data", "inputs", "workbooks")
    workbook_dir = joinpath(base_dir, splitext(basename(filename))[1])
    record_file = joinpath(workbook_dir, sheetname * ".csv")
    if !isfile(record_file)
        df = DataFrame(XLSX.readtable(filename, sheetname)...)
        mkpath(workbook_dir)
        to_csv(df, record_file)
        @info "Converted $filename $sheetname to $record_file"
    end

    return read_axis_array(record_file, num_dims)
end

function save_param(
    vals::AxisArray,
    set_names::Array,
    value_name::Symbol,
    filepath::AbstractString,
)
    # TODO: Replace this with more efficient code

    # make sure we always have the same order
    indices =
        sort!(reshape(collect(Iterators.product(AxisArrays.axes(vals)...)), length(vals)))

    # get categorical data
    result = DataFrame()
    for (i, set_name) in enumerate(set_names)
        data = []
        for index in indices
            push!(data, index[i])
        end
        result[!, set_name] = data
    end

    # get values
    data = [vals[x...] for x in indices]
    result[!, value_name] = data

    # save out
    CSV.write(filepath, result)
    return result
end

function compute_difference_one_norm(before_after_pairs)
    result = 0.0
    for (before, after) in before_after_pairs
        result += sum((
            abs(after[i...] - before[i...]) for
            i in Iterators.product(AxisArrays.axes(before)...)
        ))
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
    filename::Union{Nothing,AbstractString} = "hem.log",
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
#function initialize_param(
#    name::AbstractString,
#    indices::Dimension;
#    value = 0.0,
#    prose_name = "",
#    description = "",
#)
#    return ParamAxisArray(
#        name,
#        indices,
#        Dict(t => value for t in indices),
#        prose_name = prose_name,
#        description = description,
#    )
#end

function initialize_param(
    name::AbstractString,
    index::Dimension;
    value = 0.0,
    prose_name = "",
    description = "",
)
    return ParamAxisArray(
        name,
        (index,),
        AxisArray(
            fill!(Vector{Float64}(undef, length(index.elements)), value),
            index.elements,
        ),
        prose_name = prose_name,
        description = description,
    )
end

"""
Return an AxisArray with all values set to value, and indices formed from
Iterators.product(indices...).
"""
function initialize_param(
    name::AbstractString,
    indices...;
    value = 0.0,
    prose_name = "",
    description = "",
)
    num_dims = length(indices)
    return ParamAxisArray(
        name,
        indices,
        AxisArray(
            fill!(Array{Float64,num_dims}(undef, (length(x) for x in indices)...), value),
            (x.elements for x in indices)...,
        ),
        prose_name = prose_name,
        description = description,
    )
end

function read_dataframe(filename::AbstractString)
    ext = lowercase(splitext(filename)[2])
    if ext == ".csv"
        return read_dataframe(CSV.File, filename)
    end

    @assert false "read_dataframe does not support extension $ext"
end

function read_dataframe(::Type{CSV.File}, filename::AbstractString)
    open(filename) do io
        CSV.read(io, DataFrame)
    end
end

function to_csv(df, filename::AbstractString)
    open(filename, "w") do io
        CSV.write(io, df)
    end
end
