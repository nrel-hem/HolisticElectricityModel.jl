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
        # compare the sets TODO DT
        if ax.val != elements
            throw(
                ArgumentError(
                    "dimension elements of AxisArray ($ax) does not match passed dimension ($elements)",
                ),
            )
        end
    end
    # TODO DT: make negative test to verify that conflicting inputs are rejected.
    result =
        ParamAxisArray(name, dims, vals, prose_name = prose_name, description = description)
    @debug "Loaded $sheetname" result
    return result
end

"""
Return an AxisArray from an N-dimensional array flattened in a CSV file.

If there is one dimension then the file must have a single row with dimension names.
If there is more than one dimension then it must conform to the following format:

The header row consists of dimension names for the first N-1 dimensions followed by the last
dimension's element ids pivoted out to form data column headers.

The data rows contain dimension element ids in the first N-1 columns followed by parameter
values in the remaining columns. Each value maps to the dimension element ids listed in the
row's first N-1 columns plus the dimension element id found in that value's column.

3-dimension example:

d1_variable_name,d2_variable_name,d3_name1,d3_name2,d3_name3
d1_name1,d2_name1,1.0,1.0,1.0
d1_name2,d2_name2,1.0,1.0,1.0
"""
function read_axis_array(filename::AbstractString, num_dims::Int)
    @debug "read_axis_array" filename num_dims
    file = open(filename) do io
        CSV.File(io)
    end

    isempty(file) && error("$filename is empty")
    return read_axis_array(file, num_dims)
end

function read_axis_array(file::CSV.File, num_dims)
    if num_dims == 1
        data = [getproperty(file, x)[1] for x in file.names]
        return AxisArray(data, file.names)
    end

    index_names = Vector{Vector{Symbol}}(undef, num_dims)
    for i in 1:(num_dims - 1)
        index_names[i] = Symbol.(unique(file.columns[i]))
    end
    index_names[num_dims] = Symbol.(file.names[num_dims:end])

    data =
        AxisArray(Array{Float64, num_dims}(undef, length.(index_names)...), index_names...)
    for i in 1:(file.rows)
        indices = [Symbol(file.columns[j][i]) for j in 1:(num_dims - 1)]
        data[indices...] =
            [file.columns[x + num_dims - 1][i] for x in 1:length(index_names[end])]
    end

    return data
end

function read_record_file(::Type{DataFrame}, filename, sheetname)
    # TODO: Remove when all .xlsx files have been converted.
    base_dir = joinpath("..", "HolisticElectricityModel-Data", "inputs")
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
    base_dir = joinpath("..", "HolisticElectricityModel-Data", "inputs")
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
Returns a ParamAxisArray with all values set to value.
"""
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
            fill!(Array{Float64, num_dims}(undef, (length(x) for x in indices)...), value),
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
