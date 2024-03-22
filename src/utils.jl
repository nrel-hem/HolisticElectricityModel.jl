# Question: Also, when should I specify types?
# Answer: When necessary for dispatch, and often for clarity. 
#         Always use AbstractString, not String.
#         Iterables can be tricky--go with duck typing when you can in that case

"""
    read_set(filename, filename)

Reads the column names of csv dirpath/filename in as symbols.
"""
function read_set(
    dirpath::AbstractString,
    filename::AbstractString,
    name::AbstractString;
    prose_name = "",
    description = "",
)
    vals = read_record_file(DataFrame, dirpath, filename)
    result = Dimension(
        name,
        Symbol.(names(DataFrame(vals))),
        prose_name = prose_name,
        description = description,
    )
    @info "Loaded $filename = $result" result
    return result
end

"""
Reads parameter data from the csv file at dirpath/filename. Assumes the data has a 
single header row, and that the index elements comprise the column names.

Returns the data loaded into a ParamAxisArray.
"""
function read_param(
    name::AbstractString,
    dirpath::AbstractString,
    filename::AbstractString,
    index::Dimension;
    prose_name::AbstractString = "",
    description::AbstractString = "",
)
    vals = read_record_file(KeyedArray, dirpath, filename, 1)
    result = ParamArray(
        name,
        (index,),
        vals,
        prose_name = prose_name,
        description = description,
    )
    @debug "Loaded $filename" result
    return result
end

"""
Reads parameter data from the csv file at dirpath/filename. Assumes the data has 
a single header row, and that the column_index elements are the right-most column 
names. All other data dimension values are listed in the first length(row_indices) 
columns of the file, corresonding one-to-one and in the same order as the row_indices.

Returns the data loaded into a ParamAxisArray with dimensions 
(row_indices..., column_index).
"""
function read_param(
    name::AbstractString,
    dirpath::AbstractString,
    filename::AbstractString,
    column_index::Dimension,
    row_indices::Vector{Dimension};
    prose_name::AbstractString = "",
    description::AbstractString = "",
)
    dims = Tuple(push!(copy(row_indices), column_index))
    vals = read_record_file(KeyedArray, dirpath, filename, length(dims))
    ar_axes = AxisKeys.axiskeys(vals)
    if length(ar_axes) != length(dims)
        throw(
            ArgumentError(
                "dimension length of KeyedArray ($(length(ar_axes))) does not match passed dimensions ($(length(dims)))",
            ),
        )
    end

    for (i, ax) in enumerate(ar_axes)
        elements = dims[i].elements
        if ax != elements
            throw(
                ArgumentError(
                    "dimension elements of KeyedArray axis $i ($ax) does not match passed dimension $i's ($elements)",
                ),
            )
        end
    end
    # TODO DT: make negative test to verify that conflicting inputs are rejected.
    result =
        ParamArray(name, dims, vals, prose_name = prose_name, description = description)
    @debug "Loaded $filename" result
    return result
end

"""
Return a KeyedArray from an N-dimensional array flattened in a CSV file.

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
function read_keyed_array(filename::AbstractString, num_dims::Int)
    @debug "read_keyed_array" filename num_dims
    file = open(filename) do io
        CSV.File(io)
    end

    isempty(file) && error("$filename is empty")
    return read_keyed_array(file, num_dims)
end

function read_keyed_array(file::CSV.File, num_dims)
    if num_dims == 1
        data = [getproperty(file, x)[1] for x in file.names]
        return KeyedArray(data, (file.names,))
    end

    index_names = Vector{Vector{Symbol}}(undef, num_dims)
    for i in 1:(num_dims - 1)
        index_names[i] = Symbol.(unique(Tables.getcolumn(file, i)))
    end
    index_names[num_dims] = Symbol.(file.names[num_dims:end])

    data =
        KeyedArray(Array{Float64, num_dims}(undef, length.(index_names)...), Tuple(index_names))
    for i in 1:(file.rows)
        indices = [Symbol(Tables.getcolumn(file, j)[i]) for j in 1:(num_dims - 1)]
        data(indices...,:,:) .=
            [Tables.getcolumn(file, x + num_dims - 1)[i] for x in 1:length(index_names[end])]
    end

    return data
end

function read_record_file(::Type{DataFrame}, dirpath, filename)
    record_file = joinpath(dirpath, filename * ".csv")
    if !isfile(record_file)
        @error "Missing $record_file"
    end

    return read_dataframe(record_file)
end

function read_record_file(::Type{KeyedArray}, dirpath, filename, num_dims)
    record_file = joinpath(dirpath, filename * ".csv")
    if !isfile(record_file)
        @error "Missing $record_file"
    end

    return read_keyed_array(record_file, num_dims)
end

function save_param(
    vals::KeyedArray,
    set_names::Array,
    value_name::Symbol,
    filepath::AbstractString,
)
    # TODO: Replace this with more efficient code

    # make sure we always have the same order
    indices =
        sort!(reshape(collect(Iterators.product(AxisKeys.axiskeys(vals)...)), length(vals)))

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
    data = [vals(x...) for x in indices]
    result[!, value_name] = data

    # save out
    CSV.write(filepath, result)
    return result
end

function compute_difference_one_norm(before_after_pairs)
    result = 0.0
    for (before, after) in before_after_pairs
        result += sum((
            abs(after(i...) - before(i...)) for
            i in Iterators.product(AxisKeyes.axiskeys(before)...)
        ))
    end
    return result
end

function compute_difference_percentage_one_norm(before_after_pairs)
    result_vec = []
    for (before, after) in before_after_pairs
        result_one = sum((
            before(i...) == 0.0 ? abs(after(i...) - before(i...)) : abs(after(i...) - before(i...)) / before(i...) for
            i in Iterators.product(AxisKeys.axiskeys(before)...)
        ))
        push!(result_vec, result_one)
    end
    result = maximum(result_vec)
    return result
end

function compute_difference_percentage_maximum_one_norm(before_after_pairs)
    result_vec = []
    for (before, after) in before_after_pairs
        for i in Iterators.product(AxisKeys.axiskeys(before)...)
            result_one = 
                before(i...) == 0.0 ? abs(after(i...) - before(i...)) : abs(after(i...) - before(i...)) / before(i...)
            push!(result_vec, result_one)
        end
    end
    result = maximum(result_vec)
    return result
end

"""
Sets the log level to, e.g., Logging.Debug, Logging.Warn. Currently only 
implemented for ConsoleLogger.
"""
function set_log_level(level::Base.CoreLogging.LogLevel)
    current_logger = global_logger()
    new_logger = nothing
    if current_logger isa ConsoleLogger
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
Returns a ParamArray with all values set to value.
"""
function initialize_param(
    name::AbstractString,
    index::Dimension;
    value = 0.0,
    prose_name = "",
    description = "",
)
    return ParamArray(
        name,
        (index,),
        initialize_keyed_array(index; value=value),
        prose_name = prose_name,
        description = description,
    )
end

"""
Return a ParamArray with all values set to value, and indices formed from
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
    param = try
        ParamArray(
            name,
            indices,
            initialize_keyed_array(indices...; value=value),
            prose_name = prose_name,
            description = description,
        )
    catch e
        @info "Failed to initialize parameter $name"
        rethrow(e)
    end
    return param
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
