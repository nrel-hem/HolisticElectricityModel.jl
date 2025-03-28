abstract type AbstractData end

# Required interface functions:
# - Constructor(src::AbstractData, name; prose_name = "", description = "")
# Optional interface functions:
# - Default function exists for types with a field called name.
# - get_name()
# - Default function exists for types with a field called name.
# - get_description()
# - Default function exists for types with a field called prose_name.
# - get_prose_name()
# - get_symbol()

# HERE -- What properties should all AbstractDatas have? 
# HERE -- What is key for working with JuMP?
# HERE -- Why doesn't JuMP have better support for this kind of stuff?
# 1. Start with JuMP docs
# 2. Maybe take a quick look at Victor Zavela's packages
# 3. Other JuMP-compatible packages that may be helpful?

get_name(data::AbstractData) = data.name
get_symbol(data::AbstractData) = Symbol(get_name(data))
get_prose_name(data::AbstractData) = data.prose_name
get_description(data::AbstractData) = data.description

abstract type AbstractDimension <: AbstractData end

# TODO: Seems unused. Remove?
const DimensionKey{N} = NTuple{N, Symbol}

"""
Behaves like a Vector of Symbols with additional metadata fields.
"""
mutable struct Dimension <: AbstractDimension
    name::String
    prose_name::String
    description::String
    elements::Vector{Symbol}
end

function Dimension(
    name::AbstractString,
    elements::Vector{Symbol};
    prose_name = "",
    description = "",
)
    return Dimension(name, prose_name, description, elements)
end

get_pair(dim::Dimension) = Pair(get_symbol(dim), getproperty(dim, :elements))

@forward Dimension.elements Base.IteratorSize,
Base.IteratorEltype,
Base.size,
Base.axes,
Base.ndims,
Base.length,
Base.iterate,
Base.getindex,
Base.lastindex

Dimension(x::Dimension) = Dimension(x.name, x.prose_name, x.description, copy(x.elements))

function Dimension(
    dim::Dimension,
    name::AbstractString;
    prose_name::AbstractString = "",
    description::AbstractString = "",
)
    return Dimension(name, prose_name, description, copy(dim.elements))
end

mutable struct DimensionSet{N} <: AbstractDimension
    name::String
    prose_name::String
    description::String
    dims::NTuple{N, Dimension}
    set::Vector{NTuple{N, Symbol}}

    function DimensionSet{N}(
        name::String,
        prose_name::String,
        description::String,
        dims::NTuple{N, Dimension},
        set::Vector{NTuple{N, Symbol}},
    ) where N
        # verify that the ith element of each tuple in set is a member of the ith dim
        for i in 1:N
            vals = unique([elem[i] for elem in set])
            diff = setdiff(vals, dims[i].elements)
            @assert (length(diff) == 0) "Column $i in set contains values $diff not in $(dims[i])"
        end
        return new(name, prose_name, description, dims, set)
    end
end

function DimensionSet{N}(
    name::String,
    prose_name::String,
    description::String,
    dims::NTuple{N, Dimension},
    set::Vector{NTuple{N, String}},
) where N
    return DimensionSet{N}(
        name, 
        prose_name,
        description,
        dims,
        [Tuple(Symbol(x) for x in elem) for elem in set]
    )
end

function DimensionSet{N}(
    name::String,
    dims::NTuple{N, Dimension},
    set::Vector{NTuple{N, Any}};
    prose_name = "",
    description = "",
) where N
    return DimensionSet{N}(name, prose_name, description, dims, set)
end

DimensionSet{N}(set::DimensionSet{N}) where N = DimensionSet{N}(
    set.name,
    set.prose_name,
    set.description,
    set.dims,
    copy(set.set),
)

function DimensionSet{N}(
    set::DimensionSet{N},
    name::AbstractString;
    prose_name::AbstractString = "",
    description::AbstractString = "",
) where N
    return DimensionSet{N}(name, prose_name, description, set.dims, copy(set.set))
end

@forward DimensionSet.set Base.iterate

"""
Returns DimensionSet.set data as a Dict{Symbol,Vector} where the keys correspond to the 
ith dimension of the DimensionSet and the vector contains either Symbols (if N == 2) or
Tuple{N-1,Symbol} (if N > 2).
"""
function get_one_to_many_dict(set::DimensionSet{N}, i::Integer) where N
    @assert ((i <= N) & (i > 0)) "Invalid index for key."
    @assert N > 1 "DimensionSets with only one dimension (which shouldn't really exist) cannot be turned into Dicts"

    result = Dict()
    for tup in set
        key = tup[i]
        if !haskey(result, key)
            result[key] = []
        end
        val = Tuple(x for (j, x) in enumerate(tup) if j != i)
        @assert length(val) == N-1
        if N == 2
            push!(result[key], val[1])
        else
            push!(result[key], val)
        end
    end
    return result
end

function get_one_to_many_dict(set::DimensionSet{N}, dim_name::String) where N
    for (i, dim) in enumerate(set.dims)
        if dim.name == dim_name
            return get_one_to_many_dict(set, i)
        end
    end
    @error "Did not find $dim_name in $set."
end

function get_one_to_many_dict(set::DimensionSet{N}, dim_sym::Symbol) where N
    for (i, dim) in enumerate(set.dims)
        if get_symbol(dim) == dim_sym
            return get_one_to_many_dict(set, i)
        end
    end
    @error "Did not find $dim_sym in $set."
end

abstract type HEMParameter <: AbstractData end

"""
Behaves like a number with additional metadata fields.
"""
mutable struct ParamScalar{T <: Number} <: HEMParameter
    name::String
    prose_name::String
    description::String
    value::T
end

function ParamScalar(name::AbstractString, value::Number; prose_name = "", description = "")
    return ParamScalar(name, prose_name, description, value)
end

@forward ParamScalar.value Base.isless, 
Base.isgreater, 
Base.:+, 
Base.:*, 
Base.:-, 
Base.:/

Base.:+(x::Number, y::ParamScalar) = x + y.value
Base.:-(x::Number, y::ParamScalar) = x - y.value
Base.:*(x::Number, y::ParamScalar) = x * y.value
Base.:*(x::ParamScalar, y::ParamScalar) = x.value * y.value
Base.:/(x::Number, y::ParamScalar) = x / y.value
Base.:/(x::ParamScalar, y::ParamScalar) = x.value / y.value
Base.isless(x::Number, y::ParamScalar) = isless(x, y.value)

ParamScalar(param::ParamScalar) =
    ParamScalar(param.name, param.prose_name, param.description, param.value)

function ParamScalar(
    param::ParamScalar,
    name::AbstractString;
    prose_name::AbstractString = "",
    description::AbstractString = "",
)
    return ParamScalar(name, prose_name, description, copy(param.value))
end

update!(param::ParamScalar, value) = param.value = value

mutable struct ParamArray{N} <: HEMParameter
    name::String
    prose_name::String
    description::String
    dims::NTuple{N, Dimension}
    values::KeyedArray
end

ParamArray(param::ParamArray) = ParamArray(
    param.name,
    param.prose_name,
    param.description,
    param.dims,
    copy(param.values),
)

function ParamArray(
    param::ParamArray,
    name::AbstractString;
    prose_name::AbstractString = "",
    description::AbstractString = "",
)
    return ParamArray(name, prose_name, description, param.dims, copy(param.values))
end

# TODO: Validate vals against dims
function ParamArray(
    name::AbstractString,
    dims::NTuple{N, Dimension},
    vals::KeyedArray;
    prose_name = "",
    description = "",
) where {N}
    return ParamArray(name, prose_name, description, dims, vals)
end

function ParamArray(
    name::AbstractString,
    dims::NTuple{N,Dimension},
    vals::Array{Float64,N}; # TODO: Replace with something more general than Float64
    prose_name = "",
    description = "",
) where {N}
    return ParamArray(
        name,
        prose_name,
        description,
        dims,
        KeyedArray(
            vals;
            [get_pair(dim) for dim in dims]...
        ),
    )
end

@forward ParamArray.values Base.length, 
    Base.getindex,
    Base.setindex!,
    Base.iterate,
    Base.findmax,
    Base.fill!,
    Base.size,
    Base.axes,
    Base.eachindex,
    Base.keys,
    Base.strides,
    Base.transpose,
    Base.IndexStyle,
    AxisKeys.axiskeys

Base.pointer(A::ParamArray, i::Integer) = Base.pointer(A.values, i)

Base.stride(A::ParamArray, d::Integer) = Base.stride(A.values, d)

AxisKeys.axiskeys(A::ParamArray, d::Int) = AxisKeys.axiskeys(A.values, d)

# implement KeyedArray callable syntax for ParamArrays
(P::ParamArray)(args...) = P.values(args...)

#Base.:+(x::ParamArray, y::ParamArray) = x.values + y.values
#Base.:-(x::ParamArray, y::ParamArray) = x.values - y.values
#Base.:*(x::ParamArray, y::ParamArray) = x.values * y.values
#Base.:(*)(x::ParamArray, y::ParamArray) = x.values .* y.values
#Base.:(*)(x::Matrix, y::ParamArray) = x .* y.values

# TODO PERF: turn off fill_nan when we are confident in the code.
"""
Return an uninitialized KeyedArray from any number of Dimension values.
"""
function make_keyed_array(indices...; fill_nan = true)
    array = KeyedArray(
        Array{Float64, length(indices)}(undef, length.(indices)...);
        [get_pair(x) for x in indices]...,
    )
    fill_nan && fill!(array.data, NaN)
    return array
end

function initialize_keyed_array(indices...; value = 0.0)
    array = try
        KeyedArray(
            Array{Float64, length(indices)}(undef, length.(indices)...);
            [get_pair(x) for x in indices]...,
        )
    catch e
        @info "Error encountered in initialize_keyed_array with indices" indices length(indices) length.(indices) [get_pair(x) for x in indices]
        rethrow(e)
    end
    fill!(array.data, value)
    return array
end

# HERE
# How can we make it so that ParamArray = 1.0 works (after ParamArray has already been initialized)?
# Base.convert(ParamArray, x::Number) doesn't work, because that doesn't provide the instance of 
# ParamArray whose value you are trying to set. I need some function that has the lefthand and righthand
# sides of x = y so I can define UNKNOWNMETHOD(x::ParamArray, y::Number) = fill!(x.values.data, y) and
# UNKNOWNMETHOD(x::KeyedArray, y::Number) = fill!(x.data, y)

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
