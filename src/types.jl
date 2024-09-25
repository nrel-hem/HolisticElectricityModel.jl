# ------------------------------------------------------------------------------
# Base Data Types
# ------------------------------------------------------------------------------

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
Base.iterate

Dimension(x::Dimension) = Dimension(x.name, x.prose_name, x.description, copy(x.elements))

function Dimension(
    dim::Dimension,
    name::AbstractString;
    prose_name::AbstractString = "",
    description::AbstractString = "",
)
    return Dimension(name, prose_name, description, copy(dim.elements))
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
    dims::NTuple{N, Dimension},
    vals::Array{Float64, N}; # TODO: Replace with something more general than Float64
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

Base.:+(x::ParamArray, y::ParamArray) = x.values + y.values
Base.:-(x::ParamArray, y::ParamArray) = x.values - y.values
#Base.:*(x::ParamArray, y::ParamArray) = x.values * y.values
Base.:(*)(x::ParamArray, y::ParamArray) = x.values .* y.values
Base.:(*)(x::KeyedArray, y::ParamArray) = x .* y.values
Base.:(*)(x::ParamArray, y::KeyedArray) = x.values .* y

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

function initialize_keyed_array(indices...; value=0.0)
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
