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

const DimensionKey{N} = NTuple{N,Symbol}

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
mutable struct ParamScalar{T<:Number} <: HEMParameter
    name::String
    prose_name::String
    description::String
    value::T
end

function ParamScalar(name::AbstractString, value::Number; prose_name = "", description = "")
    return ParamScalar(name, prose_name, description, value)
end

@forward ParamScalar.value Base.isless, Base.isgreater, Base.:+, Base.:*, Base.:-, Base.:/

Base.:+(x::Number, y::ParamScalar) = x + y.value
# Base.:+(x::Ref{Float64}, y::Number) = x[] + y
Base.:-(x::Number, y::ParamScalar) = x - y.value
Base.:*(x::Number, y::ParamScalar) = x * y.value
Base.:*(x::ParamScalar, y::ParamScalar) = x.value * y.value
# Base.:*(x::Ref{Float64}, y::ParamScalar) = x[] * y.value
# Base.:*(x::Ref{Float64}, y::Number) = x[] * y
# Base.:*(x::Ref{Int64}, y::JuMP.VariableRef) = x[] * y
Base.:/(x::Number, y::ParamScalar) = x / y.value
Base.:/(x::ParamScalar, y::ParamScalar) = x.value / y.value
# Base.:/(x::Ref{Float64}, y::Number) = x[] / y
Base.isless(x::Number, y::ParamScalar) = isless(x, y.value)
# Base.zero(x::ParamScalar) = zero(x.value[])

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

mutable struct ParamAxisArray{N} <: HEMParameter
    name::String
    prose_name::String
    description::String
    dims::NTuple{N,Dimension}
    values::AxisArray
end

ParamAxisArray(param::ParamAxisArray) = ParamAxisArray(
    param.name,
    param.prose_name,
    param.description,
    param.dims,
    copy(param.values),
)

function ParamAxisArray(
    param::ParamAxisArray,
    name::AbstractString;
    prose_name::AbstractString = "",
    description::AbstractString = "",
)
    return ParamAxisArray(name, prose_name, description, param.dims, copy(param.values))
end

function ParamAxisArray(
    name::AbstractString,
    dims::NTuple{N,Dimension},
    vals::AxisArray;
    prose_name = "",
    description = "",
) where {N}
    return ParamAxisArray(name, prose_name, description, dims, vals)
end

@forward ParamAxisArray.values Base.getindex,
Base.setindex!,
Base.findmax,
Base.fill!,
Base.length,
AxisArrays.axes

# TODO PERF: turn off fill_nan when we are confident in the code.
"""
Return an uninitialized AxisArray from any number of Dimension values.
"""
function make_axis_array(indices...; fill_nan = true)
    array = AxisArray(
        Array{Float64, length(indices)}(undef, length.(indices)...),
        [getproperty(x, :elements) for x in indices]...,
    )
    fill_nan && fill!(array.data, NaN)
    return array
end

# ------------------------------------------------------------------------------
# Optimization Solvers
# ------------------------------------------------------------------------------

abstract type HEMSolver end

struct XpressSolver <: HEMSolver
    solver::Any
end

function get_new_jump_model(hem_solver::XpressSolver)
    return Model(hem_solver.solver.Optimizer)
end

struct GurobiSolver <: HEMSolver
    solver::Any
    env::Any
end

function get_new_jump_model(hem_solver::GurobiSolver)
    return Model(() -> hem_solver.solver.Optimizer(hem_solver.env))
end

struct Ipopt_Solver <: HEMSolver
    solver::Any
end

function get_new_jump_model(hem_solver::Ipopt_Solver)
    return Model(hem_solver.solver.Optimizer)
end
