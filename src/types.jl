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

"""
Behaves like a dictionary with one Dimension and additional metadata fields.
"""
struct ParamVector <: HEMParameter
    name::String
    prose_name::String
    description::String
    dim::Dimension
    values::Dict{Symbol, Float64}

    function ParamVector(
        name::AbstractString,
        prose_name::AbstractString,
        description::AbstractString,
        dim::Dimension,
        vals::Dict{Symbol, Float64},
    )
        # Check that values are defined over dim
        set_symbols = Set([sym for sym in dim])
        value_symbols = Set(keys(vals))
        invalid_symbols = setdiff(value_symbols, set_symbols)
        if !isempty(invalid_symbols)
            error(
                "Attempted ParamVector definition of $name over $name(dim) is invalid. " *
                "Values provided for symbols $invalid_symbols that are not in $dim.",
            )
        end
        missing_symbols = setdiff(set_symbols, value_symbols)
        if !isempty(missing_symbols)
            @warn "ParamVector definition of $name over $(get_name(dim)) is missing " *
                  "values for $missing_symbols"
        end

        new(name, prose_name, description, dim, vals)
    end
end

function ParamVector(
    name::AbstractString,
    dim::Dimension,
    vals::Dict{Symbol, Float64};
    prose_name = "",
    description = "",
)
    return ParamVector(name, prose_name, description, dim, vals)
end

@forward ParamVector.values Base.getindex, Base.setindex!, Base.keys, Base.findmax

ParamVector(param::ParamVector) = ParamVector(
    param.name,
    param.prose_name,
    param.description,
    param.dim, # sets are fully static, unlike parameters; therefore reuse them
    copy(param.values),
)

function ParamVector(
    param::ParamVector,
    name::AbstractString;
    prose_name::AbstractString = "",
    description::AbstractString = "",
)
    return ParamVector(name, prose_name, description, param.dim, copy(param.values))
end

Base.empty!(x::ParamVector) = empty!(x.values)
Base.merge!(x::ParamVector, vals::Dict{Symbol}) = merge!(x.values, vals)

# TODO: We currently can't make this immutable because of how it's used in regulator.jl.
"""
Behaves like a dictionary with N Dimensions and additional metadata fields.
"""
mutable struct ParamArray{N} <: HEMParameter
    name::String
    prose_name::String
    description::String
    dims::NTuple{N, Dimension}
    values::Dict{DimensionKey{N}, Float64}

    function ParamArray(
        name::AbstractString,
        prose_name::AbstractString,
        description::AbstractString,
        dims::NTuple{N, Dimension},
        vals::Dict{DimensionKey{N}, Float64},
    ) where {N}
        for (index, dim) in enumerate(dims)
            # Check that vals are defined over dim
            set_symbols = Set(dim)
            value_symbols = Set((k[index] for k in keys(vals)))
            invalid_symbols = setdiff(value_symbols, set_symbols)
            if !isempty(invalid_symbols)
                error(
                    "Attempted ParamArray definition of $name, but definition " *
                    "over dim $index, $(get_name(dim)), is invalid. " *
                    "Values provided for symbols $invalid_symbols that are not in $dim.",
                )
            end
            missing_symbols = setdiff(set_symbols, value_symbols)
            if !isempty(missing_symbols)
                @warn "ParamArray definition of $name is missing " *
                      "values for some of $(get_name(dim))'s elements: " *
                      "$missing_symbols"
            end
        end

        # Check for expected number of combinations
        n_expected = reduce(*, map(length, dims))
        n_actual = length(vals)
        @assert n_actual <= n_expected "There are at most $n_expected combinations of $dims, but $n_actual parameter values were passed"
        if n_actual < n_expected
            @warn "ParamArray definition of $name has not defined values for all" *
                  "$n_expected possible combinations of $dims"
        end

        new{N}(name, prose_name, description, dims, vals)
    end
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
    vals::Dict{DimensionKey{N}, Float64};
    prose_name = "",
    description = "",
) where {N}
    return ParamArray(name, prose_name, description, dims, vals)
end

@forward ParamArray.values Base.getindex, Base.keys, Base.setindex!

Base.empty!(x::ParamArray) = empty!(x.values)
Base.merge!(x::ParamArray, vals) = merge!(x.values, vals)

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