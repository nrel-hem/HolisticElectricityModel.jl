using Lazy: @forward

# ------------------------------------------------------------------------------
# Base Data Types
# ------------------------------------------------------------------------------

abstract type HEMSymbol end

# HERE -- What properties should all HEMSymbols have? 
# HERE -- What is key for working with JuMP?
# HERE -- Why doesn't JuMP have better support for this kind of stuff?
# 1. Start with JuMP docs
# 2. Maybe take a quick look at Victor Zavela's packages
# 3. Other JuMP-compatible packages that may be helpful?

name(hemsym::HEMSymbol) = hemsym.name

symbol(hemsym::HEMSymbol) = Symbol(symbol(name(hemsym)))

prose_name(hemsym::HEMSymbol) = hemsym.prose_name

description(hemsym::HEMSymbol) = hemsym.description

abstract type HEMSet <: HEMSymbol end

struct Set1D <: HEMSet
    name::AbstractString
    prose_name::AbstractString
    description::AbstractString
    elements::Vector{Symbol}
end

function Set1D(name::AbstractString, elements::Vector{Symbol}; prose_name = "", description = "")
    return Set1D(name, prose_name, description, elements)
end

@forward Set1D.elements Base.IteratorSize, Base.IteratorEltype, Base.size, 
    Base.axes, Base.ndims, Base.length, Base.iterate

abstract type HEMParameter <: HEMSymbol end

struct ParamScalar{T<:Number} <: HEMParameter
    name::AbstractString
    prose_name::AbstractString
    description::AbstractString
    value::T
end

function ParamScalar(name::AbstractString, value::Number; prose_name = "", 
                     description = "")
    return ParamScalar(name, prose_name, description, value)
end

@forward ParamScalar.value Base.isless, Base.isgreater, Base.:+, Base.:*, 
    Base.:-, Base.:/

Base.:+(x::Number, y::ParamScalar) = x + y.value
Base.:*(x::Number, y::ParamScalar) = x * y.value

struct ParamVector <: HEMParameter
    name::AbstractString
    prose_name::AbstractString
    description::AbstractString
    dim::Set1D
    values::Dict{Symbol, AbstractFloat}

    function ParamVector(name::AbstractString, prose_name::AbstractString, 
        description::AbstractString, dim::Set1D, values::Dict{Symbol})
        
        # Check that values are defined over dim
        set_symbols = Set([sym for sym in dim])
        value_symbols = Set([k for k in keys(values)])
        invalid_symbols = setdiff(value_symbols, set_symbols)
        if !isempty(invalid_symbols)
            error("Attempted ParamVector definition of $name over $name(dim) is invalid. "*
                "Values provided for symbols $invalid_symbols that are not in $dim.")
        end
        missing_symbols = setdiff(set_symbols, value_symbols)
        if !isempty(missing_symbols)
            @warn "ParamVector definition of $name over $name(dim) is missing "*
                "values for $missing_symbols"
        end

        new(name, prose_name, description, dim, values)
    end
end

function ParamVector(name::AbstractString, dim::Set1D, values::Dict{Symbol}; 
                     prose_name = "", description = "")
    return ParamVector(name, prose_name, description, dim, values)
end

@forward ParamVector.values Base.getindex

struct ParamArray{N} <: HEMParameter
    name::AbstractString
    prose_name::AbstractString
    description::AbstractString
    dims::NTuple{N,Set1D}
    values::Dict{NTuple{N,Symbol},AbstractFloat}

    function ParamArray(name::AbstractString, prose_name::AbstractString, 
        description::AbstractString, dims::NTuple{N,Set1D}, 
        values::Dict{NTuple{N,Symbol}}) where {N}

        for (index, dim) in enumerate(dims)
            # Check that values are defined over dim
            set_symbols = Set([sym for sym in dim])
            value_symbols = Set([k[index] for k in keys(values)])
            invalid_symbols = setdiff(value_symbols, set_symbols)
            if !isempty(invalid_symbols)
                error("Attempted ParamArray definition of $name, but definition "*
                    "over dim $index, $name(dim), is invalid. "*
                    "Values provided for symbols $invalid_symbols that are not in $dim.")
            end
            missing_symbols = setdiff(set_symbols, value_symbols)
            if !isempty(missing_symbols)
                @warn "ParamArray definition of $name is missing "*
                    "values for some of $(name(dim))'s elements: "*
                    "$missing_symbols"
            end
        end

        # Check for expected number of combinations
        n_expected = reduce(*, map(length, dims))
        n_actual = length(values)
        @assert n_actual <= n_expected "There are at most $n_expected combinations of $dims, but $n_actual parameter values were passed"
        if n_actual < n_expected
            @warn "ParamArray definition of $name has not defined values for all"*
                "$n_expected possible combinations of $dims"
        end

        new{N}(name, prose_name, description, dims, values)
    end
end

function ParamArray(name::AbstractString, dims::NTuple{N,Set1D}, 
        values::Dict{NTuple{N,Symbol}}; prose_name = "", 
        description = "") where {N}
    return ParamArray(name, prose_name, description, dims, values)
end

@forward ParamArray.values Base.getindex, Base.keys, Base.setindex!

# ETH@20210218 - For now forwarding values to copy. This came up because of 
# x_DG_before = copy(customers.x_DG_new). The copy doesn't really have the same
# name as x_DG_new, and it's just a temporary variable, so for now I am leaving 
# it as a bare Dict of values.
@forward ParamArray.values Base.copy

# ------------------------------------------------------------------------------
# Collections
# ------------------------------------------------------------------------------

# https://docs.julialang.org/en/v1/base/base/#Base.getproperty

struct AgentData{T <: HEMSet, U <: HEMParameter}
    sets::Vector{T}
    parameters::Vector{U}
end

function get_set(agent_data::AgentData, sym::Symbol)
    for a_set in agent_data.sets
        if symbol(a_set) == sym
            return a_set
        end
    end
    return nothing
end

function get_set(agent_data::AgentData, name::AbstractString)
    for a_set in agent_data.sets
        if name(a_set) == name
            return a_set
        end
    end
    return nothing
end

function get_parameter(agent_data::AgentData, sym::Symbol)
    for a_param in agent_data.parameters
        if symbol(a_param) == sym
            return a_param
        end
    end
    return nothing
end

function get_parameter(agent_data::AgentData, name::AbstractString)
    for a_param in agent_data.parameters
        if name(a_param) == name
            return a_param
        end
    end
    return nothing
end

function Base.getproperty(agent_data::AgentData, sym::Symbol)
    val = get_set(agent_data, sym)
    if not isnothing(val)
        return val
    end
    val = get_parameter(agent_data, sym)
    if not isnothing(val)
        return val
    end
    return getfield(agent_data, sym)
end

function Base.propertynames(agent_data::AgentData)
    # static properties
    result = [fn for fn in fieldnames(typeof(agent_data))]
    
    # properties defined by contents of agent_data
    for a_set in agent_data.sets
        append!(result, symbol(a_set))
    end
    for a_param in agent_data.parameters
        append!(result, symbol(a_param))
    end

    return Tuple(result)
end

# ------------------------------------------------------------------------------
# Optimization Solvers
# ------------------------------------------------------------------------------

abstract type HEMSolver end

struct XpressSolver <: HEMSolver
    solver
end

function get_new_jump_model(hem_solver::XpressSolver)
    return Model(hem_solver.solver.Optimizer)
end

struct GurobiSolver <: HEMSolver
    solver
    env
end

function get_new_jump_model(hem_solver::GurobiSolver)
    return Model(() -> hem_solver.solver.Optimizer(hem_solver.env))
end
