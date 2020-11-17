using Lazy: @forward

const SetType1D = Vector{Symbol}
const ParamType1D = Dict{Symbol, AbstractFloat}
const ParamTypeND = Dict{Tuple, AbstractFloat}

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

@forward Set1D.elements Base.iterate

# ------------------------------------------------------------------------------
# Collections
# ------------------------------------------------------------------------------

# https://docs.julialang.org/en/v1/base/base/#Base.getproperty

struct AgentData{T <: HEMSet}
    sets::Vector{T}
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

function Base.getproperty(agent_data::AgentData, sym::Symbol)
    val = get_set(agent_data, sym)
    if not isnothing(val)
        return val
    end
    return getfield(agent_data, sym)
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