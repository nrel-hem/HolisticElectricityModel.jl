isdefined(Base, :__precompile__) && __precompile__()

module HolisticElectricityModel

################################################################################
# Exports

# Meta
export HEMSolver
export XpressSolver
export GurobiSolver
export HEMData
export HEMOptions
export Agent
export AgentOptions, NullAgentOptions
export AgentAndOptions
export configure_logging

# Agents
export Regulator
export Utility
export Customers
export IPP

# Agent Options
export RegulatorOptions

# Policies
export FlatRate, TOU                                    # tariff structures
export ExcessRetailRate, ExcessMarginalCost, ExcessZero # exported DG treatment
export VerticallyIntegratedUtility, WholesaleMarket     # regulatory structures

# Solvers
export solve_equilibrium_problem

################################################################################
# Imports

using Logging
using DataFrames
using DataStructures
using Logging
using JuMP
using XLSX
using Lazy: @forward
import InfrastructureSystems
import Distributions
import CSV

const IS = InfrastructureSystems

#################################################################################

using DocStringExtensions

@template (FUNCTIONS, METHODS) = """
                                 $(TYPEDSIGNATURES)
                                 $(DOCSTRING)
                                 """

#################################################################################
# Includes

include("types.jl")
include("utils.jl")

include("agents/common.jl")
include("agents/regulator.jl")
include("agents/utility.jl")
include("agents/customers.jl")
include("agents/ipp.jl")

end # module
