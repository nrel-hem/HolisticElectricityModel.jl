isdefined(Base, :__precompile__) && __precompile__()

module HolisticElectricityModel

################################################################################
# Exports

# Meta
export HEMSolver
export XpressSolver
export GurobiSolver
export Ipopt_Solver
export HEMData
export HEMOptions
export AbstractAgent
export AgentGroup
export Agent
export AgentOptions, NullAgentOptions
export AgentAndOptions
export configure_logging
export read_dataframe
export AgentStore
export iter_agents_and_options

# Agents
export Regulator
export Utility
export CustomerGroup
export IPPGroup

# Agent Options
export RegulatorOptions
export CustomerOptions
export IPPOptions

# Policies
export FlatRate, TOU                                    # tariff structures
export ExcessRetailRate, ExcessMarginalCost, ExcessZero # exported DG treatment
export VerticallyIntegratedUtility, WholesaleMarket     # regulatory structures
export DERAdoption, SupplyChoice                        # consumer decisions
export LagrangeDecomposition, MIQP                      # ipp algorithms

# Solvers
export solve_equilibrium_problem!
export solve_agent_problem!
# export Lagrange_Sub_Investment_Retirement_Cap
# export Lagrange_Sub_Dispatch_Cap
# export Lagrange_Feasible_Cap
# export solve_agent_problem_ipp_lagrange_cap
export save_results
export welfare_calculation!
# export solve_agent_problem_ipp_energy_cap_combined
export solve_agent_problem_decomposition_by_year
export solve_agent_problem_decomposition_by_year_feasible
export solve_agent_problem_decomposition_by_year_feasible_obj
export solve_agent_problem_decomposition_by_year_master

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
import QuadGK

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
include("agents/customer_group.jl")
include("agents/ipp_group.jl")

end # module
