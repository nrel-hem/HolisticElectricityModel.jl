isdefined(Base, :__precompile__) && __precompile__()

module HolisticElectricityModel

################################################################################
# Exports

# Meta
export HEMSolver
export XpressSolver
export GurobiSolver
export IpoptSolver
export import_gurobi
export import_ipopt
export import_xpress
export HEMData
export Options
export HEMOptions
export AbstractAgent
export AgentGroup
export Agent
export AgentOptions, NullAgentOptions
export AgentAndOptions
export AgentOrOptions
export configure_logging
export read_dataframe
export AgentStore
export iter_agents_and_options
export get_file_prefix

# Agents
export Regulator
export Utility
export CustomerGroup
export IPPGroup
export GreenDeveloper
export DistributionUtility
export DERAggregator

# Agent Options
export AgentOptionsStore
export RegulatorOptions
export CustomerOptions
export GreenDeveloperOptions
export IPPOptions
export UtilityOptions
export DERAggregatorOptions

# Policies
export FlatRate, TOU                                    # tariff structures
export ExcessRetailRate, ExcessMarginalCost, ExcessZero # exported DG treatment
export VIU, WM     # regulatory structures

# Modeling Options
export NullUseCase
export DERAdoption, SupplyChoice                                   # types of customer decisions
export DERAggregation                                              # presence of aggregators
export StandalonePVOnly, SolarPlusStorageOnly, CompeteDERConfigs   # DER types
export LagrangeDecomposition, MIQP, MPPDCMER, MPPDCMERTransStorage # ipp algorithms

# HEM Functions
export run_hem
export solve_equilibrium_problem!
export solve_agent_problem!
export save_results
export welfare_calculation!
export solve_agent_problem_decomposition_by_year
export solve_agent_problem_decomposition_by_year_feasible
export solve_agent_problem_decomposition_by_year_feasible_obj
export solve_agent_problem_decomposition_by_year_master

# Data Structures
export initialize_param
export make_keyed_array

# Solvers
export import_solver_package
export get_optimizer_for_solver
export get_new_jump_model
export addsolvers_ipp!

################################################################################
# Imports

using Logging
using DataFrames
using DataStructures
using Logging
using JuMP
using XLSX
using DelimitedFiles
using Statistics
using Lazy: @forward
import AxisKeys
import AxisKeys: KeyedArray
import CSV
import Distributions
import InfrastructureSystems
import InfrastructureSystems: @assert_op
import QuadGK
import Tables
import TimerOutputs

const IS = InfrastructureSystems

#################################################################################

using DocStringExtensions

@template (FUNCTIONS, METHODS) = """
                                 $(TYPEDSIGNATURES)
                                 $(DOCSTRING)
                                 """

#################################################################################
# Includes

include("solvers.jl")
include("types.jl")
include("utils.jl")

include("agents/common.jl")
include("agents/regulator.jl")
include("agents/utility.jl")
include("agents/customer_group.jl")
include("agents/ipp_group.jl")
include("agents/green_developer.jl")
include("agents/distribution_utility.jl")
include("agents/der_aggregator.jl")
include("run_hem.jl")

end # module
