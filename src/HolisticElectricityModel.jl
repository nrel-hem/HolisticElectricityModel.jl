module HolisticElectricityModel

################################################################################
# Exports

# Meta
export HEMData
export HEMOptions
export Agent
export AgentOptions, NullAgentOptions
export AgentAndOptions

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

import JuMP
import Cbc
import DataFrames
import Logging
import XLSX
import Distributions
import CSV
import Xpress
import QuadGK

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
