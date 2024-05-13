abstract type HEMSolver end

function make_jump_model(attributes::MOI.OptimizerWithAttributes)
    model = Model()
    JuMP.set_optimizer(model, attributes)
    return model
end

struct AnySolver <: HEMSolver
    attributes::MOI.OptimizerWithAttributes
end

get_new_jump_model(hem_solver::AnySolver) = make_jump_model(hem_solver.attributes)

"""
Construct an Ipopt optimizer. Use default settings or pass an MOI.OptimizerWithAttributes
instance.

If using the default settings, call `import_ipopt()` first. Ipopt must be installed.

# Examples
```
julia> import_ipopt()
julia> hem_solver = IpoptSolver()

julia> import Ipopt
julia> import JuMP
julia> ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol" => 1e-6, "print_level" => 0)
julia> solver = IpoptSolver(ipopt)
```
"""
struct IpoptSolver <: HEMSolver
    attributes::Union{Nothing,MOI.OptimizerWithAttributes}

    IpoptSolver(attributes = nothing) = new(attributes)
end


function get_new_jump_model(hem_solver::IpoptSolver)
    isnothing(hem_solver.attributes) && return Model(_get_ipopt())
    make_jump_model(hem_solver.attributes)
end

"""
Construct an Xpress optimizer. Use default settings or pass an MOI.OptimizerWithAttributes
instance.

If using the default settings, call `import_xpress()` first. Xpress must be installed.

# Examples
```
julia> import_xpress()
julia> hem_solver = XpressSolver()

julia> import Xpress
julia> import JuMP
julia> xpress = JuMP.optimizer_with_attributes(
    Xpress.Optimizer,
    "MIPRELSTOP" => 0.01,
    "OUTPUTLOG" => 1,
)
julia> solver = XpressSolver(xpress)
```
"""
struct XpressSolver <: HEMSolver
    attributes::Union{Nothing,MOI.OptimizerWithAttributes}

    XpressSolver(attributes = nothing) = new(attributes)
end

function get_new_jump_model(hem_solver::XpressSolver)
    isnothing(hem_solver.attributes) && return Model(_get_xpress())
    make_jump_model(hem_solver.attributes)
end

"""
Construct an Gurobi optimizer. Use default settings or pass an MOI.OptimizerWithAttributes
instance.

If using the default settings, call `import_gurobi()` first. Gurobi must be installed.

# Examples
```
julia> import_gurobi()
julia> hem_solver = GurobiSolver()

julia> import Gurobi
julia> import JuMP
julia> gurobi = JuMP.optimizer_with_attributes(
    () -> Gurobi.Optimizer(Gurobi.Env()),
    "MIPFocus" => 3,
)
julia> solver = GurobiSolver(gurobi)
```
"""
struct GurobiSolver <: HEMSolver
    env::Any
    attributes::Union{Nothing,MOI.OptimizerWithAttributes}

    GurobiSolver(env, attributes = nothing) = new(GurobiSolver(env, attributes))
end

function get_new_jump_model(hem_solver::GurobiSolver)
    isnothing(hem_solver.attributes) && return Model(() -> _get_gurobi()(hem_solver.env))
    model = Model()
    make_jump_model(hem_solver.attributes)
end

function import_gurobi()
    try
        @eval Main begin
            import Gurobi
            return Gurobi.Env()
        end
    catch
        error("The Gurobi optimizer must be installed.")
    end
end

function _get_gurobi()
    @eval Main begin
        import Gurobi
        return Gurobi.Optimizer
    end
end

function import_ipopt()
    try
        @eval Main begin
            import Ipopt
        end
    catch
        error("The Ipopt optimizer must be installed.")
    end
end

function _get_ipopt()
    @eval Main begin
        import Ipopt
        return Ipopt.Optimizer
    end
end

function import_xpress()
    try
        @eval Main begin
            import Xpress
        end
    catch
        error("The Xpress optimizer must be installed.")
    end
end

function _get_xpress()
    @eval Main begin
        import Xpress
        return Xpress.Optimizer
    end
end
