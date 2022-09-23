abstract type HEMSolver end

function make_jump_model(attributes::MOI.OptimizerWithAttributes)
    model = Model()
    JuMP.set_optimizer(model, attributes)
    return model
end

struct Any_Solver <: HEMSolver
    attributes::MOI.OptimizerWithAttributes
end

get_new_jump_model(hem_solver::Any_Solver) = make_jump_model(hem_solver.attributes)

"""
Construct an Ipopt optimizer. Use default settings or pass an MOI.OptimizerWithAttributes
instance.

If using the default settings, call `import_ipopt()` first. Ipopt must be installed.

# Examples
```
julia> import_ipopt()
julia> hem_solver = Ipopt_Solver()

julia> import Ipopt
julia> import JuMP
julia> ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol" => 1e-6, "print_level" => 0)
julia> solver = Ipopt_Solver(ipopt)
```
"""
struct Ipopt_Solver <: HEMSolver
    attributes::Union{Nothing, MOI.OptimizerWithAttributes}

    Ipopt_Solver(attributes=nothing) = new(attributes)
end


function get_new_jump_model(hem_solver::Ipopt_Solver)
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
julia> hem_solver = Xpress_Solver()

julia> import Xpress
julia> import JuMP
julia> xpress = JuMP.optimizer_with_attributes(
    Xpress.Optimizer,
    "MIPRELSTOP" => 0.01,
    "OUTPUTLOG" => 1,
)
julia> solver = Xpress_Solver(xpress)
```
"""
struct Xpress_Solver <: HEMSolver
    attributes::Union{Nothing, MOI.OptimizerWithAttributes}

    Xpress_Solver(attributes=nothing) = new(attributes)
end

function get_new_jump_model(hem_solver::Xpress_Solver)
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
julia> hem_solver = Gurobi_Solver()

julia> import Gurobi
julia> import JuMP
julia> gurobi = JuMP.optimizer_with_attributes(
    () -> Gurobi.Optimizer(Gurobi.Env()),
    "MIPFocus" => 3,
)
julia> solver = Gurobi_Solver(gurobi)
```
"""
struct Gurobi_Solver <: HEMSolver
    env::Any
    attributes::Union{Nothing, MOI.OptimizerWithAttributes}

    Gurobi_Solver(env, attributes=nothing) = new(Gurobi_Solver(env, attributes))
end

function get_new_jump_model(hem_solver::Gurobi_Solver)
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
