@testset "Construct IpoptSolver default" begin
    import_ipopt()
    solver = IpoptSolver()
    model = HEM.get_new_jump_model(solver)
    @test model isa JuMP.Model
end

@testset "Construct IpoptSolver with attributes" begin
    import_ipopt()
    attr = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol" => 1e-6, "print_level" => 0)
    solver = IpoptSolver(attr)
    model = HEM.get_new_jump_model(solver)
    @test model isa JuMP.Model
end

# This requires an Xpress license and so may need to be disabled.
@testset "Construct XpressSolver default" begin
    import_xpress()
    solver = XpressSolver()
    model = HEM.get_new_jump_model(solver)
    @test model isa JuMP.Model
end

@testset "Construct XpressSolver with attributes" begin
    import_xpress()
    attr = JuMP.optimizer_with_attributes(Xpress.Optimizer, "OUTPUTLOG" => 0)
    solver = XpressSolver(attr)
    model = HEM.get_new_jump_model(solver)
    @test model isa JuMP.Model
end

@testset "Make JuMP model" begin
    attr = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol" => 1e-6, "print_level" => 0)
    model = HEM.make_jump_model(attr)
    @test model isa JuMP.Model
end
