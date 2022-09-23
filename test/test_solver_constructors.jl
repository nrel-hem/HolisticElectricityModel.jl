@testset "Construct Ipopt_Solver default" begin
    import_ipopt()
    solver = Ipopt_Solver()
    model = HEM.get_new_jump_model(solver)
    @test model isa JuMP.Model
end

@testset "Construct Ipopt_Solver with attributes" begin
    import_ipopt()
    attr = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol" => 1e-6, "print_level" => 0)
    solver = Ipopt_Solver(attr)
    model = HEM.get_new_jump_model(solver)
    @test model isa JuMP.Model
end

# This requires an Xpress license and so may need to be disabled.
@testset "Construct Xpress_Solver default" begin
    import_xpress()
    solver = Xpress_Solver()
    model = HEM.get_new_jump_model(solver)
    @test model isa JuMP.Model
end

@testset "Construct Xpress_Solver with attributes" begin
    import_xpress()
    attr = JuMP.optimizer_with_attributes(Xpress.Optimizer, "OUTPUTLOG" => 0)
    solver = Xpress_Solver(attr)
    model = HEM.get_new_jump_model(solver)
    @test model isa JuMP.Model
end

@testset "Make JuMP model" begin
    attr = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol" => 1e-6, "print_level" => 0)
    model = HEM.make_jump_model(attr)
    @test model isa JuMP.Model
end
