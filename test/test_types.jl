"""
Tests that the values stored in two KeyedArrays are equal when accessed 
individually using their named indices.
"""
function values_equal(A::KeyedArray, B::KeyedArray)
    for inds in Base.Iterators.product(axiskeys(A)...)
        if !isapprox(A(inds...), B(inds...))
            @info "With inds=$(inds), $(A(inds...)) != $(B(inds...))"
            return false
        end
    end
    return true
end

input_dir = joinpath(BASE_DIR, "test", "data")

@testset "Read Multidimensional ParamArray" begin
    # read indices
    index_h = read_set(
        input_dir,
        "index_h",
        "index_h",
        prose_name = "customer group index h",
        description = "customer groups",
    )
    index_z = read_set(
        input_dir,
        "index_z",
        "index_z",
        prose_name = "balancing area index z",
        description = "balancing area",
    )
    index_d = read_set(
        input_dir,
        "index_d",
        "index_d",
        prose_name = "representative day index d",
        description = "ReEDS representative days",
    )
    index_t = read_set(
        input_dir,
        "index_t",
        "index_t",
        prose_name = "time index t",
        description = "ReEDS representative hour within each representative day",
    )
    
    # read parameter
    d = read_param(
        "d",
        input_dir,
        "Demand",
        index_t,
        [index_h, index_z, index_d],
    )    

    # cannot load if indices listed in wrong order
    @test_throws ArgumentError read_param(
        "d",
        input_dir,
        "Demand",
        index_h,
        [index_t, index_z, index_d],
    )

    # cannot load if indices listed in wrong order
    @test_throws ArgumentError read_param(
        "d",
        input_dir,
        "Demand",
        index_t,
        [index_h, index_d, index_z],
    )

    # cannot load if an index's elements are shuffled
    index_t_alt = Dimension(index_t)
    index_t_alt.elements = circshift(index_t.elements,1)
    @test index_t_alt.elements != index_t.elements
    @test_throws ArgumentError read_param(
        "d",
        input_dir,
        "Demand",
        index_t_alt,
        [index_h, index_z, index_d],
    )

    # this data file shuffles rows, but in a way that preserves index ordering
    d_alt = read_param(
        "d",
        input_dir,
        "Demand_rows_shuffled",
        index_t,
        [index_h, index_z, index_d],
    ) 

    # the way we read files into the KeyedArray is robust to these kinds of 
    # shuffles
    values_equal(d.values, d_alt.values)

    # this data file shuffles rows in a way that changes the order some indices 
    # are read into the KeyedArray
    @test_throws ArgumentError read_param(
        "d",
        input_dir,
        "Demand_bad_shuffle",
        index_t,
        [index_h, index_z, index_d],
    )
end

@testset "ParamArray Math" begin
    # read indices
    index_h = read_set(
        input_dir,
        "index_h",
        "index_h",
        prose_name = "customer group index h",
        description = "customer groups",
    )
    index_z = read_set(
        input_dir,
        "index_z",
        "index_z",
        prose_name = "balancing area index z",
        description = "balancing area",
    )
    index_d = read_set(
        input_dir,
        "index_d",
        "index_d",
        prose_name = "representative day index d",
        description = "ReEDS representative days",
    )
    index_t = read_set(
        input_dir,
        "index_t",
        "index_t",
        prose_name = "time index t",
        description = "ReEDS representative hour within each representative day",
    )
    
    # read parameter
    d = read_param(
        "d",
        input_dir,
        "Demand",
        index_t,
        [index_h, index_z, index_d],
    )

    d_alt = read_param(
        "d",
        input_dir,
        "Demand_rows_shuffled",
        index_t,
        [index_h, index_z, index_d],
    ) 

    result = d + d_alt
    @test result isa KeyedArray
    @test values_equal(d.values .* 2.0, result)
end