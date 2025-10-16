@testset "test read_param with county level index" begin
    index_z = Dimension("index_z", [:p129, :p130, :p131])
    index_h = Dimension("index_h",
        [:Com_50001, :Com_70001, :Com_90001, :Ind_70001, :Ind_90001, :Res_10001, :Res_50001, :Res_90001])
    index_m = Dimension("index_m", [:BTMPV])
    index_d = Dimension("index_d", [:d004])
    index_t = Dimension("index_t", [:h004, :h008])
    index_z_h_map = DimensionSet{2}("index_z_h_map", "index_z_h_map", "mapping from index_z to index_h",
        (index_z, index_h),
        [
            (:p129, :Res_10001),
            (:p129, :Com_50001),
            (:p129, :Ind_90001),
            (:p130, :Res_50001),
            (:p130, :Com_90001),
            (:p130, :Ind_70001),
            (:p131, :Res_90001),
            (:p131, :Com_70001),
            (:p131, :Ind_90001),
        ]
    )

    valid_param = read_param(
        "test_param",
        joinpath(TEST_DIR, "data"),
        "test_input_param_county",
        index_t,
        [index_h, index_m, index_z, index_d]
    )

    for (z, h) in index_z_h_map, m in index_m, d in index_d, t in index_t
        @test !(ismissing(valid_param(h, m, z, d, t)))
    end

    @test_throws ArgumentError begin
        read_param(
            "test_param",
            joinpath(TEST_DIR, "data"),
            "test_input_param_county_invalid",
            index_t,
            [index_h, index_m, index_z, index_d]
        )
    end

    param_df = CSV.read(
        joinpath(TEST_DIR, "data", "test_input_param_county.csv"),
        DataFrame,
    )

    for (z, h) in index_z_h_map, m in index_m, d in index_d, t in index_t
        @test valid_param(h, m, z, d, t) == param_df[
            (param_df.index_h.==String(h)).&(param_df.index_z.==String(z)).&(param_df.index_m.==String(m)).&(param_df.index_d.==String(d)), String(t)
            ][1]
    end

end

@testset "test save_dimension" begin
    index_z = Dimension("index_z", [:p129, :p130, :p131])
    temp_dir = mktempdir()
    HEM.save_dimension(index_z, joinpath(temp_dir, "index_z.csv"))
    saved_data = read_set(
        temp_dir,
        "index_z",
        "index_z",
    )

    @test index_z.elements == saved_data.elements

end

@testset "test_read_saved_result" begin
    index_z = Dimension("index_z", [:p129, :p130, :p131])
    index_h = Dimension("index_h", [:Commercial])
    index_y = Dimension("index_y", [Symbol("2021"), Symbol("2022")])
    index_cost_type = Dimension("index_cost_type", [:Administration])

    result_cost = read_saved_result(
        "test_result_cost",
        joinpath(TEST_DIR, "data"),
        "test_result_cost",
        [index_y, index_z, index_h, index_cost_type],
        :Cost;
        column_labels = [:Year, :Zone, :CustomerType, :CostType]
    )

    @test length(result_cost.dims) == 4
    @test [dim.name for dim in result_cost.dims] == ["index_y", "index_z", "index_h", "index_cost_type"]
    @test result_cost(Symbol("2021"), :p129, :Commercial, :Administration) == 0.1111

end
 