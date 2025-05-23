test_dir = dirname(@__FILE__)
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
        joinpath(test_dir, "data"),
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
            joinpath(test_dir, "data"),
            "test_input_param_county_invalid",
            index_t,
            [index_h, index_m, index_z, index_d]
        )
    end

    param_df = CSV.read(
        joinpath(test_dir, "data", "test_input_param_county.csv"),
        DataFrame,
    )

    for (z, h) in index_z_h_map, m in index_m, d in index_d, t in index_t
        @test valid_param(h, m, z, d, t) == param_df[
            (param_df.index_h.==String(h)).&(param_df.index_z.==String(z)).&(param_df.index_m.==String(m)).&(param_df.index_d.==String(d)), String(t)
            ][1]
    end


end