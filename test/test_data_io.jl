using RandomMeasAdd

# 整个 npz（一次 SettingRun）即一个 MeasurementGroup，共 4 个被测 site。
test_index = 2

if test_index == 1
    group_path = joinpath(
        @__DIR__, "data", "aer_independent_pidx_-1",
        "aer_independent_pidx_-1_setting1_settings729_shots1024.npz",
    )
    # 缺省链式导入
    permuted_order = [1, 2, 3, 4]
    permuted_group, permuted_sites, G = import_random_group(
        group_path; permuted_order=permuted_order, is_mitigation=true,
    )
    @show permuted_sites
    @show size(permuted_group.measurements)
    @show isnothing(G)
    println(G)
end

if test_index == 2
    pair_dir = joinpath(@__DIR__, "data", "aer_pair_pidx_-1")
    run_tag = "setting1_settings729_shots1024"
    filepath1 = joinpath(pair_dir, "aer_pair_pidx_-1_exp1_$(run_tag).npz")
    filepath2 = joinpath(pair_dir, "aer_pair_pidx_-1_exp2_$(run_tag).npz")
    # 配对导入：共用同一 permuted_order，奇位 = I_1、偶位 = I_2
    permuted_order = [1, 2, 3, 4]
    group1, group2, permuted_sites, G1, G2 = import_random_pair(
        filepath1, filepath2; permuted_order=permuted_order, is_mitigation=false,
    )
    @show permuted_sites
    @show (group1.N, group1.NU)
    @show (group2.N, group2.NU)
    @show isnothing(G1) && isnothing(G2)
end


