using RandomMeasAdd

# hist 导入手工验证脚本（不在 runtests.jl 内，手动跑）：
# 整个 npz（一次 SettingRun）即一个 HistGroup，共 4 个被测 site。
test_index = 2

if test_index == 1
    group_path = joinpath(
        @__DIR__, "data", "aer_independent_pidx_-1",
        "aer_independent_pidx_-1_setting1_settings729_shots1024.npz",
    )
    # 缺省链式导入
    permuted_order = [1, 2, 3, 4]
    hgroup, htrivial, _ = import_hist_group(
        group_path; permuted_order=permuted_order,
    )
    @show (hgroup.N, hgroup.NU, hgroup.NM)
    @show !isnothing(htrivial)
end

if test_index == 2
    pair_dir = joinpath(@__DIR__, "data", "aer_pair_pidx_-1")
    run_tag = "setting1_settings729_shots1024"
    filepath1 = joinpath(pair_dir, "aer_pair_pidx_-1_exp1_$(run_tag).npz")
    filepath2 = joinpath(pair_dir, "aer_pair_pidx_-1_exp2_$(run_tag).npz")
    # 配对导入：共用同一 permuted_order，奇位 = I_1、偶位 = I_2
    permuted_order = [1, 2, 3, 4]
    h1, h2, _, _, _ = import_hist_pair(
        filepath1, filepath2; permuted_order=permuted_order,
    )
    @show (h1.N, h1.NU, h1.NM)
    @show (h2.N, h2.NU, h2.NM)
end
