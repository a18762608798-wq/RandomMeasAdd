using RandomMeasAdd

# 单独运行时提供缺省值；被入口 include 时沿用入口的设置。
if !@isdefined(pidx)
    pidx = 1
end
if !@isdefined(permuted_order)
    permuted_order = [1, 2, 3, 4]
end
if !@isdefined(run_tag)
    run_tag = "setting0_settings81_shots1024"
end

# hamming reversal（Z_T 分子）：quark 配对小数据。
# NOTE: 配对数据的列已按 I_1、I_2 交错排好，重排 frame 下奇位 = I_1、偶位 = I_2。
pair_dir = joinpath(@__DIR__, "..", "..", "..", "data", "quark_pair_pidx_$(pidx)")
filepath1 = joinpath(pair_dir, "quark_pair_pidx_$(pidx)_exp1_$(run_tag).npz")
filepath2 = joinpath(pair_dir, "quark_pair_pidx_$(pidx)_exp2_$(run_tag).npz")
@show get_reversal_hamming(
    filepath1, filepath2; permuted_order=permuted_order, is_compute_sem=true, is_show_progress=false,
)
