using RandomMeasAdd

# 单独运行时提供缺省值；被入口 include 时沿用入口的设置。
if !@isdefined(pidx_list)
    pidx_list = [1, 0, -1]
end
if !@isdefined(permuted_order)
    permuted_order = [1, 2, 3, 4]
end
if !@isdefined(run_tag)
    run_tag = "setting1_settings729_shots1024"
end

# hamming reversal（Z_T 分子）：配对数据，逐 pidx 跑。
# NOTE: 配对数据的列已按 I_1、I_2 交错排好（见 get_data/gen_pairs.py），
# 重排 frame 下奇位 = I_1、偶位 = I_2。
for pidx in pidx_list
    println("pidx = $pidx")
    pair_dir = joinpath(@__DIR__, "..", "..", "..", "data", "aer_pair_pidx_$(pidx)")
    filepath1 = joinpath(pair_dir, "aer_pair_pidx_$(pidx)_exp1_$(run_tag).npz")
    filepath2 = joinpath(pair_dir, "aer_pair_pidx_$(pidx)_exp2_$(run_tag).npz")
    @show get_reversal_hamming(
        filepath1, filepath2; permuted_order=permuted_order, is_compute_sem=true, is_show_progress=true,
    )
end
