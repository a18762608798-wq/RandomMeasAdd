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

# shadow reversal（Z_T 分子）：普通单实验 independent 数据，逐 pidx 跑。
# NOTE: 列序 [(2,), (5,), (3,), (4,)] 下奇位 {2, 3} = I_1、偶位 {5, 4} = I_2。
for pidx in pidx_list
    println("pidx = $pidx")
    group_path = joinpath(
        @__DIR__, "..", "..", "..", "data", "aer_independent_pidx_$(pidx)",
        "aer_independent_pidx_$(pidx)_$(run_tag).npz",
    )
    @show get_reversal_shadow(
        group_path; permuted_order=permuted_order, is_compute_sem=true, is_show_progress=true,
    )
end
