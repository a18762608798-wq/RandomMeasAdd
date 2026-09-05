using RandomMeasAdd

# 单独运行时提供缺省值；被入口 include 时沿用入口的设置。
if !@isdefined(pidx)
    pidx = 1
end
if !@isdefined(permuted_order)
    permuted_order = [1, 2, 3, 4]
end
if !@isdefined(run_tag)
    run_tag = "setting1_settings729_shots1024"
end

# hamming reflect：shared 数据，
# reflect_hamming 按 permuted_order 下的相邻配对算 swap，shared 数据正好对应。
scheme = "shared"
group_path = joinpath(
    @__DIR__, "..", "..", "..", "data", "aer_$(scheme)_pidx_$(pidx)",
    "aer_$(scheme)_pidx_$(pidx)_$(run_tag).npz",
)
@show get_reflect_hamming(
    group_path; permuted_order=permuted_order, is_compute_sem=true, is_show_progress=true,
)
