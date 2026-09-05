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

# hamming reflect：quark shared 小数据。
group_path = joinpath(
    @__DIR__, "..", "..", "..", "data", "quark_shared_pidx_$(pidx)",
    "quark_shared_pidx_$(pidx)_$(run_tag).npz",
)
@show get_reflect_hamming(
    group_path; permuted_order=permuted_order, is_compute_sem=true, is_show_progress=false,
)
