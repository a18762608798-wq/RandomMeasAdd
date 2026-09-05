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

# hamming purity：independent 数据，各比特独立幺正。
scheme = "independent"
group_path = joinpath(
    @__DIR__, "..", "..", "..", "data", "aer_$(scheme)_pidx_$(pidx)",
    "aer_$(scheme)_pidx_$(pidx)_$(run_tag).npz",
)
@show get_purity_hamming(
    group_path; permuted_order=permuted_order, is_compute_sem=true, is_show_progress=true,
)
