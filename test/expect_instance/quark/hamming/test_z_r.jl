using RandomMeasAdd

# 单独运行时提供缺省值；被入口 include 时沿用入口的设置。
if !@isdefined(pidx)
    pidx = 1
end
if !@isdefined(permuted_order)
    permuted_order = [1, 2, 3, 4]
end

# hamming z_r：quark shared 文件夹，返回 (vals, sems) 两个列表。
group_dir = joinpath(
    @__DIR__, "..", "..", "..", "data", "quark_shared_pidx_$(pidx)",
)
@show get_z_r_hamming(
    group_dir; permuted_order=permuted_order, is_show_progress=false,
)
