using RandomMeasAdd

# 单独运行时提供缺省值；被入口 include 时沿用入口的设置。
if !@isdefined(pidx)
    pidx = 1
end
if !@isdefined(permuted_order)
    permuted_order = [1, 2, 3, 4]
end

# shadow z_r：quark independent 文件夹，返回 (vals, sems) 两个列表。
# NOTE: 这里对permted_order的 odd 和 even 相间要求已经被在get_data的时候考虑.
group_dir = joinpath(
    @__DIR__, "..", "..", "..", "data", "quark_independent_pidx_$(pidx)",
)
@show get_z_r_shadow(
    group_dir; permuted_order=permuted_order, is_mitigation=true, is_show_progress=false,
)
