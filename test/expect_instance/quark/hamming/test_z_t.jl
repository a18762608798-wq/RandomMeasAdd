using RandomMeasAdd

# 单独运行时提供缺省值；被入口 include 时沿用入口的设置。
if !@isdefined(pidx)
    pidx = 1
end
if !@isdefined(permuted_order)
    permuted_order = [1, 2, 3, 4]
end

# hamming z_t：quark 配对文件夹（exp1/exp2 按文件名配对），
# 返回 (vals, sems) 两个列表。
# NOTE: 配对数据的列已按 I_1、I_2 交错排好，重排 frame 下奇位 = I_1、偶位 = I_2。
pair_dir = joinpath(@__DIR__, "..", "..", "..", "data", "quark_pair_pidx_$(pidx)")
@show get_z_t_hamming(
    pair_dir; permuted_order=permuted_order, is_show_progress=false,
)
