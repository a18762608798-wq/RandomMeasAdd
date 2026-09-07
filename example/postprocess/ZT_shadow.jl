# Z_t 示例：shadow 路线（aer 模拟数据，不做误差缓解）。
using RandomMeasAdd

# 文件夹输入：逐个 .npz 计算，返回 (vals, sems)；bias 不保留。
vals, sems = get_z_t_shadow(
    joinpath(@__DIR__, "../aer/data", "aer_independent_pidx_-1");
    permuted_order=[1, 2, 3, 4], # 缺省 nothing（恒等置换）
    is_mitigation=true,        
    is_compute_sem=true,         # 缺省 false（文件夹模式内部强制求 SEM）
    is_show_progress=true,       # 缺省 false
)
@show vals sems