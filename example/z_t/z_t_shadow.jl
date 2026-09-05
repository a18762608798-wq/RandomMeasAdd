# Z_t 示例：shadow 路线（quark 真机数据，开误差缓解）。
# 数据由本目录 gen_data.py 在真机上生成（independent 方案，pidx=-1）。
# 运行：julia --project=. example/z_t/z_t_shadow.jl
using RandomMeasAdd

# 文件夹输入：逐个 .npz 计算，返回 (vals, sems)；bias 不保留。
vals, sems = get_z_t_shadow(
    joinpath(@__DIR__, "data", "quark_independent_pidx_-1");
    permuted_order=[1, 2, 3, 4], # 缺省 nothing（恒等置换）
    is_mitigation=true,          # 用 npz 内 trivial 数据算 G；缺省 false
    is_compute_sem=true,         # 缺省 false（文件夹模式内部强制求 SEM）
    is_show_progress=true,       # 缺省 false
)
@show vals sems
