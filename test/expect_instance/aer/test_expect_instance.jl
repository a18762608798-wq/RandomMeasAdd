using RandomMeasAdd

# 数据集选择：pidx 取 1/0/-1；
# shadow 估计子和 purity_hamming 要求各比特独立幺正，只能用 independent 数据；
# reflect_hamming 按 permuted_order 下的相邻配对算 swap，shared 数据正好对应。
pidx = 1
pidx_list = [1, 0, -1]
permuted_order = [1, 2, 3, 4]
run_tag = "setting1_settings729_shots1024"

# 按需注释掉不需要的测试（与 src/expect_instance 一一对应）：
include("shadow/test_purity.jl")
include("hamming/test_purity.jl")
include("shadow/test_reflect.jl")
include("hamming/test_reflect.jl")
include("shadow/test_z_r.jl")
include("hamming/test_z_r.jl")
include("shadow/test_reversal.jl")
include("hamming/test_reversal.jl")
include("shadow/test_z_t.jl")
include("hamming/test_z_t.jl")
include("hamming/test_z_t.jl")
