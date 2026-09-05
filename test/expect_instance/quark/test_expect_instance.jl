using RandomMeasAdd

# quark 入口：目前只有 pidx=1 小数据。
# 数据集选择见各测试文件头；按需注释掉不需要的测试。
pidx = 1
permuted_order = [1, 2, 3, 4]
run_tag = "setting0_settings81_shots1024"

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
