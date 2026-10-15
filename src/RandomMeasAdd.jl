module RandomMeasAdd
# 载入外部依赖。
include("imports.jl")

# 数据导入
include("hist_data.jl")
include("data_io.jl")

# shadow 工具
include("shadow_utils.jl")

# robust Hamming 标定与核（hamming 路线误差缓解用）
include("robust_hamming/calibration.jl")
include("robust_hamming/kernels.jl")

# jackknife
include("jackknife/jackknife.jl")

# 算符构造
include("ops_creator.jl")

# 改进后的方法
include("modified_method.jl")

# 具体物理量估计
include("expect_instance/expect_instance.jl")

# 导出公开 API。
include("exports.jl")
end

