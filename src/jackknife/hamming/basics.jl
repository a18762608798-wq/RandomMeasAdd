# ---------------------
# hamming 基础组件：setting 层面的 leave-one-out
# ---------------------
# 说明：hamming 估计子每个 setting 独立给出一个估计值（分子分母各自
# 普通样本平均即可），只有拼成 z_r 这种比值时才需要 jackknife。
# 此时 leave-one-out 直接在 setting 层面做，不需要置换枚举。

"""
get_setting_loos(ests)

由各 setting 的估计值算 leave-one-out 的 jackknife 估计。
输入一串数（每个 setting 的估计值），输出同样长度的一串数，
其中第 i 个 = 去掉第 i 个 setting 后剩下所有设置的均值.

参数
- ests::Vector{Float64}：每个随机幺正设置下的估计值（相互独立）。

返回
- loos::Vector{Float64}：去掉第 i 个设置后剩下设置的均值。

说明
setting i 的 leave-one-out 值 = 不含 i 的设置的均值。
"""
function get_setting_loos(ests::Vector{Float64})
    n = length(ests)
    @assert n ≥ 2 "At least 2 settings are required for leave-one-out estimation."
    ssum = sum(ests)
    return [(ssum - ests[i]) / (n - 1) for i in 1:n]
end
