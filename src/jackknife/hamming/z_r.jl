# ---------------------
# hamming 具体力学量：z_r 的估计与 loos
# ---------------------
# 说明：分子（reflect）分母（子系统纯度）的 hamming 估计各自是 setting
# 平均，只有比值 z_r 需要 jackknife（见 hamming/basics.jl）。

"""
把 reflect 和 purity 的 setting 估计组合起来，算 z_r 估计子及其
jackknife 值（hamming 版本）。

参数
- reflect_ests::Vector{Float64}：每个设置下的 reflect 估计。
- odd_ests::Vector{Float64}：每个设置下其中一个纯度分区的估计。
- even_ests::Vector{Float64}：每个设置下另一个互补纯度分区的估计。

返回
- z_r_val::Float64：组合后的估计值。
- z_r_loos::Vector{Float64}：z_r 的 jackknife leave-one-out 估计。

说明
- z_r = R / sqrt((P_odd + P_even)/2)，其中 R 为 reflect 期望，
  P_odd/P_even 为两个纯度估计；估计值用三路均值组合，
  loos 用三路 setting-level loos 组合。
"""
function get_z_r_loos_hamming(
    reflect_ests::Vector{Float64},
    odd_ests::Vector{Float64},
    even_ests::Vector{Float64},
)
    n = length(reflect_ests)
    @assert length(odd_ests) == n && length(even_ests) == n "三路估计的设置数必须一致。"
    @assert n ≥ 2 "At least 2 settings are required for leave-one-out estimation."

    Z_R(R_I_val, P_I1, P_I2) = R_I_val / sqrt((P_I1 + P_I2) / 2)
    z_r_est = Z_R(mean(reflect_ests), mean(odd_ests), mean(even_ests))
    z_r_loos = Z_R.(get_setting_loos(reflect_ests), get_setting_loos(odd_ests), get_setting_loos(even_ests))

    return z_r_est, z_r_loos
end
