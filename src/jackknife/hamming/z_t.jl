# ---------------------
# hamming 具体力学量：z_t 的估计与 loos
# ---------------------
# 说明：分子（两份实验的互关联 Z_T）分母（实验一的子系统纯度）各自是
# setting 平均，只有比值 z_t 需要 jackknife（见 hamming/basics.jl）。

"""
把互关联和纯度的 setting 估计组合起来，算 z_t 估计子及其
jackknife 值（hamming 版本）。

参数
- zt_ests::Vector{Float64}：每个配对设置下的互关联 Z_T 估计。
- odd_ests::Vector{Float64}：每个设置下 I_1 分区（奇位）的纯度估计。
- even_ests::Vector{Float64}：每个设置下 I_2 分区（偶位）的纯度估计。

返回
- z_t_val::Float64：组合后的估计值。
- z_t_loos::Vector{Float64}：z_t 的 jackknife leave-one-out 估计。

说明
- z_t = Z_T / ((P_odd + P_even)/2)^(3/2)，其中 Z_T 为互关联，
  P_odd/P_even 为两个纯度估计；估计值用三路均值组合，
  loos 用三路 setting-level loos 组合。
"""
function get_z_t_loos_hamming(
    zt_ests::Vector{Float64},
    odd_ests::Vector{Float64},
    even_ests::Vector{Float64},
)
    n = length(zt_ests)
    @assert length(odd_ests) == n && length(even_ests) == n "三路估计的设置数必须一致。"
    @assert n ≥ 2 "At least 2 settings are required for leave-one-out estimation."

    Z_T_norm(Z_T_val, P_I1, P_I2) = Z_T_val / ((P_I1 + P_I2) / 2)^(3 / 2)
    z_t_est = Z_T_norm(mean(zt_ests), mean(odd_ests), mean(even_ests))
    z_t_loos = Z_T_norm.(get_setting_loos(zt_ests), get_setting_loos(odd_ests), get_setting_loos(even_ests))

    return z_t_est, z_t_loos
end
