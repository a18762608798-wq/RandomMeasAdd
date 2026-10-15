# ---------------------
# shadow z_r 的估计与 loos
# ---------------------

"""
把 reflect 和 purity 估计组合起来，算 z_r 估计子及其 jackknife 值。

参数
- shadows::Array{Shadow, 2}：算 reflect 算符期望用的 shadow
  (n_ru, n_m)。
- odd_shadows::Array{Shadow, 2}：其中一个纯度分区的 shadow。
- even_shadows::Array{Shadow, 2}：另一个互补纯度分区的 shadow。
- reflect_op::MPO：算 reflect 期望用的算符。

关键词参数
- is_show_progress::Bool=true：置换平均时是否显示进度。

返回
- z_r_val::Float64：组合后的估计值。
- z_r_loos::Vector{Float64}：z_r 的 jackknife leave-one-out 估计。

说明
- z_r = R / sqrt((P_odd + P_even)/2)，其中 R 为 reflect 期望，
  P_odd/P_even 为两个纯度估计。
"""
function get_z_r_loos_shadow(
    shadows::Array{Shadow, 2},
    odd_shadows::Array{Shadow, 2},
    even_shadows::Array{Shadow, 2},
    reflect_op::MPO,
    is_show_progress::Bool=true,
)
    n_ru, _ = size(shadows)
    @assert n_ru ≥ 3 "At least 3 random unitaries are required for z_r estimation."

    # 一阶矩（reflect）与二阶矩（纯度腿）各走加和口径，与旧组合口径一致。
    reflect_expect, reflect_loos = get_moment1_loos_shadow(
        shadows; O=reflect_op, is_show_progress=is_show_progress
    )
    odd_expect, odd_loos = get_purity_loos_shadow(
        odd_shadows; is_show_progress=is_show_progress
    )
    even_expect, even_loos = get_purity_loos_shadow(
        even_shadows; is_show_progress=is_show_progress
    )

    # 取 z_r_loos
    Z_R(R_I_val, P_I1, P_I2) = R_I_val / sqrt((P_I1 + P_I2) / 2)
    z_r_est = Z_R(reflect_expect, odd_expect, even_expect)
    z_r_loos = Z_R.(reflect_loos, odd_loos, even_loos)

    return z_r_est, z_r_loos
end

"""
向量版简便重载：把输入 reshape 成二维再转调主函数 get_z_r_loos_shadow。
"""
function get_z_r_loos_shadow(
    shadows::Array{Shadow, 1},
    odd_shadows::Array{Shadow, 1},
    even_shadows::Array{Shadow, 1},
    reflect_op::MPO,
    is_show_progress::Bool=true,
)
    return get_z_r_loos_shadow(
        reshape(shadows, length(shadows), 1),
        reshape(odd_shadows, length(odd_shadows), 1),
        reshape(even_shadows, length(even_shadows), 1),
        reflect_op,
        is_show_progress,
    )
end
