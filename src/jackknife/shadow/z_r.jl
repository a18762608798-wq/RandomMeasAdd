# ---------------------
# shadow z_r 的估计与 loos
# ---------------------

"""
把 reflect 和 purity 估计组合起来，算 z_r 估计子及其 jackknife 值。

参数
- shadows::Array{<:AbstractShadow, 2}：算 reflect 算符期望用的 shadow
  (n_ru, n_m)。
- odd_shadows::Array{<:AbstractShadow, 2}：其中一个纯度分区的 shadow。
- even_shadows::Array{<:AbstractShadow, 2}：另一个互补纯度分区的 shadow。
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
    shadows::Array{<:AbstractShadow, 2},
    odd_shadows::Array{<:AbstractShadow, 2},
    even_shadows::Array{<:AbstractShadow, 2},
    reflect_op::MPO,
    is_show_progress::Bool=true,
)
    # 预先枚举置换（和测量的笛卡尔积）
    n_ru, _ = size(shadows)
    @assert n_ru ≥ 3 "At least 3 random unitaries are required for z_r estimation."
    reflect_combs = collect(combinations(1:n_ru, 1))
    purity_combs = collect(combinations(1:n_ru, 2))

    # 每个置换对测量求平均
    reflect_comb_avgs = get_comb_avgs_shadow(
        shadows, 1; O=reflect_op, is_show_progress=is_show_progress
    )
    reflect_expect = mean(reflect_comb_avgs)
    odd_comb_avgs = get_comb_avgs_shadow(
        odd_shadows, 2; is_show_progress=is_show_progress
    )
    odd_expect = mean(odd_comb_avgs)
    even_comb_avgs = get_comb_avgs_shadow(
        even_shadows, 2; is_show_progress=is_show_progress
    )
    even_expect = mean(even_comb_avgs)

    # Loo 分组：leave-one-out。
    # 取 reflect 的 leave-one-out。
    reflect_loos = get_combs_loos_shadow(
        n_ru,
        reflect_combs,
        reflect_comb_avgs,
    )
    # 取 purity 的 leave-one-out。
    odd_loos = get_combs_loos_shadow(
        n_ru,
        purity_combs,
        odd_comb_avgs,
    )
    even_loos = get_combs_loos_shadow(
        n_ru,
        purity_combs,
        even_comb_avgs,
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
    shadows::Array{<:AbstractShadow, 1},
    odd_shadows::Array{<:AbstractShadow, 1},
    even_shadows::Array{<:AbstractShadow, 1},
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
