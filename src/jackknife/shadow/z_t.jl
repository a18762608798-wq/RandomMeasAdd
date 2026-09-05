# ---------------------
# shadow z_t 的估计与 loos
# ---------------------

"""
把互关联和纯度估计组合起来，算 z_t 估计子及其 jackknife 值。

参数
- shadows::Array{DenseShadow, 2}：算互关联 Z_T 用的全系统 shadow
  (n_ru, n_m)。
- odd_shadows::Array{DenseShadow, 2}：I_1 分区的 shadow。
- even_shadows::Array{DenseShadow, 2}：I_2 分区的 shadow。
- uT::ITensor：u_T 算符（I_1 区逐比特 σ^y）。
- tpos::Vector{Int}：做部分转置的 site 位置（I_1 区）。

关键词参数
- is_show_progress::Bool=true：置换平均时是否显示进度。

返回
- z_t_val::Float64：组合后的估计值。
- z_t_loos::Vector{Float64}：z_t 的 jackknife leave-one-out 估计。

说明
- z_t = Z_T / ((P_odd + P_even)/2)^(3/2)，其中 Z_T 为互关联，
  P_odd/P_even 为两个纯度估计；分子分母都是 k=2 置换平均，
  loos 经 get_combs_loos_shadow 对齐后组合。
"""
function get_z_t_loos_shadow(
    shadows::Array{DenseShadow, 2},
    odd_shadows::Array{DenseShadow, 2},
    even_shadows::Array{DenseShadow, 2},
    uT::ITensor,
    tpos::Vector{Int},
    is_show_progress::Bool=true,
)
    n_ru, _ = size(shadows)
    @assert n_ru ≥ 3 "At least 3 random unitaries are required for z_t estimation."
    zt_combs = collect(combinations(1:n_ru, 2))
    purity_combs = collect(combinations(1:n_ru, 2))

    _, zt_comb_avgs = get_reversal_comb_avgs_shadow(
        vec(shadows), uT, tpos; is_show_progress=is_show_progress
    )
    zt_expect = mean(zt_comb_avgs)
    odd_comb_avgs = get_comb_avgs_shadow(
        odd_shadows, 2; is_show_progress=is_show_progress
    )
    odd_expect = mean(odd_comb_avgs)
    even_comb_avgs = get_comb_avgs_shadow(
        even_shadows, 2; is_show_progress=is_show_progress
    )
    even_expect = mean(even_comb_avgs)

    # Loo 分组：leave-one-out。
    zt_loos = get_combs_loos_shadow(
        n_ru,
        zt_combs,
        zt_comb_avgs,
    )
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
    Z_T_norm(Z_T_val, P_I1, P_I2) = Z_T_val / ((P_I1 + P_I2) / 2)^(3 / 2)
    z_t_est = Z_T_norm(zt_expect, odd_expect, even_expect)
    z_t_loos = Z_T_norm.(zt_loos, odd_loos, even_loos)

    return z_t_est, z_t_loos
end

"""
向量版简便重载：把输入 reshape 成二维再转调主函数 get_z_t_loos_shadow。
"""
function get_z_t_loos_shadow(
    shadows::Array{DenseShadow, 1},
    odd_shadows::Array{DenseShadow, 1},
    even_shadows::Array{DenseShadow, 1},
    uT::ITensor,
    tpos::Vector{Int},
    is_show_progress::Bool=true,
)
    return get_z_t_loos_shadow(
        reshape(shadows, length(shadows), 1),
        reshape(odd_shadows, length(odd_shadows), 1),
        reshape(even_shadows, length(even_shadows), 1),
        uT,
        tpos,
        is_show_progress,
    )
end
