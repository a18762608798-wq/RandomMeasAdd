# ---------------------
# shadow reversal（Z_T 分子）的估计与 loos
# ---------------------

"""
由单实验 dense shadow 估计时间反演交叠 Z_T 分子，并算 jackknife 估计。

参数
- shadows::Vector{DenseShadow}：每设置的 dense shadow（已按 shots 平均），长度 n_ru。
- uT::ITensor：u_T 算符（I_1 区逐比特 σ^y）。
- tpos::Vector{Int}：做部分转置的 site 位置（I_1 区）。

关键词参数
- is_show_progress::Bool=false：是否显示进度。

返回
- θ：标量估计（无序对交叠的均值）。
- loos::Vector{Float64}：每个随机幺正的 leave-one-out jackknife 估计。
"""
function get_reversal_loos_shadow(
    shadows::Vector{DenseShadow},
    uT::ITensor,
    tpos::Vector{Int};
    is_show_progress=false,
)
    n_ru = length(shadows)
    @assert n_ru ≥ 3 "At least 3 random unitaries are required for reversal estimation."
    combs, avgs = get_reversal_comb_avgs_shadow(
        shadows, uT, tpos; is_show_progress=is_show_progress
    )
    θ = mean(avgs)
    loos = get_combs_loos_shadow(n_ru, combs, avgs)
    return θ, loos
end
