# ---------------------
# shadow purity（二阶矩）的估计与 loos
# ---------------------

"""
由 shadow 数据估计二阶矩（纯度），并对随机幺正算 jackknife 估计。

参数
- shadows::Array{<:AbstractShadow, 2}：(n_ru, n_m) 的 shadow 数据。

关键词参数
- is_compute_renyi::Bool=false：为 true 时返回 Rényi-2 熵估计
  （对平均纯度取 log2）；否则返回纯度。
- is_show_progress::Bool=true：置换平均时是否显示进度。

返回
- θ：标量估计（纯度或 Rényi-2，取决于 is_compute_renyi）。
- loos::Vector{Float64}：每个随机幺正的 leave-one-out jackknife 估计。
"""
function get_purity_loos_shadow(
    shadows::Array{<:AbstractShadow, 2}; is_compute_renyi::Bool=false, is_show_progress::Bool=true
)
    n_ru, n_m = size(shadows)
    @assert n_ru ≥ 3 "At least 3 random unitaries are required for 2-moment estimation."
    # 预先枚举置换
    combs = collect(combinations(1:n_ru, 2))

    # 每个置换对测量求平均
    comb_avgs = get_comb_avgs_shadow(shadows, 2; is_show_progress=is_show_progress)

    # 定义平均泛函
    avgfun(x) = is_compute_renyi ? (1 / (1 - 2)) * log2(mean(x)) : mean(x)

    θ = avgfun(comb_avgs)

    # jackknife loo 分组：不含幺正 i 的置换
    loos = get_combs_loos_shadow(
        n_ru,
        combs,
        comb_avgs,
    )

    return θ, loos
end

function get_purity_loos_shadow(shadows::Vector{<:AbstractShadow}; kwargs...)
    return get_purity_loos_shadow(reshape(shadows, length(shadows), 1); kwargs...)
end
