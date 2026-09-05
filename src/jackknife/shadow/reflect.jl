# ---------------------
# shadow reflect（一阶矩）的估计与 loos
# ---------------------
# 说明：一阶矩 loos 是通用的一阶矩工具，在本库中用于 reflect 期望
# （调用时 O 传 reflect_op），故归入 reflect。

"""
由 shadow 数据估计一阶矩（期望值），并算 jackknife 估计。

参数
- shadows::Array{<:AbstractShadow, 2}：(n_ru, n_m) 的 shadow 数据。

关键词参数
- O::Union{Nothing, MPO}=nothing：在 trace product 求值中用的 MPO 算符
  （估计期望值时用，本库中传 reflect_op 算 reflect）。
- is_show_progress::Bool=true：置换平均时是否显示进度。

返回
- θ：标量估计（平均期望值）。
- loos::Vector{Float64}：每个随机幺正的 leave-one-out jackknife 估计。
"""
function get_moment1_loos_shadow(
    shadows::Array{<:AbstractShadow, 2};
    O::Union{Nothing, MPO}=nothing,
    is_show_progress::Bool=true,
)
    n_ru, _ = size(shadows)
    @assert n_ru ≥ 2 "At least 2 random unitaries are required for 1-moment estimation."
    combs = collect(combinations(1:n_ru, 1))

    # 每个置换对测量求平均
    comb_avgs = get_comb_avgs_shadow(shadows, 1; O=O, is_show_progress=is_show_progress)
    θ = mean(comb_avgs)

    # jackknife loo 分组：不含幺正 i 的置换
    loos = get_combs_loos_shadow(
        n_ru,
        combs,
        comb_avgs,
    )

    return θ, loos
end

function get_moment1_loos_shadow(shadows::Vector{<:AbstractShadow}; kwargs...)
    return get_moment1_loos_shadow(reshape(shadows, length(shadows), 1); kwargs...)
end
