# ---------------------
# shadow reflect（一阶矩）的估计与 loos
# ---------------------
# 说明：一阶矩 loos 是通用的一阶矩工具，在本库中用于 reflect 期望
# （调用时 O 传 reflect_op），故归入 reflect。

"""
由 shadow 数据估计一阶矩（期望值），并算 jackknife 估计。

参数
- shadows::Array{Shadow, 2}：(n_ru, n_m) 的 shadow 数据。

关键词参数
- O::Union{Nothing, MPO}=nothing：在 trace product 求值中用的 MPO 算符
  （估计期望值时用，本库中传 reflect_op 算 reflect）。
- is_show_progress::Bool=true：置换平均时是否显示进度。

返回
- θ：标量估计（平均期望值）。
- loos::Vector{Float64}：每个随机幺正的 leave-one-out jackknife 估计。
"""
function get_moment1_loos_shadow(
    shadows::Array{Shadow, 2};
    O::Union{Nothing, MPO}=nothing,
    is_show_progress::Bool=true,
)
    n_ru, n_m = size(shadows)
    @assert n_ru ≥ 2 "At least 2 random unitaries are required for 1-moment estimation."
    # k=1：一阶矩可直接加和（Σ_r T_r），loos 逐个摘除，与旧组合口径一致。
    T_rs = Vector{Float64}(undef, n_ru)
    @showprogress desc="Moments..." enabled=is_show_progress for r in 1:n_ru
        s = 0.0
        for m in 1:n_m
            s += real(get_trace_product(shadows[r, m]; O))
        end
        T_rs[r] = s / n_m
    end
    θ = sum(T_rs) / n_ru
    loos = [(sum(T_rs) - T_rs[i]) / (n_ru - 1) for i in 1:n_ru]

    return θ, loos
end

function get_moment1_loos_shadow(shadows::Vector{Shadow}; kwargs...)
    return get_moment1_loos_shadow(reshape(shadows, length(shadows), 1); kwargs...)
end
