# ---------------------
# shadow purity（二阶矩）的估计与 loos
# ---------------------
# 说明：k=2 配对平均经恒等式
#   Σ_{a≠b}Tr(ρ̂_a ρ̂_b) = Tr(S²) − Σ_a Tr(ρ̂_a²)（S = Σ_a ρ̂_a 的按 shot 求和版）
# 化成 O(NU) 次迹乘，不再枚举 C(NU,2) 对，数学精确一致。

"""
由 shadow 数据估计二阶矩（纯度），并对随机幺正算 jackknife 估计。

参数
- shadows::Array{Shadow, 2}：(n_ru, n_m) 的 shadow 数据。

关键词参数
- is_compute_renyi::Bool=false：为 true 时返回 Rényi-2 熵估计
  （对平均纯度取 log2）；否则返回纯度。
- is_show_progress::Bool=true：置换平均时是否显示进度。

返回
- θ：标量估计（纯度或 Rényi-2，取决于 is_compute_renyi）。
- loos::Vector{Float64}：每个随机幺正的 leave-one-out jackknife 估计。
"""
function get_purity_loos_shadow(
    shadows::Array{Shadow, 2}; is_compute_renyi::Bool=false, is_show_progress::Bool=true
)
    n_ru, n_m = size(shadows)
    @assert n_ru ≥ 3 "At least 3 random unitaries are required for 2-moment estimation."
    N = shadows[1, 1].N
    ξ = shadows[1, 1].site_indices
    npair = n_ru * (n_ru - 1) * n_m^2

    # setting 内求和 S_r 与自配对和 D_r = Σ_{m,m'} Tr(ρ̂ρ̂')。
    S_rs = Vector{ITensor}(undef, n_ru)
    D_rs = Vector{Float64}(undef, n_ru)
    @showprogress desc="Setting sums..." enabled=is_show_progress for r in 1:n_ru
        acc = shadows[r, 1].rho
        for m in 2:n_m
            acc = acc + shadows[r, m].rho
        end
        S_rs[r] = acc
        d = 0.0
        for m in 1:n_m, mp in 1:n_m
            d += real(get_trace_product(shadows[r, m], shadows[r, mp]))
        end
        D_rs[r] = d
    end

    S = sum(S_rs)
    S_wrap = Shadow(S, N, ξ)
    # Tr(S²)：先自乘再求迹（trace(S) 只是归一和，不含配对信息）。
    θ_raw = (real(get_trace_product(S_wrap, S_wrap)) - sum(D_rs)) / npair
    avgfun(x) = is_compute_renyi ? (1 / (1 - 2)) * log2(x) : x
    θ = avgfun(θ_raw)

    # loos：去掉 setting i 后同式（与旧 combs 口径一致：不含 i 的配对均值）。
    loos = Vector{Float64}(undef, n_ru)
    denom = (n_ru - 1) * (n_ru - 2) * n_m^2
    @showprogress desc="Loos..." enabled=is_show_progress for i in 1:n_ru
        Si = Shadow(S - S_rs[i], N, ξ)
        loos[i] = (real(get_trace_product(Si, Si)) - (sum(D_rs) - D_rs[i])) / denom
    end

    return θ, loos
end

function get_purity_loos_shadow(shadows::Vector{Shadow}; kwargs...)
    return get_purity_loos_shadow(reshape(shadows, length(shadows), 1); kwargs...)
end
