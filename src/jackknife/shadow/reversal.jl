# ---------------------
# shadow reversal（Z_T 分子）的估计与 loos
# ---------------------

"""
由单实验自有 shadow 估计时间反演交叠 Z_T 分子（k=2 配对经恒等式
Σ_{a≠b}f(a,b) = F(S_A,S_BT) − Σ_a f(a,a) 化成 O(NU) 次缩并，不再枚举
C(NU,2) 对，数学精确一致），并算 jackknife 估计。

参数
- shadows::Vector{Shadow}：每设置的 shadow（已按 shots 平均），长度 n_ru。
- uT::ITensor：u_T 算符（I_1 区逐比特 σ^y）。
- tpos::Vector{Int}：做部分转置的 site 位置（I_1 区）。

关键词参数
- is_show_progress::Bool=false：是否显示进度。

返回
- θ：标量估计（无序对交叠的均值）。
- loos::Vector{Float64}：每个随机幺正的 leave-one-out jackknife 估计。
"""
function get_reversal_loos_shadow(
    shadows::Vector{Shadow},
    uT::ITensor,
    tpos::Vector{Int};
    is_show_progress=false,
)
    n_ru = length(shadows)
    @assert n_ru ≥ 3 "At least 3 random unitaries are required for reversal estimation."
    N = shadows[1].N
    ξ = shadows[1].site_indices
    Ud = dag(uT)
    shadows_T = [partial_transpose(sh, tpos) for sh in shadows]

    # 预配 deltas：任取一对做出 M3，按 id 把 plev 0 腿与 plev 4 腿配起来
    probe = shadows[1].rho * prime(uT, 1) * prime(shadows_T[2].rho, 2) * prime(Ud, 3)
    r0 = [ind for ind in inds(probe) if plev(ind) == 0]
    s4 = [ind for ind in inds(probe) if plev(ind) == 4]
    deltas = ITensor[]
    for r in r0
        s = only([x for x in s4 if id(x) == id(r)])
        push!(deltas, delta(r, s))
    end

    S_A = sum(sh.rho for sh in shadows)
    S_B = sum(sh.rho for sh in shadows_T)
    F = get_reversal_pair_value(Shadow(S_A, N, ξ), Shadow(S_B, N, ξ), uT, Ud, deltas)
    selfs = [get_reversal_pair_value(shadows[r], shadows_T[r], uT, Ud, deltas) for r in 1:n_ru]
    θ = (F - sum(selfs)) / (n_ru * (n_ru - 1))

    loos = Vector{Float64}(undef, n_ru)
    denom = (n_ru - 1) * (n_ru - 2)
    @showprogress desc="Reversal loos..." enabled=is_show_progress for i in 1:n_ru
        Fi = get_reversal_pair_value(
            Shadow(S_A - shadows[i].rho, N, ξ),
            Shadow(S_B - shadows_T[i].rho, N, ξ),
            uT, Ud, deltas,
        )
        loos[i] = (Fi - (sum(selfs) - selfs[i])) / denom
    end
    return θ, loos
end
