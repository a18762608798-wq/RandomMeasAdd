# ---------------------
# shadow 基础组件：置换平均 + leave-one-out
# ---------------------

"""
对随机幺正的所有 k 元置换求 trace-product 平均，再对测量指标求平均。

参数
- shadows::Array{<:AbstractShadow, 2}：(n_ru, n_m) 的 shadow 数组，
  n_ru 为随机幺正个数，n_m 为每个幺正下的测量次数。
- k::Int：置换大小（只支持 1 或 2）。

关键词参数
- O::Union{Nothing, MPO}=nothing：传给 get_trace_product 的 MPO 算符
  （仅 k == 1 时允许）。
- is_show_progress::Bool=true：是否显示进度。

返回
- Array{Float64}，长度 choose(n_ru, k)：每个 k 置换在测量笛卡尔积上
  平均后的 trace-product 均值。
"""
function get_comb_avgs_shadow(
    shadows::Array{<:AbstractShadow, 2},
    k::Int;
    O::Union{Nothing, MPO}=nothing,
    is_show_progress=false,
)
    @assert (k == 2 && isnothing(O)) || (k == 1) "O must be nothing for k=2, and k in the range of [1, 2]"

    # 预先枚举置换和测量笛卡尔积
    n_ru, n_m = size(shadows)
    combs = collect(combinations(1:n_ru, k))
    cprod = collect(CartesianIndices(ntuple(_ -> 1:n_m, k)))
    combs_num = length(combs)

    # 每个置换对测量求平均
    comb_avgs = zeros(Float64, combs_num)
    @showprogress desc="Combinations Processing..." enabled=is_show_progress @threads for pidx in eachindex(combs)
        r = combs[pidx]
        ssum = 0.0
        for m in cprod
            ssum += real(get_trace_product((shadows[r[i], m[i]] for i in 1:k)...; O))
        end
        comb_avgs[pidx] = ssum / length(cprod)
    end

    return comb_avgs
end

"""
get_combs_loos_shadow(n_ru, combs, avgs)

由置换平均算 leave-one-out 的 jackknife 估计。

参数
- n_ru::Int64：随机幺正总数。
- combs::Vector{Vector{Int64}}：幺正指标的 k 元组合向量。
- avgs::Vector{Float64}：每个组合的 trace-product 平均。

返回
- loos::Vector{Float64}：每个随机幺正的 leave-one-out jackknife 估计。

说明
幺正 i 的 leave-one-out 值，等于所有不含 i 的组合的置换平均之均值。
"""
function get_combs_loos_shadow(
    n_ru::Int64,
    combs::Vector{Vector{Int64}},
    avgs::Vector{Float64},
)
    k = length(first(combs))
    ssum = sum(avgs) # 平均之和
    incident_threads = zeros(Float64, n_ru, Threads.maxthreadid()) # 含每个随机幺正的组合之和
    # 求 leave-one-out。
    @threads for idx in eachindex(combs)
        tid = threadid()
        incident_indices = combs[idx] # 记下被关联的项。
        incident_val = avgs[idx]
        for index in incident_indices
            incident_threads[index, tid] += incident_val
        end
    end
    incidents = vec(sum(incident_threads; dims=2))
    denom = binomial(n_ru - 1, k)
    loos = (ssum .- incidents) ./ denom # 得 loos

    return loos
end

"""
get_reversal_pair_value(shadow_a, shadow_b_T, U, Ud, deltas)

算一对 shadow 的时间反演交叠值 Tr(A·U·B_T·U†)（已验证过的缩并链）。

参数
- shadow_a::DenseShadow：第一个设置的 shadow（未转置）。
- shadow_b_T::DenseShadow：第二个设置的 shadow（已在 I_1 区做部分转置）。
- U::ITensor：u_T 算符（I_1 区逐比特 σ^y）。
- Ud::ITensor：U 的 dagger（`dag(U)`，prime 层已验证）。
- deltas::Vector{ITensor}：预先配好的 δ 张量（按 id 把 plev 0 与 plev 4 腿粘起来求迹）。

返回
- val::Float64：该有序对的估计值（取实部）。

说明
缩并链：M1 = A·U(+1) → M2 = M1·Bt(+2) → M3 = M2·Ud(+3)，
最后用 deltas 缩掉首尾腿。括号内为 prime 层提升数。
"""
function get_reversal_pair_value(
    shadow_a::DenseShadow,
    shadow_b_T::DenseShadow,
    U::ITensor,
    Ud::ITensor,
    deltas::Vector{ITensor},
)
    A = shadow_a.shadow_data
    Bt = shadow_b_T.shadow_data
    M1 = A * prime(U, 1)
    M2 = M1 * prime(Bt, 2)
    M3 = M2 * prime(Ud, 3)
    acc = M3
    for d in deltas
        acc = acc * d
    end
    return real(acc[])
end

"""
get_reversal_comb_avgs_shadow(shadows, uT, tpos; is_show_progress)

对随机幺正的所有 2 元组合算时间反演交叠平均（shadow 版 Z_T 分子）。

参数
- shadows::Vector{DenseShadow}：每设置的 dense shadow（已按 shots 平均）。
- uT::ITensor：u_T 算符。
- tpos::Vector{Int}：做部分转置的 site 位置（I_1 区，重排 frame 下的奇位）。

关键词参数
- is_show_progress::Bool=false：是否显示进度。

返回
- combs：`combinations(1:n_ru, 2)` 的组合向量。
- avgs::Vector{Float64}：每个无序对的交叠平均（两 orientation 的均值）。

说明
每个无序对 {a,b} 的值 = (f(a,b) + f(b,a))/2，其中 f(a,b) 把转置放在
第二个上。转置 shadow 预先算好（NU 次），deltas 预先配好（1 次）。
"""
function get_reversal_comb_avgs_shadow(
    shadows::Vector{DenseShadow},
    uT::ITensor,
    tpos::Vector{Int};
    is_show_progress=false,
)
    n_ru = length(shadows)
    @assert n_ru ≥ 3 "At least 3 random unitaries are required for reversal estimation."
    Ud = dag(uT)
    shadows_T = [partial_transpose(sh, tpos) for sh in shadows]

    # 预配 deltas：任取一对做出 M3，按 id 把 plev 0 腿与 plev 4 腿配起来
    probe = shadows[1].shadow_data * prime(uT, 1) * prime(shadows_T[2].shadow_data, 2) * prime(Ud, 3)
    r0 = [ind for ind in inds(probe) if plev(ind) == 0]
    s4 = [ind for ind in inds(probe) if plev(ind) == 4]
    deltas = ITensor[]
    for r in r0
        s = only([x for x in s4 if id(x) == id(r)])
        push!(deltas, delta(r, s))
    end

    combs = collect(combinations(1:n_ru, 2))
    avgs = Vector{Float64}(undef, length(combs))
    @showprogress desc="Reversal pairs..." enabled=is_show_progress @threads for pidx in eachindex(combs)
        a, b = combs[pidx]
        fwd = get_reversal_pair_value(shadows[a], shadows_T[b], uT, Ud, deltas)
        bwd = get_reversal_pair_value(shadows[b], shadows_T[a], uT, Ud, deltas)
        avgs[pidx] = (fwd + bwd) / 2
    end

    return combs, avgs
end
