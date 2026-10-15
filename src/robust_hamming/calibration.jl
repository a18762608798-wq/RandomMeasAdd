# ---------------------
# Robust Hamming 标定：由 trivial（|0>）直方图估计每单元的核权重 (x̂, ŷ)。
# ---------------------
# 物理背景见笔记 RES-Robust_Hamming.md：
# 三个量共用同一套单比特核形状 ĝ(s,s') = x̂（同）/ ŷ（异），
# 差别只在“哪两个比特配对”与标定量 R̂、f̂0 怎么数：
# - purity：同比特、不同 shot 配对（k≠k'），一支 trivial 数据，δ 两条腿相同；
# - Z_R：同发内镜像对 (i, ī) 配对，一支 trivial 数据，两点通道各一个 δ；
# - Z_T：两份实验同一比特配对，两支 trivial 数据（exp1 带 u_T），每实验各一个 δ。
# 无噪声极限 R̂=2/3、δ̂=0 给出 x̂=2、ŷ=-1，退化为标准 Hamming 核。
#
# 实现直接读直方图权重（精确改写，不展开）：bit 0 ↔ 物理 0。

"""
get_f0(group, cols) -> Vector{Float64}

trivial 组指定列的 0 频率 f̂0（全部 setting × shot 平均），每列一个。
δ̂ = 2f̂0 - 1 由调用方组装。
"""
function get_f0(group::HistGroup, cols::AbstractVector{Int}=1:group.N)
    f0 = zeros(Float64, length(cols))
    total = 0
    for d in group.measurements
        total += d.nm
        for (j, c) in enumerate(cols)
            s = 0
            for (k, w) in zip(d.keys, d.counts)
                s += (hist_bit(k, c) == 0) * w
            end
            f0[j] += s
        end
    end
    return f0 ./ total
end

"""
get_R_purity(group, cols) -> Vector{Float64}

purity 配对的符合率 R̂：同列、不同 shot (k≠k') 的符合分数，
先每 setting 按 matches/(NM*(NM-1)) 算，再对 setting 平均。
"""
function get_R_purity(group::HistGroup, cols::AbstractVector{Int}=1:group.N)
    nus = group.measurements
    nu = length(nus)
    R = zeros(Float64, length(cols))
    for d in nus
        nm = d.nm
        @assert nm ≥ 2 "purity 标定要求每 setting 至少 2 shots。"
        denom = nm * (nm - 1)
        for (j, c) in enumerate(cols)
            c0 = 0
            for (k, w) in zip(d.keys, d.counts)
                c0 += (hist_bit(k, c) == 0) * w
            end
            c1 = nm - c0
            R[j] += (c0 * (c0 - 1) + c1 * (c1 - 1)) / denom
        end
    end
    return R ./ nu
end

"""
get_R_mirror(group, as, bs) -> Vector{Float64}

Z_R 配对的符合率 R̂：同发内镜像对 (as[p], bs[p]) 的符合分数，
对全部 (setting, shot) 平均。每对一个。
"""
function get_R_mirror(group::HistGroup, as::AbstractVector{Int}, bs::AbstractVector{Int})
    @assert length(as) == length(bs) "镜像对两侧数量必须一致。"
    nus = group.measurements
    np = length(as)
    R = zeros(Float64, np)
    total = 0
    for d in nus
        total += d.nm
        for (k, w) in zip(d.keys, d.counts)
            for p in 1:np
                R[p] += (hist_bit(k, as[p]) == hist_bit(k, bs[p])) * w
            end
        end
    end
    return R ./ total
end

"""
get_R_cross(group1, group2, cols) -> Vector{Float64}

Z_T 配对的符合率 R̂：两份实验同列、全部 (k,k') 配对的符合分数，
每 setting 按 (c0₁c0₂+c1₁c1₂)/NM² 算，再对 setting 平均。
两份数据的 NM 必须一致（import_hist_pair 已保证）。
plain trivial 的 R̂ 正常应在 2/3 附近；若远低于 1/2（如 ~1/3），
说明 trivial_1 混进了 u_T（|1> 态），标定无意义，此处直接告警。
"""
function get_R_cross(group1::HistGroup, group2::HistGroup, cols::AbstractVector{Int}=1:group1.N)
    nus1 = group1.measurements
    nus2 = group2.measurements
    @assert length(nus1) == length(nus2) "两份实验的设置数必须一致。"
    nu = length(nus1)
    R = zeros(Float64, length(cols))
    for (d1, d2) in zip(nus1, nus2)
        nm = d1.nm
        @assert nm == d2.nm "两份实验的每设置 shots 数必须一致。"
        for (j, c) in enumerate(cols)
            c0a = 0
            c0b = 0
            for (k, w) in zip(d1.keys, d1.counts)
                c0a += (hist_bit(k, c) == 0) * w
            end
            for (k, w) in zip(d2.keys, d2.counts)
                c0b += (hist_bit(k, c) == 0) * w
            end
            R[j] += (c0a * c0b + (nm - c0a) * (nm - c0b)) / nm^2
        end
    end
    R ./= nu
    if any(R .< 0.45)
        @warn "跨实验符合率 R̂ 远低于正常值（min=$(minimum(R))，正常 ~2/3）：trivial_1 可能混进了 u_T（|1> 态给出 ~1/3），请检查标定电路。" maxlog = 1
    end
    return R
end

"""
get_xy(R, da, db) -> (x, y)

由符合率 R̂ 与两条腿的 δ̂ 组装核权重：
ε̂ = da*db，Â = 6R̂-2-2ε̂，q̂ = 3/(Â-1-ε̂)，p̂ = 1-ε̂q̂，
x̂ = (p̂+q̂)/2，ŷ = (p̂-q̂)/2。
"""
function get_xy(R::Float64, da::Float64, db::Float64)
    eps = da * db
    A = 6R - 2 - 2eps
    denom = A - 1 - eps
    @assert !iszero(denom) "标定分母 Â-1-ε̂ 为零（R̂=$(R)，δ=($(da),$(db))）：trivial 数据退化，无从标定。"
    q = 3 / denom
    p = 1 - eps * q
    return (p + q) / 2, (p - q) / 2
end

"""
get_purity_xy(triv, cols) -> (xs, ys)

purity 标定的成套核权重：δ 两条腿相同（δ̂=2f̂0-1），R̂ 取跨 shot 符合率。
"""
function get_purity_xy(triv::HistGroup, cols::AbstractVector{Int}=1:triv.N)
    f0 = get_f0(triv, cols)
    R = get_R_purity(triv, cols)
    n = length(cols)
    xs = Vector{Float64}(undef, n)
    ys = Vector{Float64}(undef, n)
    for j in 1:n
        d = 2 * f0[j] - 1
        xs[j], ys[j] = get_xy(R[j], d, d)
    end
    return xs, ys
end

"""
get_mirror_xy(triv, as, bs) -> (xs, ys)

Z_R 标定的成套核权重：δ 按对内左右两侧各自的 f̂0 估计，R̂ 取同发镜像符合率。
"""
function get_mirror_xy(triv::HistGroup, as::AbstractVector{Int}, bs::AbstractVector{Int})
    fa = get_f0(triv, as)
    fb = get_f0(triv, bs)
    R = get_R_mirror(triv, as, bs)
    n = length(as)
    xs = Vector{Float64}(undef, n)
    ys = Vector{Float64}(undef, n)
    for p in 1:n
        xs[p], ys[p] = get_xy(R[p], 2 * fa[p] - 1, 2 * fb[p] - 1)
    end
    return xs, ys
end

"""
get_cross_xy(triv1, triv2, cols) -> (xs, ys)

Z_T 标定的成套核权重：δ₁/δ₂ 取两份实验各自的 f̂0，R̂ 取跨实验符合率。
"""
function get_cross_xy(triv1::HistGroup, triv2::HistGroup, cols::AbstractVector{Int}=1:triv1.N)
    f1 = get_f0(triv1, cols)
    f2 = get_f0(triv2, cols)
    R = get_R_cross(triv1, triv2, cols)
    n = length(cols)
    xs = Vector{Float64}(undef, n)
    ys = Vector{Float64}(undef, n)
    for j in 1:n
        xs[j], ys[j] = get_xy(R[j], 2 * f1[j] - 1, 2 * f2[j] - 1)
    end
    return xs, ys
end
