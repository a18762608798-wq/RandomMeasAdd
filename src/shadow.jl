# ---------------------
# 自有 Shadow 类型与代数：不依赖上游 RandomMeas。
# ---------------------
# 数学与上游 DenseShadow 逐项一致（Elben/Vermersch 经典影子）：
# - 单比特影子核：ρ̂_i = α·ψ'ψ† + β·δ，ψ = U†|s>，α=3/(2G-1)，β=(G-2)/(2G-1)；
# - setting 影子：经验分布 P（直方图权重直写，无需 2^N countmap 中间量）
#   经 Hamming 张量 + U/U† 变换，ρ̂ = T(P)；
# - 代数：multiply（trace 积）、trace、get_expect_shadow、partial_transpose。
# 约定：outcome 值 {1,2}（bit+1），与旧 MeasurementProbability 一致。

"""自有经典影子：单 setting 的影子算符（2N 条腿 ITensor）。"""
struct Shadow
    rho::ITensor
    N::Int
    site_indices::Vector{Index{Int64}}
    function Shadow(rho::ITensor, N::Int, site_indices::Vector{Index{Int64}})
        @assert length(site_indices) == N "site 数须为 N。"
        inds_all = inds(rho)
        @assert length(inds_all) == 2 * N "影子 ITensor 须有 2N 条腿。"
        unprimed = [i for i in inds_all if plev(i) == 0]
        primed = [i for i in inds_all if plev(i) == 1]
        @assert length(unprimed) == N && length(primed) == N "prime 层不对。"
        @assert Set(unprimed) == Set(site_indices) "unprimed 腿须匹配 site_indices。"
        @assert Set(primed) == Set(map(prime, site_indices)) "primed 腿不匹配。"
        new(rho, N, site_indices)
    end
end

"""单比特 Hamming 张量（上游 get_h_tensor 同式）：对角 (α+β)/2，非对角 β/2。"""
function get_h_tensor(s::Index, s_prime::Index, G::Float64=1.0)
    α = 3.0 / (2.0 * G - 1.0)
    β = (G - 2.0) / (2.0 * G - 1.0)
    H = ITensor(Float64, s, s_prime)
    H[s => 1, s_prime => 1] = (α + β) / 2.0
    H[s => 2, s_prime => 2] = (α + β) / 2.0
    H[s => 1, s_prime => 2] = β / 2.0
    H[s => 2, s_prime => 1] = β / 2.0
    return H
end

"""
Shadow(d, U, ξ, G) -> Shadow

由单个 setting 的直方图构造影子：经验分布 P[s] = w/NM（{1,2} 取值），
ρ = 2^N·P 经逐比特 Hamming 张量 + U/U† 变换（与上游 DenseShadow 同链）。
- d::HistData：该 setting 的直方图（bit 0 ↔ 物理 0，取值时 +1）。
- U::Vector{Matrix{ComplexF64}}：每比特测量幺正（qmeas settings 同矩阵）。
- ξ：site index 向量；G：读出缓解向量（默认全 1）。
"""
function Shadow(
    d::HistData,
    U::Vector{Matrix{ComplexF64}},
    ξ::Vector{Index{Int64}},
    G::Vector{Float64}=ones(Float64, d.n),
)
    N = d.n
    @assert length(U) == N && length(ξ) == N && length(G) == N "U/ξ/G 长度须为 N。"
    # 经验分布经稠密数组转 ITensor（只填 uniques，与上游 countmap 同值）。
    P = zeros(ComplexF64, fill(2, N)...)
    for (k, w) in zip(d.keys, d.counts)
        P[map(j -> hist_bit(k, j) + 1, 1:N)...] = w / d.nm
    end
    rho = 2^N * ITensor(P, ξ...)
    for i in 1:N
        s = ξ[i]
        rho *= get_h_tensor(s, s'', G[i])
        rho *= δ(s, s', s'')
        # uit 腿序 (s', s)，与上游 basis_transformation 一致。
        ut = ITensor(Matrix(U[i]), prime(s), s) * δ(s'', s) * δ(s, s')
        rho = mapprime(ut * rho, 2, 0)
        utd = dag(ITensor(Matrix(U[i]), prime(s), s)) * δ(s'', s)
        rho = mapprime(utd * rho, 2, 1)
    end
    return Shadow(rho, N, ξ)
end

"""
get_dense_shadows(h, ξ; G) -> Vector{Shadow}

每 setting 一个影子（已按 shots 平均），与旧 `get_dense_shadows(group; G)`
同语义。子系统用约化后的 HistGroup + 对应 ξ 子集分别调。
"""
function get_dense_shadows(
    h::HistGroup,
    ξ::Vector{Index{Int64}};
    G::Vector{Float64}=ones(Float64, h.N),
)
    @assert length(ξ) == h.N && length(G) == h.N "ξ/G 长度须为 N。"
    return map(h.measurements, 1:h.NU) do d, s
        U = [Matrix(h.settings[s, j, :, :]) for j in 1:h.N]
        Shadow(d, U, ξ, G)
    end
end

"""两个影子的 trace 积影子（上游 DenseShadow multiply 同式）。"""
function multiply(shadow1::Shadow, shadow2::Shadow)::Shadow
    @assert shadow1.N == shadow2.N "比特数不一致。"
    @assert shadow1.site_indices == shadow2.site_indices "site 不一致。"
    product = mapprime(shadow1.rho * prime(shadow2.rho), 2, 1)
    return Shadow(product, shadow1.N, shadow1.site_indices)
end

"""影子求迹 Tr[ρ̂]（上游 DenseShadow trace 同式）。"""
function trace(shadow::Shadow)
    t = copy(shadow.rho)
    for i in 1:shadow.N
        t *= δ(shadow.site_indices[i], prime(shadow.site_indices[i]))
    end
    return scalar(t)
end

"""MPO 期望 Tr[O·ρ̂]（上游 DenseShadow get_expect_shadow 同式）。"""
function get_expect_shadow(O::MPO, shadow::Shadow)
    N = shadow.N
    ξ = shadow.site_indices
    X = 1 * shadow.rho'
    for i in 1:N
        s = ξ[i]
        X *= O[i] * δ(s, s'')
    end
    return X[]
end

"""
get_trace_product(shadows...; O) -> 标量

影子乘积求迹（O==nothing）或 MPO 期望（给 O），与上游同语义。
"""
function get_trace_product(shadows::Shadow...; O::Union{Nothing,MPO}=nothing)
    result = shadows[1]
    for shadow in shadows[2:end]
        result = multiply(result, shadow)
    end
    if O === nothing
        return trace(result)
    else
        return get_expect_shadow(O, result)
    end
end

"""影子部分转置（上游 DenseShadow partial_transpose 同式，返回新对象）。"""
function partial_transpose(shadow::Shadow, subsystem::Vector{Int})::Shadow
    @assert all(i -> 1 ≤ i ≤ shadow.N, subsystem) "子系统越界。"
    @assert length(unique(subsystem)) == length(subsystem) "子系统须无重复。"
    A = shadow.rho
    for i in subsystem
        a = shadow.site_indices[i]
        A = swapind(A, a, prime(a))
    end
    return Shadow(A, shadow.N, shadow.site_indices)
end

"""
modified_get_calibration_vector(tgroup; is_show_progress) -> Vector{Float64}

由 trivial（|0>）直方图估计读出缓解向量 G（Vitale et al, PRXQ 2025，
与旧实现同式）：逐比特、逐 setting 算实测 Born 分布与理想 |0> 分布
（P_expected[s] = |U[s,1]|²，直接取 settings 矩阵第一列）的互关联，
G[i] += 3*(cross-self)/NU + 1/NU。
"""
function modified_get_calibration_vector(
    tgroup::HistGroup; is_show_progress=false,
)
    N = tgroup.N
    NU = tgroup.NU
    G_e = zeros(Float64, N)
    @showprogress dt=1 enabled=is_show_progress for i in 1:N
        red = reduce_to_subsystem(tgroup, [i])
        for (s_idx, d) in enumerate(red.measurements)
            c0 = 0
            for (k, w) in zip(d.keys, d.counts)
                c0 += (hist_bit(k, 1) == 0) * w
            end
            pm0 = c0 / d.nm
            Pm = [pm0, 1 - pm0]
            U = tgroup.settings[s_idx, i, :, :]
            Pe = [abs2(U[1, 1]), abs2(U[2, 1])]
            cross_corr = Pm[1] * Pe[1] + Pm[2] * Pe[2]
            self_corr = Pe[1]^2 + Pe[2]^2
            G_e[i] += 3 * (cross_corr - self_corr) / NU + 1 / NU
        end
    end
    return G_e
end
