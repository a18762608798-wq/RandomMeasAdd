# ---------------------
# Robust Hamming 核估计子：把上游固定核 (2/-1) 换成逐单元标定核 (x̂/ŷ)。
# ---------------------
# 缩并链与上游 get_overlap(MeasurementProbability, MeasurementProbability)
# 完全一致，唯一区别是每比特的核张量对角取 x̂、非对角取 ŷ，
# 且不再乘全局 2^N（权重已含在 x̂、ŷ 里）。

"""
由核权重构造单比特权重张量：对角 x̂，非对角 ŷ。
等值比特（无论 0/1 侧）走对角，不等走非对角，故与 0/1↔1/2 的标号约定无关。
"""
function _get_weight_tensor(s::Index, s_prime::Index, x::Float64, y::Float64)
    W = ITensor(Float64, s, s_prime)
    W[s => 1, s_prime => 1] = x
    W[s => 2, s_prime => 2] = x
    W[s => 1, s_prime => 2] = y
    W[s => 2, s_prime => 1] = y
    return W
end

"""
get_overlap_robust(data1, data2, xs, ys; apply_bias_correction=false)

单 setting 的标定核交叠：Σ_{s,s'} Π_i g_i(s_i,s'_i) P₁(s)P₂(s')，
经逐比特权重张量缩并求值（与上游同复杂度）。

bias 修正（purity 用）：对角项 K(s,s) = Πx̂ 与 s 无关，
故修正量为 (Πx̂)/(NM-1)，即 raw*NM²/(NM*(NM-1)) - Πx̂/(NM-1)。
标准核 x̂=2 时退化为上游公式（2^N 项）。
"""
function get_overlap_robust(
    data1::MeasurementData,
    data2::MeasurementData,
    xs::AbstractVector{Float64},
    ys::AbstractVector{Float64};
    apply_bias_correction::Bool=false,
)
    N = data1.N
    @assert N == data2.N "两份数据的比特数必须一致。"
    @assert length(xs) == length(ys) == N "核权重长度须等于比特数。"

    prob1 = MeasurementProbability(data1)
    prob2 = MeasurementProbability(data2)
    xi1 = prob1.site_indices
    xi2 = prob2.site_indices

    overlap = prob1.measurement_probability
    for i in 1:N
        overlap *= _get_weight_tensor(xi1[i], prime(xi2[i]), xs[i], ys[i])
    end
    raw = real(scalar(overlap * prob2.measurement_probability'))

    if apply_bias_correction
        @assert data1.NM == data2.NM "bias 修正要求两份数据的 shots 数一致。"
        NM = data1.NM
        @assert NM ≥ 2 "bias 修正要求每 setting 至少 2 shots。"
        return raw * NM^2 / (NM * (NM - 1)) - prod(xs) / (NM - 1)
    end
    return raw
end

"""
get_reflect_robust(data, as, bs, xs, ys)

单 setting 的标定核反射估计：(1/NM) Σ_m Π_p g_p(s_{a_p}, s_{b_p})，
镜像对逐发配对（Z_R 是一阶量，无 shot 配对、无 bias 修正）。
xs/ys 按对编号（非比特编号）。
"""
function get_reflect_robust(
    data::MeasurementData,
    as::AbstractVector{Int},
    bs::AbstractVector{Int},
    xs::AbstractVector{Float64},
    ys::AbstractVector{Float64},
)
    np = length(as)
    @assert length(bs) == np "镜像对两侧数量必须一致。"
    @assert length(xs) == length(ys) == np "核权重长度须等于对数。"
    res = data.measurement_results
    nm = data.NM
    ssum = 0.0
    for m in 1:nm
        prod = 1.0
        for p in 1:np
            prod *= res[m, as[p]] == res[m, bs[p]] ? xs[p] : ys[p]
        end
        ssum += prod
    end
    return ssum / nm
end
