# ---------------------
# Robust Hamming 核估计子：把上游固定核 (2/-1) 换成逐单元标定核 (x̂/ŷ)。
# ---------------------
# 直接读直方图权重（精确改写，不展开）：单 setting 的交叠
# Σ_{a,b} w_a w_b Π_i g_i(s^a_i, s^b_i) / (NM₁NM₂)，g 取 x̂（同）/ ŷ（异）；
# 标准核 x̂=2、ŷ=-1 时 Πg = 2^N·(-2)^{-D}，与上游张量缩并数学一致。
# 本文件不再用 ITensor（Phase 1 hamming 路零上游类型）。

"""
单 setting 标定核交叠：Σ_{a,b} w_a w_b Π_i g_i / (NM₁NM₂)，
g_i 取 xs[i]（同）/ ys[i]（异）。

bias 修正（purity 用）：对角项 K(s,s) = Πx̂ 与 s 无关，
故修正量为 (Πx̂)/(NM-1)，即 raw*NM²/(NM*(NM-1)) - Πx̂/(NM-1)。
标准核 x̂=2 时退化为上游公式（2^N 项）。
"""
function get_overlap_robust(
    data1::HistData,
    data2::HistData,
    xs::AbstractVector{Float64},
    ys::AbstractVector{Float64};
    apply_bias_correction::Bool=false,
)
    N = data1.n
    @assert N == data2.n "两份数据的比特数必须一致。"
    @assert length(xs) == length(ys) == N "核权重长度须等于比特数。"
    nm1, nm2 = data1.nm, data2.nm
    raw = 0.0
    for (ka, wa) in zip(data1.keys, data1.counts)
        for (kb, wb) in zip(data2.keys, data2.counts)
            prod = 1.0
            for i in 1:N
                prod *= hist_bit(ka, i) == hist_bit(kb, i) ? xs[i] : ys[i]
            end
            raw += wa * wb * prod
        end
    end
    raw /= nm1 * nm2

    if apply_bias_correction
        @assert nm1 == nm2 "bias 修正要求两份数据的 shots 数一致。"
        NM = nm1
        @assert NM ≥ 2 "bias 修正要求每 setting 至少 2 shots。"
        return raw * NM^2 / (NM * (NM - 1)) - prod(xs) / (NM - 1)
    end
    return raw
end

"""
get_overlap(d1, d2; apply_bias_correction)（HistData 派发，显式扩展上游同名函数）

单 setting 标准核交叠（x̂=2、ŷ=-1），数学同上游
`get_overlap(::MeasurementData, ::MeasurementData)`。
"""
function RandomMeas.get_overlap(
    data1::HistData,
    data2::HistData;
    apply_bias_correction::Bool=false,
)
    N = data1.n
    return get_overlap_robust(
        data1, data2, fill(2.0, N), fill(-1.0, N);
        apply_bias_correction,
    )
end

"""
get_reflect_hamming(data)

单个 setting 的 Hamming 距离反射估计：
(1/NM) Σ_u w_u 2^np (-2)^{-D(u)}，D 为相邻奇偶对内差异数，
系统按相邻奇偶 site 配成对（与旧实现同式）。
"""
function get_reflect_hamming(data::HistData)
    N = data.n
    np = N ÷ 2
    @assert 2 * np == N "reflect 要求偶数比特。"
    ssum = 0.0
    for (k, w) in zip(data.keys, data.counts)
        D = 0
        for p in 1:np
            D += hist_bit(k, 2p - 1) != hist_bit(k, 2p)
        end
        ssum += w * 2.0^np * (-2.0)^(-D)
    end
    return ssum / data.nm
end

"""
get_reflect_robust(data, as, bs, xs, ys)

单 setting 的标定核反射估计：(1/NM) Σ_u w_u Π_p g_p，
镜像对逐发配对（Z_R 是一阶量，无 shot 配对、无 bias 修正）。
xs/ys 按对编号（非比特编号）。
"""
function get_reflect_robust(
    data::HistData,
    as::AbstractVector{Int},
    bs::AbstractVector{Int},
    xs::AbstractVector{Float64},
    ys::AbstractVector{Float64},
)
    np = length(as)
    @assert length(bs) == np "镜像对两侧数量必须一致。"
    @assert length(xs) == length(ys) == np "核权重长度须等于对数。"
    ssum = 0.0
    for (k, w) in zip(data.keys, data.counts)
        prod = 1.0
        for p in 1:np
            prod *= hist_bit(k, as[p]) == hist_bit(k, bs[p]) ? xs[p] : ys[p]
        end
        ssum += w * prod
    end
    return ssum / data.nm
end

"""
get_reversal_hamming(data1, data2)（HistData 派发）

单个配对设置的互关联（Z_T 分子）：标准核交叠，不加 bias 修正
（两份实验 shots 独立，无对角自配对）。
"""
function get_reversal_hamming(data1::HistData, data2::HistData)
    return get_overlap(data1, data2; apply_bias_correction=false)
end
