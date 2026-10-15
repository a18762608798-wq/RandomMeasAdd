# ---------------------
# 自有直方图数据模型（hist-v1）：不依赖上游 RandomMeas 类型。
# ---------------------
# 约定（与 qmeas hist-v1 一致）：
# - key 为 UInt64 整数，bit j（0-based，即 (key >> j) & 1）= 第 j 列的值，
#   列按 meas_indices 展平顺序；0 ↔ 物理 0（与旧矩阵 0/1 一致，无需 2.- 变换）。
# - 每个 setting 存去重 key（升序）+ 权重 counts，权重和 == shots。
# - settings 矩阵沿用旧含义 (NU, N, 2, 2)，Phase 2 直接张量用。

"""单个 setting 的直方图数据。"""
struct HistData
    keys::Vector{UInt64}
    counts::Vector{Int64}
    nm::Int
    n::Int
    function HistData(keys, counts, nm, n)
        @assert length(keys) == length(counts) "keys/counts 长度不一致"
        @assert issorted(keys) "keys 须升序"
        @assert sum(counts) == nm "权重和 $(sum(counts)) != shots $nm"
        @assert all(>=(0), counts) "counts 须非负"
        @assert n <= 64 "n=$n 超出 UInt64 key"
        new(Vector{UInt64}(keys), Vector{Int64}(counts), nm, n)
    end
end

"""一次 SettingRun 的直方图 group（对标上游 MeasurementGroup）。"""
struct HistGroup
    measurements::Vector{HistData}
    settings::Array{ComplexF64,4}
    N::Int
    NU::Int
    NM::Int
    function HistGroup(measurements, settings, N, NU, NM)
        @assert length(measurements) == NU "setting 数不一致"
        @assert size(settings) == (NU, N, 2, 2) "settings 形状须为 (NU, N, 2, 2)"
        for d in measurements
            @assert d.n == N && d.nm == NM "HistData 的 n/nm 与 group 不一致"
        end
        new(measurements, settings, N, NU, NM)
    end
end

"""key 的第 j 列（1-based）→ 0/1，0 ↔ 物理 0。"""
hist_bit(key::UInt64, j::Int) = Int((key >> (j - 1)) & 0x1)

"""两 key 的 Hamming 距离（key 高位恒零，只数 N 位内差异）。"""
hist_dist(a::UInt64, b::UInt64) = count_ones(a ⊻ b)

"""
reduce_to_subsystem(h, cols) -> HistGroup

HistGroup 版子系统约化（显式扩展上游同名函数，多派发）：
只保留 cols 列，key 按 cols 顺序重排重编，settings 同步切片。
"""
function RandomMeas.reduce_to_subsystem(h::HistGroup, cols::AbstractVector{Int})
    @assert all(1 .<= cols .<= h.N) "列越界"
    @assert length(unique(cols)) == length(cols) "列不能重复"
    ms = map(h.measurements) do d
        agg = Dict{UInt64,Int64}()
        for (k, c) in zip(d.keys, d.counts)
            nk = UInt64(0)
            for (p, c0) in enumerate(cols)
                nk |= UInt64(hist_bit(k, c0)) << (p - 1)
            end
            agg[nk] = get(agg, nk, 0) + c
        end
        ks = sort(collect(keys(agg)))
        HistData(ks, [agg[k] for k in ks], d.nm, length(cols))
    end
    return HistGroup(ms, h.settings[:, cols, :, :], length(cols), h.NU, h.NM)
end

"""
expand_hist(h) -> Array{Int,3}

HistGroup 展开成旧格式 (NU, NM, N) 0/1 矩阵（shadow 过渡 shim 用，
行顺序按 key 升序；估计子对 shot 置换对称，顺序无关）。
"""
function expand_hist(h::HistGroup)
    mat = zeros(Int, h.NU, h.NM, h.N)
    for (s, d) in enumerate(h.measurements)
        row = 1
        for (k, c) in zip(d.keys, d.counts)
            for _ in 1:c
                for j in 1:h.N
                    mat[s, row, j] = hist_bit(k, j)
                end
                row += 1
            end
        end
        @assert row == h.NM + 1 "展开行数不一致"
    end
    return mat
end

"""
permute_hist_keys(keys, perm, N) -> Vector{UInt64}

列置换：新第 j 列 ← 旧第 perm[j] 列（1-based，Julia 侧编号）。
"""
function permute_hist_keys(keys::AbstractVector{UInt64}, perm::AbstractVector{Int}, N::Int)
    @assert sort(perm) == collect(1:N) "perm 须为 1:N 的置换"
    return map(keys) do k
        nk = UInt64(0)
        for j in 1:N
            nk |= UInt64(hist_bit(k, perm[j])) << (j - 1)
        end
        nk
    end
end
