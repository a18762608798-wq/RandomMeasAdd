using Test
using RandomMeasAdd

# Robust Hamming 标定与核的单元测试：全手算小例子，不读数据文件。
# 约定：hist key 的 bit 0 ↔ 物理 0（与旧矩阵 0/1 一致）。

"""由 0/1 行向量建单比特 HistData（keys 去重升序）。"""
function _h1(rows)
    uk = sort(unique(rows))
    return RandomMeasAdd.HistData(
        UInt64.(uk), Int64[sum(==(k), rows) for k in uk], length(rows), 1,
    )
end

"""由 0/1 矩阵（每行一个 outcome）建 HistData。"""
function _hm(mat)
    n = size(mat, 2)
    keys = UInt64[sum(row[j] << (j - 1) for j in 1:n) for row in eachrow(mat)]
    uk = sort(unique(keys))
    return RandomMeasAdd.HistData(
        uk, Int64[sum(==(k), keys) for k in uk], length(keys), n,
    )
end

"""HistGroup 壳（settings 占位，标定/核不用它）。"""
_hg(d::RandomMeasAdd.HistData) = RandomMeasAdd.HistGroup(
    [d], zeros(ComplexF64, 1, d.n, 2, 2), d.n, 1, d.nm,
)

@testset "calibration f0/R" begin
    # f̂0：列值 [0,0,0,1] → 3/4。
    g = _hg(_h1([0, 0, 0, 1]))
    @test RandomMeasAdd.get_f0(g, [1]) ≈ [0.75]

    # R̂ purity：列值 [0,0,1,1]，c0=c1=2 → 4/(4*3)=1/3。
    g2 = _hg(_h1([0, 0, 1, 1]))
    @test RandomMeasAdd.get_R_purity(g2, [1]) ≈ [1 / 3]

    # R̂ mirror：行 [0,0],[0,1],[1,1] → 2/3。
    g3 = _hg(_hm([0 0; 0 1; 1 1]))
    @test RandomMeasAdd.get_R_mirror(g3, [1], [2]) ≈ [2 / 3]

    # R̂ cross：g1 列 c0=2，g2 列 c0=1，NM=3 → (2*1+1*2)/9=4/9。
    ga = _hg(_h1([0, 0, 1]))
    gb = _hg(_h1([0, 1, 1]))
    @test RandomMeasAdd.get_R_cross(ga, gb, [1]) ≈ [4 / 9]
end

@testset "calibration get_xy" begin
    # 无噪声极限退化为标准核。
    x, y = RandomMeasAdd.get_xy(2 / 3, 0.0, 0.0)
    @test x ≈ 2.0
    @test y ≈ -1.0
    # 退化数据（全零 trivial：R̂=1，δ=1）分母为零，应报错而非吐 Inf。
    @test_throws AssertionError RandomMeasAdd.get_xy(1.0, 1.0, 1.0)
end

@testset "weighted kernels match brute force" begin
    # 暴力 ground truth： outcomes 逐个比，不经直方图也不经张量。
    res = [0 1; 1 1; 0 0; 1 0; 0 1] # (NM=5, N=2)，0/1 取值
    outs = [Vector(res[k, :]) for k in 1:size(res, 1)]
    NM = length(outs)
    std_kernel(a, b) = prod(x == y ? 2.0 : -1.0 for (x, y) in zip(a, b))
    K = [std_kernel(a, b) for a in outs, b in outs]
    brute = (sum(K) - sum(K[i, i] for i in 1:NM)) / (NM * (NM - 1))

    d = _hm(res)
    xs = [2.0, 2.0]
    ys = [-1.0, -1.0]
    @test RandomMeasAdd.get_overlap_robust(d, d, xs, ys; apply_bias_correction=true) ≈ brute atol = 1e-12
    @test RandomMeasAdd.get_overlap(d, d; apply_bias_correction=true) ≈ brute atol = 1e-12
    # 不加 bias：与 K 全平均一致。
    @test RandomMeasAdd.get_overlap(d, d; apply_bias_correction=false) ≈ sum(K) / NM^2 atol = 1e-12

    # reflect：暴力按相邻对算。
    function brute_reflect(outs)
        s = 0.0
        for o in outs
            D = (o[1] != o[2])
            s += 2.0^1 * (-2.0)^(-D)
        end
        return s / length(outs)
    end
    @test RandomMeasAdd.get_reflect_robust(d, [1], [2], [2.0], [-1.0]) ≈ brute_reflect(outs) atol = 1e-12
    @test RandomMeasAdd.get_reflect_hamming(d) ≈ brute_reflect(outs) atol = 1e-12
end
