using Test
using RandomMeasAdd

# Robust Hamming 标定与核的单元测试：全手算小例子，不读数据文件。
# 约定：{1,2} 计数中 2 ↔ 物理 0（与 data_io 的 2.-bit 一致）。

@testset "calibration f0/R" begin
    # f̂0：列 [2,2,1,2] → 3/4。
    d = RandomMeasAdd.MeasurementData([2; 2; 1; 2;;])
    g = RandomMeasAdd.MeasurementGroup([d])
    @test RandomMeasAdd.get_f0(g, [1]) ≈ [0.75]

    # R̂ purity：列 [2,2,1,1]，c0=c1=2 → matches=2+2=4，4/(4*3)=1/3。
    d2 = RandomMeasAdd.MeasurementData([2; 2; 1; 1;;])
    g2 = RandomMeasAdd.MeasurementGroup([d2])
    @test RandomMeasAdd.get_R_purity(g2, [1]) ≈ [1 / 3]

    # R̂ mirror：行 [2,2],[2,1],[1,1] → 2/3。
    d3 = RandomMeasAdd.MeasurementData([2 2; 2 1; 1 1])
    g3 = RandomMeasAdd.MeasurementGroup([d3])
    @test RandomMeasAdd.get_R_mirror(g3, [1], [2]) ≈ [2 / 3]

    # R̂ cross：g1 列 c0=2，g2 列 c0=1，NM=3 → (2*1+1*2)/9=4/9。
    a = RandomMeasAdd.MeasurementData([2; 2; 1;;])
    b = RandomMeasAdd.MeasurementData([2; 1; 1;;])
    ga = RandomMeasAdd.MeasurementGroup([a])
    gb = RandomMeasAdd.MeasurementGroup([b])
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

@testset "robust kernels match standard kernels" begin
    # 标准核 (2/-1) 下，robust 路径须与上游路径一致。
    res = [2 1; 1 1; 2 2; 1 2; 2 1]
    d = RandomMeasAdd.MeasurementData(res)
    xs = [2.0, 2.0]
    ys = [-1.0, -1.0]
    @test RandomMeasAdd.get_overlap_robust(d, d, xs, ys; apply_bias_correction=true) ≈
        RandomMeasAdd.get_overlap(d, d; apply_bias_correction=true) atol = 1e-10
    @test RandomMeasAdd.get_overlap_robust(d, d, xs, ys; apply_bias_correction=false) ≈
        RandomMeasAdd.get_overlap(d, d; apply_bias_correction=false) atol = 1e-10
    @test RandomMeasAdd.get_reflect_robust(d, [1], [2], [2.0], [-1.0]) ≈
        RandomMeasAdd.get_reflect_hamming(d) atol = 1e-10

    # bias 修正恒等式：缩并修正 == 暴力 k≠k' 求和。
    NM = size(res, 1)
    K = zeros(NM, NM)
    for k in 1:NM, kp in 1:NM
        v = 1.0
        for j in 1:2
            v *= res[k, j] == res[kp, j] ? 2.0 : -1.0
        end
        K[k, kp] = v
    end
    brute = (sum(K) - sum(K[i, i] for i in 1:NM)) / (NM * (NM - 1))
    @test RandomMeasAdd.get_overlap_robust(d, d, xs, ys; apply_bias_correction=true) ≈ brute atol = 1e-10
end
