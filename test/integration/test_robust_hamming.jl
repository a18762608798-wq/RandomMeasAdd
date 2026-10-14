using Test
using RandomMeasAdd
using Random

# Robust Hamming 集成测试：aer setting0 数据（243 settings，无噪声仿真）。
# 1. 无噪声极限：缓解版 ≈ 标准版（trivial 自带统计涨落，留 5% 相对容差）。
# 2. 合成读出噪声：主数据注入非对称翻转后，标准估计偏掉、缓解版恢复。

const _DATADIR = joinpath(@__DIR__, "..", "data")
const _PERM = [1, 2, 3, 4]
const _IND0 = joinpath(
    _DATADIR, "aer_independent_pidx_1",
    "aer_independent_pidx_1_setting0_settings243_shots1024.npz",
)
const _SH0 = joinpath(
    _DATADIR, "aer_shared_pidx_1",
    "aer_shared_pidx_1_setting0_settings243_shots1024.npz",
)
const _PAIR1 = joinpath(
    _DATADIR, "aer_pair_mit_pidx_1",
    "aer_pair_mit_pidx_1_exp1_setting0_settings243_shots1024.npz",
)
const _PAIR2 = joinpath(
    _DATADIR, "aer_pair_mit_pidx_1",
    "aer_pair_mit_pidx_1_exp2_setting0_settings243_shots1024.npz",
)

"""内存中注入非对称读出噪声：物理 0（值 2）以 p01 翻成 1，物理 1 翻成 0 以 p10。"""
function _inject_readout_noise(group; p01=0.08, p10=0.02, seed=1234)
    ng = deepcopy(group)
    rng = MersenneTwister(seed)
    for data in ng.measurements
        res = data.measurement_results
        for idx in eachindex(res)
            v = res[idx]
            r = rand(rng)
            if v == 2 && r < p01
                res[idx] = 1
            elseif v == 1 && r < p10
                res[idx] = 2
            end
        end
    end
    return ng
end

@testset "noiseless limit" begin
    group, _, _, trivial = import_random_group(
        _IND0; permuted_order=_PERM, is_mitigation=true, is_return_trivial=true,
    )
    # 核级：243-setting 标定统计下，核权重落在标准核附近。
    xs, ys = RandomMeasAdd.get_purity_xy(trivial, 1:4)
    @test all(abs.(xs .- 2.0) .< 0.25)
    @test all(abs.(ys .+ 1.0) .< 0.25)
    # 估计子级：标定涨落在 4 比特乘积结构中放大，留 15% 相对容差。
    p_std = get_purity_hamming(group; is_compute_sem=false)
    p_rob = get_purity_hamming(group; trivial, is_mitigation=true)
    @test isapprox(p_rob, p_std; rtol=0.15)

    # reflect / z_r 要求镜像对共享转角，用 shared 数据。
    sgroup, _, _, strivial = import_random_group(
        _SH0; permuted_order=_PERM, is_mitigation=true, is_return_trivial=true,
    )
    r_std = get_reflect_hamming(sgroup; is_compute_sem=false)
    r_rob = get_reflect_hamming(sgroup; trivial=strivial, is_mitigation=true)
    @test isapprox(r_rob, r_std; rtol=0.15)

    zr_std = get_z_r_hamming(sgroup)
    zr_rob = get_z_r_hamming(sgroup; trivial=strivial, is_mitigation=true)
    @test isapprox(zr_rob, zr_std; rtol=0.15)

    g1, g2, _, _, _, t1, t2 = import_random_pair(
        _PAIR1, _PAIR2; permuted_order=_PERM,
        is_mitigation=true, is_return_trivial=true,
    )
    zt_std = get_reversal_hamming(g1, g2)
    zt_rob = get_reversal_hamming(g1, g2; trivial1=t1, trivial2=t2, is_mitigation=true)
    @test isapprox(zt_rob, zt_std; rtol=0.15)

    z_std = get_z_t_hamming(g1, g2)
    z_rob = get_z_t_hamming(g1, g2; trivial1=t1, trivial2=t2, is_mitigation=true)
    @test isapprox(z_rob, z_std; rtol=0.15)
end

@testset "synthetic readout noise recovery" begin
    # 方法论要点：噪声同时进主数据和 trivial 数据（同一台设备、同 Λ），
    # 用不同随机流模拟两次独立运行。
    group, _, _, trivial = import_random_group(
        _IND0; permuted_order=_PERM, is_mitigation=true, is_return_trivial=true,
    )
    noisy = _inject_readout_noise(group; seed=1234)
    noisy_trivial = _inject_readout_noise(trivial; seed=4321)
    # 注入的是 p01=0.08/p10=0.02 非对称翻转，对应通道 δ=p10-p01=-0.06；
    # 标定应能把它认出来（f̂0 由 setting 级 Haar 涨落主导，σ(δ̂)~0.03，
    # 故只断言 4 比特均值，容差 ~2σ）。
    dhat = 2 * RandomMeasAdd.get_f0(noisy_trivial, 1:4) .- 1
    @test abs(sum(dhat) / length(dhat) + 0.06) < 0.035
    clean = get_purity_hamming(group; is_compute_sem=false)
    biased = get_purity_hamming(noisy; is_compute_sem=false)
    fixed = get_purity_hamming(noisy; trivial=noisy_trivial, is_mitigation=true)
    @test abs(biased - clean) > 1e-3 # 噪声确实打偏了标准估计
    @test abs(fixed - clean) < abs(biased - clean) # 缓解把偏差收回来

    sgroup, _, _, strivial = import_random_group(
        _SH0; permuted_order=_PERM, is_mitigation=true, is_return_trivial=true,
    )
    snoisy = _inject_readout_noise(sgroup; seed=1234)
    snoisy_trivial = _inject_readout_noise(strivial; seed=4321)
    r_clean = get_reflect_hamming(sgroup; is_compute_sem=false)
    r_biased = get_reflect_hamming(snoisy; is_compute_sem=false)
    r_fixed = get_reflect_hamming(snoisy; trivial=snoisy_trivial, is_mitigation=true)
    @test abs(r_biased - r_clean) > 1e-3
    @test abs(r_fixed - r_clean) < abs(r_biased - r_clean)

    g1, g2, _, _, _, t1, t2 = import_random_pair(
        _PAIR1, _PAIR2; permuted_order=_PERM,
        is_mitigation=true, is_return_trivial=true,
    )
    n1 = _inject_readout_noise(g1; seed=1234)
    n2 = _inject_readout_noise(g2; seed=4321)
    nt1 = _inject_readout_noise(t1; seed=111)
    nt2 = _inject_readout_noise(t2; seed=222)
    z_clean = get_reversal_hamming(g1, g2)
    z_biased = get_reversal_hamming(n1, n2)
    z_fixed = get_reversal_hamming(n1, n2; trivial1=nt1, trivial2=nt2, is_mitigation=true)
    @test abs(z_biased - z_clean) > 1e-3
    @test abs(z_fixed - z_clean) < abs(z_biased - z_clean)
end
