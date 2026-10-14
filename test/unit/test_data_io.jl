using Test
using RandomMeasAdd
using NPZ

# data_io 单元测试：NM 失配断言（临时 npz，不依赖物理有效性，只切真实数组）。

const _DATADIR = joinpath(@__DIR__, "..", "data")
const _SMALL_NPZ = joinpath(
    _DATADIR, "aer_independent_pidx_1",
    "aer_independent_pidx_1_setting0_settings243_shots1024.npz",
)

"""由真实文件切出 (NU=5, NM=nm) 的最小 npz 写进临时目录。"""
function _write_cut_npz(dir, name, nm)
    raw = npzread(_SMALL_NPZ)
    cut = Dict(
        "measurement_results" => raw["measurement_results"][1:5, 1:nm, :],
        "measurement_settings" => raw["measurement_settings"][1:5, :, :, :],
        "meas_indices" => raw["meas_indices"],
    )
    path = joinpath(dir, name)
    npzwrite(path, cut)
    return path
end

@testset "import_random_pair NM mismatch" begin
    mktempdir() do dir
        f1 = _write_cut_npz(dir, "a.npz", 100)
        f2 = _write_cut_npz(dir, "b.npz", 200)
        @test_throws AssertionError import_random_pair(
            f1, f2; permuted_order=[1, 2, 3, 4]
        )
        # 同 NM 通过。
        g1, g2, _, _, _ = import_random_pair(
            f1, f1; permuted_order=[1, 2, 3, 4]
        )
        @test g1.NM == g2.NM == 100
    end
end

@testset "import_random_group trivial return" begin
    # 默认三元组形状不变。
    out3 = import_random_group(_SMALL_NPZ; permuted_order=[1, 2, 3, 4])
    @test length(out3) == 3
    # 开关打开给四元组；该文件生成时开 mitigation，trivial 应存在且对齐。
    group, _, _, trivial = import_random_group(
        _SMALL_NPZ; permuted_order=[1, 2, 3, 4],
        is_mitigation=true, is_return_trivial=true,
    )
    @test !isnothing(trivial)
    @test (trivial.N, trivial.NU, trivial.NM) == (group.N, group.NU, group.NM)
    # 无 mitigation 时第四元为 nothing。
    _, _, _, trivial_off = import_random_group(
        _SMALL_NPZ; permuted_order=[1, 2, 3, 4], is_return_trivial=true,
    )
    @test isnothing(trivial_off)
end
