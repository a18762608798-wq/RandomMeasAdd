using Test
using RandomMeasAdd
using NPZ

# data_io 单元测试（hist-v1）：NM 失配断言 + trivial 返回形状。
# 手造最小 hist npz，不依赖真实数据文件。

const _PERM4 = [1, 2, 3, 4]

"""写最小 hist-v1 npz：单 setting，keys/counts 给定，settings 为零数组占位。"""
function _write_hist_npz(dir, name, keys_list, counts_list, nm; n=4)
    nu = length(keys_list)
    offs = [0; cumsum(length.(keys_list))]
    npzwrite(joinpath(dir, name), Dict(
        "format_version" => Int64(1),
        "measurement_hist_keys" => UInt64[k for ks in keys_list for k in ks],
        "measurement_hist_counts" => Int64[c for cs in counts_list for c in cs],
        "measurement_hist_offsets" => Int64.(offs),
        "measurement_settings" => zeros(ComplexF64, nu, n, 2, 2),
        "meas_indices" => Int64.(collect(0:n-1)),
        "num_shots" => Int64(nm),
    ))
    return joinpath(dir, name)
end

@testset "import_hist_pair NM mismatch" begin
    mktempdir() do dir
        f1 = _write_hist_npz(dir, "a.npz", [[0x0, 0x3]], [[60, 40]], 100)
        f2 = _write_hist_npz(dir, "b.npz", [[0x0, 0x3]], [[120, 80]], 200)
        @test_throws AssertionError import_hist_pair(f1, f2; permuted_order=_PERM4)
        # 同 NM 通过。
        h1, h2, _, _, _ = import_hist_pair(f1, f1; permuted_order=_PERM4)
        @test h1.NM == h2.NM == 100
        @test h1.N == h2.N == 4
    end
end

@testset "import_hist_group trivial return" begin
    mktempdir() do dir
        # 无 trivial 键：后两元为 nothing。
        f = _write_hist_npz(dir, "s.npz", [[0x0, 0x1]], [[30, 70]], 100)
        h, t, p = import_hist_group(f; permuted_order=_PERM4)
        @test (h.N, h.NU, h.NM) == (4, 1, 100)
        @test isnothing(t) && isnothing(p)
    end
end

@testset "import_hist_group shape" begin
    # hist 原生导入：group 对齐，trivial 缺席时为 nothing。
    mktempdir() do dir
        f = _write_hist_npz(dir, "s.npz", [[0x5, 0xa]], [[50, 50]], 100)
        h, t, p = import_hist_group(f; permuted_order=_PERM4)
        @test (h.N, h.NU, h.NM) == (4, 1, 100)
        @test isnothing(t) && isnothing(p)
    end
end
