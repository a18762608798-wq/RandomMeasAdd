# ---------------------
# 从 qmeas.random 输出（一次 SettingRun）导入数据。
# ---------------------
# hist-v1（当前唯一格式）：计数存 per-setting 去重 key + 权重，
# 见 qmeas src/qmeas/random/io.py counts_to_hist。旧展开格式不再支持。
# - import_hist_group / import_hist_pair：hist 原生 HistGroup（hamming 路）。
# - import_random_group / import_random_pair：legacy 上游 MeasurementGroup
#   （shadow 路过渡 shim，内部 hist→expand；Phase 2 删除）。

"""npz 须为 hist-v1，否则报错（旧展开格式请重生成）。"""
function _check_hist_format(group_data, filepath)
    @assert haskey(group_data, "measurement_hist_keys") "$filepath 不是 hist-v1：缺少 measurement_hist_keys（旧展开格式已不再支持，请用新版 qmeas 重生成）。"
    if haskey(group_data, "format_version")
        @assert Int(group_data["format_version"][]) == 1 "$filepath 的 format_version 不是 1。"
    end
end

"""由 hist 三数组建某前缀（"" / "trivial_" / "trivial_plain_"）的 per-setting 数据。缺键返回 nothing。"""
function _parse_hist_prefix(group_data, prefix, N)
    kc = prefix * "hist_keys"
    cc = prefix * "hist_counts"
    oc = prefix * "hist_offsets"
    (haskey(group_data, kc) && haskey(group_data, cc) && haskey(group_data, oc)) || return nothing
    keys = Vector{UInt64}(vec(group_data[kc]))
    counts = Vector{Int64}(vec(group_data[cc]))
    offsets = Vector{Int64}(vec(group_data[oc]))
    NU = length(offsets) - 1
    @assert NU >= 1 && offsets[1] == 0 "hist_offsets 非法"
    slices = [(keys[offsets[s]+1:offsets[s+1]], counts[offsets[s]+1:offsets[s+1]]) for s in 1:NU]
    nms = [sum(c) for (_, c) in slices]
    @assert all(==(nms[1]), nms) "各 setting 权重和不一致"
    per = [HistData(k, c, nms[1], N) for (k, c) in slices]
    return per, nms[1], NU
end

"""
由 npz dict 建 HistGroup（三件套：main/trivial/plain，缺啥给 nothing）。
返回 (hgroup, htrivial, hplain, N, NU, NM)。
"""
function _read_hist_groups(group_data, permuted_order)
    settings = ComplexF64.(group_data["measurement_settings"])
    NU_file, N = size(settings, 1), size(settings, 2)
    if isnothing(permuted_order)
        permuted_order = collect(1:N)
    end
    @assert sort(permuted_order) == collect(1:N) "permuted_order 须为 1:N 的置换"
    main, NM, NU = _parse_hist_prefix(group_data, "measurement_", N)
    @assert !isnothing(main) && NU == NU_file "主 hist 缺失或设置数与 settings 不一致"
    t = _parse_hist_prefix(group_data, "trivial_", N)
    p = _parse_hist_prefix(group_data, "trivial_plain_", N)
    if !isnothing(t)
        @assert t[3] == NU && t[2] == NM "trivial hist 的设置数/shots 与主数据不一致"
    end
    if !isnothing(p)
        @assert p[3] == NU && p[2] == NM "plain hist 的设置数/shots 与主数据不一致"
    end
    # 列置换：keys 比特重排 + settings 切片。
    permute_one(per, set) = begin
        ms = map(per) do d
            HistData(permute_hist_keys(d.keys, permuted_order, N), d.counts, d.nm, N)
        end
        HistGroup(ms, set, N, NU, NM)
    end
    hgroup = permute_one(main, settings[:, permuted_order, :, :])
    htrivial = isnothing(t) ? nothing :
        permute_one(t[1], ComplexF64.(group_data["trivial_measurement_settings"])[:, permuted_order, :, :])
    hplain = if isnothing(p)
        nothing
    else
        # plain 复用 trivial 的 settings（角度绑定相同，见 qmeas io）。
        pm, _, _ = p
        ms = map(pm) do d
            HistData(permute_hist_keys(d.keys, permuted_order, N), d.counts, d.nm, N)
        end
        HistGroup(ms, ComplexF64.(group_data["trivial_measurement_settings"])[:, permuted_order, :, :], N, NU, NM)
    end
    return hgroup, htrivial, hplain, N, NU, NM
end

"""
import_hist_group(filepath; permuted_order) -> (hgroup, htrivial, hplain)

hist 原生导入（hamming 路）：主/trivial/plain 三个 HistGroup，
缺键给 nothing。permuted_order 缺省链式（恒等置换）。
"""
function import_hist_group(filepath::String; permuted_order=nothing)
    group_data = npzread(filepath)
    _check_hist_format(group_data, filepath)
    hgroup, htrivial, hplain, _, _, _ = _read_hist_groups(group_data, permuted_order)
    return hgroup, htrivial, hplain
end

"""
import_hist_pair(filepath1, filepath2; permuted_order)
    -> (h1, h2, ht1, ht2, hp1)

hist 原生配对导入（hamming Z_T 路）：两份数据同 frame，
N/NU/NM 须一致（逐行配对前提）。缺键给 nothing。
"""
function import_hist_pair(filepath1::String, filepath2::String; permuted_order=nothing)
    h1, ht1, hp1 = import_hist_group(filepath1; permuted_order)
    h2, ht2, _ = import_hist_group(filepath2; permuted_order)
    @assert h1.N == h2.N "两份实验的比特数必须一致。"
    @assert h1.NU == h2.NU "两份实验的设置数必须一致且逐行配对。"
    @assert h1.NM == h2.NM "两份实验的每设置 shots 数必须一致。"
    return h1, h2, ht1, ht2, hp1
end

"""
从 HistGroup 展开成旧 0/1 矩阵（shadow 过渡 shim 用）。
"""
function _expand_for_legacy(hgroup::HistGroup)
    return 2 .- expand_hist(hgroup)
end

"""
从 qmeas.random 生成的 npz 文件中导入一次 SettingRun 的完整数据，
构造一个 MeasurementGroup。

输入须为 hist-v1（内部展开成旧矩阵再走原逻辑；shadow 路过渡用）。
summary json 只给人看，导入不需要它。

参数
- filepath::String：qmeas.random 生成的单个 .npz 文件路径。
- permuted_order：全部被测 site 的置换向量（Julia 从 1 开始编号），
  长度等于被测 qubit 总数。缺省为 `nothing`，表示链式（恒等置换），
  即按 npz 列顺序直接导入。

返回
- permuted_group::MeasurementGroup
- permuted_sites：重排后的 site Index 对象。
- G：校准向量（is_mitigation=false 时为 nothing）。
- is_return_trivial=true 时追加第 4 个返回 trivial_group
  （is_mitigation=false 时为 nothing，供 hamming 缓解标定用）。
- is_return_trivial_plain=true 时（须同时开 is_return_trivial）追加第 5 个
  返回 trivial_plain_group：配对实验一的不带 u_T 的 plain 标定 group，
  供 hamming Z_T 缓解用；缺键时报错（需新版 qmeas 重新生成 pair 数据）。
  非 pair 文件该键恒缺席，此时只能为 nothing（断言保证）。
"""
function import_random_group(filepath::String; permuted_order=nothing, is_mitigation=false, is_return_trivial=false, is_return_trivial_plain=false)
    @assert !is_return_trivial_plain || is_return_trivial "is_return_trivial_plain 须配合 is_return_trivial 使用。"
    @assert !is_return_trivial_plain || is_mitigation "is_return_trivial_plain 须配合 is_mitigation=true 使用（需读 trivial 数据）。"
    group_data = npzread(filepath)
    _check_hist_format(group_data, filepath)
    hgroup, htrivial, hplain, N, _, _ = _read_hist_groups(group_data, permuted_order)
    # python 从 0 开始 -> julia 从 1 开始的 qubit 编号
    qubits_jl = vec(Int64.(group_data["meas_indices"])) .+ 1
    # 全系统 site 数取最大 qubit 编号，在函数内直接构造 sites
    sites = siteinds("Qubit", maximum(qubits_jl))
    # hist 列顺序即 npz 列顺序（key bit j = 第 j 列），_read_hist_groups 内已按
    # permuted_order 重排；此处只选 site 子集，不再索引矩阵。
    order = isnothing(permuted_order) ? collect(1:N) : permuted_order
    site_indices = sites[qubits_jl]
    permuted_indices = site_indices[order]

    meas_res = _expand_for_legacy(hgroup)
    settings = ComplexF64.(group_data["measurement_settings"])
    permuted_group = MeasurementGroup(
        meas_res,
        settings[:, order, :, :],
        permuted_indices,
    )

    # trivial group
    G = nothing
    permuted_trivial_group = nothing
    if is_mitigation
        @assert !isnothing(htrivial) "npz 缺少 trivial hist：生成该文件时未测 trivial 校准数据，is_mitigation 须为 false"
        trivial_meas_res = _expand_for_legacy(htrivial)
        trivial_settings = ComplexF64.(group_data["trivial_measurement_settings"])
        permuted_trivial_group = MeasurementGroup(
            trivial_meas_res,
            trivial_settings[:, order, :, :],
            permuted_indices,
        )
        # get group
        ψ0 = MPS(permuted_indices, "1")
        G = modified_get_calibration_vector(ψ0, permuted_trivial_group)
    end
    if is_return_trivial
        if is_return_trivial_plain
            @assert !isnothing(hplain) "npz 缺少 trivial_plain hist：该 pair 文件系旧版 qmeas 生成（trivial_1 含 u_T，不可用于 hamming Z_T 标定），请用新版 qmeas 重新生成。"
            plain_meas_res = _expand_for_legacy(hplain)
            permuted_plain_group = MeasurementGroup(
                plain_meas_res,
                trivial_settings[:, order, :, :],
                permuted_indices,
            )
            return permuted_group, permuted_indices, G, permuted_trivial_group, permuted_plain_group
        end
        return permuted_group, permuted_indices, G, permuted_trivial_group
    end
    return permuted_group, permuted_indices, G
end

"""
从一对共轭配对实验（Z_T 用，见 qmeas ConjugatePair）的两个 npz 文件
导入两个 MeasurementGroup，共用同一重排 frame。

参数
- filepath1::String：实验一的 .npz 文件路径。
- filepath2::String：实验二的 .npz 文件路径。
- permuted_order：全部被测 site 的置换向量（Julia 从 1 开始编号），
  两份数据必须共用同一个。缺省为 `nothing`，表示链式（恒等置换）。

返回
- group1::MeasurementGroup：实验一的重排 group。
- group2::MeasurementGroup：实验二的重排 group，设置与 group1 逐行配对。
- permuted_sites：重排后的 site Index 对象（两份共用）。
- G1：实验一的校准向量（is_mitigation=false 时为 nothing）。
- G2：实验二的校准向量（is_mitigation=false 时为 nothing）。
- is_return_trivial=true 时追加第 6、7 个返回 (trivial_group1, trivial_group2)
  （is_mitigation=false 时为 (nothing, nothing)，供 hamming 缓解标定用）。
- is_return_trivial_plain=true 时（须同时开前两个开关）追加第 8 个返回
  plain_group1：实验一的不带 u_T 的 plain 标定 group，供 hamming Z_T
  缓解用（此时 trivial1 照传，供 cross_xy 的只有 plain1 + trivial2）。
  实验二的 plain 就是它自己的 trivial2，不另传。

说明
两份数据的被测 qubit 布局（meas_indices）必须一致，
且设置数 NU 必须相等，否则逐行互关联无意义。
"""
function import_random_pair(filepath1::String, filepath2::String; permuted_order=nothing, is_mitigation=false, is_return_trivial=false, is_return_trivial_plain=false)
    if is_return_trivial_plain
        group1, permuted_indices1, G1, trivial_group1, plain_group1 = import_random_group(
            filepath1; permuted_order, is_mitigation,
            is_return_trivial=true, is_return_trivial_plain=true,
        )
    else
        group1, permuted_indices1, G1, trivial_group1 = import_random_group(
            filepath1; permuted_order, is_mitigation, is_return_trivial=true,
        )
        plain_group1 = nothing
    end
    group2, permuted_indices2, G2, trivial_group2 = import_random_group(
        filepath2; permuted_order, is_mitigation, is_return_trivial=true
    )
    # 同一 frame、逐行配对：布局与设置数都须一致
    # （Index 的 == 比内部 id，两次导入各自构造 sites，故比 tags 与维数）
    @assert tags.(permuted_indices1) == tags.(permuted_indices2) "两份实验的 site 布局不一致：检查 meas_indices 与 permuted_order 是否相同。"
    @assert dim.(permuted_indices1) == dim.(permuted_indices2) "两份实验的 site 布局不一致：检查 meas_indices 与 permuted_order 是否相同。"
    @assert group1.N == group2.N "两份实验的比特数必须一致。"
    @assert group1.NU == group2.NU "两份实验的设置数必须一致且逐行配对。"
    @assert group1.NM == group2.NM "两份实验的每设置 shots 数必须一致。"
    if is_return_trivial
        if is_return_trivial_plain
            return group1, group2, permuted_indices1, G1, G2, trivial_group1, trivial_group2, plain_group1
        end
        return group1, group2, permuted_indices1, G1, G2, trivial_group1, trivial_group2
    end
    return group1, group2, permuted_indices1, G1, G2
end


