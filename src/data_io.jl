# ---------------------
# 从 qmeas.random 输出（一次 SettingRun）导入数据。
# ---------------------
# hist-v1（当前唯一格式）：计数存 per-setting 去重 key + 权重，
# 见 qmeas src/qmeas/random/io.py counts_to_hist。旧展开格式不再支持。
# - import_hist_group / import_hist_pair：hist 原生 HistGroup（全管线统一入口）。

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
