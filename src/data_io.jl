# ---------------------
# 从 qmeas.random 输出（一次 SettingRun）导入 MeasurementGroup，并按 permuted_order 重排。
# ---------------------

"""
从 qmeas.random 生成的 npz 文件中导入一次 SettingRun 的完整数据，
构造一个 MeasurementGroup。

npz 文件存了一次 SettingRun 下所有被测 qubit 的数据
（`measurement_results` 形状为 (num_settings, num_shots, n_meas)），
另有 `meas_indices`（展平、python 从 0 开始编号）记录被测 qubit。
本函数取全部列，qubit 编号转为 Julia 惯例（+1），
再用 `permuted_order` 对结果、设置和 site 统一重排。
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
    # python 从 0 开始 -> julia 从 1 开始的 qubit 编号
    qubits_jl = vec(Int64.(group_data["meas_indices"])) .+ 1
    # 全系统 site 数取最大 qubit 编号，在函数内直接构造 sites
    sites = siteinds("Qubit", maximum(qubits_jl))
    site_indices = sites[qubits_jl]

    meas_res = 2 .- Int64.(group_data["measurement_results"])
    settings = ComplexF64.(group_data["measurement_settings"])

    # 缺省链式：保持 npz 列顺序
    if isnothing(permuted_order)
        permuted_order = collect(1:length(site_indices))
    end
    permuted_indices = site_indices[permuted_order]
    permuted_group = MeasurementGroup(
        meas_res[:, :, permuted_order],
        settings[:, permuted_order, :, :],
        permuted_indices,
    )

    # trivial group
    G = nothing
    permuted_trivial_group = nothing
    if is_mitigation
        @assert haskey(group_data, "trivial_measurement_results") "npz 缺少 trivial_measurement_results：生成该文件时未测 trivial 校准数据，is_mitigation 须为 false"
        @assert haskey(group_data, "trivial_measurement_settings") "npz 缺少 trivial_measurement_settings：生成该文件时未测 trivial 校准数据，is_mitigation 须为 false"
        # get group
        trivial_meas_res = 2 .- Int64.(group_data["trivial_measurement_results"])
        trivial_settings = ComplexF64.(group_data["trivial_measurement_settings"])
        permuted_trivial_group = MeasurementGroup(
            trivial_meas_res[:, :, permuted_order],
            trivial_settings[:, permuted_order, :, :],
            permuted_indices,
        )
        # get G
        ψ0 = MPS(permuted_indices, "1")
        G = modified_get_calibration_vector(ψ0, permuted_trivial_group)
    end
    if is_return_trivial
        if is_return_trivial_plain
            @assert haskey(group_data, "trivial_plain_measurement_results") "npz 缺少 trivial_plain_measurement_results：该 pair 文件系旧版 qmeas 生成（trivial_1 含 u_T，不可用于 hamming Z_T 标定），请用新版 qmeas 重新生成。"
            plain_meas_res = 2 .- Int64.(group_data["trivial_plain_measurement_results"])
            permuted_plain_group = MeasurementGroup(
                plain_meas_res[:, :, permuted_order],
                trivial_settings[:, permuted_order, :, :],
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



