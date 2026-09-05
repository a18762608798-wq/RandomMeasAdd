# -------------
# reversal（时间反演交叠 Z_T 分子）的 hamming 估计
# -------------
# 物理背景（Elben 等，Sci. Adv. 6, eaaz3666 (2020)，Eq. 5–6）：
# Z_T = Tr(ρ_I u_T ρ_I^{T_1} u_T†)，u_T = ∏_{i∈I_1} σ^y_i，
# T_1 为对 I_1 区的部分转置。Z_T 用两份实验的互关联估计：
# Z_T = 2^{2n} Σ_{s,s'} (-2)^{-D[s,s']} P^{(1)}(s) P^{(2)}(s')，
# 其中实验一幺正 U^{(1)} = U_{I_1} u_T ⊗ U_{I_2}（u_T 先旋转态，不能吸收进采样），
# 实验二幺正 U^{(2)} = U_{I_1}^* ⊗ U_{I_2}（I_1 区取 Haar 部分的复共轭，I_2 区相同）。
# 两份数据的设置必须逐行配对（qmeas ConjugatePair 保证）。
# 按本项目惯例：重排 frame 下奇位 = I_1 区，偶位 = I_2 区。
"""
get_reversal_hamming(data1::MeasurementData, data2::MeasurementData)

对单个配对设置，用 Hamming 距离法算两份实验的互关联（Z_T 分子）。

参数
- data1::MeasurementData：实验一（U_{I_1} u_T ⊗ U_{I_2}）的单个设置数据。
- data2::MeasurementData：实验二（U_{I_1}^* ⊗ U_{I_2}）的同行设置数据。

返回
- zt_est::Float64：该设置下的 Z_T 估计值。

说明
直接调 get_overlap(data1, data2)（不加 bias 修正）：
两份实验的 shots 相互独立，没有对角自配对，无需修正。
"""
function get_reversal_hamming(
    data1::MeasurementData,
    data2::MeasurementData,
)
    zt_est = get_overlap(data1, data2; apply_bias_correction=false)

    return zt_est
end

"""
get_reversal_hamming(group1, group2; is_compute_sem, is_show_progress)

用 Hamming 距离法从两份重排好的测量 group 算 Z_T 分子，
对各配对设置求平均。

参数
- group1：实验一的重排 group。
- group2：实验二的重排 group，设置须与 group1 逐行配对。

关键词参数
- is_compute_sem::Bool：为 true 时计算各配对设置间的均值标准误差（SEM）。
- is_show_progress::Bool：为 true 时显示进度。

返回
- is_compute_sem == false：返回 zt_est::Float64。
- is_compute_sem == true：返回 (zt_est::Float64, sem::Float64)。

说明
本函数对每个配对设置调用
get_reversal_hamming(::MeasurementData, ::MeasurementData)，
再对各设置求平均。线性平均不需要 jackknife。
"""
function get_reversal_hamming(
    group1,
    group2;
    is_compute_sem=false,
    is_show_progress=false,
)
    # 两份数据的几何与设置数必须一致，否则互关联无意义
    @assert group1.N == group2.N "两份实验的比特数必须一致。"
    @assert group1.NU == group2.NU "两份实验的设置数必须一致且逐行配对。"
    u_num = group1.NU
    datas1 = group1.measurements
    datas2 = group2.measurements
    zt_ests = Vector{Float64}(undef, u_num)

    @showprogress desc="hamming_est..." enabled=is_show_progress @threads for u_idx = 1:u_num
        zt_ests[u_idx] = get_reversal_hamming(datas1[u_idx], datas2[u_idx])
    end

    zt_est = mean(zt_ests)

    if is_compute_sem
        sem = std(zt_ests) / sqrt(u_num)
        return zt_est, sem
    else
        return zt_est
    end
end

"""
get_reversal_hamming(filepath1, filepath2; permuted_order, is_mitigation, is_compute_sem, is_show_progress)

文件版重载：经 `import_random_pair` 导入配对数据再调核心方法。
注意 hamming 方法暂未实现误差缓解，`is_mitigation` 只能为 false。

参数
- filepath1::String：实验一的 .npz 文件路径。
- filepath2::String：实验二的 .npz 文件路径。

关键词参数
- permuted_order：全部被测 site 的置换向量，缺省 `nothing` 表示链式。
  两份实验测的是同一区间，必须共用同一 permuted_order。
- is_mitigation::Bool=false：占位输入，hamming 暂未实现误差缓解，只能为 false。
- is_compute_sem::Bool：为 true 时计算各配对设置间的均值标准误差（SEM）。
- is_show_progress::Bool：为 true 时显示进度。
"""
function get_reversal_hamming(
    filepath1::String,
    filepath2::String;
    permuted_order=nothing,
    is_mitigation=false,
    is_compute_sem=false,
    is_show_progress=false,
)
    @assert !is_mitigation "get_reversal_hamming 暂未实现误差缓解，is_mitigation 只能为 false"
    group1, group2, _, _, _ = import_random_pair(
        filepath1, filepath2; permuted_order, is_mitigation
    )
    return get_reversal_hamming(
        group1,
        group2;
        is_compute_sem, is_show_progress,
    )
end
