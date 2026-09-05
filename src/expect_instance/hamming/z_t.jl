# -------------
# z_t 的 hamming 估计
# -------------
# 物理背景（Elben 等，Sci. Adv. 6, eaaz3666 (2020)，Eq. 5）：
# 归一化时间反演不变量 z_t = Z_T / ((Tr ρ_{I_1}^2 + Tr ρ_{I_2}^2)/2)^(3/2)，
# trivial 相 +1、拓扑相 -1、对称破缺时为 0。
# 分子 Z_T 取自两份配对实验的互关联（见 reversal.jl），
# 分母纯度按论文用同一批实验数据经 Eq. 4 估计，此处取两份实验的平均。
# 布局约定：重排 frame 下 I_1、I_2 须奇偶相间（奇位 = I_1、偶位 = I_2），
# 生成数据时 meas_indices 取 [(i1_1,), (i2_1,), ...] 交错排列。
"""
get_z_t_hamming(group1, group2; is_compute_sem, is_show_progress)

用 Hamming 距离法从两份重排好的测量 group 算 z_t，
对各配对设置求平均后再组合。

参数
- group1：实验一的重排 group。
- group2：实验二的重排 group，设置须与 group1 逐行配对。

关键词参数
- is_compute_sem::Bool：为 true 时同时返回 jackknife 偏差估计和 SEM。
- is_show_progress::Bool：为 true 时显示进度。

返回
- is_compute_sem == false：返回 z_t_val::Float64。
- is_compute_sem == true：返回三元组 (z_t_val::Float64, bias_estimate::Float64, sem::Float64)。

说明
本函数对每个配对设置分别算互关联估计与两份实验平均的奇偶子系统纯度估计，
再经 get_z_t_loos_hamming 按比值组合（含 setting 层面 jackknife）。
"""
function get_z_t_hamming(
    group1,
    group2;
    is_compute_sem=false,
    is_show_progress=false,
)
    # 两份数据的几何与设置数必须一致，否则互关联无意义
    @assert group1.N == group2.N "两份实验的比特数必须一致。"
    @assert group1.NU == group2.NU "两份实验的设置数必须一致且逐行配对。"

    # 按相邻配对拆出 I_1（奇位）、I_2（偶位）两个子系统，
    # 纯度取两份实验的平均（都是同一 rho 的合法 Haar 数据）
    qubits_num = group1.N
    pairs_num = qubits_num ÷ 2
    odd_order = [2i - 1 for i in 1:pairs_num]
    even_order = [2i for i in 1:pairs_num]
    odd_group1 = reduce_to_subsystem(group1, odd_order)
    even_group1 = reduce_to_subsystem(group1, even_order)
    odd_group2 = reduce_to_subsystem(group2, odd_order)
    even_group2 = reduce_to_subsystem(group2, even_order)

    # 每个配对设置各算一路估计
    u_num = group1.NU
    datas1 = group1.measurements
    datas2 = group2.measurements
    odd_datas1 = odd_group1.measurements
    even_datas1 = even_group1.measurements
    odd_datas2 = odd_group2.measurements
    even_datas2 = even_group2.measurements
    zt_ests = Vector{Float64}(undef, u_num)
    odd_ests = Vector{Float64}(undef, u_num)
    even_ests = Vector{Float64}(undef, u_num)

    @showprogress desc="hamming_est..." enabled=is_show_progress @threads for u_idx = 1:u_num
        zt_ests[u_idx] = get_reversal_hamming(datas1[u_idx], datas2[u_idx])
        odd_ests[u_idx] = (
            get_overlap(odd_datas1[u_idx], odd_datas1[u_idx]; apply_bias_correction=true) +
            get_overlap(odd_datas2[u_idx], odd_datas2[u_idx]; apply_bias_correction=true)
        ) / 2
        even_ests[u_idx] = (
            get_overlap(even_datas1[u_idx], even_datas1[u_idx]; apply_bias_correction=true) +
            get_overlap(even_datas2[u_idx], even_datas2[u_idx]; apply_bias_correction=true)
        ) / 2
    end

    z_t_est, z_t_loos = get_z_t_loos_hamming(zt_ests, odd_ests, even_ests)

    if is_compute_sem
        variance = (u_num - 1)^2 / u_num * var(z_t_loos)
        sem = sqrt(variance)
        z_t_jack = u_num * z_t_est - (u_num - 1) * mean(z_t_loos)
        return z_t_est, z_t_est - z_t_jack, sem
    else
        return z_t_est
    end
end

"""
get_z_t_hamming(filepath1, filepath2; permuted_order, is_mitigation, is_compute_sem, is_show_progress)

文件版重载：经 `import_random_pair` 导入配对数据再调核心方法。
注意 hamming 方法暂未实现误差缓解，`is_mitigation` 只能为 false。

参数
- filepath1::String：实验一的 .npz 文件路径。
- filepath2::String：实验二的 .npz 文件路径。

关键词参数
- permuted_order：全部被测 site 的置换向量，缺省 `nothing` 表示链式。
  两份实验测的是同一区间，必须共用同一 permuted_order。
- is_mitigation::Bool=false：占位输入，hamming 暂未实现误差缓解，只能为 false。
- is_compute_sem::Bool：为 true 时同时返回 jackknife 偏差估计和 SEM。
- is_show_progress::Bool：为 true 时显示进度。

说明
filepath1 也可传文件夹：此时把其中文件名含 exp1/exp2 的 .npz
各自排序后按顺序配对（要求两两文件名除 exp 标记外完全一致），
逐对计算并返回 (vals::Vector{Float64}, sems::Vector{Float64}) 两个列表
（bias 不保留）。此时 filepath2 被忽略。
"""
function get_z_t_hamming(
    filepath1::String,
    filepath2::String="";
    permuted_order=nothing,
    is_mitigation=false,
    is_compute_sem=false,
    is_show_progress=false,
)
    @assert !is_mitigation "get_z_t_hamming 暂未实现误差缓解，is_mitigation 只能为 false"
    if isdir(filepath1)
        all_files = sort(filter(
            f -> endswith(f, ".npz"), readdir(filepath1; join=true)
        ))
        files1 = filter(f -> occursin("exp1", basename(f)), all_files)
        files2 = filter(f -> occursin("exp2", basename(f)), all_files)
        @assert !isempty(files1) && length(files1) == length(files2) "配对文件数量不一致或为空：$filepath1"
        @assert all(
            replace(basename(a), "exp1" => "exp2") == basename(b)
            for (a, b) in zip(files1, files2)
        ) "exp1/exp2 文件名除 exp 标记外必须一一对应。"
        vals = Vector{Float64}(undef, length(files1))
        sems = Vector{Float64}(undef, length(files1))
        for (i, (a, b)) in enumerate(zip(files1, files2))
            v, _, s = get_z_t_hamming(
                a, b; permuted_order, is_mitigation,
                is_compute_sem=true, is_show_progress,
            )
            vals[i], sems[i] = v, s
        end
        return vals, sems
    end
    group1, group2, _, _, _ = import_random_pair(
        filepath1, filepath2; permuted_order, is_mitigation
    )
    return get_z_t_hamming(
        group1,
        group2;
        is_compute_sem, is_show_progress,
    )
end


