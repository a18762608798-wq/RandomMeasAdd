# -------------
# z_r 的 hamming 估计
# -------------
"""
get_z_r_hamming(group; is_compute_sem, is_show_progress)

用 Hamming 距离法从重排好的测量 group 算 z_r，
对各随机幺正设置求平均。

参数
- group::MeasurementGroup：已重排好的测量 group。

关键词参数
- is_compute_sem::Bool：为 true 时同时返回 jackknife 偏差估计和 SEM。
- is_show_progress::Bool：为 true 时显示进度。

返回
- is_compute_sem == false：返回 z_r_val::Float64。
- is_compute_sem == true：返回三元组 (z_r_val::Float64, bias_estimate::Float64, sem::Float64)。

说明
本函数对每个设置分别算 reflect 估计与奇偶子系统纯度估计，
再经 get_z_r_loos_hamming 按比值组合（含 setting 层面 jackknife）。
"""
function get_z_r_hamming(
    group;
    is_compute_sem=false,
    is_show_progress=false,
)
    # 按相邻配对拆成奇偶两个子系统
    qubits_num = group.N
    pairs_num = qubits_num ÷ 2
    odd_order = [2i - 1 for i in 1:pairs_num]
    even_order = [2i for i in 1:pairs_num]
    odd_group = reduce_to_subsystem(group, odd_order)
    even_group = reduce_to_subsystem(group, even_order)

    # 每个设置各算一路估计
    u_num = group.NU
    datas = group.measurements
    odd_datas = odd_group.measurements
    even_datas = even_group.measurements
    reflect_ests = Vector{Float64}(undef, u_num)
    odd_ests = Vector{Float64}(undef, u_num)
    even_ests = Vector{Float64}(undef, u_num)

    @showprogress desc="hamming_est..." enabled=is_show_progress @threads for u_idx = 1:u_num
        reflect_ests[u_idx] = get_reflect_hamming(datas[u_idx])
        odd_data = odd_datas[u_idx]
        odd_ests[u_idx] = get_overlap(odd_data, odd_data; apply_bias_correction=true)
        even_data = even_datas[u_idx]
        even_ests[u_idx] = get_overlap(even_data, even_data; apply_bias_correction=true)
    end

    z_r_est, z_r_loos = get_z_r_loos_hamming(reflect_ests, odd_ests, even_ests)

    if is_compute_sem
        variance = (u_num - 1)^2 / u_num * var(z_r_loos)
        sem = sqrt(variance)
        z_r_jack = u_num * z_r_est - (u_num - 1) * mean(z_r_loos)
        return z_r_est, z_r_est - z_r_jack, sem
    else
        return z_r_est
    end
end

"""
get_z_r_hamming(filepath; permuted_order, is_mitigation, is_compute_sem, is_show_progress)

文件版重载：输入与 `import_random_group` 一致，先导入再调核心方法。
注意 hamming 方法暂未实现误差缓解，`is_mitigation` 只能为 false。

参数
- filepath::String：qmeas.random 生成的 .npz 文件路径。

关键词参数
- permuted_order：全部被测 site 的置换向量，缺省 `nothing` 表示链式。
- is_mitigation::Bool=false：占位输入，hamming 暂未实现误差缓解，只能为 false。
- is_compute_sem::Bool：为 true 时同时返回 jackknife 偏差估计和 SEM。
- is_show_progress::Bool：为 true 时显示进度。

说明
filepath 也可传文件夹：此时按文件名排序对其中每个 .npz 依次计算，
返回 (vals::Vector{Float64}, sems::Vector{Float64}) 两个列表
（bias 不保留）。
"""
function get_z_r_hamming(
    filepath::String;
    permuted_order=nothing,
    is_mitigation=false,
    is_compute_sem=false,
    is_show_progress=false,
)
    @assert !is_mitigation "get_z_r_hamming 暂未实现误差缓解，is_mitigation 只能为 false"
    if isdir(filepath)
        files = sort(filter(
            f -> endswith(f, ".npz"), readdir(filepath; join=true)
        ))
        @assert !isempty(files) "文件夹中没有 .npz 文件：$filepath"
        vals = Vector{Float64}(undef, length(files))
        sems = Vector{Float64}(undef, length(files))
        for (i, f) in enumerate(files)
            v, _, s = get_z_r_hamming(
                f; permuted_order, is_mitigation,
                is_compute_sem=true, is_show_progress,
            )
            vals[i], sems[i] = v, s
        end
        return vals, sems
    end
    group, _, _ = import_random_group(
        filepath; permuted_order, is_mitigation
    )
    return get_z_r_hamming(
        group;
        is_compute_sem, is_show_progress,
    )
end
