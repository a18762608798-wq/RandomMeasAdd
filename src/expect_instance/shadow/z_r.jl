# --------------------------------------
# ----------反射量 (Z_r)----------
# --------------------------------------

"""
get_z_r_shadow(permuted_group, permuted_indices, G; is_compute_sem, is_show_progress)

用经典 shadow 估计 Z_r 量（反射相关的观测量），
做法是把系统按相邻 site 配对，再由这些配对算出 Z_r。

参数
- permuted_group::MeasurementGroup：已重排好的测量 group。
- permuted_indices：重排后的 site Index 对象。
- G：已处于重排 frame 的校准权重向量（`nothing` 表示全 1）；
  子系统权重由本函数按奇偶位置索引得到，不需要调用方预拆。

关键词参数
- is_compute_sem::Bool
    为 true 时同时返回 jackknife 偏差估计和 SEM。
- is_show_progress::Bool
    为 true 时计算过程中显示进度。

返回
- is_compute_sem == false：返回 z_r_val::Float64（估计值）。
- is_compute_sem == true：返回三元组 (z_r_val::Float64, bias_estimate::Float64, sem::Float64)，
  其中 bias_estimate = z_r_val - z_r_jack（由 jackknife 值算出），sem 为标准误差。

说明
本函数：
- 把 qubit 按相邻配对拆成奇偶两个子系统，
- 为全系统和子系统分别构造 dense shadow，
- 构造相邻 swap 算符，
- 经 get_z_r_loos_shadow 算 jackknife 值，得到估计和（可选的）SEM。
"""
function get_z_r_shadow(
    permuted_group,
    permuted_indices,
    G=nothing;
    is_compute_sem=false,
    is_show_progress=false,
)
    # 取三个系统的信息
    qubits_num = length(permuted_indices)
    pairs_num = qubits_num ÷ 2
    odd_order = [2i - 1 for i in 1:pairs_num]
    even_order = [2i for i in 1:pairs_num]
    odd_group = reduce_to_subsystem(permuted_group, odd_order)
    even_group = reduce_to_subsystem(permuted_group, even_order)
    # 取每个系统的 G（G 已处于重排 frame，只需按奇偶位置索引）
    permuted_G = isnothing(G) ? ones(qubits_num) : G
    odd_G = permuted_G[odd_order]
    even_G = permuted_G[even_order]
    # 生成 shadow
    shadows = get_dense_shadows(permuted_group; G=permuted_G)
    odd_shadows = get_dense_shadows(odd_group; G=odd_G)
    even_shadows = get_dense_shadows(even_group; G=even_G)
    # 生成算符
    adjacent_swap_op = create_adjacent_swap_op(permuted_indices)

    # 算期望和 sem
    # 取 jackknife 信息
    n_ru = size(shadows, 1)
    z_r_est, z_r_loos = get_z_r_loos_shadow(
        shadows, odd_shadows, even_shadows, adjacent_swap_op, is_show_progress
    )
    # 取 sem
    if is_compute_sem
        variance = (n_ru - 1)^2 / n_ru * var(z_r_loos)
        sem = sqrt(variance)
        z_r_jack = n_ru * z_r_est - (n_ru - 1) * mean(z_r_loos)
        return z_r_est, z_r_est - z_r_jack, sem
    else
        return z_r_est
    end
end

"""
get_z_r_shadow(filepath; permuted_order, is_mitigation, is_compute_sem, is_show_progress)

文件版重载：输入与 `import_random_group` 一致，先导入再调核心方法。

参数
- filepath::String：qmeas.random 生成的单个 .npz 文件路径。

关键词参数
- permuted_order：全部被测 site 的置换向量，缺省 `nothing` 表示链式。
- is_mitigation::Bool=false：为 true 时用 trivial 数据算校准向量 G。
- is_compute_sem::Bool：为 true 时同时返回 jackknife 偏差估计和 SEM。
- is_show_progress::Bool：为 true 时显示进度。

说明
filepath 也可传文件夹：此时按文件名排序对其中每个 .npz 依次计算，
返回 (vals::Vector{Float64}, sems::Vector{Float64}) 两个列表
（bias 不保留）。
"""
function get_z_r_shadow(
    filepath::String;
    permuted_order=nothing,
    is_mitigation=false,
    is_compute_sem=false,
    is_show_progress=false,
)
    if isdir(filepath)
        files = sort(filter(
            f -> endswith(f, ".npz"), readdir(filepath; join=true)
        ))
        @assert !isempty(files) "文件夹中没有 .npz 文件：$filepath"
        vals = Vector{Float64}(undef, length(files))
        sems = Vector{Float64}(undef, length(files))
        for (i, f) in enumerate(files)
            v, _, s = get_z_r_shadow(
                f; permuted_order, is_mitigation,
                is_compute_sem=true, is_show_progress,
            )
            vals[i], sems[i] = v, s
        end
        return vals, sems
    end
    permuted_group, permuted_indices, G = import_random_group(
        filepath; permuted_order, is_mitigation
    )
    return get_z_r_shadow(
        permuted_group, permuted_indices, G;
        is_compute_sem, is_show_progress,
    )
end


