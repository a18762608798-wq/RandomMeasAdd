"""
get_purity_hamming(group; is_compute_sem, is_show_progress)

用基于重叠的方法（"hamming" 变体）从重排好的测量 group 估计纯度 Tr(ρ²)，
对各随机幺正设置求平均。

参数
- group::MeasurementGroup：已重排好的测量 group。

关键词参数
- is_compute_sem::Bool：为 true 时计算各随机幺正设置间的均值标准误差（SEM）。
- is_show_progress::Bool：为 true 时显示进度。

返回
- is_compute_sem == false：返回 purity_est::Float64。
- is_compute_sem == true：返回 (purity_est::Float64, sem::Float64)。

说明
本函数对每个随机幺正设置调用
get_overlap(data, data; apply_bias_correction=true) 算纯度，
再对各设置求平均。
"""
function get_purity_hamming(
    group;
    is_compute_sem=false,
    is_show_progress=false,
)
    u_num = group.NU
    datas = group.measurements
    purity_ests = Vector{Float64}(undef, u_num)

    @showprogress desc="hamming_est..." enabled=is_show_progress @threads for u_idx = 1:u_num
        data = datas[u_idx]
        purity_ests[u_idx] = get_overlap(data, data, apply_bias_correction=true)
    end

    purity_est = mean(purity_ests)

    if is_compute_sem
        sem = std(purity_ests) / sqrt(u_num)
        return purity_est, sem

    else
        purity_est = mean(purity_ests)
        return purity_est

    end
end

"""
get_purity_hamming(filepath; permuted_order, is_mitigation, is_compute_sem, is_show_progress)

文件版重载：输入与 `import_random_group` 一致，先导入再调核心方法。
注意 hamming 方法暂未实现误差缓解，`is_mitigation` 只能为 false。

参数
- filepath::String：qmeas.random 生成的 .npz 文件路径。

关键词参数
- permuted_order：全部被测 site 的置换向量，缺省 `nothing` 表示链式（按 npz 列顺序）。
- is_mitigation::Bool=false：占位输入，hamming 暂未实现误差缓解，只能为 false。
- is_compute_sem::Bool：为 true 时计算各随机幺正设置间的均值标准误差（SEM）。
- is_show_progress::Bool：为 true 时显示进度。
"""
function get_purity_hamming(
    filepath::String;
    permuted_order=nothing,
    is_mitigation=false,
    is_compute_sem=false,
    is_show_progress=false,
)
    @assert !is_mitigation "get_purity_hamming 暂未实现误差缓解，is_mitigation 只能为 false"
    group, _, _ = import_random_group(
        filepath; permuted_order, is_mitigation
    )
    return get_purity_hamming(
        group;
        is_compute_sem, is_show_progress,
    )
end


