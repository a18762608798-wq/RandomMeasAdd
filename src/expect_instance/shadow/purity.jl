# ----------
# purity（纯度）
# ----------
"""
get_purity_shadow(permuted_group, G; is_compute_sem, is_show_progress)

用经典 shadow 从重排好的测量 group 估计纯度 Tr(ρ²)。

参数
- permuted_group::MeasurementGroup：已重排好的测量 group。
- G：已处于重排 frame 的校准权重向量（`nothing` 表示全 1）。

关键词参数
- is_compute_sem::Bool：为 true 时同时计算均值标准误差（SEM）。
- is_show_progress::Bool：为 true 时显示进度。

返回
- is_compute_sem == false：返回 purity::Float64。
- is_compute_sem == true：返回 (purity::Float64, sem::Float64)。

说明
本函数为重排好的 group 构造 dense shadow，
再交给 modified_get_purity_shadow 做纯度估计。
"""
function get_purity_shadow(
    permuted_group,
    G=nothing;
    is_compute_sem=false,
    is_show_progress=false,
)
    n_site = permuted_group.N
    permuted_G = isnothing(G) ? ones(n_site) : G
    shadows = get_dense_shadows(permuted_group; G=permuted_G)

    if is_compute_sem
        purity, _, sem = modified_get_purity_shadow(
            shadows;
            is_compute_sem=is_compute_sem,
            is_show_progress=is_show_progress,
        )
        return purity, sem

    else
        purity = modified_get_purity_shadow(
            shadows;
            is_compute_sem=is_compute_sem,
            is_show_progress=is_show_progress,
        )
        return purity

    end
end

"""
get_purity_shadow(filepath; permuted_order, is_mitigation, is_compute_sem, is_show_progress)

文件版重载：输入与 `import_random_group` 一致，先导入再调核心方法。

参数
- filepath::String：qmeas.random 生成的 .npz 文件路径。

关键词参数
- permuted_order：全部被测 site 的置换向量，缺省 `nothing` 表示链式（按 npz 列顺序）。
- is_mitigation::Bool=false：为 true 时用 trivial 数据算校准向量 G。
- is_compute_sem::Bool：为 true 时同时计算均值标准误差（SEM）。
- is_show_progress::Bool：为 true 时显示进度。
"""
function get_purity_shadow(
    filepath::String;
    permuted_order=nothing,
    is_mitigation=false,
    is_compute_sem=false,
    is_show_progress=false,
)
    permuted_group, _, G = import_random_group(
        filepath; permuted_order, is_mitigation
    )
    return get_purity_shadow(
        permuted_group, G;
        is_compute_sem, is_show_progress,
    )
end
