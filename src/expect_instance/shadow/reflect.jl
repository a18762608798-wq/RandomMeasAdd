# -------------
# reflect 的 shadow 估计
# -------------

"""
get_reflect_shadow(permuted_group, permuted_indices, G; is_compute_sem, is_show_progress)

用经典 shadow 计算反射算符 Z_r 的期望值。

参数
- permuted_group::MeasurementGroup：已重排好的测量 group。
- permuted_indices：重排后的 site Index 对象。
- G：已处于重排 frame 的校准权重向量（`nothing` 表示全 1）。

关键词参数
- is_compute_sem::Bool：为 true 时同时计算均值标准误差（SEM）。
- is_show_progress::Bool：为 true 时计算过程中显示进度。

返回
- is_compute_sem == false：返回 real(expectation)::Float64。
- is_compute_sem == true：返回 (real(expectation)::Float64, sem::Float64)。

说明
本函数为重排后的系统构造 dense shadow，
再构造反射用的相邻 swap 算符，把期望/SEM 估计交给
modified_get_expect_shadow。
"""
function get_reflect_shadow(
    permuted_group,
    permuted_indices,
    G=nothing;
    is_compute_sem=false,
    is_show_progress=false,
)
    permuted_G = isnothing(G) ? ones(length(permuted_indices)) : G
    permuted_shadows = get_dense_shadows(permuted_group; G=permuted_G)
    adjacent_swap_op = create_adjacent_swap_op(permuted_indices)

    if is_compute_sem
        reflect_expect, sem = modified_get_expect_shadow(
            adjacent_swap_op,
            permuted_shadows;
            is_compute_sem=is_compute_sem,
            is_show_progress=is_show_progress,
        )
        return real(reflect_expect), sem
    else
        reflect_expect = modified_get_expect_shadow(
            adjacent_swap_op,
            permuted_shadows;
            is_compute_sem=is_compute_sem,
            is_show_progress=is_show_progress,
        )
        return real(reflect_expect)
    end
end

"""
get_reflect_shadow(filepath; permuted_order, is_mitigation, is_compute_sem, is_show_progress)

文件版重载：输入与 `import_random_group` 一致，先导入再调核心方法。

参数
- filepath::String：存好的 shadow/group 数据路径。

关键词参数
- permuted_order：在算 shadow 之前对 site 做的置换顺序，缺省 `nothing` 表示链式。
- is_mitigation::Bool=false：为 true 时用 trivial 数据算校准向量 G。
- is_compute_sem::Bool：为 true 时同时计算均值标准误差（SEM）。
- is_show_progress::Bool：为 true 时计算过程中显示进度。
"""
function get_reflect_shadow(
    filepath::String;
    permuted_order=nothing,
    is_mitigation=false,
    is_compute_sem=false,
    is_show_progress=false,
)
    permuted_group, permuted_indices, G = import_random_group(
        filepath; permuted_order, is_mitigation
    )
    return get_reflect_shadow(
        permuted_group, permuted_indices, G;
        is_compute_sem, is_show_progress,
    )
end
