# -------------
# reflect 的 shadow 估计
# -------------

"""
get_reflect_shadow(hgroup, ξ, G; is_compute_sem, is_show_progress)

用自有 shadow 计算反射算符的期望值。

参数
- hgroup::HistGroup：已重排好的直方图 group。
- ξ::Vector{Index{Int64}}：重排后的 site index 向量。
- G：已处于重排 frame 的校准权重向量（`nothing` 表示全 1）。

关键词参数
- is_compute_sem::Bool：为 true 时同时计算均值标准误差（SEM）。
- is_show_progress::Bool：为 true 时计算过程中显示进度。

返回
- is_compute_sem == false：返回 real(expectation)::Float64。
- is_compute_sem == true：返回 (real(expectation)::Float64, sem::Float64)。

说明
本函数为重排后的系统构造自有 shadow，再构造反射用的相邻 swap 算符，
逐 setting 求期望后平均（SEM 取设置间标准误）。
"""
function get_reflect_shadow(
    hgroup::HistGroup,
    ξ::Vector{Index{Int64}},
    G=nothing;
    is_compute_sem=false,
    is_show_progress=false,
)
    permuted_G = isnothing(G) ? ones(length(ξ)) : G
    shadows = get_dense_shadows(hgroup, ξ; G=permuted_G)
    adjacent_swap_op = create_adjacent_swap_op(ξ)

    exps = [real(get_expect_shadow(adjacent_swap_op, sh)) for sh in shadows]
    reflect_expect = mean(exps)
    if is_compute_sem
        sem = std(exps) / sqrt(length(exps))
        return reflect_expect, sem
    else
        return reflect_expect
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
    hgroup, htrivial, _ = import_hist_group(filepath; permuted_order)
    ξ = siteinds("Qubit", hgroup.N)
    @assert !is_mitigation || !isnothing(htrivial) "is_mitigation=true 要求 npz 含 trivial hist（生成时须开 mitigation）。"
    G = is_mitigation ? modified_get_calibration_vector(htrivial) : nothing
    return get_reflect_shadow(
        hgroup, ξ, G;
        is_compute_sem, is_show_progress,
    )
end
