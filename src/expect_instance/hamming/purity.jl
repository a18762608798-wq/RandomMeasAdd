"""
get_purity_hamming(group; trivial, is_mitigation, is_compute_sem, is_show_progress)

用基于重叠的方法（"hamming" 变体）从重排好的测量 group 估计纯度 Tr(ρ²)，
对各随机幺正设置求平均。

参数
- group::MeasurementGroup：已重排好的测量 group。

关键词参数
- trivial：已重排好的 trivial（|0>）group（与 group 同 frame、列对齐），
  仅 is_mitigation=true 时需要（通常经 import_random_group 的
  is_return_trivial 开关取得）；为 nothing 时表示无缓解。
- is_mitigation::Bool=false：为 true 时用 trivial 数据按 Robust Hamming
  公式（RES-Robust_Hamming.md purity 节）标定逐比特核 (x̂,ŷ)；
  无噪声极限退化为标准核，结果与 is_mitigation=false 一致。
- is_compute_sem::Bool：为 true 时计算各随机幺正设置间的均值标准误差（SEM）。
- is_show_progress::Bool：为 true 时显示进度。

返回
- is_compute_sem == false：返回 purity_est::Float64。
- is_compute_sem == true：返回 (purity_est::Float64, sem::Float64)。

说明
本函数对每个随机幺正设置调用
get_overlap(data, data; apply_bias_correction=true) 算纯度
（缓解时调 get_overlap_robust(data, data, xs, ys)），
再对各设置求平均。
"""
function get_purity_hamming(
    group;
    trivial=nothing,
    is_mitigation=false,
    is_compute_sem=false,
    is_show_progress=false,
)
    if is_mitigation
        @assert !isnothing(trivial) "is_mitigation=true 要求传入 trivial group（经 import_random_group 的 is_return_trivial 开关取得）。"
        xs, ys = get_purity_xy(trivial, 1:group.N)
    end
    u_num = group.NU
    datas = group.measurements
    purity_ests = Vector{Float64}(undef, u_num)

    @showprogress desc="hamming_est..." enabled=is_show_progress @threads for u_idx = 1:u_num
        data = datas[u_idx]
        purity_ests[u_idx] = is_mitigation ?
            get_overlap_robust(data, data, xs, ys; apply_bias_correction=true) :
            get_overlap(data, data, apply_bias_correction=true)
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

参数
- filepath::String：qmeas.random 生成的 .npz 文件路径。

关键词参数
- permuted_order：全部被测 site 的置换向量，缺省 `nothing` 表示链式（按 npz 列顺序）。
- is_mitigation::Bool=false：为 true 时用 npz 内 trivial 数据做 Robust Hamming
  缓解（要求生成该文件时开 mitigation）；为 false 时忽略 trivial 数据。
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
    group, _, _, trivial = import_random_group(
        filepath; permuted_order, is_mitigation, is_return_trivial=true
    )
    return get_purity_hamming(
        group;
        trivial, is_mitigation, is_compute_sem, is_show_progress,
    )
end


