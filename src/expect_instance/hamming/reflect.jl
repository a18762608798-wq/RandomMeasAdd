# -------------
# reflect 的 hamming 估计
# -------------
"""
get_reflect_hamming(data::MeasurementData)

对单个 MeasurementData 用 Hamming 距离法算反射期望值。

参数
- data::MeasurementData：单个测量数据集。

返回
- reflect_est::Float64：反射期望的估计值。

说明
估计子公式：R = (1/NM) * Σ_m 2^pairs * (-2)^(-hamming_dist(m))，
其中系统按相邻奇偶 site 配成对。
"""
function get_reflect_hamming(
    data::MeasurementData,
)
    # 取数据
    qubits_num = data.N
    m_num = data.NM
    results = data.measurement_results
    pairs_num = qubits_num ÷ 2
    odd_order = [2i - 1 for i in 1:pairs_num]
    even_order = [2i for i in 1:pairs_num]

    # 算 reflect 估计
    ssum = 0
    for m_idx = 1:m_num
        result = results[m_idx, :]
        odd_result = @view result[odd_order]
        even_result = @view result[even_order]
        hamming_dist = sum(odd_result .!= even_result)
        ssum += (-2.0)^(-hamming_dist)
    end
    ssum *= 2.0^pairs_num

    reflect_est = ssum / m_num

    return reflect_est
end

"""
get_reflect_hamming(group; trivial, is_mitigation, is_compute_sem, is_show_progress)

用 Hamming 距离法从重排好的测量 group 算反射期望值，
对各随机幺正设置求平均。

参数
- group::MeasurementGroup：已重排好的测量 group。

关键词参数
- trivial：已重排好的 trivial（|0>）group（与 group 同 frame、列对齐），
  仅 is_mitigation=true 时需要。
- is_mitigation::Bool=false：为 true 时用 trivial 数据按 Robust Hamming
  公式（RES-Robust_Hamming.md ZR 节）标定逐对核 (x̂,ŷ)；
  无噪声极限退化为标准核，结果与 is_mitigation=false 一致。
- is_compute_sem::Bool：为 true 时计算各随机幺正设置间的均值标准误差（SEM）。
- is_show_progress::Bool：为 true 时显示进度。

返回
- is_compute_sem == false：返回 reflect_est::Float64。
- is_compute_sem == true：返回 (reflect_est::Float64, sem::Float64)。

说明
本函数对每个随机幺正设置调用
get_reflect_hamming(::MeasurementData) 算 Hamming 距离反射估计
（缓解时调 get_reflect_robust)，再对各设置求平均。
"""
function get_reflect_hamming(
    group;
    trivial=nothing,
    is_mitigation=false,
    is_compute_sem=false,
    is_show_progress=false,
)
    # 取数据
    u_num = group.NU
    datas = group.measurements
    qubits_num = group.N
    pairs_num = qubits_num ÷ 2
    odd_order = [2i - 1 for i in 1:pairs_num]
    even_order = [2i for i in 1:pairs_num]
    if is_mitigation
        @assert !isnothing(trivial) "is_mitigation=true 要求传入 trivial group（经 import_random_group 的 is_return_trivial 开关取得）。"
        xs, ys = get_mirror_xy(trivial, odd_order, even_order)
    end
    reflect_ests = Vector{Float64}(undef, u_num)

    # 算 reflect 估计
    ssum = 0
    @showprogress desc="hamming_est..." enabled=is_show_progress @threads for u_idx = 1:u_num
        data = datas[u_idx]
        reflect_ests[u_idx] = is_mitigation ?
            get_reflect_robust(data, odd_order, even_order, xs, ys) :
            get_reflect_hamming(data)
    end

    reflect_est = mean(reflect_ests)

    if is_compute_sem
        sem = std(reflect_ests) / sqrt(u_num)
        return reflect_est, sem
    else
        reflect_est = mean(reflect_ests)
        return reflect_est
    end

end

"""
get_reflect_hamming(filepath; permuted_order, is_mitigation, is_compute_sem, is_show_progress)

文件版重载：输入与 `import_hist_group` 一致，先导入再调核心方法。

参数
- filepath::String：qmeas.random 生成的 .npz 文件路径。

关键词参数
- permuted_order：全部被测 site 的置换向量，缺省 `nothing` 表示链式。
- is_mitigation::Bool=false：为 true 时用 npz 内 trivial 数据做 Robust Hamming
  缓解（要求生成该文件时开 mitigation）；为 false 时忽略 trivial 数据。
- is_compute_sem::Bool：为 true 时计算各随机幺正设置间的均值标准误差（SEM）。
- is_show_progress::Bool：为 true 时显示进度。
"""
function get_reflect_hamming(
    filepath::String;
    permuted_order=nothing,
    is_mitigation=false,
    is_compute_sem=false,
    is_show_progress=false,
)
    group, trivial, _ = import_hist_group(filepath; permuted_order)
    return get_reflect_hamming(
        group;
        trivial, is_mitigation, is_compute_sem, is_show_progress,
    )
end

