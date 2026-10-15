# -------------
# reversal（时间反演交叠 Z_T 分子）的 shadow 估计
# -------------
# 物理背景（Elben 等，Sci. Adv. 6, eaaz3666 (2020)，Eq. 5）：
# Z_T = Tr(ρ_I u_T ρ_I^{T_1} u_T†)，u_T = ∏_{i∈I_1} σ^y_i，
# T_1 为对 I_1 区的部分转置。
# shadow 路线用两份独立 shadow 估计双线性量：
# E[Tr(ρ̂_a u_T ρ̂_b^{T_1} u_T†)] = Z_T（a≠b），
# 只需普通单实验 independent 数据，不需要共轭配对。
# 按本项目惯例：重排 frame 下奇位 = I_1 区，偶位 = I_2 区。
"""
get_reversal_shadow(hgroup, ξ, G; is_compute_sem, is_show_progress)

用自有 shadow 算时间反演交叠 Z_T 分子。

参数
- hgroup::HistGroup：已重排好的直方图 group（independent 数据）。
- ξ::Vector{Index{Int64}}：重排后的 site index 向量。
- G：已处于重排 frame 的校准权重向量（`nothing` 表示全 1）。

关键词参数
- is_compute_sem::Bool：为 true 时同时计算均值标准误差（SEM）。
- is_show_progress::Bool：为 true 时显示进度。

返回
- is_compute_sem == false：返回 zt_est::Float64。
- is_compute_sem == true：返回 (zt_est::Float64, sem::Float64)。

说明
本函数为重排后的系统构造自有 shadow，再构造 u_T 算符，
经 get_reversal_loos_shadow 算估计和 jackknife SEM
（k=2 路线，与 purity shadow 对齐）。
"""
function get_reversal_shadow(
    hgroup::HistGroup,
    ξ::Vector{Index{Int64}},
    G=nothing;
    is_compute_sem=false,
    is_show_progress=false,
)
    qubits_num = length(ξ)
    @assert iseven(qubits_num) "bitsnum must be even"
    pairs_num = qubits_num ÷ 2
    odd_order = [2i - 1 for i in 1:pairs_num]

    permuted_G = isnothing(G) ? ones(qubits_num) : G
    shadows = get_dense_shadows(hgroup, ξ; G=permuted_G)
    uT = create_unitary_part_reversal_op(
        odd_order, ξ; op_type="ITensor"
    )

    n_ru = length(shadows)
    zt_est, loos = get_reversal_loos_shadow(
        shadows, uT, odd_order; is_show_progress=is_show_progress
    )

    if is_compute_sem
        variance = (n_ru - 1)^2 / n_ru * var(loos)
        sem = sqrt(variance)
        return zt_est, sem
    else
        return zt_est
    end
end

"""
get_reversal_shadow(filepath; permuted_order, is_mitigation, is_compute_sem, is_show_progress)

文件版重载：输入与 `import_random_group` 一致，先导入再调核心方法。
注意 shadow 路线用普通单实验 independent 数据，不需要配对文件。

参数
- filepath::String：qmeas.random 生成的 .npz 文件路径。

关键词参数
- permuted_order：全部被测 site 的置换向量，缺省 `nothing` 表示链式。
- is_mitigation::Bool=false：为 true 时用 trivial 数据算校准向量 G。
- is_compute_sem::Bool：为 true 时同时计算均值标准误差（SEM）。
- is_show_progress::Bool：为 true 时显示进度。
"""
function get_reversal_shadow(
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
    return get_reversal_shadow(
        hgroup, ξ, G;
        is_compute_sem, is_show_progress,
    )
end
