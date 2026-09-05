# --------------------
# 求期望 shadow
# --------------------

"""
modified_get_expect_shadow(O, shadows; is_compute_sem=false, is_show_progress=false)

用经典 shadow 数组估计算符 O 的期望值。

参数
- O::MPO
    待估计期望的算符（MPO 形式）。
- shadows::AbstractArray{<:AbstractShadow, 2}
    (settings_num, shots) 的二维 shadow 数组，每个元素是
    对应随机幺正和 shot 的经典 shadow（AbstractShadow）。

关键词参数
- is_compute_sem::Bool=false
    为 true 时返回 (mean, sem)，sem 为各随机幺正设置间的均值标准误差。
- is_show_progress::Bool=true
    为 true 时内层循环显示进度条。

返回
- is_compute_sem == false：返回 mean_value::ComplexF64（期望均值）。
- is_compute_sem == true：返回 (mean_value::ComplexF64, sem_value::Float64)。

说明
本函数对每个 (setting, shot) 对调 get_expect_shadow 求期望，
再对 shot 和 setting 求平均，需要时算各设置间的标准误差（SEM）。
另有一维向量版重载，会 reshape 成二维再调本函数。
"""
function modified_get_expect_shadow(
    O::MPO,
    shadows::AbstractArray{<:AbstractShadow,2};
    is_compute_sem::Bool=false,
    is_show_progress::Bool=true,
)
    # 确保 shadow 数组非空且为矩阵
    @assert !isempty(shadows) "Array of shadows is empty."

    # 算所有期望值
    settings_num, shots = size(shadows)
    expect_values = Matrix{ComplexF64}(undef, settings_num, shots)
    @showprogress desc="Expectation Processing..." enabled=is_show_progress @threads for settings in
                                                                                      1:settings_num

        for shot in 1:shots
            shadow = shadows[settings, shot]
            expect_values[settings, shot] = get_expect_shadow(O, shadow)
        end
    end

    # 求均值（对每个 setting 数）
    mean_values = mean(expect_values; dims=2)
    mean_value = mean(mean_values)

    if is_compute_sem
        # 算均值标准误差（SEM）
        sem_value = std(mean_values) / sqrt(settings_num)
        return mean_value, sem_value
    else
        return mean_value
    end
end

"""
modified_get_expect_shadow(O, shadows::AbstractShadow[]; kwargs...)

一维 shadow 向量版重载，内部 reshape 成单列二维数组再调二维实现。
"""
function modified_get_expect_shadow(
    O::MPO, shadows::AbstractArray{<:AbstractShadow,1}; kwargs...
)
    return modified_get_expect_shadow(O, reshape(shadows, length(shadows), 1); kwargs...)
end

# --------------------
# 求 trace 矩
# --------------------

"""
modified_get_trace_moment(shadows, kth_moment; O=nothing, is_compute_sem=false, is_compute_renyi=false, is_show_progress=false)

用经典 shadow 估计 k 阶 trace 矩（O==nothing 时如 Tr(ρ^k)；给了 O 则为广义 trace 积）。

参数
- shadows::Array{<:AbstractShadow, 2}
    (n_ru, n_m) 的二维 shadow 数组，n_ru 为随机幺正设置数，
    n_m 为每个设置下的测量 shot 数。
- kth_moment::Int
    要估计的矩 k。

关键词参数
- O::Union{Nothing, MPO}=nothing
    插入 trace 积的可选算符。O==nothing 时估计态的纯 trace 矩。
- is_compute_sem::Bool=false
    为 true 时计算并返回 jackknife 偏差和标准误差。
- is_compute_renyi::Bool=false
    为 true 时把估计换算成 Rényi 熵形式（用 (1/(1-k)) * log2(mean)）。
- is_show_progress::Bool=true
    是否显示进度。

返回
- is_compute_sem == false：返回标量矩估计::Float64。
- is_compute_sem == true：返回 (estimate::Float64, bias::Float64, sem::Float64)。

说明
本函数转调支持向量 k 输入和协方差估计的 modified_get_trace_moments。
一维 shadows 重载会 reshape 成二维再调。
"""
function modified_get_trace_moment(
    shadows::Array{<:AbstractShadow,2},
    kth_moment::Int;
    O::Union{Nothing,MPO}=nothing,
    is_compute_sem::Bool=false,
    is_compute_renyi::Bool=false,
    is_show_progress::Bool=true,
)
    if is_compute_sem
        s, bias, cov = modified_get_trace_moments(
            shadows,
            [kth_moment];
            O=O,
            is_compute_cov=is_compute_sem,
            is_compute_renyi=is_compute_renyi,
            is_show_progress=is_show_progress,
        )
        return s[1], bias[1], sqrt(cov[1, 1])
    else
        s = modified_get_trace_moments(
            shadows,
            [kth_moment];
            O=O,
            is_compute_cov=is_compute_sem,
            is_compute_renyi=is_compute_renyi,
            is_show_progress=is_show_progress,
        )
        return s[1]
    end
end

"""
modified_get_trace_moments(shadows, k_vec; O=nothing, is_compute_cov=false, is_compute_renyi=false, is_show_progress=false)

估计 k_vec 指定的多个 trace 矩，返回估计向量，需要时返回 jackknife 偏差和协方差矩阵。

算法
- 枚举 k 个不同随机幺正设置的全部置换，以及这些设置上 shot 的笛卡尔积。
- 对每个置换算 trace 积（经 get_trace_product）再做平均。
- is_compute_renyi 时做 Rényi 换算。
- is_compute_cov 为 true 时，逐个去掉一个幺正算 jackknife 值，
  由此构造协方差矩阵和偏差修正估计。

参数
- shadows::Array{<:AbstractShadow, 2}
    (n_ru, n_m) 排列的 shadow。
- k_vec::Vector{Int}
    要算的整数矩列表。

关键词参数
- O::Union{Nothing, MPO}=nothing
- is_compute_cov::Bool=false
- is_compute_renyi::Bool=false
- is_show_progress::Bool=true

返回
- is_compute_cov == false：返回长度 nK 的 θ_est::Vector{Float64}。
- is_compute_cov == true：返回 (θ_est, bias_vec, Σ)，其中 bias_vec = θ_est - θ_jack，Σ 为协方差矩阵。
"""
function modified_get_trace_moments(
    shadows::Array{<:AbstractShadow,2},
    k_vec::Vector{Int};
    O::Union{Nothing,MPO}=nothing,
    is_compute_cov::Bool=false,
    is_compute_renyi::Bool=false,
    is_show_progress::Bool=true,
)
    n_ru, n_m = size(shadows)
    k_vec_sorted = sort(unique(k_vec)) # 只算互不相同的 k，从小到大
    nK = length(k_vec_sorted)

    # 容器
    θ_est = zeros(Float64, nK)
    jackmat = is_compute_cov ? zeros(Float64, n_ru, nK) : nothing

    # --- 辅助函数：单个 k 的估计（含可选 jackknife）----------------
    function single_k(k::Int)
        @assert !(k == 1 && is_compute_renyi) "is_compute_renyi must be false when k == 1."
        # 预先枚举置换和测量的笛卡尔积
        perms = collect(permutations(1:n_ru, k))
        cprod = collect(CartesianIndices(ntuple(_ -> 1:n_m, k)))
        n_perm = length(perms)

        # 每个置换对测量求平均
        perm_avg = zeros(Float64, n_perm)
        @showprogress desc="Permutations Processing..." enabled=is_show_progress @threads for pidx in
                                                                                           eachindex(
            perms
        )
            r = perms[pidx]

            ssum = 0.0
            for m in cprod
                ssum += real(get_trace_product((shadows[r[i], m[i]] for i in 1:k)...; O))
            end

            perm_avg[pidx] = ssum / length(cprod)
        end

        # 定义平均泛函
        avgfun(x) = is_compute_renyi ? (1/(1-k))*log2(mean(x)) : mean(x)

        θ = avgfun(perm_avg)

        if !is_compute_cov
            return θ, nothing
        end

        # jackknife 分组：不含幺正 i 的置换
        jackvals = zeros(Float64, n_ru)
        @threads for i in 1:n_ru
            s = 0.0
            count = 0
            for (idx, r) in enumerate(perms)
                if i ∉ r
                    s += perm_avg[idx]
                    count += 1
                end
            end
            μ = s / count
            jackvals[i] = is_compute_renyi ? (1 / (1 - k)) * log2(μ) : μ
        end
        return θ, jackvals
    end
    # -----------------------------------------------------------------------

    # 遍历要算的矩
    for (idx, k) in enumerate(k_vec_sorted)
        θ_est[idx], jv = single_k(k)
        if is_compute_cov
            jackmat[:, idx] = jv
        end
    end

    # 需要时构造协方差
    if is_compute_cov
        Σ = zeros(Float64, nK, nK)
        for a in 1:nK, b in a:nK           # 对称
            cov =
                (n_ru-1)^2/n_ru * dot(
                    jackmat[:, a] .- mean(jackmat[:, a]),
                    jackmat[:, b] .- mean(jackmat[:, b]),
                ) / (n_ru-1)
            Σ[a, b] = Σ[b, a] = cov
        end

        θ_jack = n_ru * θ_est .- (n_ru - 1) * vec(mean(jackmat; dims=1))

        return θ_est, θ_est - θ_jack, Σ
    else
        return θ_est
    end
end

"""
modified_get_trace_moment(shadows::Vector{<:AbstractShadow}, kth_moment::Int; kwargs...)

shadow 向量版重载，reshape 成二维再调 modified_get_trace_moment。
"""
function modified_get_trace_moment(
    shadows::Vector{<:AbstractShadow}, kth_moment::Int; kwargs...
)
    return modified_get_trace_moment(
        reshape(shadows, length(shadows), 1), kth_moment; kwargs...
    )
end

# --------------------
# 求纯度
# --------------------

# 本函数只用于 O == nothing 的情形。
"""
modified_get_purity_shadow(shadows; is_compute_sem=false, is_compute_renyi=false, is_show_progress=false)

用经典 shadow 估计底层量子态的纯度 Tr(ρ^2)。

参数
- shadows::Array{<:AbstractShadow, 2}
    (n_ru, n_m) 排列的 shadow。

关键词参数
- is_compute_sem::Bool=false
    为 true 时计算并返回 jackknife 偏差和标准误差。
- is_compute_renyi::Bool=false
    为 true 时对 jackknife 值做 Rényi 换算（如适用）。
- is_show_progress::Bool=true

返回
- is_compute_sem == false：返回 purity_estimate::Float64。
- is_compute_sem == true：返回 (purity_estimate::Float64, bias::Float64, sem::Float64)。

说明
本函数靠 calculate_purity_jackvals 拿纯度估计的 jackknife 值，
需要时再算偏差和 SEM。
"""
function modified_get_purity_shadow(
    shadows::Array{<:AbstractShadow,2};
    is_compute_sem::Bool=false,
    is_compute_renyi::Bool=false,
    is_show_progress::Bool=true,
)
    n_ru, n_m = size(shadows)

    # 遍历要算的矩
    θ_est, loos = get_purity_loos_shadow(
        shadows; is_compute_renyi=is_compute_renyi, is_show_progress=is_show_progress
    )

    # 需要时构造协方差
    if is_compute_sem
        variance = (n_ru - 1)^2 / n_ru * var(loos)
        sem = sqrt(variance)
        θ_jack = n_ru * θ_est - (n_ru - 1) * mean(loos)
        return θ_est, θ_est - θ_jack, sem
    else
        return θ_est
    end
end

"""
modified_get_purity_shadow(shadows::Vector{<:AbstractShadow}; kwargs...)

shadow 向量版重载，reshape 成二维再转调主函数。
"""
function modified_get_purity_shadow(shadows::Vector{<:AbstractShadow}; kwargs...)
    return modified_get_purity_shadow(reshape(shadows, length(shadows), 1); kwargs...)
end

# --------------------
# 校准向量 G（复刻上游，加进度条开关）
# --------------------

"""
modified_get_calibration_vector(ψ0, measurement_group; is_show_progress=false)

由校准初态 ψ0（通常取 |00…0>）和 trivial 测量数据估计校准向量 G，
用于 process shadow 求值。参考文献 Vitale et al, PRXQ 2025。

复刻自 RandomMeas 的 `get_calibration_vector`，唯一区别是上游进度条
写死常开，这里加 `is_show_progress` 开关（默认关闭）。

参数
- ψ0::MPS：校准初态。
- measurement_group::MeasurementGroup：trivial 测量数据。

关键词参数
- is_show_progress::Bool=false：为 true 时显示逐 qubit 进度条。

返回
- G_e::Vector{Float64}：每个 qubit 的校准值。
"""
function modified_get_calibration_vector(
    ψ0::MPS, measurement_group::MeasurementGroup; is_show_progress=false
)
    N = length(ψ0)
    NU = measurement_group.NU
    # G_e 存每个 qubit 的 G 估计值。
    G_e = zeros(Float64, N)
    @showprogress dt=1 enabled=is_show_progress for i in 1:N
        reduced_measurement_group = reduce_to_subsystem(measurement_group, [i])

        for r in 1:NU
            # 用实际测量数据算第 i 个 qubit 的 Born 概率
            data = reduced_measurement_group.measurements[r]
            P_measured = MeasurementProbability(data).measurement_probability
            # 用理想无噪声态 ψ0 算第 i 个 qubit 的期望 Born 概率
            P_expected = MeasurementProbability(reduce_to_subsystem(ψ0, collect(i:i)), data.measurement_setting).measurement_probability
            cross_corr = (P_measured * P_expected)[]
            self_corr = (P_expected * P_expected)[]

            G_e[i] += 3 * (cross_corr - self_corr) / NU + 1 / NU
        end
    end
    return G_e
end
