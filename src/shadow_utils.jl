# --------------------
# 把 shadow 转成 mpo
# --------------------

"""
把单个 FactorizedShadow 转成 MPO 表示。

参数
- shadow::FactorizedShadow：单个 factorized shadow 实例。

返回
- 由 shadow.shadow_data 构造的 MPO。
"""
function get_factorized_shadow_mpo(shadow::FactorizedShadow)
    return MPO(shadow.shadow_data)
end

"""
把 FactorizedShadow 矩阵转成 MPO 矩阵。

参数
- shadows::Matrix{FactorizedShadow}：(settings, shots) 的二维数组。

返回
- 同尺寸的 Matrix{MPO}，每个元素是对应 FactorizedShadow 的 MPO 形式。
"""
function get_factorized_shadow_mpo(shadows::Matrix{FactorizedShadow})
    settings_num, shots = size(shadows)
    mpo_shadows = Matrix{MPO}(undef, settings_num, shots)
    for settings in 1:settings_num
        for shot in 1:shots
            shadow = shadows[settings, shot]
            mpo_shadow = get_factorized_shadow_mpo(shadow)
            mpo_shadows[settings, shot] = mpo_shadow
        end
    end
    return mpo_shadows
end

"""
把 FactorizedShadow 向量转成单列 MPO 矩阵。

向量版简便重载：把向量 reshape 成 (N,1) 矩阵再调 Matrix{FactorizedShadow} 方法。
"""
function get_factorized_shadow_mpo(shadows::Vector{FactorizedShadow})
    return get_factorized_shadow_mpo(reshape(shadows, length(shadows), 1))
end

