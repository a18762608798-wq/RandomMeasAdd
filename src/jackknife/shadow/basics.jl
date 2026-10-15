# ---------------------
# shadow 基础组件：成对缩并值
# ---------------------

"""
get_reversal_pair_value(shadow_a, shadow_b_T, U, Ud, deltas)

算一对 shadow 的时间反演交叠值 Tr(A·U·B_T·U†)（已验证过的缩并链）。

参数
- shadow_a::Shadow：第一个设置的 shadow（未转置）。
- shadow_b_T::Shadow：第二个设置的 shadow（已在 I_1 区做部分转置）。
- U::ITensor：u_T 算符（I_1 区逐比特 σ^y）。
- Ud::ITensor：U 的 dagger（`dag(U)`，prime 层已验证）。
- deltas::Vector{ITensor}：预先配好的 δ 张量（按 id 把 plev 0 与 plev 4 腿粘起来求迹）。

返回
- val::Float64：该有序对的估计值（取实部）。

说明
缩并链：M1 = A·U(+1) → M2 = M1·Bt(+2) → M3 = M2·Ud(+3)，
最后用 deltas 缩掉首尾腿。括号内为 prime 层提升数。
"""
function get_reversal_pair_value(
    shadow_a::Shadow,
    shadow_b_T::Shadow,
    U::ITensor,
    Ud::ITensor,
    deltas::Vector{ITensor},
)
    A = shadow_a.rho
    Bt = shadow_b_T.rho
    M1 = A * prime(U, 1)
    M2 = M1 * prime(Bt, 2)
    M3 = M2 * prime(Ud, 3)
    acc = M3
    for d in deltas
        acc = acc * d
    end
    return real(acc[])
end
