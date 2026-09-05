# --------------------
# 局域 swap
# --------------------

"""
create_reflect_op(siteindices)

构造反射算符：对每个 i，把第 i 个 site 和第 (N - i + 1) 个 site 做 swap，
即沿一维链中心做空间反射。

参数
- siteindices::Vector{Index{Int64}}：全系统的 site index 向量。

返回
- reflect_op::MPO：反射的矩阵乘积算符表示。

说明
- site 数（bitsnum）必须为偶数。
- 做法是在单位 MPO 上作用 SWAP 门。
"""
function create_reflect_op(siteindices::Vector{Index{Int64}})
    bitsnum = length(siteindices)
    @assert iseven(bitsnum) "bitsnum must be even" # 约束
    pairsnum = bitsnum ÷ 2
    identity_op = MPO(siteindices, "Id")
    gates = [op("SWAP", siteindices[i], siteindices[bitsnum-i+1]) for i in 1:pairsnum]
    reflect_op = apply(gates, identity_op)
    return reflect_op
end

# 相邻 swap
"""
create_adjacent_swap_op(siteindices)

构造相邻配对 swap 算符：(1, 2)、(3, 4)、…、(N-1, N)，
即一维链上的近邻 swap。

参数
- siteindices::Vector{Index{Int64}}：全系统的 site index 向量。

返回
- adjacent_swap_op::MPO：实现相邻 swap 的矩阵乘积算符。

说明
- site 数（bitsnum）必须为偶数。
- 做法是在单位 MPO 上作用 SWAP 门。
"""
function create_adjacent_swap_op(siteindices::Vector{Index{Int64}})
    bitsnum = length(siteindices)
    @assert iseven(bitsnum) "bitsnum must be even" # 约束
    pairsnum = bitsnum ÷ 2
    identity_op = MPO(siteindices, "Id")
    gates = [op("SWAP", siteindices[2i-1], siteindices[2i]) for i in 1:pairsnum]
    adjacent_swap_op = apply(gates, identity_op)
    return adjacent_swap_op
end

# ---------------------
# 时间反演算符
# ---------------------

"""
create_unitary_part_reversal_op(part1, site_indices; op_type="MPO")

在 `part1` 指定的 site 上作用 Pauli-Y 门，构造时间反演算符的幺正部分，
即作用在所选子系统上的 ⊗_{i∈part1} Y_i。

参数
- part1::Vector{Int64}：要作用 Pauli-Y 门的 site 编号。
- site_indices::Vector{Index{Int64}}：全系统的 site index 向量。

关键词参数
- op_type::String="MPO"：输出类型，支持 "MPO"（返回 MPO）
  和 "ITensor"（把 MPO 缩成单个 ITensor）。

返回
- op_type == "MPO"：返回 unitary_part_reversal_op::MPO。
- op_type == "ITensor"：返回缩并后的 ITensor。

说明
- Y 门作用在全系统的单位 MPO 上。
- op_type 不是 "MPO"/"ITensor" 时报错。
"""
function create_unitary_part_reversal_op(
    part1::Vector{Int64}, site_indices::Vector{Index{Int64}}; op_type="MPO"
)
    qubits_num = length(site_indices)
    identity_op = MPO(site_indices, "Id")
    gates = [op("Y", site_indices[i]) for i in part1]
    unitary_part_reversal_op = apply(gates, identity_op)
    if op_type == "MPO"
        return unitary_part_reversal_op
    elseif op_type == "ITensor"
        return contract(unitary_part_reversal_op)
    else
        error("wrong output type: $op_type")
    end
end
