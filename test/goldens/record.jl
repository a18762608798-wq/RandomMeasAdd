# Golden 矩阵（Step 0/Phase 1 共享）：估计子 × 数据集。
#
# - 直接执行本文件：跑全部条目，写 test/goldens/values.toml（录制）。
# - 被 compare.jl include：只注册 THUNKS，不执行（对拍）。
# 每个 thunk 自包含（自己 import），以便对拍逐条重跑。
# 约定：全部数据集都是 4 比特，permuted_order 统一 [1,2,3,4]。
# mit 开关只在 npz 含对应 trivial 键时开启（pair_mit 的 Z_T 用 plain1+trivial2，
# 与 get_z_t_hamming 文件方法一致）。
using RandomMeasAdd
using TOML

const ROOT = joinpath(@__DIR__, "..", "..")
const DATADIR = joinpath(ROOT, "test", "data")
const PERM = [1, 2, 3, 4]

f(d, n) = joinpath(DATADIR, d, n)

# setting0 文件（243 settings aer / 81 settings quark），速度与覆盖折中。
const IND_AER = f("aer_independent_pidx_1", "aer_independent_pidx_1_setting0_settings243_shots1024.npz")
const SH_AER = f("aer_shared_pidx_1", "aer_shared_pidx_1_setting0_settings243_shots1024.npz")
const P1_AER = f("aer_pair_mit_pidx_1", "aer_pair_mit_pidx_1_exp1_setting0_settings243_shots1024.npz")
const P2_AER = f("aer_pair_mit_pidx_1", "aer_pair_mit_pidx_1_exp2_setting0_settings243_shots1024.npz")
const M1_AER = f("aer_pair_mit_pidx_-1", "aer_pair_mit_pidx_-1_exp1_setting0_settings243_shots1024.npz")
const M2_AER = f("aer_pair_mit_pidx_-1", "aer_pair_mit_pidx_-1_exp2_setting0_settings243_shots1024.npz")
const NP1_AER = f("aer_pair_pidx_1", "aer_pair_pidx_1_exp1_setting0_settings243_shots1024.npz")
const NP2_AER = f("aer_pair_pidx_1", "aer_pair_pidx_1_exp2_setting0_settings243_shots1024.npz")
const IND_Q = f("quark_dongling_pidx_1_independent", "quark_dongling_pidx_1_independent_setting0_settings81_shots1024.npz")
const SH_Q = f("quark_dongling_pidx_1_shared", "quark_dongling_pidx_1_shared_setting0_settings81_shots1024.npz")
const QP1_Q = f("quark_dongling_pidx_1_pair", "quark_dongling_pidx_1_pair_exp1_setting0_settings81_shots1024.npz")
const QP2_Q = f("quark_dongling_pidx_1_pair", "quark_dongling_pidx_1_pair_exp2_setting0_settings81_shots1024.npz")
const NQP1 = f("quark_pair_pidx_1", "quark_pair_pidx_1_exp1_setting0_settings81_shots1024.npz")
const NQP2 = f("quark_pair_pidx_1", "quark_pair_pidx_1_exp2_setting0_settings81_shots1024.npz")

const THUNKS = Dict{String,Pair{Vector{String},Function}}()

"""注册一条 golden 条目（只注册不执行）。注意 do 块固定为第一个参数。"""
function rec!(thunk, name, files)
    THUNKS[name] = files => thunk
end

"""执行 thunk，返回值元组拆成向量列表（标量包成单元素向量）。"""
function run_thunk(thunk)
    r = thunk()
    vals = r isa Tuple ? collect(r) : [r]
    return [v isa AbstractArray ? vec(Float64.(v)) : Float64[v] for v in vals]
end

# ── hamming：purity（independent）──
rec!("ham.purity.aer.std", [IND_AER]) do
    g, _, _ = import_random_group(IND_AER; permuted_order=PERM)
    get_purity_hamming(g; is_compute_sem=true)
end
rec!("ham.purity.aer.mit", [IND_AER]) do
    g, _, _, triv = import_random_group(IND_AER; permuted_order=PERM, is_mitigation=true, is_return_trivial=true)
    get_purity_hamming(g; trivial=triv, is_mitigation=true, is_compute_sem=true)
end
rec!("ham.purity.quark.std", [IND_Q]) do
    get_purity_hamming(IND_Q; permuted_order=PERM, is_compute_sem=true)
end
rec!("ham.purity.quark.mit", [IND_Q]) do
    get_purity_hamming(IND_Q; permuted_order=PERM, is_mitigation=true, is_compute_sem=true)
end

# ── hamming：reflect / z_r（shared）──
rec!("ham.reflect.aer.std", [SH_AER]) do
    s, _, _ = import_random_group(SH_AER; permuted_order=PERM)
    get_reflect_hamming(s; is_compute_sem=false)
end
rec!("ham.reflect.aer.mit", [SH_AER]) do
    s, _, _, striv = import_random_group(SH_AER; permuted_order=PERM, is_mitigation=true, is_return_trivial=true)
    get_reflect_hamming(s; trivial=striv, is_mitigation=true)
end
rec!("ham.zr.aer.std", [SH_AER]) do
    s, _, _ = import_random_group(SH_AER; permuted_order=PERM)
    get_z_r_hamming(s; is_compute_sem=true)
end
rec!("ham.zr.aer.mit", [SH_AER]) do
    s, _, _, striv = import_random_group(SH_AER; permuted_order=PERM, is_mitigation=true, is_return_trivial=true)
    get_z_r_hamming(s; trivial=striv, is_mitigation=true, is_compute_sem=true)
end
rec!("ham.reflect.quark.mit", [SH_Q]) do
    get_reflect_hamming(SH_Q; permuted_order=PERM, is_mitigation=true)
end
rec!("ham.zr.quark.mit", [SH_Q]) do
    get_z_r_hamming(SH_Q; permuted_order=PERM, is_mitigation=true, is_compute_sem=true)
end

# ── hamming：reversal / z_t（pair；mit 用 plain1+trivial2）──
rec!("ham.reversal.aer.std", [P1_AER, P2_AER]) do
    g1, g2, _, _, _ = import_random_pair(P1_AER, P2_AER; permuted_order=PERM)
    get_reversal_hamming(g1, g2)
end
rec!("ham.reversal.aer.mit", [P1_AER, P2_AER]) do
    g1, g2, _, _, _, _, t2, p1 = import_random_pair(
        P1_AER, P2_AER; permuted_order=PERM,
        is_mitigation=true, is_return_trivial=true, is_return_trivial_plain=true,
    )
    get_reversal_hamming(g1, g2; trivial1=p1, trivial2=t2, is_mitigation=true)
end
rec!("ham.zt.aer.std", [P1_AER, P2_AER]) do
    g1, g2, _, _, _ = import_random_pair(P1_AER, P2_AER; permuted_order=PERM)
    get_z_t_hamming(g1, g2; is_compute_sem=true)
end
rec!("ham.zt.aer.mit", [P1_AER, P2_AER]) do
    g1, g2, _, _, _, _, t2, p1 = import_random_pair(
        P1_AER, P2_AER; permuted_order=PERM,
        is_mitigation=true, is_return_trivial=true, is_return_trivial_plain=true,
    )
    get_z_t_hamming(g1, g2; trivial1=p1, trivial2=t2, is_mitigation=true, is_compute_sem=true)
end
# 拓扑相（pidx=-1，真值为负）
rec!("ham.zt.aer.neg.mit", [M1_AER, M2_AER]) do
    m1, m2, _, _, _, _, mt2, mp1 = import_random_pair(
        M1_AER, M2_AER; permuted_order=PERM,
        is_mitigation=true, is_return_trivial=true, is_return_trivial_plain=true,
    )
    get_z_t_hamming(m1, m2; trivial1=mp1, trivial2=mt2, is_mitigation=true, is_compute_sem=true)
end
# 无 trivial pair：std only
rec!("ham.reversal.aer.notriv", [NP1_AER, NP2_AER]) do
    n1, n2, _, _, _ = import_random_pair(NP1_AER, NP2_AER; permuted_order=PERM)
    get_reversal_hamming(n1, n2)
end
rec!("ham.zt.aer.notriv", [NP1_AER, NP2_AER]) do
    n1, n2, _, _, _ = import_random_pair(NP1_AER, NP2_AER; permuted_order=PERM)
    get_z_t_hamming(n1, n2)
end
# 真机 pair（有 trivial 无 plain：reversal 可缓解，z_t 只能 std）
rec!("ham.reversal.quark.mit", [QP1_Q, QP2_Q]) do
    q1, q2, _, _, _, qt1, qt2 = import_random_pair(
        QP1_Q, QP2_Q; permuted_order=PERM, is_mitigation=true, is_return_trivial=true,
    )
    get_reversal_hamming(q1, q2; trivial1=qt1, trivial2=qt2, is_mitigation=true)
end
rec!("ham.zt.quark.std", [QP1_Q, QP2_Q]) do
    q1, q2, _, _, _ = import_random_pair(QP1_Q, QP2_Q; permuted_order=PERM)
    get_z_t_hamming(q1, q2)
end
rec!("ham.zt.quark.notriv", [NQP1, NQP2]) do
    w1, w2, _, _, _ = import_random_pair(NQP1, NQP2; permuted_order=PERM)
    get_z_t_hamming(w1, w2)
end

# ── shadow（文件版；G 开关各一遍）──
rec!("shadow.purity.aer.std", [IND_AER]) do
    get_purity_shadow(IND_AER; permuted_order=PERM, is_compute_sem=true)
end
rec!("shadow.purity.aer.mit", [IND_AER]) do
    get_purity_shadow(IND_AER; permuted_order=PERM, is_mitigation=true, is_compute_sem=true)
end
rec!("shadow.purity.quark.mit", [IND_Q]) do
    get_purity_shadow(IND_Q; permuted_order=PERM, is_mitigation=true, is_compute_sem=true)
end
rec!("shadow.reflect.aer.mit", [SH_AER]) do
    get_reflect_shadow(SH_AER; permuted_order=PERM, is_mitigation=true, is_compute_sem=true)
end
rec!("shadow.zr.aer.mit", [IND_AER]) do
    get_z_r_shadow(IND_AER; permuted_order=PERM, is_mitigation=true, is_compute_sem=true)
end
rec!("shadow.reversal.aer.mit", [IND_AER]) do
    get_reversal_shadow(IND_AER; permuted_order=PERM, is_mitigation=true, is_compute_sem=true)
end
rec!("shadow.zt.aer.mit", [IND_AER]) do
    get_z_t_shadow(IND_AER; permuted_order=PERM, is_mitigation=true, is_compute_sem=true)
end

# ── 标定量（count 层整数数学，Phase 1 要求逐位一致）──
rec!("calib.purity_xy.aer", [IND_AER]) do
    _, _, _, triv = import_random_group(IND_AER; permuted_order=PERM, is_mitigation=true, is_return_trivial=true)
    RandomMeasAdd.get_purity_xy(triv, 1:4)
end
rec!("calib.purity_f0R.aer", [IND_AER]) do
    _, _, _, triv = import_random_group(IND_AER; permuted_order=PERM, is_mitigation=true, is_return_trivial=true)
    (RandomMeasAdd.get_f0(triv, 1:4), RandomMeasAdd.get_R_purity(triv, 1:4))
end
rec!("calib.mirror_xy.aer", [SH_AER]) do
    _, _, _, striv = import_random_group(SH_AER; permuted_order=PERM, is_mitigation=true, is_return_trivial=true)
    RandomMeasAdd.get_mirror_xy(striv, [1, 3], [2, 4])
end
rec!("calib.cross_xy.aer", [P1_AER, P2_AER]) do
    _, _, _, _, _, _, t2, p1 = import_random_pair(
        P1_AER, P2_AER; permuted_order=PERM,
        is_mitigation=true, is_return_trivial=true, is_return_trivial_plain=true,
    )
    RandomMeasAdd.get_cross_xy(p1, t2, 1:4)
end
rec!("calib.G.aer", [IND_AER]) do
    _, _, G = import_random_group(IND_AER; permuted_order=PERM, is_mitigation=true)
    (G,)
end

# 直接执行本文件 = 录制 values.toml；被 compare.jl include 时只注册不跑。
if abspath(PROGRAM_FILE) == @__FILE__
    entries = Dict{String,Any}(
        "_meta" => Dict{String,Any}(
            "commit" => strip(read(`git -C $ROOT rev-parse HEAD`, String)),
            "perm" => PERM,
        ),
    )
    for name in sort(collect(keys(THUNKS)))
        files, thunk = THUNKS[name]
        out = run_thunk(thunk)
        entries[name] = Dict{String,Any}("files" => files, "out" => out)
        println(name, " => ", out)
    end
    open(joinpath(@__DIR__, "values.toml"), "w") do io
        TOML.print(io, entries)
    end
    println("wrote ", length(entries) - 1, " entries")
end
