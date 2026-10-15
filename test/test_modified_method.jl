using RandomMeasAdd
using ITensorMPS

# 自有 shadow 手工验证脚本（不在 runtests.jl 内，手动跑）：
# hist -> Shadow -> trace 矩 / purity，与 golden 对拍。

group_path = joinpath(
    @__DIR__, "data", "aer_independent_pidx_1",
    "aer_independent_pidx_1_setting1_settings729_shots1024.npz",
)
hgroup, _, _ = import_hist_group(group_path; permuted_order=[1, 2, 3, 4]);
sub_indices = siteinds("Qubit", hgroup.N);
sub_shadows = get_dense_shadows(hgroup, sub_indices);

test_indices = [1, 2]

if 1 in test_indices
    reflect_op = create_reflect_op(sub_indices)
    @show linkdims(reflect_op)
    # NOTE: Which should be equal.
    @show modified_get_trace_moment(sub_shadows, 1; O=reflect_op, is_compute_sem=true, is_show_progress=true)
end
if 2 in test_indices
    @show modified_get_trace_moment(sub_shadows, 2; is_compute_sem=true, is_show_progress=true)
    @show modified_get_purity_shadow(sub_shadows; is_compute_sem=true, is_show_progress=true)
end
