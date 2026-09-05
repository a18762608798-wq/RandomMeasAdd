using RandomMeasAdd

# get sub shadows
group_path = joinpath(
    @__DIR__, "data", "aer_independent_pidx_1",
    "aer_independent_pidx_1_setting1_settings729_shots1024.npz",
)
sub_group, sub_indices, _ = import_random_group(group_path; permuted_order=[1, 2, 3, 4]);
sub_shadows = get_dense_shadows(sub_group);

test_indices = [1, 2]

if 1 in test_indices
    reflect_op = create_reflect_op(sub_indices)
    @show linkdims(reflect_op)
    # NOTE: Which should be equal.
    @show modified_get_expect_shadow(reflect_op, sub_shadows, is_compute_sem=true, is_show_progress=true)
    @show modified_get_trace_moment(sub_shadows, 1; O=reflect_op, is_compute_sem=true, is_show_progress=true)
end
if 2 in test_indices
    @show modified_get_trace_moment(sub_shadows, 2; is_compute_sem=true, is_show_progress=true)
    @show modified_get_purity_shadow(sub_shadows; is_compute_sem=true, is_show_progress=true)
end
