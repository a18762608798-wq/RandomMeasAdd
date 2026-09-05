using RandomMeasAdd

# settings
N = 6
site_indices = siteinds("Qubit", N)

test_indices = [1, 2, 3, 4]

if 1 in test_indices
    # reflect op
    reflect_op = create_reflect_op(site_indices)
    @show linkdims(reflect_op)
end
if 2 in test_indices
    # swap op but adjacent
    adjacent_swap_op = create_adjacent_swap_op(site_indices)
    @show linkdims(adjacent_swap_op)
end
if 3 in test_indices
    part1 = [1, 2, 3, 4]
    unitary_part_reversal_op = create_unitary_part_reversal_op(
        part1, site_indices; op_type="ITensor"
    )
    unitary_part_reversal_op = create_unitary_part_reversal_op(
        part1, site_indices; op_type="MPO"
    )
    @show linkdims(unitary_part_reversal_op)
end
