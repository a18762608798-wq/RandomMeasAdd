import asyncio
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from qmeas.random import (
    AerOptions,
    ConjugatePair,
    RandomMeasConfig,
    run_random,
)

from qmeas.models.xxz import get_initial_state

from common import DATA_DIR, N_QUBITS, SEED, setting_runs

# 时间反演配对数据（Z_T hamming 用）：列顺序与 independent 一致 [(2,), (5,), (3,), (4,)]，
# 奇位 {2, 3} = I_1、偶位 {5, 4} = I_2，故 i1_groups=(0, 2)。
# hamming 暂不支持误差缓解，此处 mitigation=False。
pair_indices = [(2,), (5,), (3,), (4,)]
pair_opts = AerOptions(
    method="matrix_product_state", device="CPU", precision="double", mitigation=False
)


def gen_pairs():
    """生成时间反演配对数据（Z_T hamming 用）。"""
    for pidx in (1, 0, -1):
        qc = get_initial_state(N_QUBITS, pidx=pidx)
        name = f"aer_pair_pidx_{pidx}"
        meas_config = RandomMeasConfig(
            qc=qc,
            setting_runs=setting_runs,
            meas_indices=pair_indices,
            ensemble="haar",
            conjugate_pair=ConjugatePair(i1_groups=(0, 2)),
            runner_opts=pair_opts,
            seed=SEED,
            output_dir=DATA_DIR / name,
            name=name,
        )
        asyncio.run(run_random(config=meas_config))
        print(f"done: {name}")


if __name__ == "__main__":
    gen_pairs()
