import asyncio
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from qmeas.random import AerOptions, RandomMeasConfig, run_random

from qmeas.models.xxz import get_initial_state

from common import DATA_DIR, N_QUBITS, SEED, setting_runs

# 两种测量方案：independent 每个比特独立幺正，shared 组内共享旋转参数。
# NOTE: independent 没有顺序要求, 但是建议考虑zr等结构奇偶绑定.
schemes = {
    "independent": [(2,), (5,), (3,), (4,)],
    "shared": [(2, 5), (3, 4)],
}
aer_opts = AerOptions(
    method="matrix_product_state", device="CPU", precision="double", mitigation=True
)


def gen_schemes():
    """生成 independent/shared 两套方案数据（shadow 与 z_r 用）。"""
    for scheme_name, meas_indices in schemes.items():
        for pidx in (1, 0, -1):
            qc = get_initial_state(N_QUBITS, pidx=pidx)
            name = f"aer_{scheme_name}_pidx_{pidx}"
            meas_config = RandomMeasConfig(
                qc=qc,
                setting_runs=setting_runs,
                meas_indices=meas_indices,
                ensemble="haar",
                runner_opts=aer_opts,
                seed=SEED,
                output_dir=DATA_DIR / name,
                name=name,
            )
            asyncio.run(run_random(config=meas_config))
            print(f"done: {name}")


if __name__ == "__main__":
    gen_schemes()
