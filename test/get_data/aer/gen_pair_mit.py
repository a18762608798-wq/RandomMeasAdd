"""生成带 trivial 标定的 pair 数据（hamming 缓解测试用，真机前先在模拟机验证）。

conjugate-pair 方案 + mitigation=True，pidx=1，243 settings x 1024 shots。
输出到 test/data/aer_pair_mit_pidx_1（gitignored，本地测试用）。
"""

import asyncio
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from qmeas.models.xxz import get_initial_state
from qmeas.random import (
    AerOptions,
    ConjugatePair,
    RandomMeasConfig,
    SettingRun,
    run_random,
)

from common import DATA_DIR, N_QUBITS, SEED

# 与 gen_pairs.py 相同的配对几何，区别只在 mitigation=True（两份实验各带 trivial）。
pair_indices = [(2,), (5,), (3,), (4,)]
pair_opts = AerOptions(
    method="matrix_product_state", device="CPU", precision="double", mitigation=True
)


def gen_pair_mit():
    """生成带标定的时间反演配对数据（Z_T hamming 缓解用）。"""
    qc = get_initial_state(N_QUBITS, pidx=1)
    name = "aer_pair_mit_pidx_1"
    meas_config = RandomMeasConfig(
        qc=qc,
        setting_runs=[SettingRun(num_settings=243, num_shots=1024)],
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
    gen_pair_mit()
