"""生成 z_t shadow 示例数据（aer 模拟，independent 方案，pidx=-1）。

参考 test/get_data/aer/；输出到本目录 data/，共两组小 SettingRun。
运行：python example/z_t_shadow/gen_data.py（需 qmeas 环境）
"""

import asyncio
from pathlib import Path

from qmeas.models.xxz import get_initial_state
from qmeas.random import AerOptions, RandomMeasConfig, SettingRun, run_random

SEED = 521
N_QUBITS = 8
DATA_DIR = Path(__file__).resolve().parent / "data"

# independent 方案各比特独立随机旋转（shadow z_t 用）；
# shadow 路线开误差缓解，npz 内含 trivial 数据。
aer_opts = AerOptions(
    method="matrix_product_state", device="CPU", precision="double", mitigation=True
)


def main():
    qc = get_initial_state(N_QUBITS, pidx=-1)
    name = "aer_independent_pidx_-1"
    config = RandomMeasConfig(
        qc=qc,
        setting_runs=[
            SettingRun(num_settings=3**5, num_shots=1024),
            SettingRun(num_settings=3**6, num_shots=1024),
        ],
        meas_indices=[(2,), (5,), (3,), (4,)],
        ensemble="haar",
        runner_opts=aer_opts,
        seed=SEED,
        output_dir=DATA_DIR / name,
        name=name,
    )
    asyncio.run(run_random(config=config))
    print(f"done: {name}")


if __name__ == "__main__":
    main()
