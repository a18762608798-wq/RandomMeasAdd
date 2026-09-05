"""生成 z_t shadow 示例数据（quark 真机，independent 方案，pidx=-1）。

参考 test/get_data/quark/；需 QUARK_TOKEN 环境变量与芯片比特映射；
输出到本目录 data/，共两组小 SettingRun。
"""

import asyncio
import os
from pathlib import Path

from qmeas.models.xxz import get_initial_state
from qmeas.random import QuarkOptions, RandomMeasConfig, SettingRun, run_random

SEED = 521
N_QUBITS = 8
DATA_DIR = Path(__file__).resolve().parent / "data"

CHIP = "Dongling"
TARGET_QUBITS: list = []  # 按实际任务填写


def main():
    qc = get_initial_state(N_QUBITS, pidx=-1)
    name = "quark_independent_pidx_-1"
    config = RandomMeasConfig(
        qc=qc,
        setting_runs=[
            SettingRun(num_settings=3**3, num_shots=1024),
            SettingRun(num_settings=3**4, num_shots=1024),
        ],
        meas_indices=[(2,), (5,), (3,), (4,)],
        ensemble="haar",
        runner_opts=QuarkOptions(
            chip=CHIP,
            token=os.environ.get("QUARK_TOKEN"),
            target_qubits=TARGET_QUBITS,
            mitigation=True,  # shadow 路线开误差缓解，npz 内含 trivial 数据
        ),
        seed=SEED,
        output_dir=DATA_DIR / name,
        name=name,
    )
    asyncio.run(run_random(config=config))
    print(f"done: {name}")


if __name__ == "__main__":
    main()
