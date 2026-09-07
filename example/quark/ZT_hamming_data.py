"""生成 z_t hamming 示例数据（quark 真机，conjugate-pair 方案，pidx=-1）。

参考 test/get_data/quark/；需 QUARK_TOKEN 环境变量与芯片比特映射；
输出到本目录 data/，共两组小 SettingRun。
"""

import asyncio
import os
from pathlib import Path

from qmeas.models.xxz import get_initial_state
from qmeas.random import (
    ConjugatePair,
    QuarkOptions,
    RandomMeasConfig,
    SettingRun,
    run_random,
)

SEED = 521
N_QUBITS = 8
DATA_DIR = Path(__file__).resolve().parent / "data"

CHIP = "Dongling"
TARGET_QUBITS: list = []  # 按实际任务填写

# 时间反演配对数据（Z_T hamming 用）：列顺序与 independent 一致 [(2,), (5,), (3,), (4,)]，
# 奇位 {2, 3} = I_1、偶位 {5, 4} = I_2，故 i1_groups=(0, 2)。
# hamming 不支持误差缓解，故 mitigation=False。
PAIR_INDICES = [(2,), (5,), (3,), (4,)]


def main():
    qc = get_initial_state(N_QUBITS, pidx=-1)
    name = "quark_pair_pidx_-1"
    config = RandomMeasConfig(
        qc=qc,
        setting_runs=[
            SettingRun(num_settings=3**3, num_shots=1024),
            SettingRun(num_settings=3**4, num_shots=1024),
        ],
        meas_indices=PAIR_INDICES,
        ensemble="haar",
        conjugate_pair=ConjugatePair(i1_groups=(0, 2)),
        runner_opts=QuarkOptions(
            chip=CHIP,
            token=os.environ.get("QUARK_TOKEN"),
            target_qubits=TARGET_QUBITS,
            mitigation=False,  # hamming 路线不支持误差缓解
        ),
        seed=SEED,
        output_dir=DATA_DIR / name,
        name=name,
    )
    asyncio.run(run_random(config=config))
    print(f"done: {name}")


if __name__ == "__main__":
    main()
