from pathlib import Path

from qmeas.random import SettingRun

# 生成脚本共享配置
SEED = 521
N_QUBITS = 8
DATA_DIR = Path(__file__).resolve().parent.parent.parent / "data"

setting_runs = [
    SettingRun(num_settings=3**5, num_shots=1024),
    SettingRun(num_settings=3**6, num_shots=1024),
]
