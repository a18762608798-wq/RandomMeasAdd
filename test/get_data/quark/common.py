import os
from pathlib import Path

from qmeas.random import QuarkOptions, SettingRun

# 生成脚本共享配置（真机）
SEED = 521
N_QUBITS = 8
DATA_DIR = Path(__file__).resolve().parent.parent.parent / "data"

# 真机按需调小设置数与 shot 数；芯片与比特映射按实际任务填写
setting_runs = [
    SettingRun(num_settings=3**4, num_shots=1024),
]

QUARK_TOKEN = os.environ.get("QUARK_TOKEN")
CHIP = "Baihua"
TARGET_QUBITS: list = []


def make_quark_opts(mitigation=False):
    """构造 QuarkOptions；token 从 QUARK_TOKEN 环境变量取。"""
    return QuarkOptions(
        chip=CHIP,
        token=QUARK_TOKEN,
        target_qubits=TARGET_QUBITS,
        mitigation=mitigation,
    )
