"""旧展开格式 npz -> hist-v1 的无损转换（一次性迁移工具，可重入）。

旧：measurement_results (NU, NM, N) uint8 0/1 展开矩阵（+ trivial/plain 同理）。
新：per-setting 去重 key（整数，bit j = 第 j 列，即 meas_indices 展平第 j 比特）
    + counts + offsets，settings/theta/phi/标量原样保留，加 format_version=1。
已是 hist-v1 的文件跳过。文件名、summary json 不动。
"""
import sys
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parent.parent / "data"


def to_hist(mat):
    """(NM, N) 0/1 矩阵 -> (keys升序, counts)。"""
    keys = mat.dot(1 << np.arange(mat.shape[1], dtype=np.uint64))
    uk, inv, cnt = np.unique(keys, return_inverse=True, return_counts=True)
    order = np.argsort(uk)
    return uk[order].astype(np.uint64), cnt[order].astype(np.int64)


def convert(path):
    d = np.load(path)
    if "format_version" in d.files and int(d["format_version"]) == 1:
        return False
    assert "measurement_results" in d.files, f"{path} 既非旧格式也非 hist-v1"
    out = {}
    n_settings = None
    for prefix in ["measurement", "trivial_measurement", "trivial_plain_measurement"]:
        mat_key = prefix + "_results"
        if mat_key not in d.files:
            continue
        mat = d[mat_key]
        assert mat.ndim == 3 and mat.dtype == np.uint8, (path, mat_key, mat.shape)
        nu, nm, _ = mat.shape
        n_settings = nu if n_settings is None else n_settings
        assert nu == n_settings
        all_k, all_c, offs = [], [], [0]
        for s in range(nu):
            k, c = to_hist(np.asarray(mat[s]))
            assert c.sum() == nm, (path, mat_key, s)
            all_k.append(k)
            all_c.append(c)
            offs.append(offs[-1] + len(k))
        stem = prefix.replace("measurement", "measurement_hist").replace("trivial_measurement_hist", "trivial_hist").replace(
            "trivial_plain_measurement_hist", "trivial_plain_hist"
        )
        # measurement -> measurement_hist_*；trivial_measurement -> trivial_hist_*
        out[stem + "_keys"] = np.concatenate(all_k)
        out[stem + "_counts"] = np.concatenate(all_c)
        out[stem + "_offsets"] = np.asarray(offs, dtype=np.int64)
    for k in d.files:
        if k in ("measurement_results", "trivial_measurement_results",
                 "trivial_plain_measurement_results"):
            continue
        out[k] = d[k]
    out["format_version"] = np.int64(1)
    np.savez(path, **out)
    return True


def main():
    files = sorted(DATA.glob("*/*.npz"))
    assert files, f"{DATA} 下没有 npz"
    n = sum(convert(p) for p in files)
    print(f"{len(files)} files, converted {n}")


if __name__ == "__main__":
    sys.exit(main())
