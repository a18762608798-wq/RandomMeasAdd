# 安装

## Julia 包

从 GitHub 安装：

```julia
julia --project=/你的项目路径 # 进入目标目录后.
using Pkg
Pkg.add(url="git@github.com:a18762608798-wq/RandomMeasAdd.git")
```

更新（`Pkg.add` 是惰性的：装过就跳过，不去云端看；要拿新版显式升级）：

```julia
using Pkg
Pkg.update("RandomMeasAdd")
```

## Python 依赖：qmeas

数据生成脚本（[`test/get_data/`](test/get_data/)、[`example/*/gen_data.py`](example/)）需要 Python 库 `qmeas`，
由 `CondaPkg.toml` 自动从 GitHub 安装（另含 `qiskit`，aer 模拟用；quark 真机另需 `quarkstudio`），
Julia 首次运行时装好，无需手动 pip。

### 用 CondaPkg 的 Python 跑脚本

想用 Julia 项目自带的 Conda 环境里的 Python（已含上述依赖）跑数据生成脚本时，在 Pkg REPL 里用 `conda run`：

```julia
julia --project=.
julia> using CondaPkg
pkg> conda run python example/z_r/gen_data.py
```

不用 REPL 的话，直接调该环境的 `python` 也行（具体路径以后端为准，可用 `CondaPkg.which("python")` 查询；
本机 pixi 后端一般是 `.CondaPkg/.pixi/envs/default/bin/python`）。
