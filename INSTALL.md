# 安装

## Julia 包

从 GitHub 安装：

```julia
# cd 进入目标目录后.
using Pkg
Pkg.activate(".")
Pkg.add(url="https://github.com/a18762608798-wq/RandomMeasAdd.git")
```

更新（`Pkg.add` 是惰性的：装过就跳过，不去云端看；要拿新版显式升级）：

```julia
using Pkg
Pkg.update("RandomMeasAdd")
```

## Python 依赖：qmeas

### python venv

参考qmeas安装流程(自动安装相关依赖):

[INSTALL](https://github.com/a18762608798-wq/qmeas/blob/master/INSTALL.md)

### CondaPkg

```julia
julia> using CondaPkg
julia> CondaPkg.resolve()
```

想用 Julia 项目自带的 Conda 环境里的 Python（已含上述依赖）跑数据生成脚本时，在
Pkg REPL 里用 `conda run`：

```julia
# cd 进入目标目录后.
using Pkg
Pkg.activate(".")
using CondaPkg
pkg> conda run python example/aer/shadow_data.py
```
