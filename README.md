# RandomMeasAdd

基于 `RandomMeas.jl` 的扩展分析工具（Jackknife、expectation 实例、shadow 工具等）。

## 安装

```julia
using Pkg
Pkg.add(url="git@github.com:a18762608798-wq/RandomMeasAdd.git")
```

本地可编辑安装：

```julia
using Pkg
Pkg.develop(path=expanduser("~/RandomMeasAdd"))
```

## 测试

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

测试数据生成脚本在 `test/get_data/`（Python，需 `qiskit` / `quarkstudio`）。

## License

Apache-2.0，见 `LICENSE`。

## 致谢 / 上游依赖协议

本包通过 Julia `Project.toml` 引用以下第三方库（未复制其源码），各库版权归各自作者所有：

- RandomMeas.jl — Apache-2.0，B. Vermersch、A. Elben 及其贡献者（`RandomMeas = "0.3.1"`）。
- Combinatorics、NPZ、ProgressMeter、Statistics — 各自采用 MIT / BSD 类宽松协议。
- CondaPkg.jl — MIT；本包另经 CondaPkg 引用 `qiskit`（Apache-2.0）、`quarkstudio`（MIT）。
