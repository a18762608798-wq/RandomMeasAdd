# RandomMeasAdd

随机测量后处理工具（Jackknife、expectation 实例、shadow 与 Hamming 估计子等）。
输入为 `qmeas` 生成的 hist-v1 `npz`；张量计算经 `ITensors.jl` 直接实现，不依赖上游测量包。

## 安装

详见 [INSTALL.md](INSTALL.md).

## License

Apache-2.0，见 `LICENSE`。

## 致谢 / 上游依赖协议

本包通过 Julia `Project.toml` 引用以下第三方库（未复制其源码），各库版权归各自作
者所有：

- ITensors.jl / ITensorMPS.jl — Apache-2.0（张量后端）。
- Combinatorics、NPZ、ProgressMeter、Statistics — 各自采用 MIT / BSD 类宽松协议。
- CondaPkg.jl — MIT；本包另经 CondaPkg 引用 `qiskit`（Apache-2.0）、`quarkstudio`（MIT）。
