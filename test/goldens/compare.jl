# Golden 对拍脚本（Step 0/Phase 1 验收）：逐条重跑 THUNKS，与 values.toml 比对。
#
# 运行（RandomMeasAdd 目录）：
#   julia --project=. test/goldens/compare.jl
# 约定：calib.* 条目要求逐位一致（整数数学）；估计子条目 rtol=1e-9。
# 注：shadow 估计子内层有 @threads 求和，run-to-run 有 ~1e-11 量级相对 nondeterminism
#（实测 shadow.zt.aer.mit[2] 两次运行差 5.7e-14/4.6e-3 ≈ 1.2e-11），故 rtol 取 1e-9；
# 真 bug（权重/归一化错）至少是 1/NU ~ 1e-3 量级，不会被 1e-9 漏掉。
# 输入文件先按 inputs.toml 验 sha256，输入变了直接报错（不比对）。
include(joinpath(@__DIR__, "record.jl"))
using TOML
using SHA

const GOLDIR = @__DIR__
const VALUES = TOML.parsefile(joinpath(GOLDIR, "values.toml"))
const INPUTS = TOML.parsefile(joinpath(GOLDIR, "inputs.toml"))

function sha256hex(path)
    return bytes2hex(open(SHA.sha256, path))
end

function check_entry(name)
    files, thunk = THUNKS[name]
    haskey(VALUES, name) || error("golden 缺条目: $name（先跑 record.jl）")
    for f in files
        rel = relpath(f, ROOT)
        haskey(INPUTS, rel) || error("输入清单缺文件: $rel")
        sha256hex(f) == INPUTS[rel]["sha256"] || error("输入变了: $rel（sha256 不一致，先确认数据）")
    end
    got = run_thunk(thunk)
    want = VALUES[name]["out"]
    length(got) == length(want) || error("$name 输出元数变了")
    # 只有 count 层整数数学（f0/R/xy）要求逐位一致；G 含张量缩并浮点求和，
    # expand 行顺序变化可带来 ulp 漂移，走容差。
    exact = occursin(r"f0R|_xy", name)
    for (i, (g, w)) in enumerate(zip(got, want))
        gv, wv = vec(g), vec(Float64.(w))
        length(gv) == length(wv) || error("$name[$i] 长度变了")
        ok = exact ? gv == wv :
            all(isapprox.(gv, wv; rtol=1e-9, atol=1e-12))
        ok || error("$name[$i] 失配:\n got  = $gv\n want = $wv")
    end
    println("PASS ", name)
end

function main()
    names = sort(collect(keys(THUNKS)))
    for name in names
        check_entry(name)
    end
    println("ALL PASS (", length(names), " entries)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
