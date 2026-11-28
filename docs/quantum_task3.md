# TASK3：直接在 main.c 中运行 Quantum 对比仿真

2026-11-28。按本次要求重新实现，流程参考 [cpicic_decoder/main.c](https://github.com/yinanli-bit/cpicic_decoder/blob/439fd79012971bc46f8aa4cdf711f5e81f94b6cb/main.c)，固定 commit 为 `439fd79012971bc46f8aa4cdf711f5e81f94b6cb`。

## 1. 文件与恢复范围

删除上一版的：

- `include/bplcgcd/quantum_simulation.h`
- `src/apps/quantum_sim_main.c`
- `src/apps/quantum_sim_options.c`
- `src/apps/quantum_sim_options.h`
- `src/simulation/quantum_simulation.c`
- `src/simulation/quantum_simulation_output.c`
- `tests/quantum_simulation_test.c`
- `tests/quantum_simulation_cli.py`
- `tools/requirements-plot.txt`

`CMakeLists.txt` 中相应的独立 executable、simulation library、复杂测试和 debug 开关也已撤销。上一版对 `systematic.c`、`local_constraints.c`、`viterbi.c`、`slva.c` 的修改已恢复到 TASK2 完成时的版本。

当前改动集中在：

| 文件 | 内容 |
| --- | --- |
| [main.c](../src/apps/main.c) | Quantum 参数、初始化、point/frame 循环、两算法调用、统计和 CSV |
| [CMakeLists.txt](../CMakeLists.txt) | executable 只链接 quantum；保留算法库测试，移除旧 CLI 的两项测试注册 |
| [run_quantum_sim.sh](../run_quantum_sim.sh) | 常用参数集中在开头，最终 exec 原 `bp_lcgcd` |
| [plot_quantum.py](../tools/plot_quantum.py) | 重新写成仅输出三张 PNG 的脚本 |
| [quantum_sim_cli.py](../tests/quantum_sim_cli.py) | 简短回归，覆盖确定性结果、点循环、共同停止和 candidate cap |
| [示例 Hx](../assets/quantum/five_qubit_hx.txt)、[示例 Hz](../assets/quantum/five_qubit_hz.txt) | 保留直接可用的 5-qubit 输入数据 |
| [README](../README.md)、本报告、[来源记录](third_party_sources.md) | 新入口、操作和验证说明 |

TASK1/TASK2 的 QNMS、SWT/partial SWT、LC-GCD 搜索、Pauli scoring、logical check 和 early stopping 没有修改。用户对 Guide 的原有修改保留。

## 2. main.c 流程

命名与变量组织参考 cpicic：`QuantumSimPara` 保存常用输入，`QuantumStatistic` 只有四个累计计数；`QuantumParaGet` 读取位置参数，`main` 直接执行仿真，`QuantumResultPrint` 直接 fprintf 一行 CSV。

根据后续要求，main.c 中旧 S0～S17 对应的经典单帧 CLI、AWGN/配对仿真、SLVA demo、bits/LLR 列表解析、旧帮助和 dispatch 均已删除；main.c 从 521 行减至 220 行。保留当前 `--quantum` 用法及 `--help/-h`，运行脚本不需要修改。算法库本身保留。

```text
main --quantum
  → QuantumParaGet
  → 读 Hx/Hz，quantum_code_init 一次
  → 初始化 QNMS workspace、BP-LCGCD workspace 一次
  → 打开 CSV，写表头
  → for each p
      channel probabilities → prior XYZ LLR
      清零两份统计与共同 frmNum
      while 未达到共同停止条件
        生成一次 physical Pauli error
        生成一次 measured syndrome
        QNMS 使用该 prior / syndrome
        BP-LCGCD 使用同一 prior / syndrome
        使用已有 quantum_recovery_check 判断两个输出
        累加 failures / BP iterations / LC candidates / BP convergence
      打印这个 p 的一行摘要
      写 QNMS、BP-LCGCD 两行 CSV 并 flush
  → 释放一次 workspace / code
```

只比较固定的 QNMS 和 BP-LCGCD，没有 decoder registry、factory、callback 或 owning report。两份 `QuantumStatistic` 各含 `errFrmNum/iterNum/searchNum/bpConvFrmNum`；帧数为共同标量，不保存逐帧历史。

logical check 复用 BP-LCGCD workspace 的离线 stabilizer form 和 scratch，无需另建 checker。该 API 不覆盖 recovery 或 posterior，因此可以依次检查 QNMS 和 LC 输出。

failure 仍为 `E XOR R` 不属于 `rowspan(S)`，不是 BP syndrome mismatch。若有正常 LC 截断但没有可用 recovery，这帧计为 logical failure，不删除样本。候选 cap 属于实验配置，正常有界搜索仍完成仿真；真实输入、数值或 I/O 错误返回错误码。

共同停止：到达 max frames，或两个 decoder 都达到 max logical failures。设为 0 时禁用 failure 目标。先达到目标的 decoder 继续使用相同帧，直到另一个也达到目标或触及帧上限。

每个 p 从 frame 0 使用相同 seed，延续 TASK1 的 stateless sampler；同一 p/seed 的输入可重复，各 p 有意共享 uniform stream。

## 3. CLI 参数

沿用原 executable，采用接近参考 main.c 的位置参数：

```text
bp_lcgcd --quantum HX HZ P_START P_STOP P_STEP MAX_FRAMES MAX_FAILURES
         BP_ITER NMS_SCALE DELTA CANDIDATE_LIMIT SEED OUTPUT.csv
```

共 13 个位置值，只有输入矩阵、p 范围、帧/错误目标、BP 参数、δ、候选 cap、seed 和输出路径。

| 参数 | 含义 |
| --- | --- |
| HX / HZ | 两个同尺寸 m×n 的 binary matrix 文件；格式为 m n 后接空白分隔的 0/1 |
| P_START / P_STOP / P_STEP | 升序 p 范围，`0<p_start≤p_stop<1`，step>0；落在网格上的 stop 包含在内 |
| MAX_FRAMES | 每个点的共同帧上限，正整数 |
| MAX_FAILURES | 每个 decoder 的 logical-failure 目标，0 禁用 |
| BP_ITER | 两算法使用相同的 BP iteration cap，正整数 |
| NMS_SCALE | 直接乘法 scale，`0<scale≤1` |
| DELTA | LC 的局部约束数，0..rank(Hsyn) |
| CANDIDATE_LIMIT | 每次 LC 完整候选预算，0 不限 |
| SEED | 固定的 unsigned seed |
| OUTPUT.csv | 普通 CSV 路径，当前运行覆盖文件并逐点写入 |

不常修改的设置集中在 main.c 顶部：state_limit=4096、SLVA slots=65536、SLVA bytes=64 MiB、`BP_POST_ON_FAILURE`、`QUANTUM_STOP_TRIVIAL`、degree-one message 默认政策。没有将它们暴露为 CLI。普通 check 的浮点 QNMS 公式保持 TASK1 语义。

## 4. 人工仿真准备与运行

准备 C11 compiler、CMake，以及码的 Hx/Hz。作图另需 Python3 和 matplotlib。输入矩阵应满足 stabilizer 对易条件，不需要准备 G、Hsys 或 logical strings。

先用保留的 `[[5,1,3]]` 示例验证；真实 QLDPC 实验替换脚本中的 HX/HZ。输出的物理长度是 n，逻辑维数是 `n-rank(S)`；允许冗余 stabilizer 行。

首次构建：

```bash
cmake -S . -B build/quantum -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build/quantum -j 4
ctest --test-dir build/quantum --output-on-failure
```

打开 `run_quantum_sim.sh`，在开头修改这几组值：

```bash
HX=assets/quantum/five_qubit_hx.txt
HZ=assets/quantum/five_qubit_hz.txt
P_START=0.04
P_STOP=0.12
P_STEP=0.04
MAX_FRAMES=1000
MAX_LOGICAL_FAILURES=100
BP_MAX_ITERATIONS=12
NMS_SCALE=0.625
DELTA=2
CANDIDATE_LIMIT=128
SEED=2026
OUTPUT=results/quantum/run.csv
```

然后只运行：

```bash
./run_quantum_sim.sh
```

脚本创建输出目录并只调用 `build/quantum/bp_lcgcd`。使用独立构建目录，保留现有 `build/release` 的 Windows/MinGW cache。TASK1/TASK2 原有教学打印保留，脚本将 stdout 重定向到 `/dev/null`；main.c 的参数摘要和每点简短结果写 stderr，因此终端只显示摘要。没有新增 debug/trace 配置层。

先做短测时把 MAX_FRAMES 改为 64、MAX_LOGICAL_FAILURES 改为 0。需要复现同一实验时使用相同矩阵、脚本参数和 seed；保存本次脚本配置与 CSV 即可。

## 5. CSV 和三张图

CSV 仅 10 个字段，每个 p 两行：

```text
p,decoder,frames,logical_failures,ler,bp_iterations_total,bp_iterations_avg,lc_searches_total,lc_searches_avg,bp_converged_frames
```

LER=logical_failures/frames；两种 average 都以全部 frames 为分母。`lc_searches_total` 使用 TASK2 的 `candidates_evaluated`，计完整 Pauli SWT 评分次数；BP 直接返回帧贡献 0，QNMS 的 LC 次数始终为 0。`bp_converged_frames` 单独记录 BP 的 syndrome 匹配数，不作为 logical success 的替代。

没有 JSON、版本 schema、runtime、invocation、histogram 或分位数输出。

安装 matplotlib 后直接画图：

```bash
python3 -m venv build/plot-venv
build/plot-venv/bin/python -m pip install matplotlib
build/plot-venv/bin/python tools/plot_quantum.py results/quantum/run.csv
```

只写到 CSV 所在目录的 `plots/` 下：

1. `ler.png`：LER vs p。
2. `bp_iterations.png`：average BP iterations/frame vs p。
3. `lc_searches.png`：average LC-GCD complete candidates/frame vs p。

均使用线性坐标，直接表示零错误或零搜索，无额外阈值/CDF/格式选项。需要改变输出目录时可用 `--output-dir PATH`。

## 6. 验证记录

单个 Quantum CLI 测试核查：TASK2 的 64 帧参考、两算法共同 frames、range 与单点逐行一致、相同 seed 的 CSV 逐字节重放、共同 failure 停止、candidate cap、帮助入口和旧命令已移除。保留原 17 项库测试；旧 `cli_smoke`、`simulation_cli` 因对应入口已删除而不再注册。

64 帧参考：p=.12、seed=2026、BP_ITER=1、NMS_SCALE=.625、δ=2、candidate limit=0：

| decoder | frames | logical failures | BP iterations total | LC searches total | BP converged |
| --- | ---: | ---: | ---: | ---: | ---: |
| QNMS | 64 | 23 | 64 | 0 | 43 |
| BP-LCGCD | 64 | 10 | 64 | 172 | 43 |

清理旧 main 后的最终结果：**Debug 18/18、Release 18/18、ASan/UBSan 18/18 全部通过**，sanitizer 开启 leak detection，无诊断。两项旧 CLI 测试注册随入口一起移除，算法库测试保留。TASK1/TASK2 算法库的 source/header 与 TASK2 完成 commit 逐文件核对，没有差异；其原有编译警告未在本轮修改。

已实际执行默认 `run_quantum_sim.sh`，终端只输出初始化摘要和三个 p 的短结果。输出 CSV 每点两行、每行 10 字段，三张 PNG 已实际生成。

| p | frames | QNMS logical failures | BP-LCGCD logical failures |
| --- | ---: | ---: | ---: |
| .04 | 1000 | 75 | 20 |
| .08 | 1000 | 136 | 58 |
| .12 | 916 | 197 | 100 |

最后一组由两个 decoder 均达到 failure 目标而停止，前两组由帧上限停止。旧版六图/JSON 示例输出已从结果目录移走；当前示例目录只包含本版 CSV 与三张 PNG。这些短样本用于检查流程和统计，不作为论文性能曲线复现。

当前停在 PAUSE C，不启动量化、多尺度或其他 decoder 开发。未 commit/push。
