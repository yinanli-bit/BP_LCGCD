# TASK1：通用 stabilizer 与浮点 flooding QNMS

开发管理以 [后续开发指南](BP_LCGCD_Next_Development_Guide.md) 为准；旧 S0～S17 保留为经典译码基线。只执行 TASK1，验收位置为 PAUSE A。

## 编码前的 implementation mapping

已实际打开工作区三篇 PDF 并阅读相关正文：

| 资料 / 公式 | 本轮函数 / 字段 |
| --- | --- |
| 2024 BP-LCOSD，III-A、III-B，Eq. (2)、(9)、(10) | `Pauli`、`QuantumCode.hx/hz`、`quantum_code_syndrome`、`hsyn=[Hz\|Hx]` |
| 2024 III-A 的 posterior SWT、III-B Eq. (13)～(15) | 确认 posterior 语义；SWT 和 partial SWT 留给 TASK2 |
| 2025，2.1、2.2 Eq. (3)、(4) | 信道概率转换为 `prior_x/y/z`，初始化 `v2c_x/y/z` |
| 2025，2.2 Eq. (6) | `qbp_commutation_llr`：`log((1+exp(-mu_P))/(exp(-mu_Q)+exp(-mu_R)))` |
| 2025，3.2 Eq. (12) | check update；工程参数 `nms_scale=1/alpha`，例如 0.625 |
| 2025，2.2 Eq. (7)～(10) | posterior XYZ、hard Pauli、真实 syndrome 检查、extrinsic V2C |
| 2026，II-A、II-B.1（印刷页 783，PDF 第 2 页） | generic stabilizer、单 scale quantum NMS；其直接乘法因子对应 `nms_scale` |
| `cpicic_decoder/parity_check_matrix.h/.c` | `Entry/PrtChkMtrx`、四方向 circular header、连续 edge block |
| `cpicic_decoder/dec.c` | `CheckNodeProcMinSum`、`Min1Min2Scaled`、`FloodingIterate` 的工程组织 |

2026 的 V2C 文字在三元更新公式后写了 `P=S[j,i]`，而该 edge 的同 Pauli 分量本应不减 C2V。本轮以 2025 Eq. (10) 明确给出的 `P in {X,Y,Z}` 为准，对全部三个分量更新，并仅在 anticommute 时减去本 edge C2V；这与 2026 的三元 V2C、commutation reduction 和 posterior 定义一致。

计划文件：`quantum.h/.c` 管理只读 code topology；`quantum_channel.h/.c` 管理概率、prior 与可重放 Pauli 样本；`qbp.h/.c` 管理可复用消息 workspace 和 decode 输出。两组测试分别核对表示/信道和消息数学。CMake 接入，README 更新当前开发入口。现有经典译码器不改动。

`S=[Hx|Hz]` 与 `Hsyn=[Hz|Hx]` 严格区分。构造时复制输入、检查二元元素和 stabilizer commutation，然后建立稀疏图；不执行 Gaussian elimination。一个 Y entry 仅建一条带 Y 标签的 edge。

## 数值约定

使用自然对数和 double。commutation reduction 通过 log-sum-exp 等价式避免直接 `exp(-mu)` 溢出。prior 输入要求有限；概率到 prior 的接口要求四个概率严格为正，所以 depolarizing prior 使用 `0<p<1`；采样接口允许端点 `p=0/1`。

正常 degree >= 2 的 check 将 `min1/min2` 初始化为 `INFINITY`，严格计算排除本 edge 后的 `min |lambda|`，不做默认或显式 `message_clip` 饱和。`message_clip` 仅用于 degree-one check：将 `min(empty set)=+infinity` 局部替换为该有限正值（0 取默认 50），再乘 `nms_scale` 并应用 syndrome sign。这一局部数值近似不影响普通 check。posterior/V2C 不限幅，若非有限则返回 numeric error。测试覆盖普通 check 的 minimum 超过 50、显式 `message_clip=1` 不影响普通 check，以及 singleton 的默认/指定值、空行和孤立列。

只有每轮消息更新后的真实 syndrome 匹配才设置 `syndrome_matched`，这不代表 logical success；TASK1 不输出 LER。相等 hard weight 按 I、X、Y、Z 顺序决胜。第一轮从 prior 开始，不增加论文以外的初始 prior 早停规则。

## 已实现的数据流和接口

```text
Hx/Hz -> quantum_code_init -> commutation validation -> Hsyn + QbpGraph
                                                       |
                                               qbp_workspace_init
                                                       |
p -> quantum_depolarizing_probabilities -> quantum_prior_llr -> QbpPrior
                                                       |
physical Pauli error -> quantum_code_syndrome -> measured syndrome
                                                       |
                                                   qbp_decode
                                                       |
                                  QbpResult: hard / post XYZ / iterations / matched
```

- `QuantumCode` 拥有 Hx、Hz、Hsyn 的副本和图。初始化/释放在 offline 阶段；调用方保持对象只读并保证寿命长于 workspace。
- `QbpEdge` 保存 row/col、edge_index、Pauli 与四方向链接；row/column header 是 circular sentinel。Y entry 只有一条 edge。
- `QbpWorkspace` 的一个连续 double buffer 保存 `5*edges + 3*n` 个数：三元 V2C、标量 C2V、commutation cache、三元 posterior；另外保存 n 个 hard Pauli 和 m 个 syndrome bit。
- `QbpPrior` 明确借用 `prior_x/y/z`；`QbpResult` 借用 workspace 内的 `post_x/y/z` 和 hard estimate，下次 decode/free 前有效。prior 不能与 workspace storage 重叠。
- `QUANTUM_OK` 表示译码正常执行，包括用完迭代仍不匹配。必须另外读取 `syndrome_matched`。输入/数值错误清空 result，workspace 可用于下一帧。

每轮依次执行所有 row check updates、所有 column posterior updates、真实 symplectic syndrome 检查，未停止时才按 Eq. (10) 更新 V2C。每条 edge 为常数次访问；包含孤立列和空行时总复杂度为 `O(edges+n+m)`，没有 dense `m*n` BP 扫描。

## 单轮手算

取一行 stabilizer `XZ`，syndrome=1、scale=0.625。两个 qubit 的 prior 为：

| qubit | prior X | prior Y | prior Z | commutation LLR |
| --- | --- | --- | --- | --- |
| 0，edge X | ln4 | ln2 | ln8 | ln((1+1/4)/(1/2+1/8))=ln2 |
| 1，edge Z | ln8 | ln4 | ln2 | ln((1+1/2)/(1/8+1/4))=ln4 |

所以 C2V 为 `(-0.625*ln4, -0.625*ln2)`。q0 的 Y/Z posterior 加第一条 C2V，X 不变；q1 的 X/Y posterior 加第二条 C2V，Z 不变。结果 hard=`YI`，其真实 syndrome 为 1，第一轮停止。测试同时把 syndrome 改成 0，核对 C2V sign 反转与 hard=`II`。

## 文件改动

| 文件 | 作用 |
| --- | --- |
| `include/bplcgcd/quantum.h`、`src/core/quantum.c` | Pauli、QuantumCode、commutation validation、Hsyn、四方向图与 sparse syndrome |
| `include/bplcgcd/quantum_channel.h`、`src/channel/quantum_channel.c` | depolarizing 概率、通用概率到 prior XYZ、可重放采样 |
| `include/bplcgcd/qbp.h`、`src/decoders/qbp.c` | workspace、commutation reduction、floating flooding QNMS、借用输出 |
| `tests/quantum_test.c` | sparse traversal、Pauli 表、4^5 syndrome oracle、信道 prior/采样 |
| `tests/qbp_test.c` | 手算、独立 dense 消息参考、min1/min2、hard/early stop、replay 与内存分配计数 |
| `CMakeLists.txt` | 接入 quantum 库、两项测试、Linux GNU/Clang 下的 allocator 包装 |
| `README.md` | 更新开发管理入口和 TASK1 当前状态 |
| 本文件、`docs/third_party_sources.md` | 公式映射、验收结果、来源/作者/版权/commit/改写记录 |

用户原有的 `docs/BP_LCGCD_Next_Development_Guide.md` 修改完整保留；旧规范、经典 decoder、GF2、systematic 和 simulation 源码未改动。

## 实际验证

环境：WSL2 Ubuntu，Linux `6.6.87.2-microsoft-standard-WSL2`，GCC 13.3.0，CMake 3.28.3。执行日期 2026-11-26（Asia/Shanghai）。

```bash
cmake -S . -B build/wsl-debug -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=/usr/bin/gcc -DBUILD_TESTING=ON
cmake --build build/wsl-debug -j 4
ctest --test-dir build/wsl-debug --output-on-failure
./build/wsl-debug/quantum_test
./build/wsl-debug/qbp_test
```

结果：**18/18 CTest PASS**，包括原有 16 项经典回归。新增 quantum 测试 6455 项检查，QNMS 测试 31704 项检查（含普通 check 不饱和、singleton 局部策略回归）。

- 全部 `4^5=1024` Pauli errors：sparse quantum syndrome、独立 dense symplectic 计算、`Hsyn*[EX;EZ]` 三者一致。这是 syndrome oracle，不是 TASK2 的 quantum candidate/ML oracle。
- Dense QNMS reference 使用概率质量直接求 commutation、逐 outgoing edge 排除邻居求 minimum、逐 incoming check 排除目的行求 extrinsic；不调用生产 commutation/min1/min2/syndrome/hard 函数。比较单轮和多轮 C2V、V2C、posterior、hard、早停次数；覆盖混合 XYZ edge、正/负/零 lambda、minimum ties、scale=0.625/1。
- hard decision 覆盖 I/X/Y/Z 和 ties；degree-one check 与孤立列有显式数值检查；不可能的空行非零 syndrome 用完预算仍不匹配；非有限 posterior 返回 numeric error，并验证下一帧可恢复。
- 固定 sampler：seed=2026，frame=7，p=0.3，前 32 个 Pauli 为 `IIIZIZXIIIIIZIYIIXIIIIIIIIXIIIXX`；相同帧跨调用逐字节重放。uniform stream 由 seed、frame_index、qubit 次序与 quantum domain 常量确定；p 仅决定分类阈值，改变 p 时复用底层 uniform stream。16384 个采样的 I/X/Y/Z 计数为 `11390/1701/1582/1711`，仅作信道分布 smoke。
- 单点 QNMS smoke：通用五 qubit stabilizer 示例，p=0.12，seed=2026，256 帧，max_iterations=12，scale=0.625；非零 syndrome 117 帧，syndrome matched 226 帧，迭代总数 698。每帧与 dense 消息参考比较，posterior/hard/iterations 在独立 workspace 逐字节重放。
- 单帧 frame=7：physical error=`IIIXI`、syndrome=`0001`、hard=`IIIII`、iterations=12、matched=0，验证未收敛输出仍可供未来后处理读取。
- Linux linker `--wrap` 实测完整帧循环：`malloc/calloc/realloc=0`、`free=0`。没有每帧建图、systematize 或 Gaussian elimination。

```bash
cmake -S . -B build/wsl-task1-asan -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=/usr/bin/gcc -DBUILD_TESTING=ON \
  -DCMAKE_C_FLAGS='-fsanitize=address,undefined -fno-omit-frame-pointer -fno-pie' \
  -DCMAKE_EXE_LINKER_FLAGS='-fsanitize=address,undefined -no-pie'
cmake --build build/wsl-task1-asan -j 4 --target quantum_test qbp_test
ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 \
  ctest --test-dir build/wsl-task1-asan -R '^(quantum|qbp)$' --output-on-failure
```

ASan/UBSan：**2/2 PASS**，开启 leak detection，无诊断。全新 sanitizer 构建显示原 `src/core/systematic.c` 的 printf `%d`/size_t 与 signedness 警告；量子新增源码没有编译警告。既有文件本轮未改动。

## PAUSE A 验收与审查重点

| 指南项目 | 证据 |
| --- | --- |
| 1 generic Hx/Hz；2 commutation validation | `quantum_code_init`，混合 Pauli commuting 示例、noncommuting 拒绝 |
| 3 四方向 sparse graph；4 Pauli labels | 正/反 row/column traversal、唯一 edge_index、Y 单边、1024 syndrome oracle |
| 5 flooding quantum NMS；6 非零 syndrome sign | 手算 `XZ` 与独立 dense 多轮参考 |
| 7 posterior XYZ；8 hard decision | 每轮消息参考、全部 Pauli hard/tie 测试 |
| 9 syndrome early stop | 一轮成功、预算耗尽、空行不可能 syndrome |
| 10 depolarizing single-frame | seed=2026 frame=7 及 256 帧单点 smoke |
| 11 prior/posterior 分离 | `QbpPrior.prior_*` 与 `QbpResult.post_*`，输入不变测试 |
| 12 无 hot-loop malloc | linker 包装实测为零，workspace/code 一次初始化重复使用 |
| 13 GCD/LC-GCD regression | 原 16 项 CTest 全通过 |

建议重点阅读三处：`quantum_code_init/build_graph` 的 Hx/Hz 到 edge Pauli 映射；`check_update` 的 syndrome sign 和排除本 edge 的 min1/min2；`posterior_update/variable_update` 的 anticommute 条件与 Eq. (10) 外信息。

TASK1 没有未解决的消息数学桥接问题。仅 degree-one check 使用有限 `message_clip` 的局部近似；普通 check 保留未饱和浮点 NMS。严格正概率 prior 范围按前述数值约定公开保留；短 smoke 不用于推断阈值或 logical error rate。没有进行大规模性能实验。当前入口为库 API/测试程序，统一仿真平台属于 TASK3。

**TASK1 已完成，在 PAUSE A 停止。TASK2、TASK3 和旧 S18 及之后任务均未开始。** 所有项目改动保留为本地未暂存状态，未执行 add/commit/push。
