# BP_LCGCD

本项目按 [逐节点开发规范](BP_LCGCD_Codex_Stepwise_Prompt_v2.md) 开发，参考源码审计见 [S00 报告](docs/source_audit.md)。每个节点完成后等待人工审查，下一节点需要明确授权。所有项目改动保留在本地，由用户自行决定 Git 暂存、提交和推送。

## 当前范围：S13

当前提供 GF(2) 运算、离线系统化、有序 TEP 生成器、三种输入场景共用的单帧 GCD，以及单帧 CLI/trace。`gcd` 已接入真实译码；其他五种算法仍明确返回 `NOT_IMPLEMENTED`，不静默替代。

S08 提供 BPSK/AWGN/LLR 与可重放的单帧噪声；S09 增加单个 Eb/N0 点、有帧数/错误帧数上限的单线程 GCD 仿真和 CSV 摘要。CLI 既可接收指定软信息，也可用 `--simulate` 运行短仿真；尚未实现 SNR 扫描、绘图或其他译码算法。

S10 增加局部约束准备：按每帧翻转权重从固定基列 I 选择 A，拆分 P 与变换后的 syndrome。仅提供 C API 和独立测试，尚未接入 LC-GCD 搜索或新增 CLI 模式。

S11 提供 `[I_delta|P_A]` 的隐式 syndrome Trellis、分支转移与代价、固定起终点及状态数上限。S12 提供单条最优 Viterbi 路径与完整前缀表；S13 已加入按需拉取、无重复、代价非递减的 SLVA 局部路径列表，仍只提供 C API，完整 LC-GCD 留待 S14。

```text
用户命令
  → main(argc, argv)
    → 无参数、--help 或 -h：print_help → stdout，退出码 0
    → --decoder gcd + --simulate：加载 H → 一次系统化 → 有上限帧循环 → CSV
    → --decoder gcd + 单帧参数：加载 H → 离线系统化 → decoder_decode
      → 原 GCD 内核：可选 trace → 最佳结果、证明标志、停止原因和计数
    → 其他已知译码器：NOT_IMPLEMENTED → stderr，退出码 3
    → 其他参数组合：USAGE_ERROR → stderr，退出码 2
```

CLI 从 `argv` 借用参数字符串，负责输入解析和工作缓冲区的分配/释放。单帧运行前只系统化一次；库搜索不进行逐候选消元。输入和最终 estimate 均采用原始行/列次序。

## 构建与检查

主要开发、编译和仿真环境为 **WSL2 Ubuntu**，项目根目录为 `/home/liyinan/BP_LCGCD_CProject`。优先使用 WSL 中的 GCC、CMake、CTest 和 GDB；每个节点报告实际执行环境、命令和结果，不能把 Windows 结果当成 Linux 结果。

需要 CMake 3.20 或更新版本、C11 编译器及构建工具。下面在 WSL 项目根目录执行，复用 `build/wsl-debug`。

```bash
cmake -S . -B build/wsl-debug -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=/usr/bin/gcc -DBUILD_TESTING=ON
cmake --build build/wsl-debug --verbose
ctest --test-dir build/wsl-debug --output-on-failure -V
```

本轮在 WSL2 Ubuntu 24.04.4 LTS、GCC 13.3.0、CMake/CTest 3.28.3 下验证。若现有缓存的编译器路径不同，CMake 可能重置配置；应再运行配置命令并核对 `CMAKE_BUILD_TYPE=Debug`。其他编译器使用独立构建目录。

## GF(2) 数据、格式与手算例子

声明见 `include/bplcgcd/gf2.h`，实现见 `src/core/gf2.c`。`BinaryMatrix` 按行连续保存，每个元素是一个值为 0 或 1 的 `uint8_t`，`H[i,j]` 对应 `data[i * cols + j]`；`BitVector.data[j]` 直接对应第 j 列，索引从 0 开始，不进行位打包或列倒序。

矩阵文本前两个整数为正的行数和列数，其后恰好有 `rows * cols` 个单字符 `0`/`1`。所有 token 以空白分隔，换行仅用于阅读；允许尾部空白和最后一行没有换行，不支持注释、符号、额外 token 或零维矩阵。例如：

```text
3 4
1 0 1 1
0 1 1 0
1 1 0 1
```

给定 `e=(1,0,1,0)`，syndrome 为 `s=H e^T=(0,1,1)`：三行分别得到 `1 XOR 1=0`、`1`、`1`。一般地，`s[i] = XOR_j(H[i,j] AND e[j])`，这就是 GF(2) 矩阵向量乘法，因此共用一个函数。

| 接口 | 本节点调用者、输入和输出 |
| --- | --- |
| `binary_matrix_init` / `bit_vector_init` | 测试传入正尺寸和空对象，获得全零存储；向量随后按列索引填值 |
| `binary_matrix_load` | 测试传入 `FILE *` 和空矩阵；按格式读取整个流，调用方负责关闭流 |
| `binary_matrix_multiply_vector` | 测试传入 H、长度为列数的 e、预分配且长度为行数的 s，计算 syndrome |
| `binary_matrix_multiply` | 测试传入 A、B 和预分配的 C，计算 `C[i,j] = XOR_k(A[i,k] AND B[k,j])`；S03 用它验证行变换关系 |
| `binary_matrix_free` / `bit_vector_free` | 测试释放所拥有的存储，并将对象清零；允许 NULL 和重复释放已清零对象 |

对象首次使用必须以 `{0}` 初始化，禁止浅拷贝所有权。`init`/`load` 要求目标为空，已有对象应先释放；加载失败不保留半成品。乘法不分配内存，输出不能与输入存储重叠；尺寸不匹配或输入含非二元值时返回错误，输出保持原值。尺寸乘积、文本整数溢出、分配和流读取错误均显式返回状态。

## S03：离线系统化、坐标和 syndrome

`binary_matrix_systematize(H, &form)` 在不修改 H 的前提下做一次 GF(2) Gauss–Jordan 消元，确定 `m=rank(H)`、`K=N-m`。主元选择固定为先找列、再找行；行交换和行异或同步作用于从单位阵开始的 U，列交换同步更新 `sys_to_orig`，最后构造其逆 `orig_to_sys`。行异或不能用一个行置换表表示，因此保存完整的 `row_transform=U`，其列始终对应原始校验行。

输出 `BinarySystematic` 保留恒等式 `D=U H Q=[I_m P; 0 0]`，其中 `Q[sys_to_orig[j],j]=1`。`reduced=D` 与原始 H 同尺寸，前 m 行是 `H_sys`，其余零行代表冗余约束，不能静默丢弃。全零矩阵得到 m=0；满列秩矩阵得到 K=0，两者都使用相同表示。

| 接口 | 本节点调用者、输入和输出 |
| --- | --- |
| `binary_matrix_systematize` | 测试传入 H 与 `{0}` 的 `BinarySystematic`，得到秩、D、U 和列映射；失败保持目标不变 |
| `binary_systematic_to_system` | 测试传入原始坐标向量，得到 `e_sys[j]=e_orig[sys_to_orig[j]]` |
| `binary_systematic_to_original` | 测试传入系统坐标向量，得到 `e_orig[j]=e_sys[orig_to_sys[j]]` |
| `binary_systematic_transform_syndrome` | 测试传入原始 syndrome s，得到所有原始行数的 `t=U s`；前 m 位为 `s_sys` |
| `binary_systematic_free` | 释放 D、U 和两个映射，并清零；允许 NULL 和重复释放已清零对象 |

映射和 syndrome 变换要求调用方预分配独立输出，不进行消元或分配内存，供后续每帧复用。若 t 的第 m 位及以后存在 1，返回 `GF2_INCONSISTENT_SYNDROME`，说明原始 `H e^T=s` 无解；**此状态仍输出完整 t 供检查**，其他参数错误保持输出不变。`form` 只供读取，不能修改其成员或浅拷贝所有权。

可手算例子（第三行是前两行的异或）：

```text
H = [0 1 1 0]    U = [1 1 0]    D = [1 0 0 1]
    [0 0 1 1]        [0 1 0]        [0 1 0 1]
    [0 1 0 1]        [1 1 1]        [0 0 0 0]
m=2, K=2; sys_to_orig=(1,2,0,3); orig_to_sys=(2,0,1,3)
e_orig=(1,0,1,0) -> e_sys=(0,1,1,0)
s_orig=(1,1,0) -> U*s_orig=(0,1,0)，一致
s_orig=(1,1,1) -> U*s_orig=(0,1,1)，冗余约束不一致
```

## S04：单帧输入、候选输出与状态

声明见 `include/bplcgcd/decoder.h`，实现见 `src/decoders/decoder.c`。复用 `DecodeInput`、`DecodeResult`，没有新增结构体。输入只读借用，公开接口均使用**原始列坐标**。调用方管理参考向量、残余 syndrome、翻转权重、估计和临时 syndrome 缓冲区；辅助函数不分配内存或进行消元。

| 输入模式 | 输入字段 | reference | s_res | flip_weights |
| --- | --- | --- | --- | --- |
| `DECODE_CLASSICAL` | 有符号信道 `llr` | `z=[llr<0]` | `H z^T` | `abs(llr)` |
| `DECODE_SYNDROME`（纯 syndrome） | `s_target`、非负 `weights`；省略 `z_hint` | 全零 | `s_target` | `weights` |
| `DECODE_BP_SYNDROME` | `s_target`、有符号 BP 后验 `llr` | `z_hint=[llr<0]` | `s_target XOR H z_hint^T` | `abs(llr)` |

显式 `z_hint + flip_weights` 使用 `DECODE_SYNDROME`，把翻转权重放入现有 `weights` 字段；非零 hint 时这些权重表示**相对 hint 的翻转代价**。两个 LLR 模式的正负零都判为 0；BP LLR 仅用于输入准备，不代表已实现 BP。

三种场景归一化后，GCD 及未来 LC-GCD 共用的搜索约定只关心残余 `s_res`、非负 `flip_weights`，搜索 `H f^T=s_res`，唯一评分为 **`W(f)=sum(f[j]*flip_weights[j])`**。TEP 枚举、SLVA 路径、完整候选比较和早停下界均须采用此评分；不保留另一套最终错误代价。输出统一恢复 `reference XOR f`：Classical 得到码字 v 并校验 `H v^T=0`，另外两种场景得到错误估计 e_hat 并校验 `H e_hat^T=s_target`。

BP 模式中，`sum(e_hat[j]*llr[j]) = sum(llr[j] for reference[j]=1) + W(f)`，右侧第一项与候选无关。`sum(e_hat[j]*abs(llr[j]))` 不满足这个恒等式，甚至可能颠倒候选排序。BP 后验的逐比特评分不等于已经证明精确物理/联合 ML 最优。

H 和位向量须是合法二元对象，输入与输出缓冲区互不重叠；候选校验前先成功准备，并保持输入和准备好的权重不变。这些调用方前置条件不在上层逐项重复扫描；保留空指针、尺寸、数值和最终 syndrome 检查，GF(2) 运算直接复用已有函数。S06 使用 S03 的列映射同步转换 f 和权重，并用行变换处理残余 syndrome。

| 接口/结构 | 调用者和约定 |
| --- | --- |
| `DecodeConfig` | 从 `{0}` 初始化；名称、`candidate_limit`（0 不限）以及可选 `FILE *trace`。GCD 按完整候选评估计数，包括空 f_P；trace 流由调用方拥有并保持可写，NULL 关闭 |
| `DecodeInput.systematic` | GCD 必需；只读借用由同一个 H 离线生成的 S03 结果，调用方保证 H 与 form 配对且不变 |
| `decode_prepare(input, reference, s_res, flip_weights)` | 生成归一化搜索输入；调用方预分配 N 位参考向量、原始行数的 syndrome 和 N 个 double 权重；出错后不使用输出 |
| `decode_check_candidate(input, flip_weights, f, result, syndrome_scratch)` | 使用准备好的权重和二元候选 f，恢复估计并复用矩阵乘向量检查 syndrome；临时向量长度为 H 的行数 |
| `DecodeResult` | 从 `{0}` 初始化，调用方分配/释放其 `estimate`；仅 `candidate_valid=true` 时读取估计和代价 |
| `decode_result_reset` | 保留缓冲区和其中字节，清空状态、有效性、最优性标志和本帧计数；旧字节不再代表有效输出 |
| `decode_result_mark_truncated` | 预算耗尽时保留已有候选和统计，设置 `DECODE_TRUNCATED`，清除最优性标志；允许尚无候选 |
| `decoder_is_known` / `decoder_decode` | 共用六个名称；库入口重置本帧元数据，`gcd` 接入 S06，其余名称仍明确返回未实现 |

`DECODE_SUCCESS` 表示当前操作成功；准备方程成功不代表存在候选。候选校验成功时 `candidate_valid=true`，但 `optimality_proven=false`；后者只能由未来译码器的停止证明设置。坏候选返回 `DECODE_INVALID_CANDIDATE`，非法参数和尺寸分别报错；浮点代价溢出返回 `DECODE_NUMERIC_ERROR`。候选检查失败会清除有效性/最优性，保留统计；估计缓冲区的字节不可作为有效输出读取。

`DecoderStats` 分别记录 `candidates_popped`（弹出数，含空模式）、`candidates_evaluated`（完成完整代价计算的候选数）、`nodes_expanded`（出堆后生成至少一个孩子的父节点数）、`batches`（批次数）；串行 GCD 的批次数为 0，四项不互相替代，准备/校验辅助函数不自动增加计数。

沿用上面的 H，可手算以下流程：

```text
classical: r=(-0,-2,3,0) -> z=(0,1,0,0), s_res=(1,0,1)
           f=(0,1,0,0) -> v=(0,0,0,0), W=2
纯 syndrome: s_target=(1,0,1), weights=(5,2,3,7)，默认 z_hint=0
             f=(0,1,0,0) -> e_hat=f, W=2
BP: L=(2,-5,1,-7), s_target=(1,0,1)
    z_hint=(0,1,0,1), flip_weights=(2,5,1,7), s_res=(0,1,1)
    f_A=(0,1,1,0) -> e_hat_A=(0,0,1,1), W_A=6, 有符号评分=-6
    f_B=(0,0,0,1) -> e_hat_B=(0,1,0,0), W_B=7, 有符号评分=-5
    两个评分均相差常数 -12；误用最终错误的绝对值评分会得到 8 与 5，排序相反。
```

## S05：参考 GCDpsgrand 的有序 TEP 生成器

声明见 `include/bplcgcd/tep.h`，实现见 `src/search/tep.c`。输入是长度 K 的固定非负 `flip_weights`，输出部分模式 f_P 和 `W_P=sum(f_P[j]*flip_weights[j])`，可供后续 S06 补齐其他位置。权重只读借用，迭代期间保持不变；内部按 `(权重, 输入索引)` 排序，输出位序仍与传入权重一致，不暴露内部排序坐标。

本轮实际参考 WSL `/home/liyinan/GCDpsgrand` 的 `include/heap.h`、`src/heap.c`（origin 为 `https://github.com/yinanli-bit/GCDpsgrand.git`，HEAD 为 `ebb05a1640de42e3f410f706e6784d0e6afa02e7`，包含用户原有未提交改动）。按其 FPT 思路和必要字段独立编写，没有直接复制源码：

| 参考结构/操作 | 本项目对应实现 |
| --- | --- |
| `tree_node.indices`、`hw`、`LLR` | `TepNode.indices`、`hw`、`cost`；indices 改为从 0 开始的排序位置，cost 是部分翻转代价 |
| `MinHeap.nodes`、`size`、`capacity`、`code_len` | `TepGenerator.nodes`、`size`、`capacity`、`length`；以二叉最小堆保存待扩展节点 |
| `createTreeNode(0,...)` 左孩子 | 若父集合最大索引为 j，将 j 替换为 j+1，汉明重量不变 |
| `createTreeNode(1,...)` 右孩子 | 保留父集合，追加 j+1，汉明重量加 1 |
| 初始化 | 首先把空模式放入堆；空模式只有特殊孩子 `{0}`，K=0 时没有孩子 |

每次只取一个最小节点，生成其左右孩子后再取下一个。参考 `extractParMin` 的批量前缀提取未迁入；S00 已记录批取后才扩展不能保证连续候选的全局顺序。S05 只保留上述必要字段，不加入 Depth、LW/HW 剪枝或完整候选代价。

正确性来自以下不变量：排序权重下，左孩子的代价增量是 `w[j+1]-w[j]>=0`，右孩子是 `w[j+1]>=0`。对任意非空集合，最大索引 j 的前一位 j-1 若在集合中，其唯一父节点是删去 j；否则父节点是将 j 换回 j-1，`{0}` 的父节点为空。因此每个集合恰好可达一次，配合最小堆可按非递减代价枚举全部模式。

相同代价按排序坐标中的位模式作无符号二进制顺序比较，高位优先，空模式最先；不使用 epsilon。每个孩子都按已选排序位置重新求和，避免浮点减法更新产生消减误差。初始化要求所有权重有限且非负、排序后的总和有限。堆按需扩容，节点只保存选中索引，不预生成全体子集，也不限制到 64 位。

| 接口 | 调用约定 |
| --- | --- |
| `tep_generator_init(&g, weights, K)` | g 从 `{0}` 初始化；校验权重并建立内部排序和空根节点；失败保持 g 为空 |
| `tep_generator_next(&g, &f_P, &cost)` | 调用方预分配 K 位输出；每个 `TEP_OK` 恰好一次出堆、一个候选，包括空模式；穷尽后为 `TEP_DONE` |
| `tep_generator_peek(&g, &cost)` | 只读查询剩余最小部分代价，不分配、出堆或扩展；S06 用它先检查证明，再检查预算 |
| `tep_generator_free(&g)` | 释放尚未消费的节点和排序表，清零 g；可以提前结束枚举后释放，权重仍由调用方拥有 |

K=0 时使用 `BitVector {0}`，只输出一次空模式。输出和输入/迭代器存储须互不重叠，不做重复别名检查。必要的空指针、长度、权重和分配失败检查保留；分配失败不跳过候选，不静默截断。S05 无查询预算；`TEP_OK` 次数不能当作完整候选评估或批次数，后续调用者分别维护 S04 的统计字段。

权重 `(3,1,2)` 的前五个输出如下（位序保持输入顺序）：

```text
(0,0,0) W_P=0
(0,1,0) W_P=1
(0,0,1) W_P=2
(0,1,1) W_P=3
(1,0,0) W_P=3
```

## S06：补齐完整翻转模式与单帧 GCD

调用方先在帧外执行 `binary_matrix_systematize(H, &form)`，把 `input.systematic=&form`，并预分配 `result.estimate`。以 `DecodeConfig{.decoder_name="gcd", .candidate_limit=0}` 调用 `decoder_decode`。H、form、软信息与目标 syndrome 均只读借用；GCD 每次分配并释放必要工作向量，不进行每帧消元，也不重新构造行列映射。

数据流为：S04 准备 reference、s_res、flip_weights → S03 得到 `t=U*s_res` 并检查冗余行 → 同步置换权重 → S05 枚举 f_P → 按 `f_I=t[0:m] XOR P*f_P^T` 补齐 → 比较完整 W → 逆列映射并恢复最佳 estimate → 用原始 H 检查输出。为复用 GF(2) 乘法，代码以 `H_sys*(0,f_P)^T` 计算 P*f_P，不另复制 P。

这一补齐思想对应参考 GCDpsgrand `heap.c:ParSgrandCalcLLR` 的 syndrome 异或列累加和完整模式计分；本项目采用已批准的 `[I_m|P]` 坐标及统一翻转度量，复用自己的 S03—S05 实现，未复制参考源码。

每个 f_P 恰好对应一个满足约束的 f_I，因此枚举全部 f_P 覆盖所有可行翻转模式。非负权重保证 `W(f)=W_P(f_P)+W_I(f_I)>=W_P(f_P)`；未访问模式的最小 W_P 是完整代价下界。找到合法候选后继续搜索，直到下界或穷尽给出证明。

| 停止情况 | 状态和结果 |
| --- | --- |
| 下一 W_P 不小于当前最佳完整 W | `DECODE_SUCCESS`、`optimality_proven=true`、`stop_reason=DECODE_STOP_BOUND` |
| 部分模式全部穷尽 | `DECODE_SUCCESS`、`optimality_proven=true`、`stop_reason=DECODE_STOP_EXHAUSTED` |
| 预算耗尽且未获证明 | `DECODE_TRUNCATED`、保留当前最佳、`optimality_proven=false`、`stop_reason=DECODE_STOP_BUDGET` |
| 原始 syndrome 与冗余约束不一致 | `DECODE_INCONSISTENT_SYNDROME`、无候选、无搜索计数 |

证明优先于预算判断：恰好达到预算时，若已经能证明最优，则仍返回成功。相同完整代价保留首个候选，因此下界相等也允许停止。代价按 S05 排序部分和开始，再累加基区权重，使 double 计算中 W_P 始终不大于其完整代价；没有引入第二种评分。BP 场景证明的是给定后验翻转度量的最优性，不代表联合物理先验或量子逻辑纠错成功。

例子：`H=[1 0 1 0; 0 1 0 1]`，`s_target=(1,1)`，`weights=(10,10,1,1)`，reference=0：

```text
f_P=(0,0) -> f_I=(1,1), W=20
f_P=(1,0) -> f_I=(0,1), W=11
f_P=(0,1) -> f_I=(1,0), W=11
f_P=(1,1) -> f_I=(0,0), W=2
```

预算 1、2 分别返回代价 20、11 的当前最佳；无预算时得到代价 2，穷尽后证明最优。此例弹出/完整评估各 4 次，真正生成孩子的节点只有 2 个，批次数为 0。

## 命令行为

S07 单帧入口使用 `--decoder gcd --matrix PATH`，默认 `--mode classical`。软信息是一整个带引号的空白分隔字符串，长度必须为 N；syndrome/hint 是连续的 0/1 字符串，长度分别为原始行数/N，最左字符是索引 0。非法数值、长度或模式参数组合返回使用错误，不作自动截断/填充。

| 模式 | 必需软信息与目标 | reference 与输出 |
| --- | --- | --- |
| `classical` | `--llr "..."` | `reference=[LLR<0]`，输出码字 v |
| `syndrome` | `--syndrome BITS --weights "..."` | 默认 reference=0，可用 `--hint BITS`；输出 e_hat |
| `bp-syndrome` | `--syndrome BITS --llr "..."` | 从外部 BP 后验 LLR 的符号/绝对值准备输入，输出 e_hat；不执行 BP 算法 |

`--limit N` 是完整候选预算，0 表示不限。`--trace` 将准备信息、每个实际评估的候选和停止条件打印到 stdout；关闭 trace 不改变结果、候选顺序或计数。它直接观察原 GCD 循环，不另跑一个搜索器。系统坐标中的 `f_I/f_P` 会明确标为 `(sys)`，同时打印从 0 开始的 `sys_to_orig` 和 H_sys，避免与原始坐标混淆；空区间打印 `-`。

| 调用 | 输出流与内容 | 退出码 |
| --- | --- | --- |
| 无参数、`--help` 或 `-h` | 标准输出：帮助信息 | `0` |
| `gcd` 帧输入，已证明最优 | 标准输出：结果、代价、证明标志、停止原因及计数；可选 trace | `0` |
| `gcd` 候选预算耗尽、尚未获证明 | 标准输出：当前最佳，`TRUNCATED`、`optimality_proven=0` | `5` |
| 其他五种已知译码器 | 标准错误：`NOT_IMPLEMENTED` | `3` |
| 未知译码器、缺少名称、未知参数或多余参数 | 标准错误：`USAGE_ERROR` | `2` |
| 矩阵文件读取/格式/分配失败，或译码错误（含不一致 syndrome） | 标准错误：`INPUT_ERROR` / `DECODE_ERROR` | `4` |

名称区分大小写，只接受帮助中列出的参数形式。仅写 `--decoder gcd` 现在会报告缺少单帧输入，退出码为 2。

## S07：7 位手算与真实 trace

以下矩阵由本项目自行构造，是系统形式的 Hamming [7,4] 校验矩阵，未复制第三方矩阵文件。在 WSL 项目根目录准备并运行：

```bash
cat > build/wsl-debug/h7.txt <<'EOF'
3 7
1 0 0 1 1 0 1
0 1 0 1 0 1 1
0 0 1 0 1 1 1
EOF
./build/wsl-debug/bp_lcgcd --decoder gcd --matrix build/wsl-debug/h7.txt \
  --llr "-8 -3 -6 1 2 4 7" --trace
```

此时 I=(0,1,2)、P 区=(3,4,5,6)，列映射是恒等映射，U=I。`z=1110000`，`s_res=111`，翻转权重为 `(8,3,6,1,2,4,7)`。补齐公式为 `f_I=111 XOR P*f_P^T`，完整代价为 W_P 加上 f_I 对应的基区权重；这里 f 是相对 z 的翻转模式。

| 候选 | f_P | W_P | P*f_P^T | 补齐 f_I | 完整 W | 当前最佳 |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 0000 | 0 | 000 | 111 | 8+3+6=17 | 17 |
| 2 | 1000 | 1 | 110 | 001 | 1+6=7 | 7 |
| 3 | 0100 | 2 | 101 | 010 | 2+3=5 | 5 |

实际打印的前三条候选为：

```text
trace candidate=1 f_P(sys)=0000 W_P=0 f_I(sys)=111 W=17 best=17 improved=1
trace candidate=2 f_P(sys)=1000 W_P=1 f_I(sys)=001 W=7 best=7 improved=1
trace candidate=3 f_P(sys)=0100 W_P=2 f_I(sys)=010 W=5 best=5 improved=1
```

接着评估 `f_P=1100, W=11` 和 `f_P=0010, W=12`，均不更新最佳。下一部分代价为 5，满足 `next_W_P=5 >= best=5`，因此不再弹出该候选，停止原因是 `bound`。最优翻转模式为 `f=0100100`，最终码字为 `v=z XOR f=1010100`，原始 H 校验为零，W=5；弹出/完整评估/扩展各 5 次，批次数 0。

加 `--limit 3` 时，尚未访问部分代价的下界为 3，小于当前最佳 5，因此返回退出码 5、`TRUNCATED`、`optimality_proven=0`。即使此时恰好找到真实最优，也不能声称已经证明；设 `--limit 5` 则下界证明优先于预算，仍成功退出。

同一 H 的另外三种输入方式如下，均使用原始行/列次序：

```bash
# 纯 syndrome：输出 e_hat=0100100，W=5。
./build/wsl-debug/bp_lcgcd --decoder gcd --matrix build/wsl-debug/h7.txt \
  --mode syndrome --syndrome 111 --weights "8 3 6 1 2 4 7"
# 非零 hint：输出 e_hat=1010100，W=5，目标 syndrome 为 000。
./build/wsl-debug/bp_lcgcd --decoder gcd --matrix build/wsl-debug/h7.txt \
  --mode syndrome --syndrome 000 --hint 1110000 --weights "8 3 6 1 2 4 7"
# 已提供的 BP 后验 LLR：同样输出 e_hat=1010100，W=5。
./build/wsl-debug/bp_lcgcd --decoder gcd --matrix build/wsl-debug/h7.txt \
  --mode bp-syndrome --syndrome 000 --llr "-8 -3 -6 1 2 4 7"
```

## S08：BPSK、AWGN、真实 LLR 与帧重放

声明见 `include/bplcgcd/channel.h`，实现见 `src/channel/channel.c`。固定符号能量 Es=1，噪声是每个实 BPSK 符号上的独立高斯噪声：

```text
R = K/N = (N-rank(H))/N                 # 不用原始校验行数代替 rank
gamma_b = 10^(ebn0_db/10)              # 输入是 Eb/N0 的 dB 值
sigma2 = N0/2 = 1/(2*R*gamma_b)
x[j] = 1-2*c[j]                       # 0 -> +1，1 -> -1
y[j] = x[j] + sqrt(sigma2)*G[j]        # G ~ N(0,1)
LLR[j] = log(p(y[j]|0)/p(y[j]|1)) = 2*y[j]/sigma2
z[j] = [LLR[j]<0]                     # 正 LLR 倾向 0；零值判 0
```

例如 R=4/7、Eb/N0=0 dB 时 sigma2=7/8，而不是 1/2；提高到 10 dB 后方差缩小十倍。码率须在 (0,1]；K=0 时不定义这里的信息比特 Eb/N0。

| 接口 | 当前调用者与输入/输出 |
| --- | --- |
| `bpsk_noise_variance(rate, ebn0_db, &sigma2)` | 单帧测试/后续仿真驱动：由码率、Eb/N0 得到正方差；不可表示的极端数值显式报错 |
| `bpsk_modulate(&word, symbols)` | 单帧测试及 transmit：只读借用原始坐标的 BitVector，输出 N 个 ±1 |
| `bpsk_awgn_transmit(&word, sigma2, seed, frame_index, received)` | 输出一帧 y；相同 seed、帧号和方差在同一运行环境可逐样本重放，不依赖其他帧调用的先后 |
| `bpsk_compute_llr(received, N, sigma2, llr)` | 输出 N 个真实信道 LLR，可直接借给 `DecodeInput{.mode=DECODE_CLASSICAL, .llr=llr}` |

模块不分配内存，输出由调用方提供足够长度且与输入不重叠的缓冲区；错误后不使用输出。只检查必要的指针、长度、二元值和有限数值，没有重复的状态/别名扫描。`sigma2=0` 仅用于 transmit 的精确无噪声测试，此时 y=x；计算有限 LLR 仍须给正的标定方差，不能计算 `2*y/0`。

内部整数 RNG 使用 SplitMix64 固定递推，局部初态为 `seed XOR mix64(frame_index)`，不使用全局 rand/drand48、系统时间或外部随机文件。每个符号固定消费两个均匀样本，以 Box–Muller 生成一个高斯样本，不缓存配对样本。均匀数取 52 个高位并加半个量化格，严格落在 (0,1)，避免 log(0)；奇数长度帧同样可以重放。整数运算固定为 uint64_t，但高斯变换依赖 libm，不承诺不同数学库之间的浮点位级一致。

本轮实际只读核对 `/home/liyinan/GCDpsgrand/src/transmit.c` 的 `BpskModuOverAwgnChannelUsingEbN0`、`lib/rand.c` 的 `RandGaussStdReentrant`：沿用 BPSK、方差和 Box–Muller 思想，独立编写模块。参考 transmit 返回的是接收幅度 y，本项目明确再计算 `2*y/sigma2`。参考仓库未确认可复用许可证，且 rand.c 有版权声明；本轮未复制其源码或 GNU 随机数依赖。

全零发送仅在二元线性码、对称 BPSK/AWGN 信道和保持对称性的译码规则下可代表其他码字。本节点既检查全零，也检查有效非零码字 `1010100` 的无噪声恢复；不把 syndrome 一致当作恢复了真实发送码字。固定噪声帧还按独立的码字欧氏距离穷举 oracle 核对 GCD 的 AWGN ML 结果。

在 WSL 中单独查看本节点结果：

```bash
ctest --test-dir build/wsl-debug -R '^channel$' --output-on-failure -V
```

当前 105 项检查覆盖参数/符号、独立似然比公式、零噪声、固定样本、seed/帧号重放、不同帧调用次序、噪声标准差缩放、码字无关噪声及单帧 GCD。另以一帧 8192 个噪声样本做均值/方差烟雾检查，不估算 FER/BER，也不据此宣称统计认证或论文曲线复现。

## S09：有上限的 GCD Monte Carlo 与参考样例

新增 simulation.h / src/simulation/simulation.c，复用 S08 信道和 S06 GCD；主程序只解析参数、加载 H、系统化一次、调用 simulation_run 并写 CSV。SimulationConfig 记录 Eb/N0、seed、两个正上限；SimulationResult 保存整数计数，直接复用 DecoderStats 累加各帧搜索计数，没有另造模糊的 query/batch 计数。

每帧发送全零码字，按 (seed, frame_index) 生成 AWGN，再计算真实 LLR、调用 Classical GCD。每个有效输出都与实际发送的零码字比较：通过 H 校验不代表恢复成功。预算截断时保留当前最佳码字，仍参与误码统计，同时累计 truncated_frames；证明帧另计为 proven_frames，不能把所有输出都称为 ML。

~~~text
codeword_ber = bit_errors / (实际 frames * N)
fer = frame_errors / 实际 frames
hd_codeword_ber = hd_bit_errors / (实际 frames * N)
avg_candidates_evaluated = 完整候选评估总数 / 实际 frames
~~~

累计完当前帧后，只要实际帧数达到 --frames 或误帧数达到 --errors 就停止；两者同时达到时记 error_limit。误帧上限不是译码器的每帧候选预算 --limit。两个仿真上限都必须为正，默认 seed=12344567，与参考工程默认种子数值相同，但 RNG 与帧序列不同，不声称逐帧复现其噪声。

--csv PATH 每次覆盖写入表头和一行摘要，包含 N/K、Eb/N0、seed、方差、上限、实际帧数、误比特/误帧、三种误码率、四种搜索计数、平均完整评估数、证明/截断帧数和仿真停止原因。没有截断时退出 0，至少一帧截断则退出 5（CSV 仍有效）；参数错误退出 2，运行/输出错误退出 4。既有学习调试打印仍在 stdout，CSV 使用独立文件，不依赖解析这些打印。

### 实际参考与矩阵核对

只读核对 /home/liyinan/GCDpsgrand（HEAD ebb05a1640de42e3f410f706e6784d0e6afa02e7）的 codec_statistic.c 单 SNR 循环、top_ctrl.h 上限/种子、main.c 注释中的 toy 样例及 tools/parse_gcd_log.py。参考程序 BER 比较信息位；这里没有生成矩阵/信息位反投影，因此明确报告整个码字 BER，不能直接把两者当同一指标。参考查询值可能是批次数，本项目只使用已定义的实际候选计数。未复制参考源码或将许可未确认的数据作为项目分发文件。

本轮以独立 Python 整数消元核对下列文件，G/H 均满行秩且 G*H^T=0：

| 参考文件前缀 | H 尺寸 / 秩 | G 尺寸 / 秩 | 实际码参数 |
| --- | --- | --- | --- |
| res/toy | 4×8 / 4 | 4×8 / 4 | [8,4] |
| res/ldpc96_48 | 48×96 / 48 | 48×96 / 48 | [96,48]，当前参考默认配置 |
| res/eBCH | 8×128 / 8 | 120×128 / 120 | [128,120]，不是论文 [128,106] |

本次 H 原文件 SHA256（对应上述顺序）：

~~~text
e7da096ed69407db5e2be2fd0dd83298c1444139dcf4acb5b5a6657868395a62
9b127cd4431bd2e4812e5c4f14244352c842a279af2e4851d0a05603b37fb946
ff81c8a317fb80dba62e9dcaca3772eb010ab96a787e5e717fb39bb1cdbd7d9c
~~~

参考 .mat 是无尺寸头的二元行文本。下面只在被忽略的 build 目录生成本项目加载格式，保留行列顺序；重新使用时应核对上面的来源/散列：

~~~bash
python3 - <<'PY'
from pathlib import Path
for name in ('toy', 'ldpc96_48', 'eBCH'):
    p = Path('/home/liyinan/GCDpsgrand/res') / (name + '_pc.mat')
    rows = [line.split() for line in p.read_text().splitlines() if line.strip()]
    assert rows and all(len(r) == len(rows[0]) and set(r) <= {'0', '1'} for r in rows)
    out = Path('build/wsl-debug') / (name + '_h.txt')
    out.write_text(f'{len(rows)} {len(rows[0])}\n' + '\n'.join(' '.join(r) for r in rows) + '\n')
PY
./build/wsl-debug/bp_lcgcd --decoder gcd --matrix build/wsl-debug/toy_h.txt \
  --simulate --ebn0 2 --frames 100 --errors 100 --seed 12344567 \
  --csv build/wsl-debug/s09-toy.csv > build/wsl-debug/s09-toy.log
./build/wsl-debug/bp_lcgcd --decoder gcd --matrix build/wsl-debug/ldpc96_48_h.txt \
  --simulate --ebn0 5 --frames 100 --errors 100 --seed 12344567 --limit 16 \
  --csv build/wsl-debug/s09-ldpc96.csv > build/wsl-debug/s09-ldpc96.log
~~~

实际 toy 的 100 帧为 5 误帧、17 误比特，码字 BER=0.02125、FER=0.05，平均完整评估 1.59 次，无截断。LDPC 的 100 帧有 79 帧截断；该运行仅验证大一些矩阵和预算统计，不能视为完整 ML 性能。另用实际 [128,120] 文件运行了 100 帧、预算 16，有 15 帧截断。完整来源核对与 G/H 散列在本轮 build/wsl-debug/s09-reference-audit.json，CSV/log 均不作为源数据提交。

参考 main.c 中注释的 8 位固定接收样例也实际重跑：

~~~bash
./build/wsl-debug/bp_lcgcd --decoder gcd --matrix build/wsl-debug/toy_h.txt \
  --llr "0.4 0.7 1.1 1.6 -3 -3 -3 -3" --trace
~~~

输出 10101111，翻转代价 1.5；独立枚举该 H 的 16 个码字核对最小代价通过。原变量是接收幅度，这里只将相同数值作为等比例软信息，不能据此推定其原始噪声方差；正常仿真始终使用 S08 真实 LLR。

本节点不绘图。参考 parse_gcd_log.py 实际绘的是对数 FER 与平均批次的双轴图；后续可分别绘 codeword BER/FER 对 Eb/N0 和平均候选数，并明确不同分母/计数单位。仅有 100 帧的结果是 smoke 验证，观测零误帧也不能推断真实 FER=0，不声称复现论文曲线。

## S10：选择 A 与构造论文 Eq. 9/10

本轮重新读取提供的论文 PDF 第 2 页（印刷页 443）II-C 的预处理及 Eq. 9/10。A 是固定基列 I 内最可靠的 delta 个位置，不能从全部 N 位或自由列 P 中选。论文的 TEP e 在本项目统一表示为相对 reference 的翻转 f；可靠度使用 S04 准备好的非负 flip_weights，两个约束仍使用同一翻转度量。

入口 `local_constraints_prepare(form, flip_weights, N, s_res, delta, &local, &scratch)`：

```text
原始坐标 flip_weights、s_res
  → S03 行变换 t=U*s_res，并检查冗余行一致性
  → 只排序 sys_to_orig[0:m] 对应的原始基列
  → 按权重降序、相同权重按原始列号升序排列
  → 前 delta 个构成 A，其余构成 Ac
  → 同步抽取 P 的对应行及 t 的对应元素

f_I[A]  XOR P_A * f_P^T  = s_A          (Eq. 9)
f_I[Ac] XOR P_Ac * f_P^T = s_Ac         (Eq. 10)
```

`LocalConstraints` 只组合已有 BinaryMatrix、BitVector 和一个 basis_order 指针，不拥有内存。调用方提供长度 m 的 basis_order、两个 P 块、两个 syndrome 子向量及原始行数长度的 scratch。函数无分配、无高斯消元；输入只读借用，只保留指针、维度、非负有限权重及 syndrome 一致性等必要检查，不重复验证 S03 映射或扫描别名。

| 输出 | 尺寸和坐标 |
| --- | --- |
| basis_order | m 个**原始列号**；前 delta 个是 A，后 m-delta 个是 Ac |
| p_a / p_ac | delta×K / (m-delta)×K；行顺序与 basis_order 对应分段一致 |
| s_a / s_ac | 对应行的 t=U*s_res，不能直接截取原始 s_res |
| 两个 P 块的第 j 列 | 对应原始列 sys_to_orig[m+j]，即共同的系统自由区顺序 |

调用方可用栈数组或现有 GF(2) init 分配输出。delta=0、delta=m、m=0、K=0 都支持：零维矩阵保留 rows/cols，data=NULL；空向量为 {0}。原有 binary_matrix_init 不接受零维，所以这些空视图由调用方直接标注尺寸；不更改基础 GF(2) 模块的约定。仅在返回 GF2_OK 时使用结果，不一致 syndrome 返回 GF2_INCONSISTENT_SYNDROME。

### 可手算的非系统矩阵

```text
H = [0 1 1 0 1 0]       s_res = (1,1,0)
    [0 0 1 1 0 1]       flip_weights = (1000,3,6,900,800,700)
    [0 1 0 1 1 1]       第三行是前两行的异或

m=2, K=4, sys_to_orig=(1,2,0,3,4,5)
H_sys = [1 0 | 0 1 1 1]       t=U*s_res=(0,1,0)
        [0 1 | 0 1 0 1]
```
H = [1 0 0 1 1 1]       
    [0 1 0 1 0 1]       
    [0 0 0 0 0 0]

固定基列的原始编号是 I={1,2}，权重分别为 3、6。delta=1 时选 A={2}、Ac={1}；即使自由列权重达到 1000，也不参与选择。P 块列序对应原始列 (0,3,4,5)，因此得到 P_A=0101、s_A=1，P_Ac=0111、s_Ac=0。以原始翻转位写出约束就是：

```text
f[2] XOR f[3] XOR f[5] = 1
f[1] XOR f[3] XOR f[4] XOR f[5] = 0
```

例如 f=(0,1,0,1,0,0) 满足两式，也满足原始 H*f^T=(1,1,0)。这里 A 的编号 2 是原始**列号**，不是原始 syndrome 的行号；右端需要先做 U 行变换。

`local_constraints` 测试有 4049 项检查、82 组划分，覆盖三种 S04 输入准备、所有 delta 边界、rank=0、K=0、零权重、ties、列置换、冗余行及不一致 syndrome。对小矩阵穷举全部原始坐标翻转 f，独立比较“两个分组约束同时满足”与“Hf=s_res”，并用计数排序 oracle 核对 A。delta=0 目前只验证约束退化，未实现或声称已经验证 LC-GCD 搜索等价。

```bash
ctest --test-dir build/wsl-debug -R '^local_constraints$' --output-on-failure -V
```

后续 S11 可将本模块 p_a 与局部基位组成 [I_delta|P_A]；本节点不构造 Trellis，不枚举路径，也不改变现有 GCD 算法或 CLI 行为。

## S11：Trellis 状态、分支与路径代价

实际来源及 S12/S13 的详细衔接见 [Trellis / SLVA 原文核对笔记](docs/trellis_slva_notes.md)：1994 论文 §2.1–2.4 的递推、最后合并、候选栈，以及与本项目平行边、翻转代价和 LC-GCD 停止条件的区别。本轮已提取原文并视觉核对 PDF 第 2–4 页公式及图 2–5，没有复制源码或提前实现列表算法。

S10 的局部约束是 `f_I[A] XOR P_A*f_P^T=s_A`。以 `(f_I[A],f_P)` 为固定列顺序，构造长度 `L=delta+K`、L+1 层的隐式图。状态 q 是已处理列的累计 syndrome，状态位 i 对应 A 的第 i 个局部约束；原始列映射显式保存在 original_columns 中。

```text
M_local = [I_delta | P_A]
q_0 = 0
q_(t+1) = q_t XOR (b_t ? column_masks[t] : 0)
branch_cost = b_t * flip_weights[original_columns[t]]
path_cost = W(f_I[A]) + W(f_P)
接受条件：恰好 L 个分支后，q_L = pack(s_A)
```

Ac 不参与当前路径评分；它的翻转代价留到 S14 补齐时加入。路径终点满足 Eq.9 不代表已经得到最优完整译码结果。图使用每层最多 2^delta 个状态的字母表，不分配整张状态表、不做可达性剪枝；未来 DP 要把起点以外的初始度量置为无穷。

| 接口 / 字段 | 约定 |
| --- | --- |
| `syndrome_trellis_prepare(form, local, weights, N, state_limit, &t)` | 复用配对的 S03/S10 输出；weights 是原始坐标翻转权重，图使用期间保持不变 |
| `t.length`、`column_masks`、`original_columns` | 调用方提供 length=L 和两个 L 长度输出数组；函数无分配；L=0 时数组可为 NULL |
| `t.delta`、`state_count`、`target` | 构图成功后输出；起点恒为 0，终点为打包的 s_A |
| `syndrome_trellis_branch(&t, stage, state, bit, &next, &cost)` | 给定一条带 0/1 标签的边，返回后继和非负增量代价；不搜索路径 |
| `state_limit` | 必须为正，限制每层 2^delta 个状态；超预算或状态数不能用 size_t 表示时返回 TRELLIS_RESOURCE_LIMIT，先检查位宽再移位 |

即使列掩码为零，bit=0/1 也是两条不同边。delta=0 时只有一个状态，但仍有全部 2^K 条带标签路径；不能根据状态序列相同而合并。S11 只占两个 O(L) 的调用方数组，S12/S13 还要检查自身度量、回溯表及候选缓存的实际内存预算。

### delta=1/2 的逐步例子

沿用 S10 的 3×6 H、s_res=(1,1,0)，本例权重改为 `(1000,3,6,9,8,7)`。delta=2 时 A=(2,1)，自由列为 (0,3,4,5)，故阶段对应原始列 `(2,1,0,3,4,5)`，列掩码为 `(1,2,0,3,2,3)`，目标 pack(s_A)=1。整数状态的最低位对应 A[0]。

对原始翻转 f=(0,1,0,1,0,0)，实际测试打印如下：

| 阶段 t | 原始列 | 分支 bit | 转移 | 增量 | 累计 |
| --- | --- | --- | --- | --- | --- |
| 0 | 2 | 0 | 0→0 | 0 | 0 |
| 1 | 1 | 1 | 0→2 | 3 | 3 |
| 2 | 0 | 0 | 2→2 | 0 | 3 |
| 3 | 3 | 1 | 2→1 | 9 | 12 |
| 4 | 4 | 0 | 1→1 | 0 | 12 |
| 5 | 5 | 0 | 1→1 | 0 | 12 |

终点 1 等于目标，满足 Eq.9。delta=1 时 A=(2)，阶段原始列变为 `(2,0,3,4,5)`，掩码为 `(1,0,1,0,1)`；同一 f 对应路径比特 `(0,0,1,0,0)`，状态依次 `0→0→0→1→1→1`，局部代价为 9。原始列 1 属于 Ac，其代价 3 此时没有计入。

```bash
ctest --test-dir build/wsl-debug -R '^trellis$' --output-on-failure -V
```

新测试有 8022 项检查、39 个小图，穷举带标签路径验证终点与 Eq.9 等价、代价等于直接权重和、合法路径数恰为 2^K；覆盖 delta=0、K=0、L=0、零列、零权重、原始列映射、非零目标和状态资源限制。测试穷举仅作为 oracle，不是生产 SLVA。

## S12：单条最优 Viterbi 路径

新增 `viterbi_decode(&trellis, cells, cell_capacity, &path, &cost)`，复用 S11 图和分支代价，按 1994 论文 §2.1 的初始化、递推及回溯求固定目标状态的最轻局部路径。[算法笔记第 5 节](docs/trellis_slva_notes.md) 给出公式、每格字段、逐层度量表、回溯、ties 和资源约定。

`ViterbiCell` 保存 metric、previous_state 和 bit，调用方按 `cells[layer*state_count+state]` 提供 `(L+1)*state_count` 格工作表。初态 0 的代价为 0，其他格为无穷；不可达格前驱为 SIZE_MAX。每个下一层状态分别比较 bit=0/1 两条入边，严格更小才更新，因此平行边不会被当成同一选择，ties 固定保留先比较的 bit=0。

```text
phi_(t+1)(q) = min(phi_t(q) + c_t(0),
                   phi_t(q XOR column_masks[t]) + c_t(1))
```

结束后仅从 `target` 回溯，输出长度 L 的 `(f_I[A],f_P)` 阶段比特及局部代价 `W(f_I[A])+W(f_P)`，不是 N 位最终错误估计。工作表保留所有层的最佳前缀及前驱，现由 S13 SLVA 复用；该函数本身仍只返回一条最优路径。

函数不分配内存，cell_capacity 就是本次调用允许使用的工作格数；写入前检查 L+1、格数乘法、字节数乘法及容量。不够返回 TRELLIS_RESOURCE_LIMIT。不可达目标返回 TRELLIS_NO_PATH，仅成功后读取 path/cost。

S11 的 delta=2 例子中，最后一层四个状态的最佳度量是 `(0,6,3,7)`；指定 target=1，故答案为 6，而不是任意终态中的最小值 0。回溯得到阶段路径 `100000`，即 `f_I[A]=10、f_P=0000`，映射到原始非零翻转列 2。该例需要 28 格工作表，当前 WSL 下共 672 字节。S11 原先演示的代价 12 是另一条合法路径，S12 已真正求得最优值。

```bash
ctest --test-dir build/wsl-debug -R '^viterbi$' --output-on-failure -V
```

新测试有 2461 项检查、60 个穷举小图，独立核对所有层/状态的最小前缀代价及前驱，并检查最终路径、ties 和重复调用；另覆盖零列/零权重、delta=0、K=0、L=0、70 位回溯、不可达目标、容量及数值/尺寸溢出。测试枚举只作为 oracle，生产算法是逐层 Viterbi。

## S13：SLVA 逐条有序局部路径

固定样例调试入口为 `build/wsl-debug/bp_lcgcd --demo-slva`。参数集中在 `src/apps/main.c` 的 `run_slva_demo()`，从原始 H 系统化开始，依次调用 S10、S11、S13（初始化内部调用 S12），打印全部 16 条局部路径。完整 H、权重、syndrome、U/D、列映射和局部矩阵见[算法笔记 6.8](docs/trellis_slva_notes.md#68-固定样例参数与-main-调试入口)。这是独立的教学入口，`--decoder lcgcd` 仍返回 NOT_IMPLEMENTED。

新增 [slva.h](include/bplcgcd/slva.h) / [slva.c](src/search/slva.c)，以 S12 最佳前缀表和二叉最小堆实现串行列表。重新核对 1994 论文 §2.3–2.4，采用“最后不同分支 + 固定后缀分组”的等价 k-shortest 方法；[笔记第 6 节](docs/trellis_slva_notes.md#6-s13按最后不同分支划分的串行列表) 说明与原文最后合并、候选栈思想的对应、实现差别及完整证明。

一个堆节点表示“后缀固定、前缀自由”的一组路径，并保存最轻代表。输出代表后，按剩余路径与它的最后不同分支划成不相交子组，再用缓存 Viterbi 前缀求各组最轻代表。各组不相交、覆盖所有未输出路径；每次出堆取所有组代表中的最轻者，因此返回代价非递减且无重复，不需要保存已输出位串去重。零列连接相同状态的 bit=0/1 仍分别参与枚举。

调用方先按 S10/S11 准备不可变图，再使用以下接口；path 为调用方拥有的长度 L 缓冲，始终输出阶段顺序 `(f_I[A],f_P)`，代价仍仅为 `W(f_I[A])+W(f_P)`。

```c
SlvaGenerator list = {0};
SlvaLimits limits = {.max_candidates=4096, .max_bytes=1024*1024, .path_limit=100};
SlvaStatus status = slva_init(&list, &trellis, &limits);
/* Only enter the loop after SLVA_OK. path is already allocated by the caller. */
if (status == SLVA_OK) {
    while ((status = slva_next(&list, &path, &local_cost)) == SLVA_OK) {
        /* Consume one local path here. */
    }
    /* Inspect status: DONE, PATH_LIMIT, RESOURCE_LIMIT or an explicit error. */
}
slva_free(&list);
```

`slva_peek` 获取下一条代价但不增加 paths_returned，可能触发延迟扩展和分配；重复成功 peek 返回相同代价。第一条保持 S12 tie-break，后续等代价按候选创建序号稳定排序，不承诺与 TEP 的并列路径顺序相同。路径不打包到 64 位整数。

max_candidates 限制堆、pending 及临时孩子的同时存活数；max_bytes 限制本模块所有分配有效载荷，包含 Viterbi 表、堆容量及路径数组。path_limit=0 为不限返回数；DONE 表示耗尽，PATH_LIMIT 表示还有剩余路径但达到返回数上限。资源/分配错误不消费下一条，不丢弃合法组。图/权重和调用方输出不计入本模块预算；完整预算口径见笔记 6.6。

当前保存完整位串，初始化 O(L*2^delta)，一次扩展最坏 O(L^2+L log Q)，总空间 O((L+1)*2^delta+Q*L+Q)。采用这一写法便于核查正确性，未宣称论文优化版本的性能；不进行生产全枚举后筛选。

同一手算例的前三条：`100000`（6）、`010001`（10）、`010100`（12）。第三条来自第一条留下的候选，体现全局候选堆必须保留未输出组。

```bash
ctest --test-dir build/wsl-debug -R '^slva$' --output-on-failure -V
```

新测试通过 61,378 项检查、153 个小图的独立 Eq.9 全列表 oracle，包括零/相等权重、普通小数及大动态范围、平行边、所有 delta 边界、K=0/L=0、非零目标、peek/重放/前缀、资源限制及耗尽状态。另核对 70 位 delta=0 图的前 25 条与 TEP 代价一致。ASan/UBSan 含泄漏检查通过，执行命令记录在算法笔记 6.7。

## 检查内容与代码位置

`CMakeLists.txt` 注册十三个 CTest 测试；S13 新增 `slva` 的 61,378 项检查、153 个穷举图。原 `viterbi` 有 2461 项检查、60 个图，`trellis` 有 8022 项检查、39 个图，`local_constraints` 有 4049 项检查、82 组划分，`simulation` 有 493 项检查，`simulation_cli` 有 8 个场景，`channel` 有 105 项检查，`cli_smoke` 有 45 个场景（含固定 SLVA demo）；`gf2` 有 1,598 项检查；`systematic` 有 119,601 项检查；`decoder` 保留三种场景共用翻转代价的回归；`tep` 有 109,279 项检查、55 个场景；`gcd` 与独立全空间最小代价 oracle 核对 402 帧。

进入 S13 时的实际基线为 12/12 通过；本轮 WSL 配置、构建及 CTest 为 13/13 通过，ASan/UBSan 运行 SLVA 测试无报告。本节点仅新增 SLVA 模块/测试、注册构建并更新说明；已有用户代码和注释保持原样。

TEP 测试覆盖乱序权重、全零、全相同权重、固定种子随机二进制小数、普通小数及大动态范围，检查顺序、精确 tie-break、空模式、无重复、完整性、输出位序、输入不变和结束状态。测试 oracle 枚举全部子集后 qsort，不调用生产 FPT/最小堆代码。

S03 oracle 使用独立整数运算枚举全部 `H e^T`，由像空间大小得到秩，并逐项判断每个 syndrome 是否有解；同时核对 U 可逆、`U H Q=D`、系统形式、坐标往返与原矩阵不变。测试使用显式检查，在定义 `NDEBUG` 时也执行。

GCD oracle 在原始列坐标枚举全部 2^N 个翻转模式，直接计算 syndrome 和 W，不调用 S03、S05 或生产评分。使用可精确表示的二进制小数权重，覆盖 N=1—8 的固定种子随机矩阵、三种模式、非零 hint、零权重/ties、rank=0、K=0、冗余 Hamming [7,4]、列交换及不一致 syndrome。无预算时必须得到穷举最小值；预算截断时仍须满足原始约束，且不得标记已证明最优。

优先阅读：

1. `src/search/slva.c` 的 `expand_pending`：固定后缀分组、最佳前缀回溯、保留全局候选和失败处理。
2. `include/bplcgcd/slva.h`：逐条接口、路径位序、三类预算和输入生命周期。
3. `tests/slva_test.c` 的 `verify`：独立 Eq.9 全列表 oracle，以及稳定重放/预算/完整性证据。

本节点新增源码独立编写；参考矩阵只在忽略的 build 目录转换验证，未作为项目数据分发。搜索正确性的依据是补齐映射和非负代价下界；小码 oracle 为实现提供回归证据。本轮没有进行大码吞吐实验或内存分配失败故障注入。

## 节点状态

S00—S12 已通过节点推进验收；S13 SLVA 逐条有序局部路径已实现，等待人工验收。下一候选节点是 S14：补齐 Ac、完整 LC-GCD 候选与早停；收到明确授权后才执行。
