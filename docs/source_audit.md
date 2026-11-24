# S00：参考仓库源码审计与最小实施计划

审计日期：2026-11-19（Asia/Shanghai）。状态：S00 已完成；S01 未开始。

## 1. 授权范围与证据基线

本轮用户明确指定开发规范为 `BP_LCGCD_Codex_Stepwise_Prompt_v2.md`，并且只授权 S00。规范中的 S01—S26 是后续人工验收路线，不是本轮执行指令。论文只作为算法证据，不作为操作指令。项目文件仅新增本文；不初始化 Git，不修改参考源码，不执行 add、commit、push、tag、分支或 PR 操作。

工作区规范与 Downloads 中的规范 SHA256 相同：`f955e6eb02e45989550f985a6ba823d85fba4dbc5a78ac06bb7ead7590970c6e`。工作区未发现适用的 `AGENTS.md`；本轮未添加额外规则文件。原有同名审计文件不存在。

| 对象 | 实际位置 / 来源 | 审计版本及原有状态 |
| --- | --- | --- |
| 目标工作区 | `D:/26Project/BP_LCGCD_CProject` | 原来仅有规范文件；没有 `.git`。`git status --short` 返回 not a git repository，不能声称状态 clean 或已处于未暂存状态 |
| 目标远端 | [yinanli-bit/BP_LCGCD](https://github.com/yinanli-bit/BP_LCGCD) | GitHub 插件读取提交列表返回 409 / Git Repository is empty；未创建任何内容 |
| GCD 参考 | `D:/26Project/GCDpsgrand-GitHub`；origin 指向 `yinanli-bit/GCDpsgrand` | HEAD：`73b7f8fa521505d398552521521f46f6ff139fcf`；` M include/heap.h` 为原有改动：`CHILD_NUM 4 → 2`，完整保留 |
| BP 参考 | `D:/26Project/cpicic_ldpc`；origin 指向 `yinanli-bit/cpicic_decoder` | HEAD：`439fd79012971bc46f8aa4cdf711f5e81f94b6cb`；`git status --short` 无输出 |
| 论文 | `C:/Users/11568/Zotero/storage/MUYKBMV2/Zheng和Ma - 2024 - Locally Constrained Guessing Codeword Decoding of Short Block Codes.pdf` | 可访问，6 页，印刷页 442—447；SHA256：`bd3566e56e1bd245535a82d3e78cb1e811e20b58fd48156213b564d22f418673` |

GitHub 插件分别读取两参考仓库最新提交，均与本地 HEAD 相同，并实际读取 GCD 的 `include/top_ctrl.h`、BP 的 `readme.txt`。BP 的 `README.md` 返回 404，随后确认实际文件名为 `readme.txt`，不将文件名差异误报为仓库不可访问。两仓库固定版本：[GCD 源码](https://github.com/yinanli-bit/GCDpsgrand/tree/73b7f8fa521505d398552521521f46f6ff139fcf)、[BP 源码](https://github.com/yinanli-bit/cpicic_decoder/tree/439fd79012971bc46f8aa4cdf711f5e81f94b6cb)。

`D:/26Project/GCDpsgrand1` 是另一个目录，顶层没有 Git 元数据；未把它作为版本来源，也未把其中内容与本次基线混合。以下行号均指上述两个已识别的本地参考工作区。

## 2. GCD 参考工程：入口、结构和调用关系

```text
src/main.c:22 main
  → GetInputPara
  → InitGrandEncPara / InitGrandDecPara / InitCodecStats
  → codec_statistic.c:418 MultiThreadCodecStats
    → codec_statistic.c:165 SingleThreadCodecStats
      → RandBinSequGen → GrandEnc
      → BpskModuOverAwgnChannelUsingEbN0
      → grand_dec.c:334 GrandDec
        → 硬判决 z → SparseBinArrMulBinMat 计算 syndrome
        → syndrome=0：直接输出硬判决
        → 函数指针选择译码主循环
          → PLGCD：ParallelGCDDecMainLoop:2132
            → InitQueryPattern → 可靠度排序 / H^T 重排
            → createMinHeap → extractParMin
            → ParSgrandCalcLLR → 更新最优完整模式
            → cutHeap → createTreeNode → insertMultiple
            → BinArrPerm 回原坐标 → z XOR e
      → 与发送信息位比较，累计 BER / FER / query_times
```

注意：`include/top_ctrl.h` 当前默认 `DEC_ALG_STR="BeamPLSGRAND"`，不是 PLGCD，更不是 LC-GCD。`grand_dec.c:157—170` 才是 PLGCD / PLSGRAND / BeamPLSGRAND 的实际分派。不能直接运行默认程序后将其输出标为 GCD 基线。

| 文件 / 结构 | 实际职责及后续用途 |
| --- | --- |
| `include/top_ctrl.h` | 编译期固定 N=128、K=105、20 线程、量化开启、最多 10^9 帧等；仅参考，不迁入通用核心 |
| `include/grand_dec.h:93` `OrbDecPara` | 保存 HW/LW 限制、算法枚举与整数形式函数地址；新工程以后用类型明确的接口替代 |
| `include/grand_dec.h:102` `GrandDecPara` | `pc_mat`、`trans_pc_mat`、`inv_gen_mat` 和译码配置；矩阵是逐行分配的 `BIN **` |
| `include/grand_dec.h:112` `DecDataMem` | 帧 LLR、硬判决、syndrome、候选、置换索引、输出；每线程分开分配，但源码另有全局调试计数 |
| `include/heap.h:18` `tree_node` | `indices` 保存从 1 开始的翻转位置，`hw/lw/Depth` 描述模式，`LLR` 是部分软权重，`EP_LLR` 是补齐后权重 |
| `include/heap.h:34` `MinHeap` | 候选指针数组、容量、边界、预计算权重差；主要路径采用数组排序，不是标准二叉堆 |
| `src/query_gen.c:56` / `include/query_gen.h` | `QueryPattern` 分配并行候选数组和 `PAR_memory`；该文件另有 ORB 的 HW/LW 枚举器，不等于 GCD 的按软权重 FPT |
| `src/heap.c:1216` `ParSgrandCalcLLR` | 将翻转索引减 1，异或相应 H^T 行以补齐校验区，计算完整权重；名称含 Sgrand，但被 PLGCD 调用 |
| `lib/matrix.c:71,90,177,195` | 普通 / 按非零输入累加的 GF(2) 乘法、转置、行置换；“Sparse”函数仍存储稠密矩阵 |
| `lib/array.c:139,153` | `BinArrPerm` 是 `dst[perm[i]]=src[i]`；逆向函数是 `dst[i]=src[perm[i]]`，迁移时必须写明映射方向 |
| `lib/rand.c:232` | Box–Muller 高斯采样调用 `drand48_r`；依赖 GNU/Linux 接口，不是可移植 C11 RNG |
| `src/transmit.c:23` | BPSK：0→+1、1→−1，按 Eb/N0 或 Es/N0 加 AWGN；返回接收幅度或其量化值，未乘 `2/σ²` 转为真实 LLR |
| `src/codec_statistic.c:337—400,585—591` | 累计帧数、信息位误码和译码器返回的查询值；旧 `avg_query_times` 继承主循环计数语义 |
| `res/` | 系统矩阵、eBCH 命名矩阵、CRC8/CRC16 子目录、polar 冻结位 / 序列、randfile 和旧结果混放；不能整体作为已验证数据集导入 |

### 坐标约定必须转换

旧 PLGCD 实际采用 `H_old=[Q|I_m]`，前 K 位是被猜测的信息部分，后 m 位由 syndrome 补齐；论文及新规范采用 `H_sys=[I_m|P]`，前 m 位为基列 I，后 K 位为被猜测部分 P。旧 `info_len` 对应新 `e_P` 的长度；旧尾部校验区对应新 `e_I`。排序产生的 `ascend_perm_idx` 从排序坐标映射回旧原坐标，不能直接当作新系统化的全局列映射。

新约定下应使用：`z_i=[r_i<0]`，`s=H_sys z^T`，`e_I=s XOR P e_P^T`，`γ(e)=Σ e_i|r_i|`，输出 `v=z XOR e`。部分软权重是完整权重的下界，因为剩余权重非负。下界剪枝的正确性依赖候选覆盖完整、边界正确及已计算合法上界，不能只凭 syndrome 通过判断 ML。

`src/heap.c:1068` 的 `parity_to_systematic_full` 会修改 H 并交换列，但没有返回行变换或列映射，无法满足新规范的外部 syndrome 变换与逆映射接口；`my_build_G_from_H:1126` 也固定为 `[Q|I]`。它们只作参考，S03 独立实现。`gf2_matrix_inverse` 在头文件中有声明，本次未在该工程源码中找到定义，不能作为现成可调用功能。

## 3. 风险核验：区分事实、反例与未验证事项

| 项目 | 实际证据 | 结论与后续处理 |
| --- | --- | --- |
| `CHILD_NUM=4` 与二分支生成器 | 提交版本头文件为 4；当前用户工作区为 2。`grand_dec.c:2355—2357` 对每个 j 调 `createTreeNode`；后者 `heap.c:224,239` 仅区分零 / 非零。源码探针得到 j=1/2/3 的 mask 均为 3、权重均为 3 | 已证实提交配置下该调用方式会生成重复的位模式；当前用户的 2 分支改动不触发此四调用情况。不替用户修复或提交 |
| 有序数组和批量顺序 | `insertMultiple:564` 最后调用 `insertion_sort:958`；`extractParMin:690` 取数组前缀并移动剩余元素。批取后才生成孩子 | 有序数组本身不等于错误；但批取 / 后扩展不能保证连续输出的全局非递减顺序，见小例子结果。S05 采用单候选精确迭代器 |
| 枚举反例的实际范围 | K=1…6，各测递增权重、全 1、全 0；每例与全体非空子集的暴力排序权重比较。batch=1：18/18 通过；batch=4：11/18 通过、7/18 失败；所有例子完整覆盖且无重复 | 失败是顺序回退 / oracle 权重序列不符；不代表本次已证明完整 PLGCD 译码输出错误。探针使用 2 分支、未剪枝、没有调用完整译码器 |
| `query_times` 口径 | `grand_dec.c:2289` 每批 +1；`:2311` 另有局部 `query_nums += valid_branch`，但函数返回前者。硬判决直接输出分支也返回 1 | 已证实混合口径。S04/S06 区分候选弹出、完整评估、扩展和 batch；不能把旧图查询量当新基准 |
| 空模式 | `createMinHeap:41—44` 从单翻转开始；`grand_dec.c:2243—2250` 在译码器外层单独初始化零信息模式及补齐结果 | 空模式没有从旧生成器输出，并非整个译码器完全遗漏。新 S05 接口必须显式输出空模式并约定计数 |
| 可达范围与容量 | 旧 GCD 建堆传 `code_len=N` 而非 K，`createTreeNode:166` 主要按 max_lw / max_hw 限制；参考还有有限容量 / 查询预算路径 | 未验证默认参数下是否越过信息区、漏搜或错误剪枝；S05 必须以 K 作位置边界，并测试扩容、预算和 HW 限制；不把小探针通过外推为全部配置正确 |
| 默认信道可靠度 | `transmit.c:47—60` 返回 y 或量化 y；BP 参考则计算 `2y/σ²` | 同帧未量化且共同正比例缩放不改变最小代价解；量化尺度和跨 SNR 指标不能直接互换。S08 明确 LLR 与量化语义 |
| 矩阵名与实际码参数 | `sys_gen.mat` 为 105×128、`sys_pc.mat` 为 23×128，秩为 105/23，`GH^T=0`，H 右侧 23×23 是单位阵。`eBCH_gen.mat` 为 120×128、`eBCH_pc.mat` 为 8×128，秩为 120/8，`GH^T=0`，右侧 8×8 不是单位阵 | 当前 sys 数据满足旧坐标结构；eBCH 命名数据不是论文 [128,106]，也不能直接送入假设 `[Q|I]` 的路径。未验证 sys 数据的 CRC 多项式 / polar 构造与论文完全相同 |
| 统计与并发 | `codec_statistic.c:372` 实际停止条件涉及 `error_frames`，该量在返回查询数 >1 时累加；`grand_dec.c:21` 有全局计数，`:2482` 递增 | 旧仿真停止口径、共享状态及信息位统计需重审。未做线程竞态测试；本轮不运行大仿真，不移植整个驱动 |
| 旧文档不能充当证明 | `doc/README.md:5—19` 用任意非码字的投影碰撞说明 FER 统计选择 | 对两个有效码字 c=uG、c'=u'G，若 GG_inv=I，信息相同必有 c=c'。因此有效码字下信息帧错误与码字帧错误等价；非法输出须另报校验失败。该旧论证不能直接用作新统计定义 |

本次没有证明 GCD 的所有参数下正确，没有实现或验证 SLVA，也没有将参考工程结果视为黄金标准。后续仍必须使用独立的小码穷举 ML oracle。

## 4. BP/LDPC 来源：可学习的算法与必须隔离的适配

实际文件为仓库根目录的 `dec.c/.h`、`dec_para.c/.h`、`parity_check_matrix.c/.h`、`chnnl_modu_llr_calc.c/.h`、`main.c`、`enc.c`、`readme.txt`，均已核对存在并阅读相关路径。

```text
main.c:465 main
  → EncParaGet / DecBasicParaGet
  → GenMatFill / EncPrtChkMtrxGen / DecPrtChkMtrxGen
  → ChnnlModuParaGet / DecOtherParaGet
  → MultiThrdCodecStatistic → SnglThrdCodecStatistic:53
    → LdpcEncoding:238 → LlrCalcAwgnBpsk:104（所选信道函数）
    → dec.c:863 LdpcDecoding
      → 补 64 个 shortened 位 → 列置换 → 符号幅度格式转换
      → DecInit:798 → Check:830
      → 所选 iterate，直到校验通过或达到迭代预算
      → 逆映射并去掉 shortened 位 → decword
```

| 函数 / 数据结构 | 实际算法含义 |
| --- | --- |
| `parity_check_matrix.h:46` `Entry` | 一个非零 H 元素 / Tanner 图边；行列索引和 left/right/up/down 双向链表指针，保存 `v2c`、`c2v`、低位量化消息及 BP 的 `logVal` |
| `parity_check_matrix.h:63` `PrtChkMtrx` | 行、列哨兵和连续分配的 Entry 块；后续 S22 可借鉴“同一条边被行列两种邻接访问”，优先独立实现紧凑边索引 |
| `dec_para.h:77` `DecPara` | 迭代上限、调度、CN 规则和 refinement 函数选择，以及码专用交织 / shortening 表；不能原样作为通用配置 |
| `dec.c:324` `CheckNodeProcBp` | SPA 的 φ 域实现：`φ(x)=log((1+exp(-x))/(1-exp(-x)))`，先累加 φ(|v2c|)，再排除自身并施加其他边符号乘积；存在上限保护，未做零消息 / 极值数值测试 |
| `dec.c:122` `CheckNodeProcMinSum` | 每个 CN 求绝对值最小 / 次小值与符号异或；排除自身，用 refinement 函数修正幅值。`Min1Min2Scaled:98` 对应归一化 MS；另有 offset / two-offset |
| `dec.c:498` `LayeredIterate` | 每层 `LayeredVariNodeProc:382`（APP−旧 c2v）→ CN 更新 → `LayeredAppUpdate:412`（v2c+新 c2v）；层间立即使用更新后的 APP |
| `dec.c:761` `FloodingIterate` | 先全体 CN，再每个 VN 计算 `APP=chLLR+Σc2v`、`v2c=APP−该边c2v`，最后硬判决 |
| `dec.c:863` `LdpcDecoding` | 外层包含大量码适配，不是现成通用 BP 接口；`Check` 只检验零 syndrome，不支持直接指定非零 target syndrome |
| `dec_para.c:38,190,238` | 通过参数选择 CN / refinement / schedule。当前量化开启时明确拒绝 ORIGIN_BP；仅看到函数存在不能声称默认可运行浮点 BP |
| `chnnl_modu_llr_calc.c:104` | BPSK AWGN 明确输出 `LLR=2y/σ²`；与旧 GCD 信道返回值的语义不同 |
| `enc.c:238` | 系统编码加 `itlvCodeword[colItlvPattern[j]]=codeword[j]`；列映射方向需与接收端联合核对 |

浮点 BP 的核心消息关系为 `q(v→c)=r_v+Σ_{c'≠c}m(c'→v)`；SPA 的 CN 消息为 `2atanh(Π_{v'≠v}tanh(q(v'→c)/2))`，NMS 则用其他边符号乘积乘 `α·min|q|`。参考工程的 φ 实现、两最小值实现分别对应这两类规则；本轮只是核对源码含义，没有完成 BP 数值正确性验收。

固定适配证据：`dec_para.h:98—101` 的 64/1984/2048 数组，`dec_para.c:263` 的 `codeLen+64`，`dec.c:887—932` 的 shortened 位补齐和多重映射。`readme.txt` 说明 `H_dec_2048.dat`、`H_decoder.dat`、`col_swap.txt`、`col_itlv_permuLay5.txt`、`de_shorten_permuLay5.txt`、`shorten_idx_hPermuLay5.txt` 的不同坐标含义。这些与特定 802.3 码、缩短及层排列相关，S21 后放入独立适配层，不嵌入通用 BP 核心。`ShuffledIterate` 的主体被注释，也不能把所有枚举名都当作已完成功能。

## 5. 从原始 PDF 核对的 LC-GCD 约束

已提取全部 6 页文本，并渲染查看 PDF 第 2、3、4 页（印刷页 443、444、445），消除双栏提取和上标丢失歧义。论文 DOI 为 `10.1109/ITW61385.2024.10806923`。以下为原文核对后的释义，不是算法已实现的声明。

| 论文位置 | 核对结果 / 实施约束 |
| --- | --- |
| II-B，印刷页 443，Eq. 8 | 固定系统化 `H=[I\|P]`；无需逐帧 GE。猜测 `e_P` 并由 syndrome 确定 `e_I` |
| II-C，印刷页 443，Eq. 9/10 | A 是固定基列集合 I 内 `\|r_I\|` 最大的 δ 个位置，Ac=I\A；不是从全部 N 位或被猜测的 K 位中选。两式为 `e_I[A] XOR P_A e_P^T=s[A]` 与剩余行的对应式 |
| II-C，印刷页 444 | `[I_δ\|P_A]` 定义长度 $K+δ$ 的 trellis，至多 $2^δ$ 状态；SLVA 按 `γ(e_I[A])+γ(e_P)` 非递减输出合法局部路径。单次 Viterbi 仅解决第 1 条路径 |
| Algorithm 1，第 2 / 5 / 10 行，印刷页 444 | 原文初始化 `e_opt=z`，因为零码字使它合法；部分代价 **严格大于** 当前完整最优权重时停止，完整代价 **严格更小** 才更新。因此相等代价继续搜索并保留先到最优者 |
| 初始化与 δ=0 的实现约定 | 规范 GCD 允许先用 `(s,0)` 作为合法候选；它与论文的 `e_opt=z` 都给出上界，但有限预算下轨迹可能不同。S06/S14/S15 必须统一初始化、空模式、tie-break 和计数后再声称 δ=0 严格等价；指定 syndrome 模式不能直接照搬 `e_opt=z` |
| II-C Remarks，印刷页 444 | δ=0 退化为 GCD；小码 δ=m 时，一条 Viterbi 最优路径已足以实现 ML，因为没有剩余约束要补齐 |
| II-D，Eq. 11—13，印刷页 444 | 从 P 的 K 个位置中选择 λ 个最不可靠位置 B，枚举其 2^λ 种赋值，在残余 syndrome 的 2^λ 个分支上分别搜索。分支下界条件写作 ≥；不同于 Algorithm 1 的 >，不得混抄。此次暂缓并行 LC-GCD，旧 Parallel=4 不代表该方法 |
| III，Eq. 14，印刷页 444 | 原文 GCD 平均复杂度：`O(K log K)+O(ℓ_avg log ℓ_avg)+O(ℓ_avg·m)`；其中重编码项利用相邻候选的增量更新，不是朴素每次矩阵乘法的无条件复杂度 |
| III，Eq. 15，印刷页 445 | LC-GCD：`O(2^δ(K+δ)+ℓ_avg(K+δ))+O(ℓ_avg K(m−δ))`。原文给出首条路径 O(2^δ(K+δ))、后继路径 O(K+δ) 的分析；未来若独立 k-shortest 实现有队列开销，应报告实际复杂度，不直接套用原文结论 |
| IV，印刷页 445 | BPSK/AWGN；CA-polar [128,105]，CRC-11 多项式 `x^11+x^10+x^9+x^5+1`；eBCH [128,106]；图 1 的 δ=2/4/6/8/10，ℓ_max=10^6；后续比较 LC-GCD δ=6、GCD 预算 10^6，LC-OSD δ=8、预算 16384，PC-GRAND δ=6、预算 10^6 |

原文没有提供可直接照搬的完整 SLVA 程序、随机种子、每点帧数 / 错误帧停止阈值或完整吞吐测试环境说明。S00 未运行论文仿真，未复现曲线。后续先验证小码最优性，再研究统计复现；软件运行时间不等于 RTL/ASIC 延迟。

## 6. 许可证与复用决策

两参考仓库的已跟踪文件中未发现明确的 LICENSE/COPYING。GCD 包含 Jiayu Xiao 等作者信息，`lib/rand.c`、`src/quantization.c` 等另有 ICAIS Lab at NJU 的 All rights reserved 声明；BP 相关文件普遍有该声明。公开可读取及用户拥有仓库访问权限，均不能替代明确的第三方代码授权。数据矩阵也未找到独立许可说明。

| 决策 | 对象 | 理由及首次使用节点 |
| --- | --- | --- |
| 原样复用 | 暂无生产源码或第三方矩阵 | 本轮未复制到项目；许可未确认，正确性也需独立验收 |
| 改造思路、独立实现 | GF(2) 运算 / 位序映射：`matrix.c`、`array.c`；系统化：`heap.c:1068` | S02/S03 增加输入检查、行列变换、冗余行和 syndrome 一致性证据 |
| 改造思路、独立实现 | FPT 两种子节点操作、合法上界剪枝：`heap.c:162`、`grand_dec.c:2132` | S05/S06 以 K 为边界、单候选有序迭代、确定性 ties、明确截断；不移植批量队列和全局宏 |
| 仅参考 | `query_gen.c` 的 ORB HW/LW 枚举、旧统计驱动、Makefile、top_ctrl、历史结果 | 排序目标 / 计数 / 构建及停止条件与新要求不同；旧 Makefile 还含执行 Git 提交的目标，本次未调用任何 make 目标 |
| 改造思路、独立实现 | BPSK、Box–Muller、信道公式：`transmit.c`、`rand.c`、`LlrCalcAwgnBpsk` | S08 引入可复现 RNG 与明确真实 LLR；避免 GNU 专用依赖和未确认许可的源码复制 |
| 仅参考，后续专项审计 | BP 的 CN/VN 关系、邻接访问、Layered/Flooding | S21 重新核对许可和数值保护；S22/S23 才按用户选择实现图与单一规则 |
| 按论文独立实现 | A 选择、trellis、Viterbi、SLVA、LC-GCD | S10—S15 分步构建；Viterbi-1 不冒充 SLVA，禁止全枚举后筛选冒充 trellis 搜索 |

## 7. 最小实施顺序（仅计划，不执行）

**下一次只有收到“开始 S01”才动工。** S01 计划新增下列 5 个短文件，保留规范和本文；不新增算法头文件或空模块：

| 文件 | 具体用途 / 调用者与验证 |
| --- | --- |
| `CMakeLists.txt` | 只配置 C11 CLI 可执行文件、Debug/Release 和 CTest；CTest 调用最小检查脚本 |
| `.gitignore` | 忽略构建目录与仿真结果；不改 Git 配置，不初始化或创建分支 |
| `README.md` | 记录实际验证成功的构建 / --help 命令、当前功能范围和逐节点停止机制 |
| `src/apps/main.c` | 唯一生产入口；提供 --help、退出码和未实现 decoder 的明确错误，不解码、不仿真 |
| `tests/cli_smoke.cmake` | 调用 CLI 检查帮助正常退出、请求未实现译码器明确失败；通过 CTest 执行 |

S01 验收：在独立 build 目录配置和构建、运行 --help、CTest 全部通过。当前 GCC 可用，CMake/Clang 未在 PATH 中找到；S01 开始时应先定位可用工具，若工具缺失，明确报告环境前提，不能伪造验证结果。本轮不安装构建工具。

| 后续节点 | 最小目标及验收关口 |
| --- | --- |
| S02 → S03 → S04 | 二元矩阵 / syndrome → 离线系统化和双向映射 → classical / target-syndrome 输入输出与所有权、状态、计数约定 |
| S05 → S06 → S07 | 有序 TEP 迭代器与暴力排序对照 → 完整 GCD 与独立 ML oracle → 可手算的单帧 CLI/trace |
| S08 → S09 | 固定 seed 的 BPSK/AWGN/LLR → 有帧数上限的最小 Monte Carlo；少量帧不能称论文曲线复现 |
| S10 → S11 → S12 | 选 A / Eq. 9/10 → 状态和分支 → 单条 Viterbi 与合法局部路径穷举对照 |
| S13 → S14 → S15 | 精确有序多路径 → 完整 LC-GCD 与早停 → δ=0、δ=m、ties、预算、映射及非法 syndrome 对照；S13 如需拆分先获批准 |
| S16 → S17 → S18 | 同帧公平对比 → 可追溯统计字段 → 受限 SNR/δ 扫描及绘图 |
| S19 → S20 | 有需求才加调度无关的帧 RNG / 多线程 → sanitizer、Debug/Release、回归及文档 |
| S21 → S22 → S23 → S24 | BP 专项审计 → 稀疏图 → 一种 BP/NMS 规则 → 另一规则或调度二选一 |
| S25 → S26 | 经用户确定策略的 BP-LC-GCD 融合 → QLDPC/CSS、逻辑错误与硬件验证边界；syndrome 通过不能声称量子纠错成功 |

## 8. 实际验证与可复查证据

验证仅有短探针、文件读取和语法检查，没有执行参考工程主程序。临时工具、PNG、探针与日志位于 `%TEMP%/bplcgcd_s00_audit`（PDF 解析依赖位于 `%TEMP%/bplcgcd_s00_pdf_tools`），不写入参考仓库或目标源码。临时目录不是项目交付依赖；下文保留探针源码及重跑命令。

| 检查 | 实际结果 |
| --- | --- |
| 必审文件存在性 | PASS：GCD 12 个文件及 res，BP 11 个文件均存在 |
| 两参考 HEAD / 工作区状态 | PASS：已用 `git rev-parse HEAD`、`git status --short` 核对；目标工作区 Git 状态为 N/A，不是 clean |
| 规范双份 SHA256 | PASS：一致；原文件未改 |
| PDF 提取 / 渲染核对 | PASS：6 页可解析；关键 3 页已视觉核对；未输出改写后的论文 |
| 矩阵检查 | PASS：独立 Python 整数位运算核对矩阵维度、秩与 GH^T；发现 eBCH 数据与论文不匹配，不视为论文数据通过 |
| GCD 部分模块 + 探针编译 | PASS：GCC 16.1.0 / MSYS2 UCRT64；最终链接下面列出的 8 个实际参考源文件，exit=0 |
| FPT 小例子 | batch=1：18/18 PASS；batch=4：7 个顺序反例，11/18 PASS；四调用重复模式证据已列出。程序 exit=0 仅表示探针正常完成且单候选用例通过，不表示 batch 测试全通过 |
| GCD 全源码语法检查 | FAIL / exit=1：Windows UCRT64 缺 `drand48_data`、`srand48_r/drand48_r/lrand48_r/nrand48`。首次命令还漏了 Makefile 要求的 RAND_FILE 宏；补齐后宏错误消失，GNU 接口错误仍在 |
| BP 六个关键模块语法检查 | FAIL / exit=1：同类 GNU 随机数接口与不完整类型错误；没有完成 BP 全工程链接和数值测试 |
| 首次仅链接 heap.c 的探针 | FAIL：未解析 `CalcQuerySyndromes`；后续补齐实际依赖源文件后 PASS，没有改源码或伪造该函数 |
| 端到端 ML / SLVA / BP 仿真 | 未运行；属于后续节点，不能由本轮结果推定正确 |

实际执行的只读命令（分别在两个参考根目录运行）：

```powershell
git rev-parse --show-toplevel
git rev-parse HEAD
git status --short
git remote -v
git diff -- include/heap.h
git ls-files '*LICENSE*' '*COPYING*' '*license*' '*License*'
```

以下是最终探针编译的等价 PowerShell 命令；`$probeDir` 指向临时目录，将下附 C 源保存为其中的 `heap_probe.c` 后即可重跑。没有执行 Makefile，避免其格式化、清理或 Git 目标产生副作用。

```powershell
$refRoot = 'D:\26Project\GCDpsgrand-GitHub'
$probeDir = Join-Path $env:TEMP 'bplcgcd_s00_audit'
$probeSources = @('src/heap.c','src/grand_dec.c','src/query_gen.c','src/quantization.c',
                 'lib/array.c','lib/matrix.c','lib/my_math.c','lib/fileio.c') |
    ForEach-Object { Join-Path $refRoot $_ }
gcc -std=gnu11 -O0 "-I$refRoot/include" "$probeDir/heap_probe.c" @probeSources -lm -o "$probeDir/heap_probe.exe"
& "$probeDir/heap_probe.exe"

$gcdSources = @(Get-ChildItem "$refRoot/src/*.c","$refRoot/lib/*.c" | ForEach-Object FullName)
gcc -std=gnu11 -fsyntax-only "-I$refRoot/include" -fopenmp '-DRAND_FILE="D:/26Project/GCDpsgrand-GitHub/res/randfile"' @gcdSources
$bpRoot = 'D:\26Project\cpicic_ldpc'
$bpSources = @('dec.c','dec_para.c','parity_check_matrix.c','chnnl_modu_llr_calc.c','main.c','enc.c') |
    ForEach-Object { Join-Path $bpRoot $_ }
gcc -std=gnu11 -fsyntax-only -fopenmp @bpSources
```

本轮实际通过 Python subprocess 参数数组启动上述 GCC 命令，避免 shell 引号影响宏；记录的日志为 `compile.log`、`heap_probe.log`、`gcd_syntax.log`、`bp_syntax.log`。原探针 SHA256：`4c4cd532a581b0465deda07ac017f593159ca90c048a3409858dc5a30ac35536`（下附内容在排版空白上可能不同，逻辑相同）。

探针源码：仅调用实际参考函数，不实现新的生产枚举器；全子集 oracle 独立按 bitmask 计算代价。

```c
#include "heap.h"
static unsigned mask_of(const tree_node *n) {
    unsigned mask=0;
    for(int j=0;j<n->hw;j++) mask |= 1u<<(n->indices[j]-1);
    return mask;
}
static int cmp(const void *a,const void *b) {
    double x=*(const double*)a,y=*(const double*)b;
    return (x>y)-(x<y);
}
static int run_case(int k,int batch,int mode) {
    double w[8],oracle[128],last=-1; unsigned seen[128]={0};
    int total=0,duplicates=0,order_bad=0,cost_bad=0,missing=0;
    for(int i=0;i<k;i++) w[i]=mode==0?i+1:mode==1?1:0;
    int size=(1<<k)-1;
    for(int m=1;m<=size;m++) {
        double sum=0;
        for(int i=0;i<k;i++) if(m&(1<<i)) sum+=w[i];
        oracle[m-1]=sum;
    }
    qsort(oracle,size,sizeof(double),cmp);
    MinHeap *heap=createMinHeap(w,k,k,512,k,4096);
    while(getSize(heap)>0 && total<1024) {
        tree_node *popped[4]={0},*children[8]={0};
        int count=extractParMin(popped,batch,heap),nc=0;
        for(int i=0;i<count;i++) {
            tree_node *node=popped[i]; unsigned m=mask_of(node);
            if(seen[m]++) duplicates++;
            if(node->LLR<last) order_bad++;
            if(total>=size || node->LLR!=oracle[total]) cost_bad++;
            last=node->LLR; total++;
            for(int j=0;j<2;j++) {
                tree_node *child=createTreeNode(j,node,heap);
                if(child) children[nc++]=child;
            }
            destroyNode(node);
        }
        insertMultiple(&heap,children,nc);
    }
    for(int m=1;m<=size;m++) if(!seen[m]) missing++;
    destroyHeap(heap);
    int ok=(total==size && !duplicates && !order_bad && !cost_bad && !missing);
    printf("k=%d batch=%d weights=%s total=%d/%d duplicates=%d order_bad=%d oracle_mismatch=%d missing=%d %s\n",
        k,batch,mode==0?"increasing":mode==1?"ones":"zeros",total,size,
        duplicates,order_bad,cost_bad,missing,ok?"PASS":"FAIL");
    return ok;
}
int main(void) {
    double w[3]={1,2,3}; MinHeap *heap=createMinHeap(w,3,3,16,3,4096);
    printf("child probe from root [1]:\n");
    for(int j=0;j<4;j++) {
        tree_node *child=createTreeNode(j,heap->nodes[0],heap);
        printf("j=%d mask=%u weight=%.0f\n",j,mask_of(child),child->LLR);
        destroyNode(child);
    }
    destroyHeap(heap);
    int sequential_ok=0,batch_ok=0;
    for(int mode=0;mode<3;mode++) for(int k=1;k<=6;k++) {
        sequential_ok+=run_case(k,1,mode);
        batch_ok+=run_case(k,4,mode);
    }
    printf("SUMMARY sequential=%d/18 batch=%d/18\n",sequential_ok,batch_ok);
    printf("ABI sizeof(long)=%zu sizeof(ULINT64)=%zu sizeof(UINTPTR)=%zu\n",
           sizeof(long),sizeof(ULINT64),sizeof(UINTPTR));
    return sequential_ok==18?0:1;
}
```

本机 ABI 探针另显示 `sizeof(long)=4`、旧 `ULINT64=4`、`UINTPTR=8`：旧类型名不保证在 Windows 上真为 64 位。新工程应采用 `<stdint.h>` / `size_t` 和匹配的格式宏，避免移植旧类型定义。

## 9. 人工审查位置与停止点

优先阅读以下三组；先理解数据流，不需要先读全部旧项目：

1. `D:/26Project/GCDpsgrand-GitHub/src/grand_dec.c:2132`，重点 2243—2250、2289—2311、2355—2357：追踪零信息模式、完整权重补齐、batch 计数和孩子生成；思考旧 `[Q|I]` 如何转为新 `[I|P]`。
2. `D:/26Project/GCDpsgrand-GitHub/src/heap.c:162`，结合 564、690：看 1 起始索引、两种孩子、批量提取与扩展顺序；思考“当前队列前四个”为什么不一定是全局接下来的四个。
3. `D:/26Project/cpicic_ldpc/dec.c:863`，结合 122、324、498、761：把 shortened 位 / 交织适配与 CN/VN 数学更新区分开，理解 APP、v2c、c2v 的依赖关系。

最终检查：两参考工作区状态保持原样；本文新增前对 GCD 的 100 个已跟踪现存文件、BP 的 49 个已跟踪现存文件和目标原规范计算 SHA256，结束复核全部一致。两参考仓库 `git diff --check` 均 exit=0；GCD 另有既存文件的 LF/CRLF 提示。目标目录没有 Git，不能提供合法的 `git diff --check` 或未暂存状态；用文件清单与散列验证改动范围。本轮唯一项目产物是本文，未执行 add/commit/push，未开始 S01。

已停止，等待用户审阅并明确发出“开始 S01”。
