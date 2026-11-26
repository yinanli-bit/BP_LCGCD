# 第三方来源与改写记录

## TASK1：cpicic_decoder

- 仓库：[yinanli-bit/cpicic_decoder](https://github.com/yinanli-bit/cpicic_decoder)。本轮实际读取并固定在 commit `439fd79012971bc46f8aa4cdf711f5e81f94b6cb`。
- 原版权：`Copyright (C) ICAIS Lab at NJU. 2020 - 2021. All rights reserved`。
- 原作者：Yangcan Zhou，`yczhou@smail.nju.edu.cn`。
- 原文件声明日期：2020-12-10。仓库没有单独的 LICENSE；保留源文件版权/作者，并在改写位置标明 adapted 和本记录的路径。

| 实际参考源码 | 本项目位置 | 保留的设计与改动 |
| --- | --- | --- |
| [parity_check_matrix.h，Entry / PrtChkMtrx](https://github.com/yinanli-bit/cpicic_decoder/blob/439fd79012971bc46f8aa4cdf711f5e81f94b6cb/parity_check_matrix.h) | `include/bplcgcd/quantum.h` | 保留 left/right/up/down、row/column headers 和连续 edge block；改为 size_t 坐标、Pauli label、稳定 edge_index；移除拓扑内的动态消息 |
| [parity_check_matrix.c，PrtChkMtrxAllocate / EntryInsert](https://github.com/yinanli-bit/cpicic_decoder/blob/439fd79012971bc46f8aa4cdf711f5e81f94b6cb/parity_check_matrix.c) | `src/core/quantum.c` 的 `build_graph` | 保留 circular header 与双向插入；由已排序 Hx/Hz 扫描直接尾插，省去任意位置查找和重复 entry 处理；增加 Y 单边标签 |
| [dec.c，CheckNodeProcMinSum / Min1Min2Scaled](https://github.com/yinanli-bit/cpicic_decoder/blob/439fd79012971bc46f8aa4cdf711f5e81f94b6cb/dec.c) | `src/decoders/qbp.c` 的 `check_update` | 保留 sign parity、min1/min2、minimum edge 和两次 row traversal；增加 measured syndrome sign，用 quantum commutation LLR 替代 classical V2C；消息位于 workspace，直接乘 nms_scale |
| 同上，FloodingIterate / DecInit | `qbp_decode`、`posterior_update`、`variable_update` | 参考先完成所有 checks 再按列更新的组织；所有概率语义改用论文的 quantum XYZ 公式，未移植 classical 二元 VN 公式 |

未迁入固定码长、交织、shortening、整数消息、quantization 或 layered 分支。新增概率转换、Pauli 采样、commutation 验证与测试独立编写；SplitMix64 使用当前项目已有的固定整数混合公式，量子样本另外加入 domain 常量。

## 论文

实际资料为本地 `docs/` 的三篇 PDF，逐项公式到函数映射见 [TASK1 记录](quantum_task1.md)。2024 论文负责 stabilizer/symplectic 与 soft output 语义；2025 论文负责完整 XYZ quantum NMS 公式；2026 论文负责目标 single-scale 前端架构及直接乘法参数含义。没有从 classical 实现推导或替代 quantum 消息公式。

原 GCD/LC-GCD 的参考审计继续见 [source_audit.md](source_audit.md) 和 [trellis_slva_notes.md](trellis_slva_notes.md)，本轮未修改它们及对应生产代码。
