# Trellis 与 List Viterbi：原文核对及本项目约定

本记录覆盖 S11–S17。第 1–4 节保留原文对应关系，第 5–6 节记录 Viterbi/SLVA，第 7 节记录完整 LC-GCD，第 8 节记录系统验证和精度边界，第 9 节记录同帧比较，第 10 节记录正式统计口径。

## 1. 实际阅读的来源

- [Seshadri 与 Sundberg，1994 原文](<Seshadri和Sundberg - 1994 - List Viterbi decoding algorithms with applications.pdf>)：11 页；本轮提取全文，并视觉核对 PDF 第 2–4 页的公式和图。重点是 §2.1–2.4（印刷页 314–317），不是应用部分的卷积码参数。
- 文件 SHA256：`930ec4279283ea08b6b1f7633400e513e0cf2d1fdff6665c9986825bace9ec75`。
- [Zheng 与 Ma，2024 原文](<Zheng和Ma - 2024 - Locally Constrained Guessing Codeword Decoding of Short Block Codes.pdf>)：II-C、Eq. 9/10 及 Algorithm 1 指定本项目局部约束、部分软权重和完整候选的停止条件。

| 1994 原文位置 | 核对内容 | 本项目对应 |
| --- | --- | --- |
| §2.1，Eq. 1–4，图 2 | 固定起终点、分支增量度量、最优前缀度量和前驱回溯 | S11 图模型；S12 保存逐层度量及前驱 |
| §2.2，Eq. 5–8，图 3 | 每个状态保留 L 个前缀及其前驱排名的 parallel LVA | 未采用；S13 按需输出，不预先为每个状态保留整个列表 |
| §2.3，Eq. 9、图 4 | 第二条路径的偏离、最后合并与已有局部最优前缀 | S13 生成替代路径的依据 |
| §2.3，图 5、Serial Algorithm | 第三条及后续路径、合并计数、候选栈复用 | S13 需保留候选，不重复寻找已输出的前缀排名 |
| §2.4 | 保存中间度量/路径信息的计算与存储代价 | 后续按实际数据结构核算资源，不直接照抄原文复杂度 |

## 2. 本项目的 syndrome Trellis

S10 已得到 `f_I[A] XOR P_A*f_P^T=s_A`。定义 `L=delta+K`、`M_local=[I_delta|P_A]`，第 t 个局部比特记 b_t。位序固定为 `(f_I[A],f_P)`，不进行另一次可靠度排序。

状态 q_t 是前 t 列已经累计出的 delta 位 syndrome。状态位 i 对应 A 的第 i 个局部约束，而不是原始 H 的第 i 行。列掩码 h_t 将矩阵第 t 列按相同位序打包：

```text
q_0 = 0
q_(t+1) = q_t XOR (b_t ? h_t : 0)
c_t(b_t) = b_t * flip_weights[original_columns[t]]
W_local = sum_t c_t = W(f_I[A]) + W(f_P)
接受条件：走完恰好 L 条边，q_L = pack(s_A)
```

前 delta 个 original_columns 来自 S10 的 basis_order，后 K 个来自 S03 的 sys_to_orig[m+j]。所有权重仍为原始坐标的翻转代价，排除 Ac；Ac 的代价留到 S14 完成完整候选时加入。

S11 表示完整的每层 `2^delta` 状态字母表，按需计算两条带 0/1 标签的出边，不分配完整状态/边表，不做可达性裁剪。第 0 层只有状态 0 是有效起点；后续 DP 必须把其他起点及不可达状态的度量初始化为无穷，不能因为 branch 接口可计算一条边就认定该前缀可达。

1994 的固定起终点常写作状态 1 到状态 1；本项目为 0 到 pack(s_A)，这是层内状态编号/终点的适配，不要求把非零目标 syndrome 改成零。图有 L+1 层；状态上限由调用方给定，构图前检查 `2^delta` 和机器位宽。S12/S13 另须检查实际 `(L+1)*states` 度量/回溯表及候选缓存的大小，S11 的状态上限不代表这些表已经分配。

## 3. Serial LVA 的递推含义与适配

§2.3 假设先由 VA 得到全局最佳路径，并保留每层各状态的局部最优前缀、其度量及前驱。第二条全局路径可按与最佳路径最后合并的位置分组，比较两类来源：较早已经合并的次优前缀继续走最佳分支；在当前层首次完成最后合并的替代入边。

设最佳路径在层 t 的状态为 q_t、选中的入边为 e*_t，d 为已知最佳前缀度量，g 为沿该路径计算的第二前缀度量。将原文 Eq. 9 改写成带边标签的含义为：

```text
g_t = min(g_(t-1) + c(e*_t),
          min_{e=(u->q_t), e!=e*_t} [d_(t-1)(u) + c(e)])
```

这只是第二条候选的两类来源，不是完整 k-best 实现。需要记录最后合并层及替代分支，才能由局部前缀和已知最佳后缀回溯出路径。

图 5 将第三条候选分为：其他合并位置仍有未使用的局部第二前缀；第二条路径所用合并位置还需要该前缀的下一排名。Serial Algorithm 用 C_j 追踪各合并位置已用排名，查找下一未用前缀并加上最佳后缀度量，比较后更新选中位置的计数。文末候选栈方案保留此前仍有竞争力的候选，并与针对上一输出路径新产生的替代候选比较；“栈”不能按普通 LIFO 出栈来代替代价排序。

本项目有三个必须显式处理的差异：

1. **平行边。** 当列掩码为 0，b=0/1 有相同前驱和后继；delta=0 时所有边都有此性质。原文 Eq. 9 排除最佳前驱状态的写法不能直接照搬，必须排除已用分支而保留同一前驱的另一比特分支。前驱记录需包含比特；候选身份也不能只由状态序列确定。
2. **相同代价。** 零权重及 ties 下，状态、度量都相同的不同翻转序列仍是不同路径。S12/S13 必须约定稳定次序，并按带标签路径区分最后不同分支/合并位置及重复候选；论文的严格优劣示意不能代替此边界验证。
3. **外层停止。** 1994 应用可在外码判定正确后停止。LC-GCD 中满足局部约束只是获得一条局部候选，仍需补齐 Ac、比较完整 W，并用部分代价下界证明停止；不能在第一条局部合法路径停止。

后续优先沿 §2.3 的最后合并与候选复用思路实现；若采用等价有序 k-shortest 实现，需要说明正确性、区别和实际复杂度。不能用全 2^K 枚举后筛选替代 SLVA，也不能只对首个 Viterbi 回溯结果改一位就宣称得到了第二、第三条全局路径。

## 4. 本轮已验证与后续验收

S11 对小图穷举带标签的位路径，仅作为测试 oracle：终点接受当且仅当 Eq. 9 成立，累积代价等于直接翻转权重和，每个 f_P 恰有一个 f_I[A]，因此合法路径数为 2^K。覆盖 delta=0、K=0、L=0、零列/零权重、非零目标、原始列映射及资源上限。

S12 已完成第一条路径与独立最小代价 oracle 的核对，详见第 5 节。S13 已完成列表的有序、唯一、完整和 ties 验证，详见第 6 节；S14 的完整候选和早停见第 7 节。S12 的单路径结果不能代替列表或完整译码验证。

## 5. S12：单条最优 Viterbi 的递推与回溯

本轮重新读取 1994 原文 §2.1 并视觉核对 PDF 第 2 页（印刷页 314）的 Eq. 1–4。实现仍使用 S11 的二元 XOR 图，不照搬全连接图的所有前驱循环。新增 `viterbi_decode` 只返回一条最轻局部路径，全部前缀表保留给调用方检查及后续 SLVA 使用。

### 5.1 每格存什么

```c
typedef struct {
    double metric;
    size_t previous_state;
    uint8_t bit;
} ViterbiCell;
```

令 S=state_count、L=trellis.length。`cells[t*S+q]` 对应第 t 层、状态 q，存储从初态到该点的最小前缀代价 phi_t(q)，以及实现该最小值的前驱状态和分支比特。bit 必须独立保存：零掩码时，两条不同分支的前驱和后继都相同。

输出 path 的长度为 L，位序是 `(f_I[A],f_P)`，通过 original_columns 对应原始位置；它不是 N 位完整翻转或最终错误估计。输入图/权重不变，工作表、path 和 cost 均由调用方提供；函数不分配内存。

### 5.2 初始化与递推，对应 Eq. 1/2

```text
phi_0(0) = 0
phi_0(q != 0) = +infinity

给定第 t 列掩码 h_t 和下一层目标状态 q：
  bit=0 的前驱 p0=q
  bit=1 的前驱 p1=q XOR h_t

phi_(t+1)(q) = min(phi_t(p0) + c_t(0),
                   phi_t(p1) + c_t(1))
```

代码按目标状态收集两条入边，并复用 `syndrome_trellis_branch` 计算分支代价。不可达前缀直接跳过；有限候选仅在严格更小时覆盖当前格，并同步记录前驱和 bit。没有额外的最终错误评分，代价始终为 W(f_I[A])+W(f_P)，不计 Ac。

每格固定先比较 bit=0，再比较 bit=1；代价相等保留先到者，不使用 epsilon。即使掩码和边权都是零，仍比较两个不同标签，只在该格保留一个幸存前缀，而没有从图中删除另一条边。测试以整数和可精确表示的二进制小数权重，按最后一位优先为 0 的位串序独立核对 ties；实现不把路径打包到 64 位整数，已验证 70 位路径。

最优性依据是归纳不变量：每个有限格保存到该状态的最轻前缀。任一路径的最后分支必是上面两条之一，其更早前缀可用对应状态的最佳前缀代替；逐层取小即可得到目标状态的最轻路径。这证明的是给定局部图的单路径最优性，不是完整 LC-GCD，也不是有序 k-best 列表。

### 5.3 终止与回溯，对应 Eq. 3/4

```text
局部最小代价 = phi_L(target)
q = target
for t = L, L-1, ..., 1:
    path[t-1] = cells[t*S+q].bit
    q = cells[t*S+q].previous_state
```

只从指定 target 回溯，不在最后一层任意选最小状态。初始/不可达格的前驱标为 SIZE_MAX；成功回溯的链由递推不变量保证到达初态 0。目标不可达返回 TRELLIS_NO_PATH；正常 `[I_delta|P_A]` 图因含恒等块总有合法路径，但仍显式检查终点可达性。仅在 TRELLIS_OK 后读取 path/cost。

### 5.4 逐层手算

采用 S11 的 delta=2 例子：阶段原始列 `(2,1,0,3,4,5)`，掩码 `(1,2,0,3,2,3)`，阶段权重 `(6,3,1000,9,8,7)`，target=1。

| 层 t | phi(q=0) | phi(q=1) | phi(q=2) | phi(q=3) |
| --- | --- | --- | --- | --- |
| 0 | 0 | inf | inf | inf |
| 1 | 0 | 6 | inf | inf |
| 2 | 0 | 6 | 3 | 9 |
| 3 | 0 | 6 | 3 | 9 |
| 4 | 0 | 6 | 3 | 9 |
| 5 | 0 | 6 | 3 | 9 |
| 6 | 0 | 6 | 3 | 7 |

虽然最后一层状态 0 的代价为 0，目标是状态 1，所以答案为 6。回溯状态为 `q6=1,q5=1,q4=1,q3=1,q2=1,q1=1,q0=0`，恢复 `path=100000`，即 `f_I[A]=10、f_P=0000`；原始非零翻转位置为列 2。S11 展示的代价 12 只是另一条合法路径，不是最优路径。

第 4 层、状态 3 有代价同为 9 的两种入边；本实现保留 bit=0，体现了 ties 规则。不可达格也保留在表中，例如第 1 层的状态 2/3 均为 inf。

### 5.5 资源预算与实际验证

工作表需要 `(L+1)*S` 个 ViterbiCell，实际字节数为该值乘 `sizeof(ViterbiCell)`。调用方传入已分配数组及 cell_capacity，容量同时作为该次调用的工作表预算。函数在写入前检查 L+1、格数乘法、字节数乘法及容量，超限返回 TRELLIS_RESOURCE_LIMIT；不靠 malloc 失败控制规模，也不引入隐藏缓存。path 与图数组是调用方另外拥有的存储。

时间为 O(L*S)，回溯 O(L)，完整前缀表空间 O((L+1)*S)。例如本例需要 28 格；当前 WSL GCC ABI 下每格 24 字节，工作表共 672 字节。S11 的状态数上限不能替代这里的表容量检查；S13 的候选缓存另设预算，见 6.6。

`viterbi` 测试通过 2461 项检查、60 个穷举小图：独立按 Eq.9 的向量奇偶计算每个前缀状态，核对所有层/状态的最小度量、前驱状态/比特、最终回溯及重复调用。覆盖 delta=0、K=0、L=0、零列、零权重/ties、非零目标、列映射、容量/整数溢出、数值溢出及不可达目标；不可达目标使用额外连通性探针，不把它误称为正常恒等块图的现象。

## 6. S13：按最后不同分支划分的串行列表

### 6.1 与原文的对应及差别

本轮重新核对 1994 原文 §2.3（印刷页 315–317）、Fig.4/5、Serial Algorithm 和 §2.4。论文先用 VA 保存各时刻/状态的局部最佳前缀，再按最后合并位置寻找后续路径，保存尚未输出的候选，以免每次重算整个列表。印刷页 317 还讨论了前向 VA 后接后向候选树搜索的关系。

本实现采用上述思想的**不相交固定后缀分组**：缓存 S12 最佳前缀，用最后不同的分支定义子组，用全局最小堆保留每组的最轻代表。它是 S13 允许的等价有序 k-shortest 生成器；并非逐行翻译原文的合并计数 C_j，也没有宣称达到原文具体实现的复杂度。正确性依据在 6.3，不依赖给算法贴 SLVA 名称。

输入仍是 S11 的固定起终点图，阶段位序仍是 `(f_I[A],f_P)`。目标仍为 Eq.9，返回值只有 `W_local=W(f_I[A])+W(f_P)`；不补齐 Ac、不计算完整候选、不判断完整 LC-GCD 早停。相同状态序列上的不同分支标签必须保留：零列可能使 bit=0/1 的两个前驱相同，仍是两条不同路径。

### 6.2 一个候选表示什么

记 `D(x,r)` 为所有合法路径中、阶段 `[r,L)` 与 x 相同的集合，前缀 `[0,r)` 自由。堆节点保存该集合的最轻代表 x、`free_length=r`、完整路径位串、局部代价和稳定序号。初始化根组 `r=L`，包含所有合法路径；代表就是 S12 得到的第一条最优路径。

输出 x 后，`D(x,r)` 中剩余的每一条 y 必有唯一的最大不同位置 j，且 j<r。按 j 分成以下子组：

```text
阶段 >j：保持 x 的后缀
阶段  j：取 1-x[j]
阶段 <j：自由，但必须从初态到达该分支的前驱
子组 free_length = j
```

从目标状态沿 x 后退即可得到每个 j 的后继状态 `q_(j+1)`。XOR 自反，因此替代分支要求的前驱为：

```text
p_alt = q_(j+1) XOR ((1-x[j]) ? column_masks[j] : 0)
```

若 `phi_j(p_alt)` 不可达，该子组为空，跳过。否则从 S12 的 `cells[j*S+p_alt]` 回溯出最佳前缀，接上替代分支和固定后缀，得到该子组最轻代表。后续还能在该代表的自由前缀内继续分组，因而允许多处偏离，不局限于只翻转整条路径的一位。

所有代表放进同一个二叉最小堆：父节点下标 `(i-1)/2`，左右孩子 `2*i+1、2*i+2`。每次出堆选择全局最小代表；其他组持续保留。这里复用了本项目 S05 的普通堆布局；S05 的 `TepNode` 保存排序位置集合，不能直接表达固定后缀和 Viterbi 回溯边界，因此没有强行共用该专用节点，也没有重构成泛型堆。

### 6.3 无重复、完整和有序的依据

1. **分组不相交且完整。** 对 `D(x,r)\{x}` 的每条路径，最大的不同位置唯一；按该位置划分的子组两两不交，且并集正好是剩余路径。根组覆盖整个可行集合；每次只替换被输出代表所在的组，故全局各组始终不相交，并覆盖所有未输出路径。无须保存全部历史位串或另设去重表。
2. **代表最轻。** 固定后缀唯一决定连接它的前驱状态；到该前驱的最轻前缀已经由 S12 求得。拼接它可得该组最轻完整路径，且从初态 0 到固定目标，所以每条返回路径满足局部 syndrome。
3. **出堆是全局下一条。** 每条未输出路径都属于某组，且不轻于该组代表；堆顶是所有代表中最轻者，因此也是所有未输出路径中最轻者。替换成子组后保持相同不变量，于是代价非递减。在没有返回数/资源截断时，有限图的合法带标签路径全部恰好输出一次。

这给出了算法不变量及归纳证明。小码穷举是具体 C 实现的回归证据，不将有限测试当作任意规模正确性的证明。若 max_candidates/max_bytes 不够，返回明确错误，不丢弃某组后继续冒充精确列表。

### 6.4 手算前三条及候选复用

继续使用 5.4 的图：阶段权重 `(6,3,1000,9,8,7)`、掩码 `(1,2,0,3,2,3)`、target=1。

| 输出/扩展 | 剩余候选代表（位串，代价，free_length） |
| --- | --- |
| 第 1 条 `100000`，代价 6；根 r=6 | 按 j=5、4、3、2 生成 `010001,10,5`；`110010,17,4`；`010100,12,3`；`101000,1006,2`。j=1/0 的替代前缀不可达 |
| 第 2 条取全局堆顶 `010001`，代价 10；r=5 | 原堆中的代价 12、17、1006 保留；新增 `000011,15,4`；`100101,22,3`；`011001,1010,2` |
| 第 3 条取全局堆顶 `010100`，代价 12 | 来自第一条留下的候选，不是只从第二条的新孩子中选最轻者 |

三条路径分别表示 `(f_I[A],f_P)=(10,0000)、(01,0001)、(01,0100)`。前三条测试会实际打印出来。这里原始权重仍按 original_columns 取值，阶段位串不是原始列顺序的完整错误估计。

### 6.5 调用、ties 与浮点度量

```text
slva_init：分配受限工作区，调用一次 viterbi_decode，保存第一条代表
slva_next：必要时扩展上一条的 pending 组 → 出堆 → 写入一条路径/代价
slva_peek：必要时做相同扩展 → 只读取下一代价，不增加 paths_returned
slva_free：释放前缀表、所有堆中路径和 pending 路径，清零迭代器
```

扩展发生在请求下一条/下一代价时，首次 next 不预先生成偏离候选。输出后只检查自由前缀是否还有非空子组；有则保留为 pending，没有则释放代表。这也使已输出全部路径时能直接报告 DONE。多次成功 peek 得到相同的下一代价；peek 可能分配内存和返回资源错误，不是 const 操作。

第一条沿用 S12 的入边 bit=0 优先。后续按 `(cost,sequence)` 比较，sequence 按候选创建顺序递增，每次分组按 j 从大到小生成。等代价时保证确定、可重放；不承诺整个列表都是全局字典序或 S05 的 tie-break。零权重的不同路径照常输出，不按状态/代价去重。

每个代表的代价在**阶段顺序**重新求和，与 S12 一致；不使用浮点减法的差分度量，也不把前缀/后缀按另一种括号顺序相加。不用 epsilon 把不同代价强行视为相等。固定后缀的逐项非负加法对前缀度量单调，因此更小的缓存前缀拼接后不会变得更大；浮点舍入造成的并列也遵循上述稳定序号。

### 6.6 资源边界及真实复杂度

`SlvaLimits` 必须显式给出正的 max_candidates 和 max_bytes；`path_limit=0` 表示不限制返回条数，仍受资源上限约束。初始化对象必须为 `{0}`，再次初始化前先 free；输入图及其数组、原始翻转权重在迭代期间只读借用并保持有效。每次 next 的输出由调用方提供 L 字节缓冲，不能与迭代器/图的内部缓冲别名。

- max_candidates 计算堆中的路径、pending 路径以及本次扩展全部临时孩子的同时存活数。这是扩展峰值限制，而不是仅限制下一步堆大小。
- max_bytes 计算本模块拥有的分配有效载荷：完整 Viterbi 表、堆容量乘节点大小、所有活路径的 L 字节数组，包括临时孩子。包含闲置堆容量；不包括借用的图/权重、调用方输出、malloc 管理信息及 realloc 实现内部的瞬时拷贝。
- 扩展先预检计数/字节预算，再构造全部子组，成功后才改变逻辑候选集。资源/分配失败不消费下一条，不写 next/peek 输出；释放已构造的临时孩子。堆容量可能已增长，但保留的代表、pending、稳定序号和已返回数不变。可释放并用更大预算重建迭代器来重放；未提供热调预算的额外接口。
- `SLVA_PATH_LIMIT` 只表示存在剩余路径但达到返回上限；`SLVA_DONE` 表示真正耗尽。无路径图初始化成功为空列表。L=0 时，目标为 0 有一条空路径，成本 0；之后 DONE。达到返回上限的同时恰好耗尽也返回 DONE。

设 S=2^delta，Q 为堆和临时候选规模。初始化调用 S12，时间 O(L*S)、空间 O((L+1)*S)。一次扩展最多产生 L 个孩子；当前实现每个孩子复制/回溯/重算 O(L)，并用 O(log Q) 入堆，因此最坏每次请求为 `O(L^2 + L log Q)`，候选内存 `O(Q*L+Q)`。若不发生扩展，next 为 O(L+log Q)，已准备好的 peek 为 O(1)。这是便于审阅的完整位串实现，**没有宣称论文优化版本的每条 O(L)**；尚未做大码性能实验。生产算法没有先枚举 2^K 再过滤，只有测试 oracle 使用全枚举。

### 6.7 本轮实际验证

WSL2 Ubuntu 24.04.4 / GCC 13.3.0 / CMake、CTest 3.28.3：`slva` 通过 61,378 项检查、153 个独立 Eq.9 穷举小图。逐排名比较代价；并列时核对合法性、无重复、全覆盖及再次初始化的逐位相同顺序，而不假定 oracle 与候选堆采用同一个 tie-break。

覆盖全零/相同/乱序/二进制小数及普通小数、大动态范围权重，非零目标、零列平行边、delta=0、K=0、L=0、多次 peek/next、前缀限制及耗尽边界、候选/内存限制、尺寸/数值溢出和不可达探针。另验证 70 位路径前 25 条的确定性、唯一性及与 S05 TEP 的代价一致。ASan/UBSan 含泄漏检查运行同一测试，无报告；全项目 CTest 13/13 通过。未做操作系统分配失败故障注入。

一次清理重编译曾出现 Make 时钟偏差提示，核对发现少数 build 产物时间戳超前约 40 秒。再次检查时系统时间已追平，无需改时间戳；最终再次清理重编译无该提示，CTest 仍为 13/13。未修改系统时钟或源文件时间戳。

```bash
cmake -S . -B build/wsl-debug -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=/usr/bin/gcc -DBUILD_TESTING=ON
cmake --build build/wsl-debug -j 4
cmake --build build/wsl-debug --clean-first -j 4
ctest --test-dir build/wsl-debug -R '^slva$' --output-on-failure -V
ctest --test-dir build/wsl-debug --output-on-failure
gcc -std=c11 -g -O1 -Wall -Wextra -Wpedantic -Werror \
  -fsanitize=address,undefined -fno-omit-frame-pointer -Iinclude \
  tests/slva_test.c src/search/slva.c src/search/viterbi.c \
  src/search/trellis.c src/search/local_constraints.c src/search/tep.c \
  src/core/gf2.c src/core/systematic.c -lm -o build/wsl-debug/slva_test_sanitized
ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 build/wsl-debug/slva_test_sanitized
```

### 6.8 固定样例参数与 main 调试入口

`src/apps/main.c` 的 `run_slva_demo()` 集中保存下列输入，执行 `build/wsl-debug/bp_lcgcd --demo-slva` 即可进入。它是纯 Syndrome 场景：参考向量全零，所以 `s_target=s_res`；不需要输入 LLR。所有下标从 0 开始，打印位串的最左位是数组第 0 项。

```text
H（3 行、6 列；第三行是前两行 XOR）：
0 1 1 0 1 0
0 0 1 1 0 1
0 1 0 1 1 1

flip_weights（原始列顺序） = (1000,3,6,9,8,7)
z_hint = (0,0,0,0,0,0)
s_target = s_res = (1,1,0)
N=6，m=rank(H)=2，K=N-m=4，delta=2，L=delta+K=6

S03：
sys_to_orig = (1,2,0,3,4,5)
orig_to_sys = (2,0,1,3,4,5)
U =                 D = U*H*Q =
1 1 0               1 0 0 1 1 1
0 1 0               0 1 0 1 0 1
1 1 1               0 0 0 0 0 0
U*s_res = (0,1,0)

S10（basis_order 按基列权重降序）：
I（系统化时的原始基列顺序） = (1,2)
basis_order = A = (2,1)，Ac 为空
P_A =               s_A =
0 1 0 1             1
0 1 1 1             0
P_Ac 是 0×4 矩阵，s_Ac 长度为 0

S11：M_local = [I_2 | P_A] =
1 0 0 1 0 1
0 1 0 1 1 1
original_columns = (2,1,0,3,4,5)
column_masks     = (1,2,0,3,2,3)
阶段权重          = (6,3,1000,9,8,7)
state_count=4；start=0；target=1
```

`s_A=(1,0)` 打包为 `1*2^0+0*2^1=1`，不要把打印位串 `10` 当成通常书写的二进制整数 2。这里 delta=m，A 是全部基列，但 A 的顺序与 S03 的基列顺序不同。第一条阶段路径 `100000` 对应原始翻转列 2，即原始坐标 `f=001000`；`H*f^T=(1,1,0)`，代价为原始权重第 2 项 6。

实际调用链如下；没有额外调用一次 Viterbi 来重复计算第一条路径：

```text
main (--demo-slva)
  → run_slva_demo
    → binary_matrix_systematize                         S03
    → local_constraints_prepare                        S10
    → syndrome_trellis_prepare                         S11
    → slva_init → viterbi_decode                        S13 初始化 / S12
    → slva_next → ready → expand_pending（按需）          S13 列表
    → slva_free / binary_systematic_free
```

demo 默认 max_candidates=128、max_bytes=65536、path_limit=0，列出全部 `2^K=16` 条后以 DONE 结束。权重、syndrome、预算和输出循环都在一个函数内，便于自己插入 printf；输出的是阶段路径和局部代价，没有补齐 Ac 或接入 S14。该函数的数组形状固定匹配此样例；若改 H 的秩、N 或 delta，应同时调整 K、缓冲区和 Ac 的形状。

可在 WSL 项目根目录启动 `gdb --args build/wsl-debug/bp_lcgcd --demo-slva`，在 `run_slva_demo`、`local_constraints_prepare`、`syndrome_trellis_prepare`、`viterbi_decode`、`expand_pending` 设置断点。前两条 next 之间的延迟扩展及全局候选保留，适合与 6.4 的表对照。

## 7. S14：完整 LC-GCD、Ac 补齐与最优性证明

### 7.1 原文及共用流程

本轮读取 Zheng/Ma 2024 原文 II-C、Eq.9/10、Algorithm 1（PDF 第 2–3 页，印刷页 443–444）。以本项目的翻转变量 f 表示原文的 TEP e：搜索约束是 `H*f^T=s_res`，最终 estimate 为 `reference XOR f`。三种输入模式仍复用 decode_prepare，共用现有 gcd_decode 流程；GCD 使用 TEP，LC-GCD 使用 S10→S11→SLVA，恢复、最佳候选和原始 H 校验共用。

Algorithm 1 第 2 行用 z 初始化候选，针对经典场景的零码字有效；外部 syndrome 场景的 reference 未必满足目标。本实现统一从尚无候选开始，用第一条已完成 Eq.10 的有效候选建立上界，不凭空假设一个可行解。每帧借用已配对的 S03 结果，不重新消元。新增代码独立编写，直接复用现有 GF(2)、S10 和 S13 接口。

### 7.2 从局部路径到完整翻转

SLVA 返回阶段位串 `(f_I[A],f_P)`，局部代价 `W_local=W_A+W_P`。先将前 delta 位按 basis_order 的原始列索引映射回系统坐标，将后 K 位放到系统向量 `[m,N)`。再调用已有矩阵乘向量计算：

```text
f_I[Ac] = s_Ac XOR P_Ac*f_P^T                 Eq.10
W(f) = W_local + sum_(i in Ac) f_i*flip_weights_i
```

Ac 的行序与 S10 的 basis_order 后缀一致。K=0 时乘积为零；Ac 为空时没有追加项。每个完成候选均回到原始列坐标，验证 `H*f^T=s_res`，然后严格比较完整 W，只在更小时保存 `reference XOR f`。最后再次用原始 H 校验返回 estimate：Classical 目标为零，其余模式为 s_target。不会把最终 estimate 的非负错误权重当作翻转代价。

完整代价从 SLVA 返回值开始，按 Ac 顺序只追加非负项，保证每个候选的 `W>=W_local`；避免重新按不同顺序求和导致浮点下界关系倒转。当前采用 double 的严格比较，无 epsilon。S15 的严格 oracle 和极端动态范围限制见第 8 节；S16 的同帧仿真比较尚未开展。

### 7.3 早停与截断

SLVA 顺序非递减，因此所有未访问候选的完整代价至少为下一条局部代价。`next_W_local >= best_W` 时可以证明当前 best 最优；等号允许停止，因为接口保留首个最轻候选，没有要求枚举所有并列解。这里对应 Algorithm 1 第 5–7 行，原文的 `>` 放宽为 `>=` 不影响最小代价。

刚返回路径的 `last_W_local` 也是所有后续路径的下界。若 `last_W_local>=best_W`，可直接证明，不调用会分配偏离候选的 slva_peek。特别地 delta=m 时 Ac 为空，首条的局部/完整代价相同，Viterbi-1 就足够；即使 max_candidates=1 也无需生成第二条。

| 结束条件 | status / stop_reason / 证明 |
| --- | --- |
| 上述下界成立 | SUCCESS / BOUND / true |
| 列表确已耗尽 | SUCCESS / EXHAUSTED / true |
| 达到完整候选数上限且尚无证明 | TRUNCATED / BUDGET / false，保留最佳 |
| Trellis/SLVA 资源预算不足 | TRUNCATED / RESOURCE / false，保留已有最佳；初始化失败时 candidate_valid=false |

分配器失败和数值错误仍返回各自错误，不能当作证明。顺序为已知下界→peek 的耗尽/下界→候选预算；peek 可能先遇到资源上限，该情况下只报告 RESOURCE，不声称已经检查过未知的下一条代价。

### 7.4 手算：首条合法路径还不够

```text
H = 1 0 1 0       weights=(10,10,1,1), reference=0000, s_res=11
    0 1 0 1       m=2, K=2, delta=1
A=(0), Ac=(1), P_A=(1,0), P_Ac=(0,1), s_A=s_Ac=(1)
```

| 查询 | 局部 `(f_A,f_P)` | W_local | 补齐 f_Ac | 原始 f | 完整 W |
| --- | --- | --- | --- | --- | --- |
| 1 | `010`，即 `(0,10)` | 1 | `1 XOR 0=1` | `0110` | 11 |
| 2 | `011`，即 `(0,11)` | 2 | `1 XOR 1=0` | `0011` | 2 |

第一条虽合法，但所有未访问路径下界尚小于 11；limit=1 时返回 TRUNCATED/BUDGET，保留代价 11。第二条得到完整代价 2，与 last_W_local=2 相等，所以所有后续完整代价都不小于 2，立即证明最优。limit=2 仍返回 SUCCESS，不把“正好达到上限”误判为截断。原 6 位 `100000` 例在 delta=1 时补出的 Ac 为零，首条代价 6 即可证明，是另一种情形。

### 7.5 接口、计数及实际验证

`DecodeConfig.delta` 默认为 0，LC-GCD 校验 delta<=rank(H)。state_limit/slva_max_candidates/slva_max_bytes 为 0 时分别采用 4096 状态、65536 活候选、64 MiB SLVA 有效载荷。后者仍采用 6.6 的口径，不包含共用译码缓冲区、局部 P 块和借用的 H。candidate_limit 是完成完整代价计算的候选数；SLVA 的 path_limit 固定为 0，避免重复设置返回预算。

现有 candidates_popped/evaluated 对应实际取出/完成的候选；nodes_expanded 对 LC-GCD 计成功生成孩子的组，peek 可触发它；batches 仍为 0。未开展 S17 的额外统计系统。CLI 单帧支持 `--decoder lcgcd --delta N`、三种输入模式、--limit 和 --trace；仿真入口仍仅支持 GCD。

WSL2 Ubuntu 24.04.4 / GCC 13.3.0 / CMake、CTest 3.28.3：新增 lcgcd 测试通过 520 项检查、12 个固定小例的独立 ML oracle 帧；包括 Eq.10、三种模式/非零 hint、候选与资源截断、无候选初始化截断、非一致 syndrome、delta=m 一条路径证明。全项目 14/14 CTest 通过；ASan/UBSan 含泄漏检查运行 lcgcd_test 无报告。S15 的系统随机对照尚未开展。

S14 开始时基线为 12/13；cli_smoke 把前三条 trace 当作连续字符串，额外堆打印打断后误报。当前只提取 trace 记录核对其顺序，错误场景允许 stdout 上的调试信息，仍检查退出码和 stderr 错误内容；成功场景仍要求 stderr 为空。保留用户的打印、注释和原测试的结果/计数断言。

```bash
cmake -S . -B build/wsl-debug -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=/usr/bin/gcc -DBUILD_TESTING=ON
cmake --build build/wsl-debug -j 4
ctest --test-dir build/wsl-debug --output-on-failure
gcc -std=c11 -g -O1 -Wall -Wextra -Wpedantic -fsanitize=address,undefined \
  -fno-omit-frame-pointer -Iinclude tests/lcgcd_test.c src/decoders/decoder.c \
  src/search/slva.c src/search/viterbi.c src/search/trellis.c \
  src/search/local_constraints.c src/search/tep.c src/core/gf2.c src/core/systematic.c \
  -lm -o build/wsl-debug/lcgcd_test_sanitized
ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 build/wsl-debug/lcgcd_test_sanitized
```

## 8. S15：系统验证、delta=0 回归及数值边界

### 8.1 独立 oracle 与覆盖范围

新增 `tests/lcgcd_oracle_test.c`，只从原始 H 和原始坐标翻转 f 枚举全部 `2^N` 模式，直接计算 `H*(reference XOR f)` 及权重和。oracle 不调用生产 GF(2)、系统化、局部约束、TEP、SLVA 或候选补齐代码。秩另由原始 H 的像空间大小独立核对。每组输入检查所有 delta=0..m；候选上限依次为 1、2、全部 `2^K`、不限，重复数值跳过。

严格组使用整数/四分之一权重，测试规模内可被 double 精确求和，比较容差为零。分别测试 Classical、纯 Syndrome、带 hint 的 Syndrome 和 BP 后验输入；有证明时必须等于独立 ML 最小值，预算截断时必须保持原始 syndrome、最佳代价随预算不增、证明为 false。delta=m 要求一条完整候选即可成功；m=0 按 delta=0 退化处理。并列最小值允许不同 delta 选不同解，同一配置重放必须逐位、逐状态及计数相同。

覆盖 64 个全部二进制 2×3 矩阵的四个 syndrome、固定种子 20261124 的 N=1..10 随机矩阵及配对列置换、附加 XOR 冗余行和零行、故意不一致 syndrome、K=0 恒等矩阵、带冗余行的自构 Hamming [7,4]、零/相同/乱序权重。检查 H、权重、LLR、目标/hint、系统化矩阵、行变换及双向列映射均未被修改。

### 8.2 实际发现并修复：delta=0 有限预算不等价

```text
H = 1 0 1 1       weights=(10,10,1,1), reference=0000
    0 1 0 1       target syndrome=10, candidate_limit=2
```

修复前实际 CLI 证据：

| 算法 | 第 2 个 f_P | 返回原始 estimate | 代价 | 结果 |
| --- | --- | --- | --- | --- |
| GCD | `10` | `0010` | 1 | SUCCESS，证明最优，退出 0 |
| LC-GCD(delta=0) | `01` | `1000` | 10 | TRUNCATED/BUDGET，退出 5 |

无局部约束时两个枚举器遍历同一集合，但 TEP 与 SLVA 的等代价顺序不同，有限预算下会影响返回值。这与 delta=0 精确退化为 GCD 的约定不符。新增测试在 fixture=1、delta=0、limit=2 首先报告状态不一致；修复前的断言日志和 CLI trace 保存在 build/wsl-debug 的 `s15-oracle-before.log`、`s15-delta0-*-before.log`。

唯一的生产逻辑修复是：共用内核仅在名称为 lcgcd 且 delta>0 时使用局部约束/SLVA；delta=0 直接走原 TEP 分支。因而包括并列解、停止原因和候选/扩展计数在内均与 GCD 一致，delta=0 不使用 SLVA 资源上限。S13 的独立单状态 SLVA 仍保留，不修改其并列规则。测试实际核对了 4780 组 delta=0/GCD 结果。

### 8.3 数值诊断：通过误差界不等于严格 ML

四组普通小数/大动态范围探针单独标记，使用 long double 参考（本轮 WSL 的 LDBL_MANT_DIG=64）。记正数和的尺度为 max(|a|,|b|)，诊断容差为 `4*(N+1)*DBL_EPSILON*scale`；该容差仅用于数值探针，生产比较和严格 oracle 仍无 epsilon。除检查此界外，单独累计返回解与较高精度 ML 不同的次数。

实测 164 个数值探针中，104 个报告代价与 long double 原始顺序求和不同，16 个已返回 SUCCESS 的解在较高精度下并非最轻。最大代价偏差约 0.200195312，最大较高精度 ML 超额约 0.100585938，相对超额约 1.00585937e-17；不是 16 个严格 ML 通过样例。long double 本身也有舍入，这些数值为该参考下的测量。

可复现的反例：

```text
H = 1 0 0        weights=(1e16,0.2,0.1), reference=000, s_target=11
    0 1 1
GCD 和 LC-GCD(delta=2) 实际返回 f=110，double 报告 cost=1e16。
另一合法 f=101，按输入 double 值的精确有理数权重和更轻约 0.1。
```

在 1e16 附近 double 的相邻值间隔为 2，追加 0.1 或 0.2 均可能不改变累计值，比较因此成为并列。用 Python Fraction 对实际 CLI 输出独立核对，精确超额为 `3602879701896397/36028797018963968`。这说明当前 optimality_proven 针对实现的浮点评分，不能推广成任意动态范围下精确实数总分的 ML 保证。此次记录边界，没有引入跨模块高精度评分重构，也未开展大量性能实验。反例日志为 `build/wsl-debug/s15-numeric-*.log` 及 `s15-numeric-evidence.txt`。

### 8.4 实际结果与复现命令

| 项目 | 本轮实际结果 |
| --- | --- |
| WSL 环境 | Ubuntu 24.04.4、GCC 13.3.0、CMake/CTest 3.28.3 |
| 严格组 | 12852 个 LC 配置：10936 个证明最优、518 个正确预算截断、1398 个正确不一致 syndrome，无严格断言失败 |
| delta 分类 | delta=0：4780；delta=m>0：4680；中间值：3556，包含数值探针 |
| 并列重放 | 1574 组逐位/状态/计数相同 |
| 总规模 | 353 个 H/目标组合，13016 个 LC 配置，19370 次译码调用，323442 项检查 |
| 数值组 | 164 个按声明误差界检查；16 个较高精度 ML 差异单独记录 |
| 工程回归 | CTest 15/15；ASan/UBSan 含泄漏检查无报告 |

新测试仅在传入日志路径时重定向它自己的 stdout，每次译码截断重写该文件，保留最近调用的用户调试输出，避免 CTest 保存大量重复打印。摘要和失败参数写 stderr；失败会打印 H、软值、模式、delta、预算、目标、hint 和固定种子。生产打印及注释保持原样。

```bash
cmake -S . -B build/wsl-debug -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=/usr/bin/gcc -DBUILD_TESTING=ON
cmake --build build/wsl-debug --clean-first -j 4
ctest --test-dir build/wsl-debug -R '^lcgcd_oracle$' --output-on-failure -V
ctest --test-dir build/wsl-debug --output-on-failure
gcc -std=c11 -g -O1 -Wall -Wextra -Wpedantic -fsanitize=address,undefined \
  -fno-omit-frame-pointer -Iinclude tests/lcgcd_oracle_test.c src/decoders/decoder.c \
  src/search/slva.c src/search/viterbi.c src/search/trellis.c src/search/local_constraints.c \
  src/search/tep.c src/core/gf2.c src/core/systematic.c -lm -o build/wsl-debug/lcgcd_oracle_sanitized
ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 \
  build/wsl-debug/lcgcd_oracle_sanitized build/wsl-debug/lcgcd_oracle_sanitized_debug.log
```

一次增量构建因 WSL 产物时间戳超前而跳过重新链接，仍执行旧二进制；采用 clean-first 重编译后复核通过。未修改系统时钟或源文件时间戳。既有 systematic.c 调试格式/有符号比较警告仍保留。

## 9. S16：单线程同帧公平比较

### 9.1 共享输入和预算

新增 `simulation_compare(h, form, lc_config, config, results, paired_trace)`。直接复用两个 SimulationResult，results[0] 为 GCD，results[1] 为 LC-GCD。lc_config 的 candidate_limit 同时用于两算法，delta/SLVA 资源设置用于 LC-GCD。输入 H、配对的离线系统化结果以及 Eb/N0、seed、帧上限均只提供一份，不为两算法分别预处理 H。

共用的单线程 run_frames 每帧只调用一次 bpsk_awgn_transmit 和 bpsk_compute_llr，发送同一个全零码字；两算法借用同一个 y/LLR 缓冲和同一个 H/form。采用 S08 的现有 RNG，样本由 `(seed, Eb/N0, frame_index)` 决定，R=(N-rank(H))/N、sigma2 均共用。帧号从 0 开始，只有两算法都完成当前帧并返回有效候选才递增。资源不足且没有候选时返回错误，双方都不累计这一个未完成配对帧。

完整候选评估是预算单位；SLVA peek/前缀构造不计为完整候选。popped/evaluated/expanded/batches 沿用现有定义，未添加 S17 的分位数、运行时间或新统计系统。delta=0 继续直接走 GCD/TEP，保证并列顺序、逐帧输出、证明/停止状态及计数一致。新调度和 trace 代码独立编写，复用已有信道/译码模块。

### 9.2 共同停止和输出

每帧两边都完成后，再累计各自 FER/码字 BER、证明/截断帧数及 DecoderStats。达到 max_frames，或者任一方达到 max_frame_errors，整个配对运行一起停止。因此两边统计使用同一帧集合；在配对结果中 error_limit_hit 表示共同运行因误帧上限停止，即使另一方尚未达到该上限也为 true。

CLI 使用 `--decoder lcgcd --compare --delta N`；--compare 与 --simulate 互斥，前者要求 Classical 模式。--limit 为双方每帧相同的完整候选上限；任一方存在截断帧，CLI 返回 5，CSV 仍有效。单算法 --simulate 及其 24 列 CSV 保持原口径。配对 CSV 使用原有数值行格式，仅加 decoder、delta 两列，输出一个表头和两行，GCD 的 delta 记 0。

--trace 输出 paired_config、每帧两条 paired_input 和两条 paired_result。input 的 tx 是发送比特，noise 定义为实际浮点接收值 y 减去 BPSK 符号 x；y/LLR/noise 用 `%a` 十六进制浮点打印，能往返还原 double，避免十进制截断掩盖差异。它记录的是实现的实际加性增量，而不是重新调用 RNG 生成的噪声。原有 printf 仍在 stdout，可按 `^paired_` 筛选。

### 9.3 实际三帧和误帧上限样例

本轮 WSL 中运行：

```bash
cat > build/wsl-debug/s16-h7.txt <<'EOF'
3 7
1 0 0 1 1 0 1
0 1 0 1 0 1 1
0 0 1 0 1 1 1
EOF
build/wsl-debug/bp_lcgcd --decoder lcgcd --matrix build/wsl-debug/s16-h7.txt \
  --compare --delta 1 --ebn0 0 --frames 3 --errors 3 --seed 2026 \
  --csv build/wsl-debug/s16-three-frames.csv --trace > build/wsl-debug/s16-three-frames.log
rg '^paired_|^comparison ' build/wsl-debug/s16-three-frames.log
```

| frame_index | GCD 输出 / 完整候选数 | LC-GCD 输出 / 完整候选数 |
| --- | --- | --- |
| 0 | 0000000 / 1 | 0000000 / 1 |
| 1 | 0000000 / 1 | 0000000 / 1 |
| 2 | 0100101 / 7 | 0100101 / 3 |

N=7、K=4、sigma2=0.875；每帧两条 input 的 tx/noise/y/LLR 完全一致。双方都是 3 帧、1 误帧、3 误比特；完整候选总数分别为 9 和 5。frame=2 虽两边一致且证明最优，但发送的是全零码字，所以它仍是误帧。三帧计数仅为功能演示，不据此给出性能提升结论。

另实际运行相同 H/seed、Eb/N0=-3、delta=1、limit=1、frames=100、errors=2：两边均处理 3 帧，GCD 为 2 误帧/2 截断帧，LC-GCD 为 1 误帧/1 截断帧，双方共同 stop_reason=error_limit，CLI 退出 5。CSV/log 为 `build/wsl-debug/s16-shared-stop.csv`、`s16-shared-stop.log`。

### 9.4 验证与边界

WSL2 Ubuntu 24.04.4、GCC 13.3.0、CMake/CTest 3.28.3：simulation_pair 通过 9136 项检查，独立重放 99 对帧。分别核对输入 hex 浮点往返后的字节、共同停止时刻、直接译码重放的逐帧估计及汇总、独立原始 H 的欧氏距离 ML oracle。覆盖 delta=0..m、无限/单候选预算、CSV 与完整 trace 重放、原单算法兼容、非对称误帧上限、无有效候选时不部分计数。simulation_cli 增为 15 个场景；全项目 CTest 16/16，ASan/UBSan 含泄漏检查无报告。

```bash
cmake -S . -B build/wsl-debug -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=/usr/bin/gcc -DBUILD_TESTING=ON
cmake --build build/wsl-debug -j 4
ctest --test-dir build/wsl-debug -R '^simulation(_cli|_pair)?$' --output-on-failure -V
ctest --test-dir build/wsl-debug --output-on-failure
gcc -std=c11 -g -O1 -Wall -Wextra -Wpedantic -fsanitize=address,undefined \
  -fno-omit-frame-pointer -Iinclude tests/simulation_pair_test.c src/simulation/simulation.c \
  src/channel/channel.c src/decoders/decoder.c src/search/slva.c src/search/viterbi.c \
  src/search/trellis.c src/search/local_constraints.c src/search/tep.c src/core/gf2.c \
  src/core/systematic.c -lm -o build/wsl-debug/simulation_pair_sanitized
ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 \
  build/wsl-debug/simulation_pair_sanitized build/wsl-debug/simulation_pair_sanitized_debug.log
```

没有开展 S17 的正式统计系统、SNR 扫描、多线程或大量性能实验；本次只使用全零发送。S15 的极端动态范围精度边界仍适用；固定 seed 的浮点逐字节重放限定于同一运行环境，不承诺跨数学库一致。

## 10. S17：版本化统计、CSV/JSON 与测量范围

### 10.1 字段来源及单位

本轮复用 DecoderStats/SimulationResult，新增实际计数和汇总，不改变搜索顺序、评分或预算。S16 记录的 24/26 列未版本化 CSV 在本节点迁移为 schema_version=1，单算法和配对均使用同一份 40 列表头。JSON 是含 1/2 个记录的数组，字段名与 CSV 一一对应，整数按十进制 uint64_t 写出，query_limit_hit 为 JSON boolean / CSV 0 或 1。

| 字段 | 口径与代码位置 |
| --- | --- |
| partial_patterns_generated | 成功提交给枚举队列的部分模式代表数，包含根。TEP 在 push 内逐次加一；SLVA 使用 next_sequence，即最轻合法根路径与已提交子组代表总数。分配失败时回滚的临时孩子不计入 |
| full_candidates_evaluated | 现有 DecoderStats.candidates_evaluated；补齐完整 f、计算完整代价并通过必要校验后加一 |
| patterns_popped / nodes_expanded | 现有弹出数及实际生成孩子的父节点数，分别报告，不称为完整评估数或批次数 |
| query_batches | 现有 batches；当前算法串行，无批量查询操作，实际为 0 |
| trellis_states_visited | 成功 Viterbi 初始化中，层 0..L-1 的有限前缀状态展开次数；从保留的 cells 中逐格计数。排除终点层、不可达格及 SLVA 回溯读表；GCD/delta=0 为实际的 0 |
| avg_queries | full_candidates_evaluated / 实际计入 frames |
| p50_queries / p95_queries / p99_queries | 每帧完整候选数的 nearest-rank 分位数：排序后取第 ceil(p*F) 项，不插值，F 为实际帧数 |
| early_stop_reason | 整个记录只有一种逐帧停止原因时为该原因，否则为 mixed；与仿真循环停止原因区分 |
| query_limit_hit / query_limit_hit_frames | 是否存在/多少帧真正因候选预算停止。不是“评估数等于上限”；达到上限同时获得证明时不计入 |
| resource_limit_hit_frames | stop_reason=RESOURCE 的有效截断帧，与查询预算截断分别计数 |
| stop_bound_frames / stop_exhausted_frames | 下界证明/穷尽证明的帧数 |
| fer / ber / hd_ber | frame_errors/F、bit_errors/(F*N)、hd_bit_errors/(F*N)；BER 为整个码字比特口径 |
| runtime_seconds / runtime_kind | 计入帧的 decoder_decode 调用进程 CPU 时间总和，秒；kind 固定为 process_cpu_decoder_calls |

其余字段为算法名/delta、seed/Eb/N0、N/K/rate/sigma2、运行上限、实际帧数、误帧/比特数、证明/截断帧数和 run_stop_reason。后者为 frame_limit 或 error_limit，沿用 S16 的共同停止规则。单帧文本及 paired_stats 也报告正式搜索字段；旧文本的 candidates_evaluated 是 full_candidates_evaluated 的同一计数，没有维护两套评分或重复计数器。

### 10.2 分位数、停止状态和缺失值

每个算法内部临时分配 max_frames 个 uint64_t 保存逐帧完整候选数，完成运行后仅排序实际 F 项，写回三个标量并释放缓冲；配对运行分别计算。额外容量为每算法 `8*max_frames` 字节，不在输出结构中引入动态所有权。对配置导致的字节数溢出和分配失败保留必要检查；错误运行不作为完成记录导出。

有效完成/截断帧满足：停止原因计数之和=frames，proven_frames=bound+exhausted，truncated_frames=budget+resource。生成数不小于弹出数，当前 GCD/LC-GCD 正常记录的弹出数等于完整评估数。批次数可以为 0，而完整候选数大于 0。Trellis 分配容量 `2^delta*(L+1)` 也不等于实际前向可达状态展开数。

runtime 使用 C11 clock() 包围每次 decoder_decode，排除信道生成、离线 H 系统化、paired 输入/结果输出、分位数排序和 CSV/JSON 写出；包含每帧准备、分配、搜索、校验、统计取值以及译码函数内部的用户 printf。配对顺序为 GCD 后 LC-GCD。调试输出和缓存状态会影响结果，不能解释为 ASIC 延迟，也不据此声称吞吐提升。

若任一次 clock() 不可用或发生回绕，该记录的 runtime_available=false，CSV 写 NA、JSON 写 null，避免把部分测量当作完整时间。零时间也可能是计时分辨率造成的有效测量。尚未测量的 FLOPs、硬件周期和量子逻辑错误率不出现在当前 schema 中，不填虚假的 0。调用方只在仿真返回 DECODE_SUCCESS 后导出完成记录，截断情况在帧级字段中保留。

runtime 不要求重放一致；其他样本、查询统计和结果仍要求一致。测试只在文本重放比较时排除 runtime_seconds，另检查时间非负有限和格式互通。生产搜索仍使用严格 double 比较，S15 的极端动态范围限制未改变。

### 10.3 实际样例与验证

在 WSL 中运行同一三帧样例并同时导出两种格式：

```bash
./build/wsl-debug/bp_lcgcd --decoder lcgcd --matrix build/wsl-debug/s16-h7.txt \
  --compare --delta 1 --ebn0 0 --frames 3 --errors 3 --seed 2026 \
  --csv build/wsl-debug/s17-example.csv --json build/wsl-debug/s17-example.json
```

| 指标 | GCD | LC-GCD(delta=1) |
| --- | --- | --- |
| partial_patterns_generated | 18 | 11 |
| full_candidates_evaluated | 9 | 5 |
| query_batches / trellis_states_visited | 0 / 0 | 0 / 27 |
| avg_queries | 3 | 5/3 |
| p50 / p95 / p99 | 1 / 7 / 7 | 1 / 3 / 3 |
| early_stop_reason / query_limit_hit | bound / false | bound / false |

字段均来自同一次运行的结果；CSV/JSON 的 runtime 也一致，但重新运行的时间可以不同。--json 可单独使用，也可与 --csv 同时使用；两文件不可使用相同路径。JSON 做分组换行以便人工查看。实际运行只覆盖短样例，不是性能曲线复现。

WSL2 Ubuntu 24.04.4、GCC 13.3.0、CMake/CTest 3.28.3：全项目 16/16 CTest 通过。simulation 有 518 项检查，包括手算生成数 2/3、Trellis 展开数 5/11、达到预算但已证明最优及缺失时间的 null。simulation_pair 有 10972 项检查、111 对独立重放帧，用经验 CDF 独立计算分位数，并覆盖资源截断原因。simulation_cli 为 19 个场景，逐字段核对 40 个 CSV/JSON 字段、停止计数及 schema，ASan/UBSan 含泄漏检查无报告。新改动生产源码通过 -Wall -Wextra -Wpedantic -Werror 语法检查；保留既有系统化调试警告和用户排版，未开展 S18 扫描。

```bash
cmake -S . -B build/wsl-debug -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=/usr/bin/gcc -DBUILD_TESTING=ON
cmake --build build/wsl-debug -j 4
ctest --test-dir build/wsl-debug -R '^simulation(_cli|_pair)?$' --output-on-failure -V
ctest --test-dir build/wsl-debug --output-on-failure
```
