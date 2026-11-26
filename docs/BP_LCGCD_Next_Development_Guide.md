# BP_LCGCD 后续开发指南：通用 QLDPC / Quantum BP / BP-LCGCD / 仿真平台

> **本文件替代 `BP_LCGCD_Codex_Stepwise_Prompt_v2.md` 的后续 step 规划。**
>
> 后续开发不再按旧 S18、S19、S21……顺序推进，而是按真正需要人工停顿的系统节点推进。
>
> 本项目的目标是：
>
> - 保留并继续验证已有 GCD / LC-GCD；
> - 尽可能复用 `cpicic_decoder` 的稀疏矩阵、BP/NMS、量化和仿真工程经验；
> - **从第一版 BP 开始就实现通用 QLDPC 的 quantum BP/NMS，而不是先写 classical BP 再改量子版本；**
> - 构建通用 stabilizer/QLDPC 的 BP-LCGCD；
> - 建立能够同时支持 classical 和 quantum 曲线的完整仿真平台；
> - 最后在稳定的 float reference 上继续做量化与硬件导向优化。
>
> **当前不绑定 surface code、toric code、BB code、HGP code 等任何具体码族。**
> 所有算法核心首先服务于通用 stabilizer / QLDPC 表示。

---

# 0. 强制参考资料与来源优先级

Codex 在开始相关模块编码前，**必须实际打开并阅读工作区 `./docs/` 中的三篇关键论文**，不能仅根据文件名、摘要、既有聊天总结或自身记忆实现公式。
两个参考github工程仓库：
https://github.com/yinanli-bit/GCDpsgrand.git 作为GCD算法的参考
https://github.com/yinanli-bit/cpicic_decoder.git 作为BP/NMS算法的参考

三篇关键参考文献为：

1. **Jifan Liang, Qianfan Wang, Lvzhou Li and Xiao Ma,  
   “The BP-LCOSD Algorithm for Toric Codes”, 2024.**

2. **Jifan Liang et al.,  
   “A low-complexity BP-OSD algorithm for quantum LDPC codes”, 2025.**

3. **Qianfan Wang et al.,  
   “BP-LCGCD: A Gaussian-Elimination-Free and High-Performance Decoder for Surface Codes”, 2026.**

同时必须参考：

4. **`cpicic_decoder` 仓库**
   - 重点用于：
     - 四方向稀疏 parity-check / Tanner graph 数据结构；
     - BP/NMS 的工程组织；
     - min1/min2 check update；
     - flooding / layered schedule；
     - message workspace；
     - quantization / saturation；
     - simulation point sweep；
     - frame/error stopping；
     - 统计与多线程组织。
   - **不能用 `cpicic_decoder` 的 classical BP 公式替代上述量子论文的 quantum BP 公式。**

5. **当前 `BP_LCGCD` 已验证代码**
   - GF(2)
   - systematic form
   - GCD
   - local constraints
   - Trellis / Viterbi / SLVA
   - LC-GCD
   - exhaustive oracle
   - classical channel / simulation

---

## 0.1 来源优先级

发生差异时遵守：

```text
量子算法公式：
三篇论文正文 / 算法 / equation
        >
现有项目说明
        >
cpicic_decoder classical 实现
        >
Codex 自己的推断
```

```text
工程数据结构 / 内存组织 / simulation pattern：
cpicic_decoder
        +
当前 BP_LCGCD 工程结构
```

### 禁止

如果论文没有给出某个关键算法桥接公式，Codex 不得：

```text
“按常见做法”
“应该就是”
“为了统一接口”
```

自行定义新的 quantum metric、early stopping 或 candidate ordering。

必须：

1. 标出文献未明确的点；
2. 查三篇论文对应章节；
3. 检查现有 LC-GCD 数学条件；
4. 写清楚采用的等价转换及依据；
5. 若仍无法确认，则在当前 PAUSE 节点停止，不继续编码。

---

# 1. 三篇论文如何贯通

这三篇工作不是三套互相独立的方法，而是一条非常清晰的演进线。

---

## 1.1 2024 BP-LCOSD：重点参考 quantum 后处理表示

这篇论文最值得本项目借鉴的是：

```text
quantum BP posterior
        ↓
Soft Weight Tuple (SWT)
        ↓
X/Z partial SWT
        ↓
binary reliability
        ↓
MREIS / local constraints
        ↓
SLVA
        ↓
re-encoding
```

### Soft Weight Tuple

对 qubit `j`：

\[
w_j=(w_{j,I},w_{j,X},w_{j,Y},w_{j,Z}),
\]

其中：

\[
w_{j,P}=-\log p_{j,P}.
\]

hard decision 为 soft weight 最小的 Pauli。

### Partial SWT

论文将 Pauli SWT 转换成 X/Z binary partial SWT：

\[
w^X=(w^I\boxplus w^Z,\;w^X\boxplus w^Y),
\]

\[
w^Z=(w^I\boxplus w^X,\;w^Z\boxplus w^Y),
\]

其中：

\[
w_\alpha\boxplus w_\beta
=
-\log\left(
e^{-w_\alpha}+e^{-w_\beta}
\right).
\]

binary reliability 由对应 partial SWT 两个分量的绝对差得到。

### 对本项目的意义

这部分非常重要，因为它说明：

> quantum Pauli posterior 不能简单变成 `2n` 个彼此独立的 `abs(llr)`。

特别是 Y error 同时具有 X/Z binary component，其概率不能未经证明写成：

```text
cost(Y) = cost(X-bit) + cost(Z-bit)
```

2024 文献提供了从 Pauli soft information 到 binary X/Z partial reliability 的正规转换。

### 但不要直接移植它的 GE-based LCOSD

该论文的：

```text
MREIS
partial Gaussian elimination
H'_local
H'_enc
```

是 BP-LCOSD 路径。

我们的目标算法是 2026 的 **GE-free BP-LCGCD**。

因此这篇论文在本项目中主要提供：

- Pauli SWT；
- X/Z partial SWT；
- reliability；
- symplectic syndrome；
- quantum candidate soft metric；
- local-constraint post-processing 的语义；
- simulation 组织参考。

不要重新引入在线 GE。

---

## 1.2 2025 low-complexity BP-OSD：重点参考通用 QLDPC BP/NMS 与复杂度

这篇论文把方法明确放到了通用 QLDPC 上。

重点阅读：

```text
Section 2.2 BP decoding
Eq. (3)–(10)

Section 3.2 NMS
Eq. (12)

Section 3.3.2 quantum LCOSD
Section 3.3.3 early stopping
Eq. (21)–(26)

simulation examples
```

对本项目最关键的内容：

- 每个 qubit 的 X/Y/Z prior；
- syndrome-aware check update；
- Pauli commute/anticommute message；
- posterior X/Y/Z message；
- hard Pauli decision；
- syndrome matching termination；
- NMS normalization；
- average BP iterations；
- quantum soft weight；
- logical error rate；
- post-processing search count；
- trivial / DAI stopping思想。

### 重要

论文中的 NMS 参数写法与 2026 论文不统一。

2025 论文使用形式：

\[
\frac{1}{\alpha}\min|\lambda|.
\]

示例中：

```text
alpha = 1.6
```

有效乘法 scale 为：

\[
1/1.6=0.625.
\]

工程实现不要直接使用含糊的 `alpha`。

统一定义：

```c
double nms_scale;
```

并让 check update 直接乘：

```text
nms_scale
```

例如：

```text
nms_scale = 0.625
```

---

## 1.3 2026 BP-LCGCD：目标架构以此为主

这篇论文定义了我们最终最接近的目标：

```text
Quantum NMS
      ↓
posterior reliability
      ↓
GE-free LC-GCD
      ↓
BP-LCGCD
```

它明确将 GE-based OSD / LCOSD 替换成固定 basis 的 LC-GCD。

同时提供：

- Phase I：single-scale BP；
- BP 成功时可直接返回；
- Phase II：多个 NMS scale；
- posterior reliability selection；
- LC-GCD；
- multi-scale BP-LCGCD；
- LC-GCD invocation probability；
- complexity analysis。

### 本项目优先级

第一版只实现：

```text
single nms_scale
+
single quantum BP/NMS
+
single quantum LC-GCD post-processing
```

multi-scale / list-parameterized BP-LCGCD 后置。

---

# 2. 重新调整后的开发顺序

旧版本：

```text
TASK 1 quantum BP
TASK 2 multi-point simulation
TASK 3 quantum LC-GCD
```

不再采用。

原因是：

> 如果在完整 quantum LC-GCD / BP-LCGCD 数据语义尚未确定前，就把 simulation runner、CSV schema、decoder registry 和 sweep 大规模固化，后面必然再次修改输入、metric、failure definition 和统计字段。

因此新顺序为：

```text
TASK 1
通用 QLDPC 表示 + Sparse Pauli Tanner Graph + Quantum NMS
        ↓
PAUSE A

TASK 2
Quantum LC-GCD + BP-LCGCD + Logical Success
        ↓
PAUSE B

TASK 3
完整 Classical / Quantum Simulation Platform
        ↓
PAUSE C

之后：
Quantization
Multi-scale BP-LCGCD
Layered / parallel
Hardware-oriented optimization
```

---

## 2.1 注意：仿真不是到 TASK 3 才第一次出现

TASK 1 / TASK 2 内仍然必须存在：

```text
single frame
single parameter point
small deterministic test
short Monte Carlo smoke test
```

用于验证算法。

TASK 3 的含义是：

> 在 end-to-end algorithm semantics 已稳定之后，再正式建立多点 sweep、统一统计、绘图和大规模 Monte Carlo 平台。

因此不是：

```text
先完全不仿真
最后才仿真
```

而是：

```text
前两阶段：最小验证仿真
TASK 3：正式仿真平台
```

---

# 3. TASK 1：通用 QLDPC + Quantum BP/NMS

这是下一阶段第一优先级。

---

## 3.1 通用 stabilizer representation

使用：

\[
S=[H_X\mid H_Z].
\]

其中：

```text
Hx: m × n
Hz: m × n
```

Pauli stabilizer entry：

| Hx | Hz | Pauli |
|---:|---:|---|
| 0 | 0 | I |
| 1 | 0 | X |
| 0 | 1 | Z |
| 1 | 1 | Y |

定义：

```c
typedef enum {
    PAULI_I = 0,
    PAULI_X,
    PAULI_Y,
    PAULI_Z
} Pauli;
```

不要增加：

```text
SurfaceCode
ToricCode
BBCode
HGPCode
```

等码族类型。

---

## 3.2 syndrome

error：

\[
E=(E_X,E_Z).
\]

syndrome：

\[
s(E)
=
E_XH_Z^T+E_ZH_X^T.
\]

为了复用现有 GF(2) decoder，可以构造：

\[
H_{\rm syn}=[H_Z\mid H_X].
\]

则：

\[
H_{\rm syn}
\begin{bmatrix}
E_X\\
E_Z
\end{bmatrix}
=s.
\]

代码中必须明确区分：

```text
stabilizer matrix S=[Hx|Hz]
syndrome matrix   Hsyn=[Hz|Hx]
```

不要都命名为 `H`。

---

## 3.3 stabilizer correctness

初始化 QuantumCode 时检查：

\[
H_XH_Z^T+H_ZH_X^T=0.
\]

这是必要数学检查，必须保留。

它不是“过度防御”。

---

# 4. Sparse Pauli Tanner Graph

必须重点参考 `cpicic_decoder` 的：

```text
PrtChkMtrx
Entry
rowHeader
colHeader
left
right
up
down
```

这一设计是后续工程的核心资产。

---

## 4.1 一个非 I stabilizer entry = 一个 sparse edge

建议概念：

```c
typedef struct QbpEdge {
    size_t row;
    size_t col;
    size_t edge_index;

    uint8_t pauli;

    struct QbpEdge *left;
    struct QbpEdge *right;
    struct QbpEdge *up;
    struct QbpEdge *down;
} QbpEdge;
```

### 必须保留的特点

- row traversal：沿 `right`
- column traversal：沿 `down`
- O(number_of_edges) BP iteration
- 不扫 dense `m*n`

---

## 4.2 topology 与 message workspace 分离

`cpicic_decoder` 把 message 放在 Entry 内。

本项目建议保留 topology，但将动态消息主要放入 reusable workspace。

理由：

- paired decoders；
- multi-scale BP；
- float/quantized coexist；
- 后续并行；
- 同一 code topology 不应被一次 decode state 污染。

例如：

```c
typedef struct {
    double *v2c_x;
    double *v2c_y;
    double *v2c_z;

    double *c2v;

    double *post_x;
    double *post_y;
    double *post_z;
} QbpWorkspace;
```

`edge_index` 负责从 pointer graph 定位 message array。

---

# 5. Quantum NMS：严格按论文

第一版只做：

```text
floating-point
flooding schedule
quantum NMS
```

不要同时做 exact BP、layered、quantization、multi-scale。

---

## 5.1 物理信道先验：channel model 与 prior LLR 必须分层

Quantum BP/NMS 的输入软信息来自**物理 Pauli 信道先验**，不是 BP 后验，也不是 LC-GCD 的 SWT。

每个 qubit `i` 的物理信道先验统一写成：

\[
\Pr(E_i=I)=p_i^I,\qquad
\Pr(E_i=X)=p_i^X,\qquad
\Pr(E_i=Y)=p_i^Y,\qquad
\Pr(E_i=Z)=p_i^Z,
\]

并满足：

\[
p_i^I+p_i^X+p_i^Y+p_i^Z=1.
\]

Quantum BP 使用三元 prior LLR：

\[
\Lambda_i^P
=
\log\frac{p_i^I}{p_i^P},
\qquad
P\in\{X,Y,Z\}.
\]

因此必须明确区分：

```text
physical channel probabilities
        ↓
prior XYZ LLR Lambda
        ↓
Quantum BP/NMS
        ↓
posterior XYZ LLR mu
```

### 第一版必须支持：depolarizing channel

第一版 quantum simulation 采用最常见的 depolarizing channel：

\[
p_i^I=1-p,
\qquad
p_i^X=p_i^Y=p_i^Z=\frac{p}{3}.
\]

于是：

\[
\Lambda_i^X
=
\Lambda_i^Y
=
\Lambda_i^Z
=
\log\frac{1-p}{p/3}.
\]

其中 `p` 是 physical error rate，是 quantum simulation 的主要横轴之一。

### channel 与 BP kernel 必须解耦

算法 kernel 不应写死 depolarizing。

推荐接口允许二选一：

```text
(A) caller 提供 pI/pX/pY/pZ
(B) caller 直接提供 prior XYZ LLR
```

channel 模块负责：

```text
physical model
    ↓
Pauli prior probabilities
    ↓
prior XYZ LLR
```

BP kernel 只负责：

```text
prior XYZ LLR + measured syndrome + Tanner graph
    ↓
posterior XYZ LLR
```

以后增加 biased Pauli channel、twirled channel 或其他 Pauli channel 时，不应修改 QNMS 核心公式。

### 第一版不要做的事情

不要在 QNMS 内部根据 `p` 临时生成先验；`p` 属于 channel/simulation 层。

不要把：

```text
prior XYZ LLR
```

和：

```text
posterior XYZ LLR
```

使用同一个无语义名称的 `llr` 数组长期混用。

建议字段/变量名明确表达：

```c
prior_x / prior_y / prior_z
post_x  / post_y  / post_z
```

或等价明确命名。

---

## 5.2 V2C

每条 edge：

```text
mu_X
mu_Y
mu_Z
```

初始化：

\[
\mu_{i\to j}^P=\Lambda_i^P.
\]

---

## 5.3 Pauli commutation reliability

check node 只关心：

> qubit error 与该 edge stabilizer Pauli commute 还是 anticommute。

根据论文公式计算：

\[
\lambda_{i\to j}.
\]

不要在 check update 内部临时拼公式。

单独封装：

```c
double qbp_commutation_llr(
    const double msg_xyz[3],
    Pauli stabilizer_pauli);
```

实现必须直接对应论文公式。

---

## 5.4 C2V NMS

统一工程形式：

\[
m_{j\to i}
=
(-1)^{s_j}
\cdot
nms\_scale
\cdot
\prod_{i'\ne i}\operatorname{sgn}(\lambda_{i'\to j})
\cdot
\min_{i'\ne i}|\lambda_{i'\to j}|.
\]

### `nms_scale`

统一：

```text
0 < nms_scale <= 1
```

例如：

```text
0.625
```

禁止算法核心出现语义不明确的：

```text
alpha
```

---

## 5.5 min1/min2

这部分应重点复用 `cpicic_decoder`。

每个 check 一次扫描：

```text
sign product
min1
min2
min_index
```

然后每条 outgoing edge O(1) 选择：

```text
edge == min_index ? min2 : min1
```

不要为每条 outgoing edge 再扫描整行。

---

## 5.6 Posterior

对每个 qubit：

\[
\mu_i^P
=
\Lambda_i^P
+
\sum_{
j:\langle P,S_{j,i}\rangle=1
}
m_{j\to i}.
\]

hard decision：

```text
if all XYZ posterior >= 0:
    I
else:
    argmin(X,Y,Z)
```

每轮 hard decision 后计算真实 quantum syndrome：

\[
s(\hat E).
\]

匹配 measured syndrome 即可 BP early stop。

---

## 5.7 TASK 1 输出

至少输出：

```text
hard Pauli estimate
posterior XYZ
iterations used
syndrome matched
```

后处理必须能够直接拿到 posterior soft output。

---

## 5.8 TASK 1 暂时不做

```text
LC-GCD
multi-scale
quantization
layered schedule
OpenMP
surface-specific optimization
soft syndrome
measurement error
circuit-level model
```

---

# PAUSE A：第一次人工审查

完成以下内容后必须停止：

1. generic `Hx/Hz`；
2. commutation validation；
3. sparse four-direction graph；
4. Pauli label 正确；
5. flooding quantum NMS；
6. nonzero syndrome sign；
7. posterior XYZ；
8. hard Pauli decision；
9. syndrome early stop；
10. depolarizing single-frame test；
11. channel prior 与 posterior XYZ LLR 在 API/变量命名上明确分离；
12. no per-edge/per-iteration malloc；
13. 原 GCD/LC-GCD tests 全通过。

### 测试重点

不要堆大量无意义 invalid-input 测试。

重点：

```text
one-iteration hand calculation
row/column sparse traversal
Pauli anticommute table
nonzero syndrome
posterior
hard decision
replay
```

---

# 6. TASK 2：Quantum LC-GCD + BP-LCGCD

TASK 2 的目标是在做正式 sweep 之前，先真正打通：

```text
physical Pauli error
        ↓
measured syndrome
        ↓
quantum NMS
        ↓
posterior soft information
        ↓
quantum LC-GCD
        ↓
recovery
        ↓
logical success/failure
```

---

# 7. Quantum soft information：Prior、Posterior、SWT、Partial SWT 必须严格区分

这是 TASK 2 最重要的设计点，也是后续最容易被 Codex 混淆的地方。

整个 soft-information 数据流固定为：

```text
物理 Pauli 信道
pI, pX, pY, pZ
        ↓
prior XYZ LLR Lambda
        ↓
Quantum NMS/BP
        ↓
posterior XYZ LLR mu
        ↓
posterior Pauli SWT
(wI,wX,wY,wZ)
        ↓
X/Z partial SWT
        ↓
binary hard decision + binary reliability
        ↓
LC-GCD local search / ordering
        ↓
reconstruct complete Pauli candidate
        ↓
posterior Pauli SWT full scoring
```

### 核心语义约束

```text
prior XYZ LLR     = BP 输入
posterior XYZ LLR = BP 输出
posterior SWT     = 由 BP 后验构造，供 quantum 后处理使用
partial SWT       = 从 posterior SWT 边缘化得到的 X/Z binary soft information
```

**禁止使用 channel prior SWT 代替 BP posterior SWT 做 BP-LCGCD 后处理。**

BP-LCGCD 的意义正是：

> BP 先利用 Tanner graph 和 measured syndrome 更新信道先验，再将 posterior reliability 提供给 LC-GCD。

如果 post-processing 仍使用原始 channel prior，则不再是论文中的 BP-LCGCD 语义。

---

## 7.1 Prior XYZ LLR 与 Posterior XYZ LLR

### Prior

物理信道给出：

\[
p_i^I,p_i^X,p_i^Y,p_i^Z.
\]

构造：

\[
\Lambda_i^P
=
\log\frac{p_i^I}{p_i^P},
\qquad P\in\{X,Y,Z\}.
\]

这是 BP 初始化时使用的 **a priori information**。

### Posterior

经过 Quantum NMS/BP 后：

\[
\mu_i^P
=
\Lambda_i^P
+
\sum_{j\in J(i),\langle P,S_{j,i}\rangle=1}
\mu_{j\to i},
\qquad P\in\{X,Y,Z\}.
\]

得到：

\[
\mu_i=(\mu_i^X,\mu_i^Y,\mu_i^Z).
\]

这是 **a posteriori soft information**。

hard decision 由 posterior 决定，而不是由 prior 决定。

---

## 7.2 Posterior SWT：它与 BP posterior LLR 是同一后验信息的不同表示

2024 BP-LCOSD 使用四元 Soft Weight Tuple：

\[
w_i
=
(w_{i,I},w_{i,X},w_{i,Y},w_{i,Z}),
\]

其中：

\[
w_{i,P}=-\log q_i(P),
\]

`q_i(P)` 表示 BP 迭代后 qubit `i` 为 Pauli `P` 的 posterior probability。

因此：

> SWT 不是另一套独立 soft information；它是 BP posterior probability / posterior LLR 的 Pauli 四元 soft-weight 表示。

### LLR 与 SWT 的相对关系

若 posterior LLR 采用本项目与论文一致的定义：

\[
\mu_i^P
=
\log\frac{q_i(I)}{q_i(P)},
\]

则：

\[
w_{i,P}-w_{i,I}=\mu_i^P.
\]

因此可在不显式恢复归一化 posterior probability 的情况下，先构造一个相对 SWT：

\[
\tilde w_i
=
(0,\mu_i^X,\mu_i^Y,\mu_i^Z).
\]

由于同一 qubit 的四个 soft weight 同时加减同一个常数不会改变 Pauli 相对概率，可以令：

\[
c_i
=
\min\{0,\mu_i^X,\mu_i^Y,\mu_i^Z\},
\]

并构造非负相对 SWT：

\[
w_{i,I}=-c_i,
\]

\[
w_{i,X}=\mu_i^X-c_i,
\]

\[
w_{i,Y}=\mu_i^Y-c_i,
\]

\[
w_{i,Z}=\mu_i^Z-c_i.
\]

这样：

```text
minimum SWT component = 0
```

且最可信 Pauli 的 soft weight 为 0。

### 工程要求

建议提供明确转换函数：

```c
qbp_posterior_to_swt(...)
```

不要让 QNMS 核心同时维护：

```text
posterior probability
posterior LLR
SWT
```

三套重复状态。

float reference 中应选一个主表示，例如：

```text
posterior XYZ LLR
```

需要 quantum post-processing 时再转换为 SWT。

### 必须测试

使用少量手工 posterior 值验证：

```text
hard decision from posterior XYZ LLR
==
hard decision from minimum SWT
```

并验证所有 Pauli pair 的相对 weight difference 与 posterior LLR 定义一致。

---

## 7.3 X/Z partial SWT：从 Pauli posterior 做概率边缘化

LC-GCD 的 GF(2) 搜索变量最终写成：

\[
(E_X|E_Z),
\]

但一个 binary X/Z component 对应多个 Pauli 状态。

因此不能直接取：

```text
X-bit weight = wX
Z-bit weight = wZ
```

必须严格参考 2024 BP-LCOSD 的 partial SWT。

定义 probability-domain soft-weight addition：

\[
w_\alpha\boxplus w_\beta
=
-\log\left(
 e^{-w_\alpha}+e^{-w_\beta}
\right).
\]

### X component

`E_X=0` 对应：

```text
I or Z
```

所以：

\[
w^X_0=w^I\boxplus w^Z.
\]

`E_X=1` 对应：

```text
X or Y
```

所以：

\[
w^X_1=w^X\boxplus w^Y.
\]

即：

\[
w^X=(w^I\boxplus w^Z,\;w^X\boxplus w^Y).
\]

### Z component

`E_Z=0` 对应：

```text
I or X
```

所以：

\[
w^Z_0=w^I\boxplus w^X.
\]

`E_Z=1` 对应：

```text
Z or Y
```

所以：

\[
w^Z_1=w^Z\boxplus w^Y.
\]

即：

\[
w^Z=(w^I\boxplus w^X,\;w^Z\boxplus w^Y).
\]

这一步本质上是在 posterior Pauli distribution 上做正确的 binary marginalization。

---

## 7.4 Binary hard decision 与 reliability

对每个 binary component：

```text
hard bit = lighter partial SWT component
```

reliability：

\[
r_X=|w^X_0-w^X_1|,
\]

\[
r_Z=|w^Z_0-w^Z_1|.
\]

这些 reliability 可以用于：

```text
fixed-basis local reliability ordering
A selection
local binary search ordering
```

但是必须记住：

> partial SWT / binary reliability 是搜索组织信息，不自动等价于完整 Pauli pattern 的最终 likelihood metric。

---

## 7.5 完整候选必须回到 Pauli SWT 评分

对一个完整 Pauli candidate：

\[
Q=(Q_1,\ldots,Q_n),
\qquad Q_i\in\{I,X,Y,Z\},
\]

完整 posterior soft weight：

\[
W_{\rm Pauli}(Q)
=
\sum_i w_{i,Q_i}.
\]

最终候选选择应基于这类 Pauli posterior soft weight，而不是简单使用：

\[
\sum r_X E_X + \sum r_Z E_Z.
\]

后者一般不能正确表达 Y error 的联合 posterior likelihood。

---

## 7.6 禁止直接使用 `2n independent bit cost`

禁止：

```text
weights[0:n]   = arbitrary X reliability
weights[n:2n]  = arbitrary Z reliability
candidate cost = independent binary sum
```

然后直接声称这是论文中的 quantum BP-LCGCD。

必须分别核对：

- local binary candidate ordering；
- partial candidate metric；
- complete Pauli soft weight；
- Y coupling；
- candidate ranking；
- early-stop lower bound。

如果 binary local metric 只用于生成候选，而不能构成完整 Pauli metric 的严格 lower bound，则：

```text
可以继续枚举
可以 candidate-limit 截断
但不能设置 optimality_proven=true
```

直到数学 bound 被文献或严格推导证明。

---

## 7.7 Prior SWT 与 Posterior SWT 禁止混用

这是 Codex 的硬约束。

### 错误做法

```text
depolarizing p
    ↓
channel prior SWT
    ↓
LC-GCD
```

然后把算法称为 BP-LCGCD。

### 正确做法

```text
depolarizing p
    ↓
channel prior probabilities
    ↓
prior XYZ LLR
    ↓
Quantum NMS/BP + syndrome
    ↓
posterior XYZ LLR
    ↓
posterior SWT
    ↓
partial SWT / reliability
    ↓
Quantum LC-GCD
```

如果 BP 被跳过，那么应明确称为：

```text
prior-based quantum LC-GCD
```

或其他准确名字，而不是 BP-LCGCD。

---

# 8. 如何复用现有 LC-GCD

现有 LC-GCD 最值得保留的是：

```text
offline systematic form
fixed basis
local A selection
Trellis
Viterbi
SLVA
A^c reconstruction
syndrome verification
```

但 quantum adaptation 必须拆开看。

---

## 8.1 可以直接复用的部分

只要 working binary matrix 已正确构造成：

\[
H_{\rm syn}=[H_Z|H_X],
\]

以下 GF(2) 结构可以复用：

```text
rank
column mapping
systematic form
local linear constraint
Trellis state transition
binary candidate reconstruction
syndrome consistency
```

---

## 8.2 必须重新证明/测试的部分

不能直接假设：

```text
binary partial cost
full quantum Pauli cost
early stopping
```

与 classical 完全相同。

特别是需要证明：

\[
W_{\rm partial}
\le
W_{\rm any\ full\ candidate\ in\ remaining\ group}
\]

才允许：

```text
optimality_proven = true
```

如果当前 quantum metric 下没有这样的严格 lower bound：

```text
可以继续有序枚举
可以有 candidate_limit
但不能声称 bound-proven ML
```

---

# 9. 2024 / 2025 LCOSD 与 2026 LC-GCD 的关系

Codex 必须明确：

```text
2024/2025:
MREIS + GE + LCOSD
```

```text
2026:
fixed basis + GE-free LC-GCD
```

因此：

### 复用 2024/2025

- quantum SWT；
- X/Z partial SWT；
- reliability；
- quantum soft weight；
- local constrained search 的语义；
- simulation metrics；
- early-stop ideas 作为对照。

### 不复用

- 在线 GE；
- 每帧构造 MREIS 后进行 GE；
- GE-based H'_enc 作为生产 decoder。

---

# 10. BP-LCGCD 第一版

第一版只做：

```text
one nms_scale
```

流程：

```text
QNMS
 ↓
hard Pauli estimate
 ↓
posterior SWT / partial reliability
 ↓
quantum LC-GCD
 ↓
candidate Pauli patterns
 ↓
full posterior soft-weight selection
```

---

## 10.1 BP 成功后是否继续 post-processing

三篇论文在算法目标上存在策略差异，不要把它误认为公式冲突。

2026 BP-LCGCD 为降低平均复杂度：

```text
Phase I BP succeeds
→ direct return
```

2025 的一个研究配置为了减少 BP mis-correction：

```text
always post-process
```

因此接口允许：

```c
typedef enum {
    BP_POST_ON_FAILURE,
    BP_POST_ALWAYS
} BpPostPolicy;
```

### 默认

目标 BP-LCGCD 第一版默认：

```text
BP_POST_ON_FAILURE
```

与 2026 低延迟架构保持一致。

`BP_POST_ALWAYS` 作为后续性能对照。

---

## 10.2 Multi-scale 后置

不要在 TASK 2 第一版做：

```text
alpha list
top-L posterior
multiple LC-GCD
```

先保证 single-scale end-to-end 正确。

multi-scale 属于后续性能增强。

---

# 11. Logical success

量子仿真不能把：

```text
syndrome matched
```

当成：

```text
decoding success
```

真实 physical error：

\[
E.
\]

decoder recovery：

\[
R.
\]

residual：

\[
D=E\oplus R.
\]

若：

\[
s(D)=0,
\]

仍可能是：

```text
nontrivial logical operator
```

第一版 logical success 判断：

\[
D\in\operatorname{rowspan}(S).
\]

即 residual 属于 stabilizer group。

---

## 11.1 为什么第一版不需要码族专用 logical strings

对 generic stabilizer code：

```text
residual in stabilizer row span
```

已经能判断 recovery 与 physical error 是否 stabilizer-equivalent。

因此第一版 LER 可以不依赖 surface-code logical path。

以后若需要分类：

```text
logical X
logical Z
logical Y
```

再增加 logical operator basis。

---

# 12. TASK 2 quantum oracle

对非常小的 `n`：

```text
enumerate all 4^n Pauli patterns
```

独立计算：

```text
syndrome
posterior soft weight
stabilizer equivalence
```

用于验证：

- posterior metric；
- candidate reconstruction；
- best candidate；
- degeneracy；
- logical success；
- early stop。

不要让 oracle 调用生产 LC-GCD 的：

```text
systematic
local constraints
SLVA
```

---

# PAUSE B：第二次人工审查

以下全部通过后，才开始正式 multi-point simulation：

1. channel prior → prior XYZ LLR 正确；
2. prior XYZ LLR → QNMS posterior XYZ LLR 正确；
3. BP posterior → posterior SWT 正确，并验证 SWT 来自 posterior 而不是 channel prior；
4. partial SWT 正确；
5. `[Hx|Hz] / [Hz|Hx]` mapping 正确；
6. quantum LC-GCD candidate syndrome 正确；
7. Y error metric 正确；
8. quantum full candidate cost 正确；
9. early bound 有严格依据，或明确关闭 proof；
10. BP-LCGCD single-scale end-to-end 可运行；
11. syndrome success / logical success 分离；
12. 4^n small oracle；
13. 能构造 degenerate case；
14. 能构造 syndrome-correct but logical-fail case；
15. 原 classical GCD/LC-GCD tests 不退化。

---

# 13. TASK 3：完整仿真平台

到 TASK 3 时，decoder semantics 已基本稳定。

这时正式学习 `cpicic_decoder` 的 simulation organization。

重点读取：

```text
main.c
chnnl_modu_llr_calc.*
dec_para.*
dec.*
parity_check_matrix.*
quantization.*
top_ctrl.h
test / Makefile
```

---

# 14. Simulation runner 不再写死两个 decoder

现有：

```text
GCD vs LC-GCD
```

paired simulation 应保留兼容 wrapper。

新的 core runner：

```text
one physical frame
      ↓
decoder[0]
decoder[1]
...
decoder[D-1]
```

所有 decoder 共享：

```text
same physical noise/error
same syndrome
same channel realization
```

---

# 15. Classical 和 Quantum 仿真域分开

不要强行统一物理横轴。

---

## 15.1 Classical

```text
BPSK
AWGN
Eb/N0
```

主指标：

```text
BER
FER
HD BER
avg queries
runtime
```

### 全零码字

线性码 + 对称 BPSK/AWGN 下继续使用 all-zero codeword 是正常 reference 做法。

不需要为了画曲线优先实现 random information encoder。

---

## 15.2 Quantum

第一版：

```text
depolarizing channel
```

每个 simulation point 必须由 channel 层先生成 `pI/pX/pY/pZ`，再转换为 prior XYZ LLR；BP 输出的 posterior XYZ LLR 与 posterior SWT 属于 decoder 输出/中间量，不能写回 channel prior 配置。

横轴：

\[
p.
\]

\[
P(I)=1-p,
\qquad
P(X)=P(Y)=P(Z)=p/3.
\]

主指标：

```text
Logical Error Rate
syndrome match rate
BP convergence rate
average BP iterations
average LC-GCD searches
LC-GCD invocation rate
runtime
```

不要把 Quantum 横轴写成 Eb/N0。

---

# 16. 三篇论文对 simulation 的共同启示

三篇论文的 simulation 重点高度一致：

```text
physical error-rate sweep
logical error rate curve
BP iteration complexity
post-processing search complexity
threshold/performance comparison
```

因此我们的 quantum simulation 不应只记录：

```text
final LER
```

还要记录：

```text
avg iterations
post-processing invoked frames
avg/full candidate searches
candidate/query distribution
```

这对未来硬件实现非常重要。

---

## 16.1 2024 BP-LCOSD 可参考的仿真组织

论文使用：

```text
depolarizing / twirled channel
logical error rate
multiple code sizes
max BP iterations
max TEP number
large Monte Carlo
```

其中：

```text
δ
BP max iterations
max TEP
```

必须成为 simulation metadata。

---

## 16.2 2025 low-complexity BP-OSD

重点对比：

```text
LER
average BP iterations
average number of searches
NMS parameter
```

因此 `SimulationResult` 需要明确支持：

```text
bp_iterations_total
bp_converged_frames
postprocess_frames
postprocess_candidates
```

---

## 16.3 2026 BP-LCGCD

未来 multi-scale 时还要统计：

```text
Phase I success rate
Phase II invocation rate
number of BP front-end runs
number of LC-GCD invocations
```

第一版可先预留统计语义，不必一次性全部实现。

---

# 17. Multi-point sweep

支持：

```text
start / stop / step
```

和显式 point list。

例如：

### classical

```text
Eb/N0:
0.0
0.25
0.5
...
```

### quantum

```text
p:
0.01
0.02
...
```

runner：

```text
for point:
    initialize channel parameters

    for frame:
        generate one physical sample

        for decoder:
            decode same sample
            accumulate
```

---

# 18. 每点停止条件

支持：

```text
max_frames
max_failures
```

### classical

failure：

```text
frame error
```

### quantum

failure：

```text
logical failure
```

不要用：

```text
BP syndrome mismatch
```

代替 logical failure。

它只是单独统计项。

---

# 19. Seed / replay

必须保持当前项目的可重放哲学。

建议随机样本由：

```text
global_seed
simulation_point_index / value
frame_index
domain
```

确定。

是否复用同一 random stream across points 需要固定规则并记录在文档中。

不要使用依赖 thread scheduling 的全局 mutable RNG。

多线程后也必须保持 frame-level reproducibility。

---

# 20. CSV / JSON

保持 long-form。

建议关键维度：

```text
domain
decoder
point_type
point_value

n
k
m

frames
failures

ber
fer
ler

bp_iterations_avg
bp_convergence_rate

postprocess_frames
avg_candidates
p50/p95/p99 candidates

runtime

delta
nms_scale
candidate_limit
max_bp_iterations
seed
```

不适用字段：

```text
CSV = NA
JSON = null
```

不要伪造 0。

---

# 21. Plotting

decoder library 只产生数据。

独立 Python tool：

```text
tools/plot_classical.py
tools/plot_quantum.py
```

### classical

```text
BER vs Eb/N0
FER vs Eb/N0
avg queries vs Eb/N0
runtime vs Eb/N0
```

### quantum

```text
LER vs p
avg BP iterations vs p
avg LC-GCD searches vs p
LC-GCD invocation rate vs p
```

---

# PAUSE C：仿真平台验收

正式大量 Monte Carlo 前确认：

1. multi-point works；
2. multiple decoders use same frame；
3. fixed seed replay；
4. classical old short simulation still matches；
5. quantum depolarizing sample replay；
6. LER correct；
7. BP iteration stats correct；
8. postprocess invocation/search stats correct；
9. CSV/JSON aligned；
10. plot script directly reads output；
11. no graph rebuild per frame；
12. no Gaussian elimination per frame；
13. no large malloc/free in inner simulation loop。

---

# 22. 量化：PAUSE C 后开始

量化确定要做，但必须建立在正确 float reference 上。

顺序：

```text
float QNMS
    ↓
float BP-LCGCD
    ↓
quantized QNMS
    ↓
quantized BP-LCGCD
```

---

## 22.1 复用 cpicic_decoder

重点复用：

```text
QuanChLlrFltToInt
SatuApp
SatuV2C
QuanV2CLowBitInt
SatuC2V
SignMagVnu
SignMagApu
```

但量子版本消息为：

```text
XYZ V2C triplet
scalar C2V
XYZ posterior
```

因此需要按 quantum message semantics 改造。

---

## 22.2 不优先量化 exact BP

`cpicic_decoder` 本身也没有直接支持 quantized ORIGIN_BP。

第一版硬件路线优先：

```text
quantized quantum NMS
```

而不是：

```text
phi/tanh exact BP LUT
```

---

# 23. Codex 编码约束：避免过度防御

后续 Codex 必须减少与算法无关的 defensive branches。

---

## 23.1 必须检查

只保留影响：

```text
算法数学正确性
内存安全
public API
资源预算
simulation correctness
```

的检查。

例如：

- public pointer；
- dimension；
- `Hx/Hz` compatibility；
- stabilizer commutation；
- syndrome length；
- allocation；
- message finite；
- resource limit；
- `nms_scale`；
- systematic mapping；
- logical membership prerequisite。

---

## 23.2 不要在 hot loop 重复检查 invariant

如果 constructor 已保证：

```text
edge row valid
edge col valid
pauli valid
workspace length valid
```

BP 每条 edge 不再检查。

Debug build 可使用：

```c
assert(...)
```

---

## 23.3 hot path 禁止动态分配

以下区域禁止 per-item/per-iteration malloc：

```text
BP iteration
check update
variable update
frame loop
candidate scoring
```

workspace 在 decoder 初始化时申请并重复利用。

---

# 24. Offline / online 分离

## Offline / code-level

```text
load Hx/Hz
validate stabilizer
build sparse graph
build Hsyn
systematize Hsyn
precompute mapping
allocate workspace
```

## Per frame

```text
sample Pauli error
compute syndrome
build prior
QNMS
post-process if needed
logical check
statistics
```

禁止每帧：

```text
sparse graph build
Gaussian elimination
systematic form build
large topology allocation
```

---

# 25. 不要过早抽象

不要一次设计：

```text
generic semiring BP
20 channels
10 schedules
plugin registry
arbitrary decoder factory
allocator abstraction
```

当前只需要：

```text
classical AWGN
quantum depolarizing
QNMS flooding
LC-GCD
BP-LCGCD
```

接口保持能扩展即可。

---

# 26. 不做无关重构

每个 TASK 只修改必要文件。

禁止：

```text
顺便重写 GF2
顺便统一全仓库 style
顺便改所有命名
顺便拆全部 decoder
```

已有 GCD/LC-GCD 有大量 oracle。

保持稳定比“结构更漂亮”更重要。

---

# 27. 测试哲学

测试优先验证数学性质，而不是 exhaust invalid-input combinations。

---

## Sparse graph

验证：

```text
every non-I stabilizer entry appears once
row traversal complete
column traversal complete
left/right consistency
up/down consistency
Pauli label correct
```

---

## QNMS

验证：

```text
one iteration hand calculation
commutation LLR
min1/min2
syndrome sign
posterior XYZ
hard decision
early stop
```

---

## Quantum LC-GCD

验证：

```text
partial SWT
binary reliability
syndrome consistency
Pauli candidate reconstruction
full Pauli metric
logical equivalence
4^n oracle
```

---

# 28. Attribution

从 `cpicic_decoder` 复制或明显改写的代码：

- 保留原 copyright；
- 保留原 author；
- 增加 adapted note；
- 记录源 repo / commit；
- 记录改动。

建议：

```text
docs/third_party_sources.md
```

---

# 29. Codex 每次执行前必须做什么

每个 TASK 开始前：

1. 打开当前仓库相关 source/header/test；
2. 打开 `./docs/` 对应三篇论文；
3. 阅读本 TASK 指定 section/equation；
4. 打开 `cpicic_decoder` 对应实现；
5. 写一个简短 implementation mapping：
   - paper equation → function
   - paper symbol → struct field
   - cpicic code → reused/adapted component
6. 再写代码。

不能：

```text
先写完代码
再去论文找一个看起来相近的公式
```

---

# 30. 每个 TASK 完成时 Codex 的报告格式

只报告真正有用的信息：

```text
1. changed files
2. algorithm mapping
3. reused cpicic_decoder parts
4. paper sections/equations used
5. key data structures
6. tests/oracles
7. known unresolved mathematical point
8. performance-sensitive allocation/loop notes
```

不要写大量“考虑了 30 种理论上的 NULL 情况”。

---

# 31. 最终开发路线摘要

```text
CURRENT:
Classical GCD / LC-GCD validated
        │
        ▼
TASK 1:
Generic Stabilizer / QLDPC
Sparse four-direction Pauli Tanner Graph
Floating QNMS
        │
        ▼
PAUSE A
        │
        ▼
TASK 2:
Posterior SWT
X/Z partial SWT
Quantum LC-GCD
Single-scale BP-LCGCD
Logical equivalence
4^n oracle
        │
        ▼
PAUSE B
        │
        ▼
TASK 3:
Generic multi-decoder simulation
Classical Eb/N0 sweep
Quantum p sweep
LER / iterations / searches
CSV / JSON / plotting
        │
        ▼
PAUSE C
        │
        ├── Quantized QNMS
        ├── Quantized BP-LCGCD
        ├── Multi-scale BP-LCGCD
        ├── Layered / parallel
        └── Hardware-oriented implementation
```

---

# 32. 当前最应该交给 Codex 的下一项任务

**只执行 TASK 1。**

目标：

```text
Generic Hx/Hz
+ stabilizer validation
+ cpicic-style four-direction sparse Pauli graph
+ floating flooding quantum NMS
+ posterior XYZ
+ syndrome early stop
```

要求：

```text
必须实际读取三篇论文；
必须实际参考 cpicic_decoder；
不实现 classical BP 作为中间版本；
不实现 LC-GCD quantum post-processing；
不实现 quantization；
不实现 multi-point sweep；
不实现 multi-scale；
不做 surface-specific code；
不进行无关重构；
不过度 defensive；
不在 hot loop malloc/free。
```

完成后在 **PAUSE A** 停止，进行人工审查。

---

# 33. 项目最终判断标准

后续代码是否值得保留，优先看：

1. 公式是否严格来自参考文献；
2. quantum Pauli 语义是否正确；
3. 是否保持 generic QLDPC；
4. 是否真正复用了成熟的 sparse / simulation / quantization 工程思想；
5. 是否避免 per-frame GE；
6. 是否避免 hot-path allocation；
7. GCD/LC-GCD regression 是否保持；
8. LER 是否按 stabilizer equivalence 统计；
9. float reference 是否足够清晰可作为量化 golden model；
10. 未来硬件实现是否能清楚找到：
    - message memory；
    - check unit；
    - variable unit；
    - min1/min2；
    - candidate search；
    - latency / iteration / search count。

**不要以代码量、抽象层数或 defensive check 数量衡量完成度。**
