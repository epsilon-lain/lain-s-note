# 第三批正文首读：队列37–54

2026-11-01。固定上游 `openai/math@4eaee87548250b2682affda7d05d947c05554bae`；全部721篇仍在同一个队列中，按17学科轮转，不用当前与π的相关性决定是否阅读。

本批18项均已补齐正文与附录：238个源文本文件、85,437行。累计54/721篇，667篇未开始；另有30个目录外版本文件夹仍保留在补充队列。累计676个列入记录的源文本文件均通过精确Git blob SHA和实际行数核对。数字见[结构化汇总](round_03_summary.json)、[完整台账](LEDGER.md)和[来源核对](source_audit.json)。

这些数字是正文首读的覆盖，不是独立验证的定理数量。全库首轮尚未结束；本批继续保留逐篇构造机制与缺口，不据此结束全库阅读或提前筛掉远领域。

## 逐篇记录

| 队列位置 | 稿件 | 本批原始笔记 |
|---|---|---|
| 37 | Brennan's conjecture and sharp inverse-square integral means | [batch_20](batch_20/notes.md) |
| 38 | Symplectic Balls in Symmetric Polar Products | [batch_20](batch_20/notes.md) |
| 39 | The Factor-Two Hardness Threshold for Vertex Cover | [batch_20](batch_20/notes.md) |
| 40 | A smooth three-torus diffeomorphism with simple Lebesgue spectrum | [batch_21](batch_21/NOTES.md) |
| 41 | A counterexample to Hadwiger's conjecture | [batch_22](batch_22/NOTES.md) |
| 42 | A Complete Local Domain without a Small Cohen–Macaulay Module | [batch_21](batch_21/NOTES.md) |
| 43 | A Linear Clock for Random Walk on Tree-Weighted Planar Maps | [batch_23](batch_23/NOTES.md) |
| 44 | Rigidity of the Turing degrees | [batch_21](batch_21/NOTES.md) |
| 45 | An infinite finitely presented periodic group | [batch_23](batch_23/NOTES.md) |
| 46 | Spacetime Penrose inequalities: enclosing area, charge, and rigidity | [batch_24](batch_24/NOTES.md) |
| 47 | An irrational-trace counterexample to reduced Baum–Connes | [batch_25](batch_25/NOTES.md) |
| 48 | A marked tensor obstruction to four-dimensional disk embedding | [batch_31](batch_31/NOTES.md) |
| 49 | Lipschitz Equivalent Separable Banach Spaces Need Not Be Linearly Isomorphic | [batch_25](batch_25/NOTES.md) |
| 50 | An integral scalar curvature bound for real simplicial volume | [batch_27](batch_27/NOTES.md) |
| 51 | Nonuniqueness for the periodic hard-sphere Boltzmann equation | [batch_28](batch_28/NOTES.md) |
| 52 | The two-primary Birch–Swinnerton-Dyer formula in Selmer corank at most one | [batch_29](batch_29/NOTES.md) |
| 53 | Abelian covers, Gale correspondences, and the Hodge conjecture for powers | [batch_30](batch_30/NOTES.md) |
| 54 | A strict inverse-first-power bound for univalent functions | [batch_26](batch_26/NOTES.md) |

## 长稿如何补齐

初始部分记录保留当时实际读过的区间；新阅读分支只报告自己读过的段落。最终台账按同一SHA和行数合并，不把下载缓存算成阅读，不重复累计交叠行数。

- 第41项：初读正文至08及参数附录；[continuation_13](continuation_13.md)补齐09–13，最终8065/8065行。
- 第46项：49个文本文件、23468行；[带电与纯电多端延拓](continuation_10.md)、[四维解析构造](continuation_11.md)、[局部AdS](continuation_12.md)、[三维原数据刚性](continuation_14.md)、[四维刚性与恢复](continuation_15.md)，与batch_24的三维数值构造一起组成完整闭包。
- 第50项：[continuation_18](continuation_18.md)补齐degree与transfer两节；最终6764/6764行。
- 第52项：[continuation_16](continuation_16.md)读08–11，[continuation_17](continuation_17.md)读12–15，连同初读01–07、根文件和文献表，最终9546/9546行。

## 本批新增的构造线索

以下是按原有前提登记的研究接口。每篇完整的来源锚点、所需条件和缺少的算术桥梁见原始记录及[IDEAS.md](IDEAS.md)。它们没有替代仍在进行的全库首读。

### 严格配对与等号分类可以分两步处理

第37项用每个阈值上的累计计数优势作同秩配对；第54项进一步利用临界值无并列，使配对严格提高同一个目标纤维上的权重。即使逐点增益没有统一正下界，正质量上的严格增益仍能把整体等号逼成退化；随后以独立的调和障碍排除退化。对π项目，候选规格是先找到保留同一个非零整数证人的严格改进，再单独分类等号。这里尚无满足全部秩、分母和占用约束的这种改进。

### 容量分层保留二次占用的交叉项

第42项在有限层、整数秩容量0≤n_j≤N下使用恒等式（j从0起）：

$$N\sum_j(2j+1)n_j-\left(\sum_jn_j\right)^2=\sum_jn_j(N-n_j)+2\sum_{j>\ell}n_j(N-n_\ell)\ge0.$$

这给出明确的二次下界机制，并显示跨层项不能随意删去。原文的n_j来自共同晶格商的层秩；要用于当前列高度问题，仍须先构造同一候选子式的真实层分解、容量和高度对应。不能把任意列分组直接称作这些层。

### 在每个有限尺度保持精确非零证据

第43项的单位周期cochain加上单值梯度修正后，仍保留同一个有符号配对恒等于1；解析能量估计随后作用于这同一个有限对象。另一个有限图论展开把同一行列式系数转成森林平方和。它提示先固定整数/同调证据，再在其核内做优化。当前π候选行尚无这样的精确配对、图关联表示或统一高度估计。

### 平均证书需要共同对象和共同可观察量

第39项用条件Jensen把平均量送回同一个兼容full witness；第41项区分正质量分支与逆多项式稀有相遇；第46项在等号分析中先证明活动极小切割共享相交点的切平面，才能排除有符号平均的相消。第48项则保留完整理想、幂零方向和整个square-free cube的数据。对共同基搜索，这些分别提出兼容性、概率尺度、相消和精确代数关系的要求；平均满足预算仍不自动产生一个满足全部整数约束的基。

第50项的transfer把耦合闭合条件写成一个共同整数关联矩阵的核，再舍入同一个有限向量；这保留了联合约束与无偏性。它提醒我们，分别舍入各行或各块可能破坏共同证书。其几何关联结构尚未在π候选矩阵上建立，相关秩、分母高度和非加性占用成本仍须重新证明。

### 二次等号后还要读到四阶和实际余项

第46项的AdS部分先确定非负二次缺陷的完整四维零空间；真正的degree-one零方向随后用四阶正平方与局部切向障碍排除，径向方向由精确守恒量处理。它还先把渐近质量通量写成有界的精确Green泛函，再由已证明的加权余项推出实际展开。该定理针对固定种子和给定分支，参数区间可依赖二者；没有统一邻域或任意种子的分支存在性结论。

### 积分性、单位性和归一化要分别记账

第52项初读明确分开height-one处的积分性与闭点处的单位原始性，并保留有限/奇异局部条件切换的共同尺度。剩余各节的转接条件和退化处理也已逐段补齐，详见两份续读。它们可以启发如何检查同一个整数证人的消失阶和归一化，但没有给π的候选矩阵建立这些局部算术结构。

本批其余领域的机制也按篇保存，包括弱信号的一阶信息与二阶损耗、有限局部规则对任意伪造中间态的处理、可分Banach空间构造的平方预算、低复杂度符号恢复、Boltzmann多尺度分支分离，以及Hodge稿的代数对应。暂时缺少接口没有成为跳读理由。

## 台账检查与范围

`reading.py`现要求`source_mapping_complete`与`source_closure_complete`同时明确为真，且全部列入区间实际覆盖、未读清单为空，才能接受完整标签。新增了“声明的行全读完但include映射未完成”的拒绝检查，四类不支持的完整标签均被拒绝。第一批12条旧记录已有原始scope_note明确说明主文件与递归正文闭包；这次只补上显式mapping标志，没有增加阅读区间或改变其状态。

SHA、行数、区间与闭包检查是来源和记录一致性检查，无法机械证明理解或数学正确性。标准TeX字形基础设施不算未取得的数学附录；外部引用、伴随稿的有效性、Lean或可执行有限证书并未因读完正文而自动核验。图说明和可读图源码均在各篇范围中注明；第38项独立示意PDF另有视觉读取记录。

没有据本批阅读证明π坏逼近或非坏逼近。仍需对同一个非零、已清分母的整数行列式同时建立严格小于1的解析界。下一断点是第55项 *Petty’s projection-volume conjecture in dimensions at least four*；全部未开始稿件继续保留。
