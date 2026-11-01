# 近期文献与下一步接口（2026-11-01）

接续 [真实周期结构](../pi-period-structure-2026-11-01/README.md)。检索重点为 2025—2026 年；以下依据作者稿或出版社正文的指定部分做定向阅读，不视为全文独立验证。方法迁移尚未完成，未得到 π 坏逼近性结论。

## 1. 优先：结构化扰动

Khiar, Mainar, Peña, Royo-Amondarain, *Accurate determinant computation of Cauchy-polynomial matrices*。
[出版社正文](https://doi.org/10.1007/s11075-025-02188-5)。2025-07-29 在线发表，2026 年卷期。
阅读：§2、§6 Theorem 7、§7–8 方法说明。

原文提供特定 Cauchy 型行列式的乘积分解，以及依赖节点、极点相对间隔的结构化敏感性界。
迁移建议：先处理现有 Gram 核的 Cauchy 块，考察实际节点误差。目标是对归一化对数行列式的影响，而非任意矩阵扰动的稳定性。
缺口：导数重节点、横向类型和复节点需另证；原定理不能直接套到完整 Gram 矩阵。真实算术抵消与浮点误差抵消必须区分。

## 2. 算术与解析延拓共同优化

Calegari, Dimitrov, Tang, *Arithmetic holonomy bounds and effective Diophantine approximation*。
[作者稿](https://math.uchicago.edu/~fcale/papers/ICM.pdf)；[arXiv](https://arxiv.org/abs/2510.04156)；[2026-07 出版记录](https://authors.library.caltech.edu/records/meecf-0zp15)。
初稿 2025-10-05。阅读：§2.1 Theorem 2.1、§2.2 Proposition 2.4 及其限制、§3 的 π 应用段落。

原文以系数分母类型和解析延拓几何约束 holonomic 函数空间，并给出有效逼近估计。
迁移建议：试找生成函数，把已清分母的系数族、共同周期处的抵消、可延拓区域一起编码，再检查是否满足其微分方程和独立性条件。
缺口：现有增长矩阵尚未构成所需生成函数；原文不提供 π 的坏逼近性结论，也不能把有效指数界当作端点正下界。

## 3. 插值存在与范数控制的分离

Chalmoukis, Dayan, *Simply interpolating and Carleson sequences for Hardy spaces in the polydisc*。
[出版社正文](https://doi.org/10.1007/s12220-025-02155-5)，2025-08-18。
阅读：Theorem 1.1、§2、§4 修改后的 Szegő 核及 §6 构造概要。

原文在多圆盘中给出能简单插值、却不满足相应 Carleson 控制的序列，并研究归一化核的 Gram 矩阵。
迁移建议：将我们有限核矩阵的归一化、上下谱界及其随 H 的增长分别记录，避免用满秩替代定量控制。
缺口：原文是特定函数空间中的点值序列，现有模型还有增长的导数阶和权重，不能直接宣称其反例覆盖本构造。

## 4. 无穷维方向的具体入口

Gallo, Zweck, Latushkin, *Numerical Fredholm Determinants for Matrix-Valued Kernels on the Real Line*。
[出版社正文](https://doi.org/10.1007/s00020-026-02864-8)，2026-08-23。
阅读：§3–4 定义与展开、§5–6 截断和求积误差的框架。

原文处理矩阵值核的 Fredholm 行列式；在规定的算子类别及衰减、正则性条件下控制截断和离散化误差。
迁移建议：探索能否把现有离散核整理为有控制的算子极限，并保留清分母尺度和归一化对数行列式。
缺口：普通增长 Gram 矩阵不自动等于单位算子加迹类扰动。需先证明正确的算子类别、核衰减和极限关系。

## 工作顺序

1. 先在完整连续导数块上保留 Cauchy 乘积结构，推导实际节点误差下的对数行列式变化。
2. 再处理横向类型耦合，检查上述结构还能保留多少。
3. 生成函数与算术 holonomy 作为另一条路线；无穷维极限须等核条件得到证明后推进。

这些是基于文献的研究建议，并非已经验证的迁移定理。当前记录没有宣称任何新 π 逼近纪录。
