# 整数波谱模型的独立审计与有限观测工具

日期：2026-11-02。性质：严格框架、独立证明与可复现实验；不主张这些一般谱/格对应具有文献首创性，不主张已构造 π 的新纪录或无限成功子序列。

## 1. 对象、定义域与零模

令 \(\mathbb T^2=(\mathbb R/\mathbb Z)^2\)，Hilbert 空间为复值 \(L^2(\mathbb T^2)\)，正交基
\[
e_{p,q}(x,y)=\exp(2\pi i(px+qy)),\qquad p,q\in\mathbb Z.
\]
定义 \(D_x=-i\partial_x/(2\pi),D_y=-i\partial_y/(2\pi)\)。固定实数 \(\alpha\) 与尺度 \(s>0\)，
\[
H_{\alpha,s}=s^{-1}D_y^2+s(D_x-\alpha D_y)^2.
\]
其 Fourier 实现是对角自伴算子，定义域等价于 Sobolev \(H^2\)：矩阵
\[
\begin{pmatrix}s&-s\alpha\\-s\alpha&s\alpha^2+s^{-1}\end{pmatrix}
\]
正定、行列式为 1，故固定尺度下椭圆、紧预解。谱为
\[
\lambda_{p,q}(\alpha,s)=q^2/s+s(p-\alpha q)^2.
\]
常数模 \((p,q)=(0,0)\) 总为零；此外 \(q=0\) 不代表有理逼近的合法分母。为整个构造固定使用不变子空间
\[
\mathcal H_+=\overline{\operatorname{span}}\{e_{p,q}:p\in\mathbb Z,q\in\mathbb N\}.
\]
这个选择属于数学观测器的约定，不把零模当作成功逼近。\(q<0\) 只重复同一有理方向。

整数指标来自周期边界条件/圆环拓扑；本模型并未从一个不含任何离散约束的对象中推导整数的本体起源。

## 2. 谱与完全约分后的逼近质量

令 \(d=\gcd(p,q)>0\)，\(P=p/d,Q=q/d\)。则
\[
\lambda_{p,q}=d^2\lambda_{P,Q},\qquad
Q^2|\alpha-P/Q|=Q|P-\alpha Q|\le\lambda_{p,q}/(2d^2).
\]
第二式直接由 AM–GM 得到。固定非零误差模式有
\[
\inf_{s>0}\lambda_{p,q}=2q|p-\alpha q|,
\qquad s_*=q/|p-\alpha q|.
\]
若 \(p=\alpha q\)，下确界是 0，由 \(s\to\infty\) 达到极限，有限尺度并不有零特征值。

因此，对无理 \(\alpha\)，构造尺度与模式使 \(\lambda\to0\)，等价于构造既约有理逼近质量 \(Q^2|\alpha-P/Q|\to0\)。如果模式的分母有界，有限个最近分数有正质量下界，所以成功序列的约分分母必趋于无穷。

这正是二维 unimodular 格
\[
\Lambda_{\alpha,s}=\{(q/\sqrt s,\sqrt s(p-\alpha q)):p,q\in\mathbb Z\}
\]
最短向量问题的谱表达，与已有 Dani correspondence 同属一个结构。它是可用的几何/谱接口，不是绕过逼近难点的自动定理。

## 3. 严格小于 1 的低谱锁定唯一既约有理方向

**定理。** 固定 \(\alpha\in\mathbb R,s>0\)。如果两个合法模式满足
\[
\lambda_{p_1,q_1}<1,\quad\lambda_{p_2,q_2}<1,\quad q_1,q_2>0,
\]
则 \(p_1q_2-p_2q_1=0\)。因此所有小于 1 的模式都约分为同一个有理数 \(P/Q\)，并恰为其若干整数谐波 \((kP,kQ)\)，\(k\ge1\)。若这个低谱非空，其最低模式必为既约 \((P,Q)\)。

**证明。** 取二维实向量
\[
v_j=(q_j/\sqrt s,\sqrt s(p_j-\alpha q_j)).
\]
有 \(\|v_j\|^2=\lambda_j\)，且
\[
\det(v_1,v_2)=q_1p_2-q_2p_1\in\mathbb Z.
\]
面积界给出
\[
|q_1p_2-q_2p_1|\le\|v_1\|\|v_2\|=\sqrt{\lambda_1\lambda_2}<1.
\]
所以这个整数为零。共线整数对的既约化唯一，而 \(q_j>0\) 固定方向符号。若模式公因子为 \(d>1\)，除去公因子仍在合法空间且能量除以 \(d^2\)，故最低模式必既约。证毕。

**边界不能改为 \(\le1\)。** 取 \(\alpha=1/2,s=2\)，模式 \((0,1),(1,1)\) 都有能量 1，却是不同方向。

这个定理给观测器一个可操作含义：一旦用独立误差证书确认谱已低于 1，所有低峰都属于一个唯一既约分数；多个峰可能只是谐波，不能当作多个逼近发现。

## 4. 完全有限、无漏项的谱观测

输入仅为一个已证明包含目标的有理区间 \([a,b]\)、有理尺度 \(s>0\)、有理阈值 \(L>0\)。任何 \(\alpha\in[a,b]\) 的 \(\lambda_{p,q}\le L\) 必满足
\[
1\le q\le\lfloor\sqrt{sL}\rfloor,\qquad
p\in[aq-\sqrt{L/s},bq+\sqrt{L/s}].
\]
故一个无限维谱的阈值以下部分可以完全、有限枚举，而无需目标小数或连分数输入。

对每对整数 \((p,q)\)，残差区间为
\[
[p-bq,p-aq].
\]
令 \(r_-\) 是这个区间到 0 的距离，\(r_+\) 是两个端点绝对值的较大者。则精确能量包络
\[
q^2/s+sr_-^2\le\lambda_{p,q}(\alpha,s)\le q^2/s+sr_+^2.
\]
下界大于阈值者可排除；上界不超过阈值者为确定峰；其余必须保留为未决，不能用中点数值替代。

**独立控制实验。** `verify.py` 用 Machin 恒等式和反正切交错级数，产生包含 π 的有理区间。设置 \(s=4,000,000,L=7/1000\)，完整分母上界为 167，总计检查 668 对整数；唯一确定模式为 \((355,113)\)，无未决项，
\[
0.00682697815894<\lambda_{355,113}<0.00682697815896,
\]
从谱直接认证质量小于 \(0.00341349\)。355/113 是已知基准，实验的结果是确认观测器的截断完整性、实际约分与误差核验，不是新的逼近纪录。

同脚本另外做了 1671 个低谱两模式检查和一个精确方差恒等式检查，全部通过。`checks.json` 保存有理数原值。

## 5. 可供新构造调用的试探波接口

令有限波
\[
\psi=\sum_j c_j e_{p_j,q_j},\qquad q_j\ge1,\quad\sum_j|c_j|^2=1.
\]
其 Rayleigh 能量是
\[
\mathcal E_s(\alpha)=\sum_j|c_j|^2\lambda_{p_j,q_j}
=A/s+s(B-2\alpha C+\alpha^2A),
\]
其中
\[
A=\sum|c_j|^2q_j^2,\quad
B=\sum|c_j|^2p_j^2,\quad
C=\sum|c_j|^2p_jq_j.
\]
如 \(c_j\) 为 Gaussian 有理数，归一化之前计算后这些量仍可用有理数完全表达。它是关于 \(\alpha\) 的凸二次函数，所以区间上的最大值是两个端点值的较大者。

**构造充分条件。** 只要一个来自目标特殊结构的显式波态族与有理尺度满足
\[
\max(\mathcal E_{s_m}(a_m),\mathcal E_{s_m}(b_m))\to0,
\qquad\alpha\in[a_m,b_m],
\]
其 Fourier 支持中必能抽取能量趋于 0 的模式，约分后即得目标质量趋零。这个判据允许用微分方程、积分恒等式等统一估计三个矩，而不必先输入目标的连分数位。它本身不提供这样一族波态。

**重要审计：纯对角观测没有干涉收益。** 在这个算子下，能量是各模式能量的非负加权平均，复系数的相位全部消失。简单叠加更多波或改变相位，不能让 Rayleigh 能量低于支持中的最低模式。

进一步有精确恒等式，设 \(w_j=|c_j|^2\)，
\[
A(B-2\alpha C+\alpha^2A)
=A^2(\alpha-C/A)^2+(AB-C^2),
\]
\[
AB-C^2=\sum_{i<j}w_iw_j(p_iq_j-p_jq_i)^2.
\]
右侧第二项由不同有理方向的整数行列式决定，量化混合波态的方向散布成本。要取得小 Rayleigh 能量，构造必须同时控制目标失配与方向散布，不能只声称“相消”。若新增耦合算子或非线性动力学，需要重新证明其输出能量仍控制原始 \(H\)，才能服务逼近目标。

## 6. 扩展到复数时的两个不同算术模型

### 6.1 整数共同分母：\((p_1+ip_2)/q\)

目标 \(z=a+ib\)，\(p_1,p_2,q\in\mathbb Z,q>0\)。这是两个实数的同时逼近，用三维环面。体积保持的谱族应为
\[
\lambda_{p_1,p_2,q}(t)=e^{-4t}q^2+e^{2t}|p_1+ip_2-qz|^2.
\]
设 \(\delta=|p_1+ip_2-qz|>0\)，直接微分得到
\[
\inf_t\lambda=3\,2^{-2/3}(q^{1/2}\delta)^{4/3}
=3\,2^{-2/3}\bigl(q^{3/2}|z-(p_1+ip_2)/q|\bigr)^{4/3}.
\]
因此这个模型的自然同时逼近指数是 \(3/2\)，不能无声沿用实数的指数 2。

若硬用 \(e^{-2t}q^2+e^{2t}\delta^2\)，它仍是合法正定谱族，也确有下确界 \(2q^2|z-p/q|\)，但其三维格流不保持体积；这是更强的目标，不是自然 Dirichlet 尺度。

### 6.2 Gaussian 整数分母：\(p/q,\ p,q\in\mathbb Z[i]\)

把复分母也允许为 Gaussian 整数，四个实整数指标对应四维环面。谱族
\[
\lambda_{p,q}(t)=e^{-2t}|q|^2+e^{2t}|p-zq|^2,
\qquad q\ne0,
\]
对应实四维体积保持的格流，且
\[
\inf_t\lambda=2|q|^2|z-p/q|.
\]
这里指数 2 自然，但分母的类型已变；必须用 Gaussian gcd 与单位 \(\{\pm1,\pm i\}\) 做约分，不能把这个分母成本直接替换回实有理分母。

### 6.3 为什么一个复指数函数不够

\(\exp(2\pi iz)\) 在 \(\operatorname{Im}z\ne0\) 时模长不为 1；其整数次幂回到 1 的条件不是 \(z\in\mathbb Q+i\mathbb Q\)。若想用闭合观测检测两个分量的有理性，必须使用两个独立单位相位 \((e^{2\pi ia},e^{2\pi ib})\)，共同闭合阶才编码共同分母。另一种合法做法是上一节的 Gaussian 格，而不是混用这两个模型。

## 7. 对根目录观测器的独立代码审阅

已逐行读取 `integer-wave-observation-2026-11-02/spectral_certificate.py`。其有理平方根下取整、区间残差包络、有限高度上界、primitive 与谐波区分均正确。它用 \(L=1/10\) 得到 \(q\le632\)，primitive 低峰唯一 355/113，而完整低谱是其 1、2、3 次谐波，与上述低谱唯一方向定理一致。本目录 verifier 为独立实现，用较窄阈值作交叉核查。

## 8. 文献归属

- Dmitry Kleinbock, *Flows on homogeneous spaces and Diophantine approximation on manifolds*, https://arxiv.org/abs/math/9810036 ，作者版本 https://people.brandeis.edu/~kleinboc/Pub/matrices.pdf 。谱/格流与坏逼近性的对应属于已有 Dani / Kleinbock 理论，不宣称为本项目新发现。
- Alexander Gorodnik, *Lectures on Diophantine approximation and Dynamics*, TIFR 2010，https://www.math.uzh.ch/gorodnik/tifr/index.html ，Lecture I https://www.math.uzh.ch/gorodnik/tifr/lecture1.pdf 。课程页明确把实向量逼近性质编码为格空间流轨道性质。

本笔记中的有限区间证书、低谱唯一方向证明、复数模型指数比较与代码审计均已在此自足推导。尚无由 π 的特殊结构产生低 Rayleigh 波态的渐近构造。
