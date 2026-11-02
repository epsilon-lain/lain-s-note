# 最小分母、周期小值与有理中心的真实成本

2026-11-02。这一轮继续研究同一个一维、双频率、完整源行包的原始列行列式。得到的主要推进是：**系数中的大公因子确实能改变结论，最小清分母后的实际周期值趋于零。** 之前用较大逐行阶乘清分母得到的增长结论仍然正确；此前 (n=4,ldots,10) 的有限表格不能推出最小清分母后也永远增长。

本文没有证明或否定 π 的坏逼近性，也没有独立认证原稿的关键高维插值定理。以下结论均限定在明确的稀疏一维构造；“独立复核”指另一代理重新检查论证，不等于外部同行评审。

## 1. 公因子来自哪里

记

\[
P_n(Z)=\det\bigl([z^s]e^{hz}(jZ+z)^a\bigr),\quad
j,h\in\{0,1\},\quad 0\le s,a<n,
\qquad A_n=\prod_{a=0}^{n-1}a!.
\]

将相应 Hankel 矩拆成两个端点的有符号积分，每个端点配置都带有阶乘因子。交叉部分的系数在 Laguerre 概率平均下仍是整数，这一步由 Schur 多项式的精确阶乘公式证明。最平衡的端点配置给出所有系数的公共因子：

\[
G_n=A_{\lfloor n/2\rfloor}^2A_{\lceil n/2\rceil}^2,
\qquad D'_n=A_n^2/G_n,
\qquad D'_nP_n\in\mathbb Z[Z].
\]

偶数 (n) 时，最高项恰为 \((-1)^{n/2}/D'_n\)，因此 (D'_n) 就是最小系数分母 (d_n)。奇数情形不需要猜测完整公因子：已证明的上下界只差一个 (O(n^2)) 因子，足以确定渐近式。

结合上一轮的实际周期双边界，得到

\[
\boxed{\log\bigl(d_n|P_n(2\pi i)|\bigr)
=-\tfrac12n^2\log n+
\left[\tfrac34+\log(\pi/\sqrt2)\right]n^2+O(n\log n).}
\]

所以最小清分母后的值趋于零。进一步用精确有理数上界和奇偶归纳证明：**对每个 (n\ge23)，都有 (0<|D'_nP_n(2\pi i)|<1)**。这不是由数值拟合得出的阈值。

详见[端点公因子证明](denominator-reflection-route/ENDPOINT_CONTENT_AND_DENOMINATORS.md)、[显式阈值](minimal-denominator-general/FINITE_THRESHOLD.md)和[独立整性审查](minimal-denominator-general/structure/INDEPENDENT_CONTENT_AUDIT.md)。前置双边界见[上一轮 Hankel 证明](../pi-period-mechanism-2026-11-02/raw-root-mechanism/HANKEL_PERIOD_MAGNITUDE.md)。

## 2. 近周期的零点不只是几个例子

对每个固定非零整数 (k)，证明了

\[
\limsup_{n\to\infty}n^2
\operatorname{dist}\!\left(2\pi ik,\{\rho\ne0:P_n(\rho)=0\}\right)
\le\frac{e^3(2\pi|k|)^3}{16}.
\]

证明使用反射结构给出的次数上界、精确最低项及根乘积不等式。还有统一推论：若 (K_n=o(n^{2/3}))，前 \(\pm K_n\) 个非零指数周期都能被不同零点一致逼近。这里的零点是代数数，尚未证明它们必在虚轴上、为单根，或达到最优收敛速度；更没有把它们变成有理逼近。

详见[一般近周期根定理](near-root-general/GENERAL_NEAR_PERIOD_ROOTS.md)和[独立复核](minimal-denominator-general/INDEPENDENT_ROOT_AUDIT.md)。这为“局部凹陷”提供了精确零点结构，但不把它解释为观测误差。

## 3. 真正留下的有理分母成本

偶数 (n\ge4) 时，令 (R_n=d_nP_n\)，它是最高项系数为 \(\pm1\) 的整系数偶多项式，次数 \(N_n=n^2/2\)。对既约分数 (p/q)，将 \(2p/q=a/b\) 约分；只要 (b\ge2)，就有精确等式

\[
\boxed{\operatorname{den}_{\rm reduced}R_n(2ip/q)=b^{N_n}.}
\]

因此，最优清除**系数分母**以后，**代入有理中心产生的分母**仍完整存在。例如 (n=10,p/q=355/113) 时，它恰好是 \(113^{50}\)。小于 1 的周期值本身不是非零整数，不能直接产生矛盾。

本轮把一个直接比较走到了结尾：选择偶数 (n\sim cb\)，其中 \(c>e^{3/2}\pi^2/2\)，可以严格得到一个 \(\exp[-O(q^2\log q)]\) 量级的误差下界。它远弱于坏逼近性要求的统一 \(c_0/q^2\)。这只确定当前比较法给出的界，没有证明所有更精细的比较都不可能。

详见[有理特化与导数成本](rational-comparison-next/RATIONAL_SPECIALIZATION_COST.md)及[独立审查](rational-comparison-next/INDEPENDENT_RATIONAL_COMPARISON_REVIEW.md)。

## 4. 增加整数频率能否帮忙

把频率 (0,1) 改成 (0,b)，有精确恒等式 \(P_{n,b}(Z)=P_n(bZ)\)。我们保留了相应高度预算和系数公因子，而没有把频率放大当成免费操作。

利用已有的 Breuer–Duits 中心极限定理，可以证明：每个固定整数 (b>1) 在充分大维数下，其最充分约分后的实际周期值都比 (b=1) 更大；两者仍可能同时趋于零。另一方面，端点分母上界保证 \(b_n\le n^{1/2-\varepsilon}\) 时周期值仍趋于零。这两条结论没有解决频率与维数共同变化时的最优选择。

详见[频率放大与算术内容](frequency-dilation-next/FREQUENCY_DILATION_CONTENT.md)及[独立文献适用性复核](minimal-denominator-general/INDEPENDENT_FREQUENCY_REVIEW.md)。外部输入明确标注为 [Breuer–Duits，arXiv:1309.6224v2，Theorem 2.6](https://arxiv.org/pdf/1309.6224v2)，不把其中心极限定理算作本轮的新证明。

## 验证与下一处缺口

附六个仅依赖 Python 标准库的验证器及精确结果：端点拆分、独立公因子检查、阈值、一般根距离、频率内容、有理特化。通用结论来自文中的证明，有限计算用于检查恒等式和证明所需的有限前提。具体命令、范围和文件 SHA-256 记录在 [CHECKS.json](CHECKS.json)。

接下来值得检验的是：能否在原稿关键高维参数中保留类似的公共因子，同时使有理特化的有效次数成本显著下降。单纯把这一维族的 (n) 继续增大，当前比较仍支付 \(b^{n^2/2}\)，不足以闭合坏逼近性的证明。
