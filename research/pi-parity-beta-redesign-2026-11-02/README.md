# 继续改造：去掉多余对称，先积分再计算分母

2026-11-02。接续 `febb7f302263c04f63096c85c5df1f47d19d9813`。目标保持为构造更好的实际既约分数。本轮完成了一个更自由的整数原函数族、两种更便宜的截断规则，以及保留余项空间分布的误差证书。

## 1. 新构造：只在实际积分端点安排零点

原来要求原函数分子 \(R\) 为奇多项式，使积分分子为偶多项式。这也把 \(x=1\) 的零点复制到了积分区间以外的 \(x=-1\)。现在允许任意整数多项式 \(R\)，仅规定 \(R(0)=0\)、\(\deg R\le2c-1\)，并定义

\[
N=(1+x^2)R'-2(c-1)xR+\alpha(1+x^2)^{c-1}.
\]

仍然严格有

\[
\int_0^1\frac{N(x)}{(1+x^2)^c}dx
=\frac{R(1)}{2^{c-1}}+\frac\alpha4\pi.
\]

允许奇次分子不会自动引入对数。精确地说，整个有理系数分子像空间是

\[
\boxed{\deg N\le2c,\qquad[x^{2c-1}]N=0.}
\]

特别地，\(\deg N\le2c-2\) 时自动只含有理与 \(\pi\) 周期：奇数部分为 \(xP(x^2)\)，\(\deg P\le c-2\)；换元后积分只有 \((1+t)^{-s}\)、\(s\ge2\)，没有对数项。

新的固定预算族取

\[
\boxed{c=4m+1,\quad N=x^{2m}(1-x)^{2m}H(x).}
\]

完整空间允许 \(\deg H\le4m+2\)，附单个条件 \(H_{4m+1}=2mH_{4m+2}\)；直接取 \(\deg H\le4m\) 则不需这个条件。完整整数核的秩为 \(4m+2\)，旧奇性构造在相同极点阶数、相同真实端点阶数下的秩为 \(m+2\)。额外自由度可以用来同时选择积分相消与周期成本。

## 2. 奇周期成本没有增加，解析预算有所改善

对整数原函数端点核的周期像 \(\Lambda_m\)，每个奇素数 \(p\) 的局部像严格为

\[
\Lambda_m\otimes\mathbb Z_p
=p^{\lfloor\log_p(2m)\rfloor}\mathbb Z_p\times\mathbb Z_p.
\]

所以强制奇周期成本仍恰为 \(\operatorname{odd}\operatorname{lcm}(1,\ldots,2m)\)。证明仍是两个端点的 Taylor 整性，但充分性改用普通多项式 CRT，不再要求奇性。

另一方面，初等不等式

\[
\frac{x(1-x)}{(1+x^2)^2}<\frac2{11}\qquad(0\le x\le1)
\]

给出新候选的严格质量上界。令 \(h=\|H\|_1\)，\(p/q=-R(1)/(2^{c-3}\alpha)\) 完全约分，则

\[
\boxed{q^2|\pi-p/q|<\frac{22}{7}\frac{q^2h}{|\alpha|}
\left(\frac4{121}\right)^m.}
\]

由此形成新的联合构造预算 \(q^2h/|\alpha|=o((121/4)^m)\)。它是新坐标下的候选逐项上界。把旧分子原样嵌入时会产生 \((1+x)^{2m}\) 因子，可能把 \(h\) 放大 \(4^m\)；因此没有把上述解析因子当作同一旧候选的免费改善。

证明及全部坐标：[PARITY_FREE_ENDPOINT_FAMILY.md](new-kernel-design/PARITY_FREE_ENDPOINT_FAMILY.md)。

## 3. 有界构造检验：较低阶恢复已知基准

使用预先规定的有理采样点和五组成本权重，在 \(m=1,2,3\) 的完整饱和整数核中生成候选。生成阶段不输入 \(\pi\) 的小数、连分数或目标分数。事后对实际质量作精确有理认证。

\(m=1,c=5\) 时生成了

\[
\alpha=339,\quad R(1)=-4260,\qquad \frac pq=\frac{355}{113}.
\]

这是有限实现基准，不能据此声称旧数学族必须提高到 \(c=9\) 才能表示它。此处的作用是核查去掉奇性后，较低阶搜索确实能利用新的自由度。

全部 24 个保留候选均经独立标准库程序核对原函数恒等式、端点因子、实际 gcd、质量区间及上面的范数上界。三个有限周期像还通过二阶子式核对：

\[
\Lambda_1=\mathbb Z^2,\quad
\Lambda_2=3\mathbb Z\times\mathbb Z,\quad
\Lambda_3=15\mathbb Z\times2\mathbb Z.
\]

采样选择不是积分上界，保留的 LLL 基行也不是该格中的全局最优候选。

## 4. 原截断构造不必固定到 \(13m\)

原中心截断 \(S_K\) 满足精确一阶微分残量

\[
(1+t)S'_K+cS_K=\gamma_K(2t-1)^K,
\quad\gamma_K=\frac{c2^c(-1)^K\binom{c+K}{K}}{3^{c+K}}.
\]

只要 \(2K-c+3>0\)，就有保留空间分布的点态余项界

\[
\boxed{|S_K(t)-(1+t)^{-c}|\le
\frac{|\gamma_K|}{2K-c+3}\frac{|2t-1|^{K+1}}{1+t}.}
\]

它避免把整个区间的误差都按最坏点收费。对同一组 \(K=13m\) 分数，\(m=6\) 的有理质量上界约缩小 \(1.57\times10^7\) 倍，系数和既约分母没有变化；这是证书改善而不是新候选的质量改善。

两个截断矩还有正系数三项递推。连同额外权重 \(1-t\) 的同类递推，严格得到

\[
0<\mu_1(K)<\mu_0(K),\qquad \alpha_H<0
\quad(m\ge1,\ K\ge0).
\]

因此不需为非零周期保留 \(13m\) 的保护规模。采用更短的 \(K=2m\)，上述余项证书已经足够证明通常误差收敛。独立复核给出的一个上界为

\[
\left|\pi-\frac{p_m}{q_m}\right|
\le\frac43(6m+1)(2m+1)(7/8)^m\to0.
\]

这使以后可以直接按实际分母和误差成本选择截断阶数。原有两矩分母组合也被写成“边界项＋一个最高阶残量矩”的精确恒等式，可用于安排共同因子；没有把恒等式本身当成额外约分已经发生。

详见 [TELESCOPING_AND_PROFILED_REMAINDER.md](moment-denominator/TELESCOPING_AND_PROFILED_REMAINDER.md) 及 [独立审查](denominator-review/INDEPENDENT_DENOMINATOR_REVIEW.md)。

## 5. 先做 Beta 积分，可以消掉单独的增长因子

另一个可复用的规则是在 \(t=1\) 展开：

\[
S_K^+(t)=2^{-c}\sum_{k=0}^K\binom{c+k-1}{k}
\left(\frac{1-t}{2}\right)^k.
\]

它的每一项为正。先将每项与端点权重积分，所得有理矩的共同分母满足

\[
\boxed{L_K^+\mid2^{2m+1}\operatorname{odd}\operatorname{lcm}
(1,\ldots,6m+2K+3).}
\]

表面上的 \(2^{-k}\) 恰被 Beta 积分的 \(2^k\) 抵消。这消去了中心展开分母包络中单独的增长因子 \(3^{4m+K+1}\)，而且二进分母指数不随 \(K\) 增加。它描述的是统一分母包络，实际既约分母仍按 gcd 计算。

这种正项结构还给出单向误差、递减的项比和直接的有理尾界。固定 \(K=4m\) 有全规模通常收敛证明。这个简单端点版本的有界对照表现较差，因此保留其分母抵消机制，不扩大该版本的扫描。

详见 [POSITIVE_ENDPOINT_BETA_RULE.md](bernstein-design/POSITIVE_ENDPOINT_BETA_RULE.md)。

## 6. 一个不含 \(3^K\) 系数账的有理基对照

利用恒等式

\[
1+x^2=(1+x)\left(1-\frac{x(1-x)}{1+x}\right),
\]

还测试了有理截断

\[
S_K^{\rm rat}(x)=\sum_{k=0}^K\binom{c+k-1}{k}
\frac{x^k(1-x)^k}{(1+x)^{c+k}}.
\]

它的系数为整数。与新端点权重积分后，前三个矩可以写成 \(a_j+b_j\ln2\)，其中 \(b_j\in\mathbb Z\)，而有理部分的共同分母整除

\[
2^{c+K-1}\operatorname{lcm}(1,\ldots,c+K-1).
\]

选二次乘子的系数为两个周期向量 \(a,b\) 的叉积，就能同时抵消截断积分的有理与 \(\ln2\) 部分。最终 \(\pi\) 积分仍符合第 1 节的无对数次数预算。该机制去掉了单独的 \(3^K\) 系数分母，但叉积的系数高度也必须付费。

固定二次乘子、提高极点以消掉对数矩、以及保留增长自由度并加入周期成本特征，这几组有界对照都没有给出值得扩大扫描的候选。保留其显式规则、分母公式与数据，以免把同一成本转移重复当成收益。详见 [RATIONAL_BASIS_TRUNCATION.md](constructive-blocks/RATIONAL_BASIS_TRUNCATION.md)。

线性乘子还存在一个精确的成本检查：不限制系数大小时，它本就能编码任意有理数；强制 gcd 的收益必须与重建这些系数的高度一起核算。公式和有限整数审计见 [LINEAR_PERIOD_MAP_AND_GCD.md](constructive-blocks/LINEAR_PERIOD_MAP_AND_GCD.md)。

## 7. 当前继续构造的对象

主构造改用第 1 节的完整 \(R\) 空间，使自由乘子不必额外在 \(-1\) 消失；同时把第 4、5 节的正项递推、积分后分母抵消和点态余项保留为选择工具。

三项收益需要继续一起使用：合法的周期空间更大，解析预算更小，截断系数的分母可在积分后重新计算。本轮尚未把它们接成一条满足最终二次质量目标的统一序列；接续研究仍围绕实际约分分母与积分相消，不转回旧素数池排除。

## 复现

证明及输出复核使用 Python 标准库；完整饱和整数核生成还需要 SymPy。保存的候选可直接独立复核，无需重跑 LLL：

```bash
python new-kernel-design/verify_parity_free_output.py
python moment-denominator/profile_remainder.py
python bernstein-design/endpoint_truncation.py
python denominator-review/verify_denominator_independent.py
```

生成器 `new-kernel-design/full_parity_search.py` 只覆盖文中声明的三个规模。来源、结果范围和检查摘要见 [SOURCES.md](SOURCES.md)、[CHECKS.json](CHECKS.json)。
