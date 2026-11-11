# 观测约分因子的确定性分布与平均性质

2026-11-11。接续仓库的有限观测 gcd 模型。以下结论对**每个固定无理数**成立；取目标为 π 时，不需要假设 π 的数位或观测结果随机。这些是经典无理旋转与原始格点条带定理在当前模型中的推论，不宣称原创，也不证明 π 的坏逼近性。

## 1. 对象与频率

固定无理数 α∈ℝ。对正整数分辨率 B，定义

\[
N_B=\lfloor B\alpha+1/2\rfloor,\quad
g_B=\gcd(N_B,B),\quad
Q_B=B/g_B,\quad
E_B=B\alpha-N_B,\quad
G(B)=g_B^2/B.
\]

其中 N_B/B 是四舍五入读数，Q_B 是该读数的真实既约分母，E_B∈(−1/2,1/2) 是乘回分辨率后的带符号残差。α 无理，故不存在半整数取整歧义，且 E_B≠0。gcd 对负整数按绝对值定义，gcd(0,B)=B。

对集合 S⊂ℕ，若极限存在，定义自然密度

\[
\operatorname{dens}(S)
=\lim_{X\to\infty}\frac1X\#\{B\in S:B\le X\}.
\]

本文的“频率”“平均”均指上述确定性整数采样极限，不是对目标 α 取概率。

## 2. 固定因子整除的频率

对每个固定正整数 d，有精确等价关系

\[
d\mid g_B
\iff B=dq\text{，且存在 }p\in\mathbb Z
\text{ 使 }|q\alpha-p|<1/(2d).
\tag{1}
\]

证明：左侧意味着 B=dq、N_B=dp，而四舍五入条件恰为 |dqα−dp|<1/2。反向同样保证 N_B=dp，因此 d 整除分子和分母。由于 1/(2d)≤1/2，该 p 若存在则唯一，条件也可写成 ||qα||<1/(2d)。

无理旋转序列 {qα} 在 [0,1) 等分布。条件 ||qα||<1/(2d) 所占长度为 1/d，而 q≤X/d，因此

\[
\boxed{\operatorname{dens}\{B:d\mid g_B\}=d^{-2}.}
\tag{2}
\]

有限容斥进一步给出：对任意有限素数集合 S 和指定非负整数 a_p，

\[
\boxed{
\operatorname{dens}\{B:v_p(g_B)=a_p\ (p\in S)\}
=\prod_{p\in S}(1-p^{-2})p^{-2a_p}.}
\tag{3}
\]

这里 v_p 是 p 进赋值。证明只对有限多个附加整除条件容斥，并用 (2)；没有把有限容斥未经控制地延伸到无限多个素数。

## 3. 原始格点条带定理及其对 round 的适配

记 ζ(s)=Σ_{n≥1}n^−s（s>1）为 Riemann ζ 函数。所需的更强输入是 Estermann 的经典定理：对每个无理斜率 θ 和固定宽度 a>0，令

\[
F_\theta(Y,a)
=\#\{(q,p)\in\mathbb N\times\mathbb Z:
q\le Y,\ 0<q\theta-p<a,\ \gcd(p,q)=1\}.
\]

则

\[
\boxed{F_\theta(Y,a)/Y\longrightarrow a/\zeta(2).}
\tag{4}
\]

端点取严格或非严格不影响极限：一个固定端点方程 qθ−p=t 至多有一个整数解，因为两个不同解相减会使 θ 为有理数。Estermann 原始论文公式 (2) 证明了这种边沿通过原点的下侧条带版本；现代作者论文另陈述了任意位置条带版本。原论文特别处理了无限整除和的尾部，不能仅从 (2) 的固定 d 频率直接交换无限 Möbius 和与极限。

固定正整数 d，事件 g_B=d 精确对应

\[
B=dq,\quad N_B=dp,\quad \gcd(p,q)=1,
\quad |q\alpha-p|<1/(2d).
\tag{5}
\]

中心条带的下半侧直接用 (4)，宽度为 1/(2d)。上半侧将 p 换成 −p、斜率换成 −α，也成为同样宽度的下侧条带；原始性保持不变。两半相加，且横向上限为 X/d，得到

\[
\boxed{\operatorname{dens}\{B:g_B=d\}
=\frac1{\zeta(2)d^2}.}
\tag{6}
\]

这一步没有把 floor 互素结论直接套到 round，而是使用了 round 本身对应的中心条带。

因为 ζ(2)=π²/6，g_B=1 的密度为 6/π²。该常数来自原始格点的互素结构，对每个无理目标相同，并非目标恰为 π 时的特有性质。

## 4. 固定 gcd 与残差的联合分布

令 I⊂(−1/2,1/2) 为任意固定区间，长度记为 |I|。由 (5)，条件 E_B∈I 变成 qα−p∈I/d。

若 I 在零的同一侧，可将相应计数写成 (4) 中两个宽度的计数之差；正负两侧通过反射互相转换。若 I 跨越零，先将它分成两段。于是

\[
\boxed{\operatorname{dens}\{B:g_B=d,\ E_B\in I\}
=\frac{|I|}{\zeta(2)d^2}.}
\tag{7}
\]

因此，在渐近经验分布中，固定 gcd 取值与残差独立，而残差在长度为 1 的区间 (−1/2,1/2) 上均匀。这里是极限频率意义的独立性，没有断言相邻 B 的事件独立，也没有给出随 d 或区间宽度变化的统一误差估计。

## 5. 完整尾部与平方根指标的密度极限

对每个固定正整数 R，(6) 的有限求和给出

\[
\boxed{\operatorname{dens}\{B:g_B>R\}
=1-\frac1{\zeta(2)}\sum_{d=1}^{R}d^{-2}
\sim\frac1{\zeta(2)R}.}
\tag{8}
\]

ζ(2) 级数总质量为 1，因此这也给出 gcd 分布的尾部紧性。随后可以证明

\[
\boxed{G(B)=g_B^2/B\longrightarrow0
\quad\text{依自然密度}.}
\tag{9}
\]

其含义是：任取 η>0，集合 {B:G(B)>η} 的自然密度为零。证明如下：固定 R，当 B>R²/η 时，G(B)>η 强迫 g_B>R。因此该集合的上密度至多为 (8) 的右侧；令 R→∞ 即得零。

这个结果对坏逼近数、Liouville 数以及 π 同样成立。它与稀疏子列上 G(B) 很大或趋于无穷兼容。判定坏逼近性所需的 sup G 或 limsup G，不能由此处的密度极限替换。

## 6. 平均 gcd 发散，但平均保留分母比例收敛

对每个固定 R，由非负性与 (6)，

\[
\liminf_{X\to\infty}\frac1X\sum_{B\le X}g_B
\ge\frac1{\zeta(2)}\sum_{d=1}^{R}\frac1d.
\]

令 R 任意增大，调和和发散，因此

\[
\boxed{\frac1X\sum_{B\le X}g_B\longrightarrow+\infty.}
\tag{10}
\]

证明只用了有限截断的下界，不需要交换无界函数的平均与无限求和。这里没有给出对所有无理 α 统一的发散速度；平均 gcd 发散也不表明 g_B/√B 无界。

另一方面，Q_B/B=1/g_B 是有界量。对任意固定 s≥0，(6)、(8) 给出

\[
\boxed{
\lim_{X\to\infty}\frac1X\sum_{B\le X}
\left(\frac{Q_B}{B}\right)^s
=\frac{\zeta(s+2)}{\zeta(2)}.}
\tag{11}
\]

证明：先在 g_B≤R 上有限求和。剩余平均的绝对值至多为 {g_B>R} 的比例，因为 0<g_B^{−s}≤1。先取 X→∞，再由 (8) 令 R→∞，即可得到 ζ(s+2)/ζ(2)。这也明确给出了交换极限所需的尾部控制。

特别地，实际分母相对原分辨率的平均比例满足

\[
\boxed{\frac1X\sum_{B\le X}\frac{Q_B}{B}
\longrightarrow\frac{\zeta(3)}{\zeta(2)}
\approx0.73076.}
\tag{12}
\]

相对分母保留比例有极限，而其倒数 g_B 的平均发散，两者并不矛盾。这里不将固定 gcd 分布擅自推广为任何正阶无界矩的收敛公式。

## 7. 对 π 问题的作用与边界

本文提供确定性的总体统计基准：小素因子频率、固定 gcd 的完整分布、残差联合分布以及有界分母比例矩，都不需要 π 特有的随机性假设。

主文研究的临界问题仍是 g_B/√B 是否统一有界。密度为零的集合仍可以承载全部决定性尖峰；固定 R 的尾部公式也不能直接令 R=√B 来获得缩小靶的结论。若要研究 π 的坏逼近性，仍须对这种随 B 增长的阈值建立额外的 π 专属控制。

## 来源

- Theodor Estermann, *On the Number of Primitive Lattice Points in a Parallelogram*, Canadian Journal of Mathematics **5** (1953), 456–459，公式 (2) 与完整证明。[原始论文 PDF](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/03A2BF45B5DEA6F84F1BB7C97E890BDB/S0008414X00044151a.pdf/on-the-number-of-primitive-lattice-points-in-a-parallelogram.pdf)，[DOI 10.4153/CJM-1953-051-6](https://doi.org/10.4153/CJM-1953-051-6)。
- Yitwah Cheung, Arek Goetz and Anthony Quas, *Piecewise Isometries, Uniform Distribution and 3 log 2 − π²/8*，§4，Theorem 9 的条带表述及后续等分布证明。[作者 PDF](https://web.uvic.ca/~aquas/papers/paper38.pdf)。

本附录以以上条带定理为外部输入，round 观测的换元与后续推论均在正文展开。没有将有限实验当作这些渐近定理的证明，未运行 Lean，也未声称 π 的无限平方根界已得到认证。
