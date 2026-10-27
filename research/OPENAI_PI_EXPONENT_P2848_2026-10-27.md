# OpenAI 关于 π 的无理指数，与波面积的联系

2026-10-27。本文记录外部结果及其在既有模型中的条件推论；不宣称本文证明了 π 是或不是坏逼近数。波面积的初等比较未另行进行 Lean 核验。

## 1. 原稿说了什么

OpenAI 第 017 组原稿 [The irrationality exponent of pi is 2](https://github.com/openai/math/blob/main/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/paper.pdf) 主张 π 的无理指数 \(\mu(\pi)=2\)，更精确地：对每个实数 \(\varepsilon>0\)，存在阈值 \(Q(\varepsilon)\)，对所有整数 \(a\) 和 \(q\ge Q(\varepsilon)\)，

\[
\left|\pi-\frac aq\right|\ge\frac{1}{q^{2+\varepsilon}}.
\]

论文说阈值依赖 \(\varepsilon\)，其方法没有有效地给出阈值。它明确指出，这不是一个对所有 \(q\) 成立的 \(C/q^2\) 下界。

[官方 Lean 说明](https://github.com/openai/math/blob/main/lean/docs/017.md) 列有该主结论的形式化。[MathVet 的独立报告](https://github.com/mathvet/mathvet/blob/main/reviews/openai-math/STATUS.md) 记载 PiExponent 的构建、依赖闭包检查及 Comparator 检查通过；这是对该项目记录版本与形式化语句的报告。我没有独立重跑这些检查，也没有完成论文与全部 Lean 依赖的人工审稿。Lean challenge 文件里的 sorry 是待解题目占位，与另存的证明入口不是同一文件。论文另提 Flint–Hills 级数收敛；官方说明该结论在这项选定的 Comparator 语句范围外。

## 2. 波面积的推论

按已有的[长方形积分笔记](DIGIT_RECTANGLE_INTEGRAL_2026-10-27.md)，设 \(x=\pi-3\in(0,1)\)，\(q\ge2\)，\(b=q+1\)。设 \(d_j\) 为 \(x\) 的 \(b\) 进制第 \(j\) 位，定义

\[
E_b=\sum_{j=2}^{\infty}\frac{|d_j-d_1|}{b^j},\qquad
D_q=\min_{a\in\mathbb Z}\left|x-\frac aq\right|.
\]

其中 \(E_b\) 是先前模型里数位矩形波与最佳水平线之间的非负面积，最佳线的高度是 \(d_1/q\)。旧笔记给出存在只依赖 \(x\) 的 \(\lambda_x>0\)，使充分大的 \(q\) 满足

\[
\lambda_x E_{q+1}\le D_q\le E_{q+1}.
\]

例如 \(D_q\le E_{q+1}\) 直接来自对第一位数对应横线使用三角不等式；反向比较利用最邻近分子对应的横线、第一位横线的最优性以及数位尾项的界。

若接受原稿的主结论，就对每个 \(\varepsilon>0\) 与所有足够大的 \(q\) 有

\[
q^2E_{q+1}\ \ge\ q^2D_q\ \ge\ q^{-\varepsilon}.
\]

由于 \(x=\pi-3\) 只做整数平移，所有分母 \(q\) 的最小距离与 π 相同。这说明规范化面积不可能在无界的分母上，比任何固定负幂 \(q^{-\varepsilon}\) 还小；但是这个下界仍趋于零，也不知道从多大 \(q\) 开始成立。

进一步，应用两侧比较和连分数中的标准逼近定理，对任意固定无理数 \(x\in(0,1)\)：

\[
\mu(x)=2
\quad\Longleftrightarrow\quad
\forall\varepsilon>0\ \exists Q_\varepsilon\ 
\forall q\ge Q_\varepsilon:\quad
q^2E_{q+1}\ge q^{-\varepsilon}.
\]

右推左时，先取 \(\varepsilon/2\)，由 \(D_q\ge\lambda_xE_{q+1}\) 得 \(D_q\ge\lambda_x q^{-2-\varepsilon/2}\)；当 \(q\) 充分大时 \(\lambda_x q^{\varepsilon/2}\ge1\)，所以 \(D_q\ge q^{-2-\varepsilon}\)。结合所有无理数都满足 \(\mu(x)\ge2\)，得到等价。这个等价是把无理指数用既有面积表达，并不额外证明 π 的无理指数。

相对地，π 是坏逼近数等价于：存在 \(c>0\)，对所有足够大的整数 \(q\)，\(q^2E_{q+1}\ge c\)。这里需要固定常数；原稿指数估计不能直接提供它。

## 3. 为何指数等于二仍可能没有固定下界

举一个明确的无理数：

\[
\alpha=[0;1,2,3,4,\ldots].
\]

它的第 \(n\) 个部分商为 \(a_n=n\)，收敛分数的分母满足 \(q_n\ge n!\)。由标准公式

\[
\mu(\alpha)=2+\limsup_{n\to\infty}
\frac{\log a_{n+1}}{\log q_n}=2.
\]

另一方面，连分数逼近给出
\[
q_n^2\left|\alpha-\frac{p_n}{q_n}\right|
<\frac1{a_{n+1}}=\frac1{n+1}\longrightarrow0.
\]
所以 \(\alpha\) 不具备 \(c/q^2\) 的固定正下界。结合面积比较，它的 \(q_n^2 E_{q_n+1}\) 同样存在趋零子列。由此可见，无理指数为二与坏逼近是两个不同条件；这个例子不能判定 π 的情况。

## 4. 原稿证明的端点障碍

论文的[参数选择部分](https://github.com/openai/math/blob/main/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/conclusion.tex) 对固定 \(\nu>2\) 选择有理数 \(t\)，使

\[
\frac12<t<1-\frac1\nu.
\]

若把 \(\nu\) 直接取为 \(2\)，两端都变成 \(1/2\)，区间为空。因此该处的具体参数构造不能直接推到指数 \(2\) 的端点。这个观察不排除其他方法、更强估计或不同的构造能处理端点；也不表示直接设 \(\nu=2\) 就能自动得到任意常数 \(c\) 的坏逼近结论。

原稿以解析插值和行列式上下界组织反证，不是数位波模型。把其估计翻译为面积不等式只说明两种表述之间的数学联系，不是识别两套证明。

## 5. 待研究

怎样在 \(q^2E_{q+1}\) 上控制次幂级的趋零，或者反过来构造趋零子列，仍是独立的难点。原稿的非有效阈值无法直接检验某个具体分母，也不能证明 π 的数位“随机”或“无规律”。有限数值实验只能寻找结构，不能决定无穷多个分母的行为。


## 6. 将指数下界迁移到进制空间

取 \(x=\pi-3\in(0,1)\)。对实数进制 \(r>1\)，令

\[
F(r)=\min_{n\in\mathbb Z}\left|x-\frac{n}{r-1}\right|
=\frac{\operatorname{dist}((r-1)x,\mathbb Z)}{r-1}.
\]

因为取整格距离函数连续，\(F\) 在 \(r>1\) 连续。它表示目标值到一组水平数位线的最小**数值距离**；它不是第 2 节的非负波面积 \(E_r\)。对于 \(r\ge8\)，使最小值成立的整数 \(n\) 至少为 1，且小于 \(r-1\)，所以是可用的数位高度。这里不需要假设目标 \(x\) 的实进制贪心数位波对 \(r\) 连续。

对整数 \(n\ge1\)，无限同数位展开的值为 \(n/(r-1)\)。它与目标 \(x\) 相等时，交点为

\[
r_n=1+\frac{n}{x},\qquad
\mathcal R_x=\left\{r_n:n=1,2,\ldots\right\}.
\]

\(F(r_n)=0\)。因此对所有充分大的实数进制成立的正下界不可能存在；以下只讨论整数进制的采样值。

设整数 \(q\ge7\)，进制 \(r=q+1\)，并记从该进制点到最近交点的轴向距离为

\[
d_q=\operatorname{dist}(q+1,\mathcal R_x).
\]

对每个有效的横线高度 \(n\)，有精确恒等式

\[
\left|x-\frac nq\right|
=\frac{x}{q}\left|(q+1)-r_n\right|.
\]

在两边分别取最小值（最近 \(n\) 确实是有效数位）便得到

\[
D_q=F(q+1)=\frac{x}{q}d_q,\qquad
q^2D_q=xq\,d_q.
\]

所以，接受第 1 节所述主结论后，对每个固定 \(b>0\)，存在依赖 \(b\) 的阈值 \(Q(b)\)，使得所有 \(q\ge Q(b)\) 满足

\[
d_q\ge\frac{1}{x\,q^{1+b}}.
\]

而 π 是坏逼近数等价于存在与 \(q\) 无关的 \(C>0\)，使所有足够大的整数 \(q\) 满足

\[
d_q\ge\frac{C}{xq},
\quad\text{等价于}\quad
\inf_{q\ge7}q\,d_q>0.
\]

两者的差别恰好是把规范化间距 \(q d_q\) 的下界从 \(q^{-b}/x\) 加强为固定正数。可以把 \(q^{-b}\) 写成会随 \(q\) 变化的系数，但不能把固定 \(b\) 的结果直接应用到 \(b=b(q)\)；阈值 \(Q(b)\) 没有给出关于 \(b\to0\) 的一致控制。

例如，\(q=7\) 时最近高度为 \(n=1\)，交点 \(r_1\approx8.0625133059\)，离整数进制 8 约 \(0.0625133059\)；\(q=113\) 时最近高度为 \(n=16\)，交点 \(r_{16}\approx114.0002128949\)，离整数进制 114 约 \(0.0002128949\)。后一个对应 \(355/113\)；这些近似数值只用来说明公式，不决定无穷情形。


## 7. 整个无理指数结论的交点序列表达

仍取 \(x=\pi-3\)、\(r_n=1+n/x\)，令

\[
\Delta_n=\operatorname{dist}(r_n,\mathbb Z)
=\left\|\frac{n}{x}\right\|,\qquad n\ge1,
\]

其中 \(\|t\|\) 是 \(t\) 到最近整数的距离。\(\Delta_n\) 是第 \(n\) 条同数位横线与目标相交的位置，距离最近整数进制的距离。

对任意非零无理数 \(x\)，倒数不改变无理指数；整数平移也不改变无理指数。一个直接的理由是：当既约有理数 \(a/q\to x\ne0\) 时，倒数 \(q/a\to1/x\)，其分母 \(a\asymp q\)，且

\[
\left|\frac1x-\frac qa\right|
=\frac{|x-a/q|}{|x|\cdot|a/q|}.
\]

于是两端误差只差趋于固定正数的因子，分母也只差固定比例；反向取倒数同理。常数可以用任意小的指数余量吸收。因此，OpenAI 原稿的 \(\mu(\pi)=2\) 等价于 \(\mu(1/x)=2\)。

据此，主结论在该交点列中的完整表达是：

1. 对每个固定 \(\varepsilon>0\)，存在 \(N(\varepsilon)\)，使所有 \(n\ge N(\varepsilon)\) 满足
   \[
   \Delta_n\ge n^{-1-\varepsilon}.
   \]
2. 存在无穷多个 \(n\) 使 \(\Delta_n<1/n\)。这是对无理数 \(1/x\) 的 Dirichlet 逼近结论，用抽屉原理可得；不依赖 OpenAI 的估计。

两者合起来等价于

\[
\boxed{\displaystyle
\limsup_{n\to\infty}\frac{-\log\Delta_n}{\log n}=1.}
\]

这里的上极限指交点接近整数进制的**最优幂指数**，不是说每个 \(n\) 的距离都大约等于 \(1/n\)。该上极限式与 \(\mu(\pi)=2\) 互相等价；它是已有定理的进制空间改写，而不是新的 π 逼近定理。所说的主结论及其形式化与独立检查状态见第 1 节。

π 是坏逼近数则等价于存在统一的 \(c>0\)，使所有 \(n\ge1\) 满足

\[
\Delta_n\ge\frac{c}{n}
\qquad\Longleftrightarrow\qquad
\inf_{n\ge1}n\Delta_n>0.
\]

理由是坏逼近性在整数平移和取非零倒数下保持，而它对 \(1/x\) 恰好等价于 \(n\|n/x\|\) 的正下界。或者直接令 \(q\) 为离 \(n/x\) 最近的整数，则 \(q\sim n/x\)，且 \(\Delta_n=(q/x)|x-n/q|\)；反向也利用近似时的 \(n\sim xq\)。有限多个小指标不影响正下界，因为 \(x\) 无理。

指数等于 2 仍允许 \(n\Delta_n\) 沿子列趋于零；它只排除该量在无限多指标上以某个固定负幂的速度趋零。比如 \(1/\log n\) 是说明逻辑差距的可能尺度，不是对 π 的实际断言。\(\varepsilon\) 的阈值可以依赖 \(\varepsilon\)，因此不能令 \(\varepsilon=\varepsilon(n)\) 而直接推得统一常数。

数值例子：\(n=16\) 时 \(r_{16}\approx114.0002128948967\)，\(\Delta_{16}\approx0.0002128948967\)，故 \(16\Delta_{16}\approx0.00340631835\)。单个样本不决定上述上极限或下确界。


## 8. 坏逼近等价于阈值的统一增长界

沿用第 7 节的 \(\Delta_n=\|n/x\|\)，此处假定 \(x\in(0,1)\) 无理，且 \(\mu(x)=2\)。对每个 \(\varepsilon>0\)，定义**最小有效阈值**

\[
N_x(\varepsilon)=\min\left\{N\ge2:
\forall n\ge N,\quad \Delta_n\ge n^{-1-\varepsilon}\right\}.
\]

因为 \(\mu(1/x)=\mu(x)=2\)，每个这样的阈值都有限。它是进制空间的阈值，不直接等同于 OpenAI 原文中分母的 \(Q(\nu)\)。

记
\[
\ell_x=\liminf_{n\to\infty}n\Delta_n.
\]
Dirichlet 逼近保证 \(0\le\ell_x\le1\)。实际上有以下严格恒等式：

\[
\boxed{\displaystyle
\limsup_{\varepsilon\downarrow0}
 \varepsilon\log N_x(\varepsilon)
 =-\log\ell_x}
\]
其中 \(-\log0\) 按 \(+\infty\) 解释。于是 \(x\) 是坏逼近数，当且仅当存在 \(K,\varepsilon_0>0\)，使
\[
N_x(\varepsilon)\le \exp(K/\varepsilon)
\qquad(0<\varepsilon<\varepsilon_0).
\]

**证明。** 令 \(A_n=n\Delta_n\)。若 \(\ell_x>0\)，任取 \(0<t<\ell_x\)，充分大的 \(n\) 满足 \(A_n\ge t\)；只要 \(n\ge t^{-1/\varepsilon}\)，就有 \(n^{-\varepsilon}\le t\)。因此
\[
N_x(\varepsilon)\le
\max\!\left(n_0,\left\lceil t^{-1/\varepsilon}\right\rceil\right),
\quad
\limsup_{\varepsilon\downarrow0}\varepsilon\log N_x(\varepsilon)
\le-\log t.
\]
令 \(t\uparrow\ell_x\) 得到上界。

如果 \(\ell_x<1\)，任取 \(\ell_x<t<1\)，存在无穷多个 \(n\) 满足 \(A_n<t\)。对这些 \(n\) 取
\[
\varepsilon_n=\frac{-\log t}{\log n},
\]
则 \(A_n<n^{-\varepsilon_n}\)，所以 \(N_x(\varepsilon_n)>n\)，从而
\[
\varepsilon_n\log N_x(\varepsilon_n)>-\log t.
\]
令 \(t\downarrow\ell_x\) 得到下界；若 \(\ell_x=0\)，下界可以任意大。若 \(\ell_x=1\)，上界为 0，且左侧非负，故也成立。有限多个小 \(n\) 的 \(\Delta_n\) 均为正，因此 \(\ell_x>0\) 等价于坏逼近。

不使用上面的精确恒等式，也可以直接看到关键量词：若对所有小 \(\varepsilon\) 有 \(N_x(\varepsilon)\le e^{K/\varepsilon}\)，对充分大的 \(n\) 取 \(\varepsilon=K/\log n\)，则 \(n\ge N_x(\varepsilon)\) 并有
\[
n\Delta_n\ge n^{-\varepsilon}=e^{-K}.
\]
反过来，若 \(n\Delta_n\ge c>0\)，选择任何 \(K>\log(1/c)\)，对足够小的 \(\varepsilon\) 即有 \(N_x(\varepsilon)\le e^{K/\varepsilon}\)；向上取整只影响常数余量。

**对 π 的含义。** OpenAI 主结论经第 7 节的变量转换，给出每个 \(N_{\pi-3}(\varepsilon)\) 有限，但没有给出这个统一增长界。原文明确说其分母阈值没有有效给出。原文参数证明对固定 \(\nu=2+\varepsilon\) 选取 \(1/2<t<1-1/\nu\)，可用区间宽度为 \(\varepsilon/[2(2+\varepsilon)]\)，在端点趋零。这是该参数构造的障碍，并不证明其他改进不可能。

仅知 \(\mu=2\) 不足以推出统一增长界。确切反例：取 \(z=[0;1,2,3,\ldots]\) 且 \(x=1/(1+z)\in(0,1)\)，则进制交点满足 \(\Delta_n=\|nz\|\)。由连分数公式 \(\mu(z)=2\)，但其收敛分数分母 \(q_k\) 满足 \(q_k\|q_kz\|<1/(k+1)\)，所以 \(\ell_x=0\)。因此所有单独的阈值仍然有限，但 \(\limsup_{\varepsilon\downarrow0}\varepsilon\log N_x(\varepsilon)=\infty\)。

## 9. 对全部非平方整数的平方根给出固定正下界

令 \(D\ge2\) 为非平方整数，\(m=\lfloor\sqrt D\rfloor\)，\(s=D-m^2\ge1\)，目标小数部分为
\[
x_D=\sqrt D-m\in(0,1).
\]
其横线交点仍是 \(r_n=1+n/x_D\)。由于
\[
\alpha=\frac1{x_D}=\frac{\sqrt D+m}{s},\qquad
\alpha'=\frac{m-\sqrt D}{s},
\]
有 \(\Delta_n=\|n\alpha\|\)。取 \(k\) 为最接近 \(n\alpha\) 的整数。平方根无理，故以下整数非零：
\[
J=s k^2-2mkn-n^2
 =s(k-n\alpha)(k-n\alpha').
\]
利用 \(|k-n\alpha|\le1/2\) 和 \(\alpha-\alpha'=2\sqrt D/s\)，得到
\[
1\le |J|
\le \left(2\sqrt D\,n+\frac s2\right)\Delta_n,
\]
故
\[
\boxed{\displaystyle
\Delta_n\ge \frac{1}{2\sqrt D\,n+s/2}
\ge \frac{1}{(2\sqrt D+s/2)n}.}
\]

这给全部非平方整数 \(D\) 一个明确的统一常数；它依赖于目标 \(D\)，与 \(n\) 无关。特别地，\(D=2\) 时 \(m=s=1\)，
\[
r_n=1+n(\sqrt2+1),\quad
\Delta_n=\|n\sqrt2\|
\ge\frac1{(2\sqrt2+1/2)n}>\frac1{4n}.
\]
因此平方根的交点序列满足坏逼近，并且直接给出
\[
N_{x_D}(\varepsilon)\le
\max\!\left(2,\left\lceil(2\sqrt D+s/2)^{1/\varepsilon}\right\rceil\right).
\]

这一步使用了平方根满足整数二次方程形成的非零整数 \(J\)。对 π 的交点，不能仅把它写成这种平方根形式就得到同一个整数证书；还需要新的算术信息。Python 以 70 位十进制精度检查了 138,000 组非平方 \(D\in[2,150]\)、\(1\le n\le1000\) 的公式，数值检查不代替上面的证明。以上转换与估计不宣称文献首创，也未做 Lean 核验。


## 10. 原论文整组参数在指数 2 的障碍

这一节只分析原文第 4 节的参数不等式；它不证明其他论证或其他形式的行列式估计不能处理端点。按[原文参数引理](https://github.com/openai/math/blob/main/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/conclusion.tex)，参数满足
\[
0<\theta<A<B<1,\qquad
B<1/C,\quad B<C\theta,\quad
\nu(A-\theta)>1-\theta.
\]
前两条 \(B\) 上界相乘得到 \(B^2<\theta\)，所以 \(A^2<\theta\)。若令 \(\nu=2\)，最后一条却要求
\[
A>\frac{1+\theta}{2}\ge\sqrt\theta,
\]
这与 \(A<\sqrt\theta\) 矛盾。因此并非只是原文某一个参数选取区间在 \(\nu=2\) 变空：**整组现有严格条件在该端点不可行**。

还能定量限制该系统的余量。令 \(\nu=2+\varepsilon\)，原文的正余量为
\[
g=\nu\bigl(A(1-\eta)-\theta\bigr)-(1-\theta)>0.
\]
由 \(A<\sqrt\theta\) 和 \(\eta>0\)，放 \(s=\sqrt\theta\) 得
\[
g<(1-s)\bigl((1+\varepsilon)s-1\bigr)
\le\frac{\varepsilon^2}{4(1+\varepsilon)}.
\]
最后一步是配方：右边与中间表达式之差等于
\[
(1+\varepsilon)
\left(s-\frac{2+\varepsilon}{2(1+\varepsilon)}\right)^2\ge0.
\]
因此该比较系统中的正余量至多为 \(O(\varepsilon^2)\)，无法在 \(\varepsilon\downarrow0\) 时维持固定正值。这本身不提供第 8 节的 \(N_x(\varepsilon)\) 下界，也不说明指数式的阈值控制不可能；要处理坏逼近，需要新的量化控制或改变现有比较条件。


## 11. 阅读实际插值与行列式构造后的研究方向（2026-10-28）

阅读对象：[interpolation.tex](https://github.com/openai/math/blob/main/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/interpolation.tex)、[determinant.tex](https://github.com/openai/math/blob/main/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/determinant.tex) 以及原文参数部分。本节区分原文构造、可以直接推导的细化和仍未证明的目标；不宣称完成坏逼近证明。

### 11.1 实际构造及进制交点的接口

原文同时选取若干有理逼近 \(a_i/q_i\)。在 \(K\) 个中心处，以加权多项式为列、指定的 Taylor 系数为行，构造矩阵。插值满射保证存在使用所有行的非零方形子式。清除分母后，该子式变为非零高斯整数，所以模至少为 1；这是算术下界的来源。

解析上，把中心 \(2\mathrm i j a_i/q_i\) 平移到精确周期 \(2\pi\mathrm i j\)，利用 \(e^{2\pi\mathrm i j}=1\)。这给出的不是独立的统计样本，而是同一指数函数周期上的多点算术联系。

与进制模型的接口是精确的。令 \(x=\pi-3\)，选择同数位高度 \(n_i\)，令 \(q_i+1\) 为其交点 \(r_{n_i}=1+n_i/x\) 最邻近的整数进制，设 \(a_i=3q_i+n_i\)。则
\[
\frac{a_i}{q_i}-\pi
=\frac{x}{q_i}\bigl(r_{n_i}-(q_i+1)\bigr).
\]
因此原文中心的横向位移恰好为
\[
2\mathrm i j\left(\frac{a_i}{q_i}-\pi\right)
=\frac{2\mathrm i jx}{q_i}\bigl(r_{n_i}-(q_i+1)\bigr).
\]
可以把若干进制交点直接带入原构造；构造本身仍依赖指数函数的精确周期，不能只依赖数位波的外观。

### 11.2 原文的两种解析收益

平移后的行按横向多重指标 \(a\) 分组。每组内部进一步展开成 Taylor 级数；若同一组的两行选到同一 Taylor 次数，该项的两条系数向量相同，行列式严格为零。若组内有 \(N_a\) 行，非零项的次数必须互异，次数之和至少为 \(\binom{N_a}{2}\)。原文据此得到形如
\[
\exp\left(-c\sum_aN_a^2\right),\qquad c=(\log2)/4
\]
的压小因素。它是精确的行列式抵消，不是数位随机性的假设。

另一项收益来自逼近误差的幂次。最终比较分成两种情况：低横向指标多，就用上面的组内抵消；高指标多，就用逼近误差的幂次。原文对每种情况各舍弃另一项非负收益。

### 11.3 保留坏逼近问题真正关心的小系数

若要反证坏逼近失败，假设存在逼近序列使
\[
\kappa_i=q_i^2\left|\pi-\frac{a_i}{q_i}\right|\longrightarrow0.
\]
在进制交点附近，\(\kappa_i=xq_i\Delta_{n_i}\sim n_i\Delta_{n_i}\)。不能把它默认改成固定的 \(q_i^{-\varepsilon}\)：无理指数等于 2 时，沿任何 \(\kappa_i\to0\) 的无界分母序列，都有
\[
\frac{\log(1/\kappa_i)}{\log q_i}\longrightarrow0.
\]

不过原文的精确行平移恒等式可以保留这些系数。令 \(w_i=\lceil\log q_i\rceil\)，并取 \(0<\kappa_i\le1\)。在原文行平移中，行的原多重指标为 \(\beta\)，选取的横向指标为 \(a\ge\beta\)，Taylor 截断尾项的次数多重指标为 \(d\)，且 \(0\le d\le a-\beta\)。原误差常数按 \(\nu=2\) 代入记为 \(E_{\rm tr}(2)\)。逐项重复原文估计，标量因子 \(\xi\) 满足
\[
|\xi|\le
 \exp\{-2w(a-\beta)+H E_{\rm tr}(2)\}
 \prod_i\kappa_i^{\,a_i-\beta_i-d_i}.
\]
理由是精确位移满足
\[
|2\mathrm i j(a_i/q_i-\pi)|\le2K\kappa_iq_i^{-2},
\]
而尾项预算、组合系数界及 Cauchy 界与原证明相同。这只是行平移引理在实际误差系数下的细化，不是把原文的端点比较定理直接宣布成立。

这说明应研究“小系数收益”能否在指数 2 处抵消剩余损失，而不是在一开始把它转换为一个固定的正指数余量。该问题仍未闭合。

### 11.4 两个可严格保留的平均量

原文行集合是
\[
v_0s+w\beta/\theta<H.
\]
固定维数 \(m\)、权重及中心，令 \(H\to\infty\)。单纯形上的一阶矩直接给出
\[
\bar b_H=\frac{1}{M_HH}\sum_{\mathrm{rows}}w\beta
\longrightarrow\frac{m\theta}{m+2}.
\]
这比原文比较时使用的 \(\bar b_H\le\theta\) 更精确。每一个尺度坐标在该 \(m+1\) 维单纯形上的均值为预算的 \(1/(m+2)\)；其中 \(m\) 个横向坐标各带比例 \(\theta\)。有限枚举检查了 \(m=1,2,3\) 的收敛方向；该公式的证明来自固定数据下的 Riemann 和，不能与 \(m\to\infty\) 的极限任意交换。

对选中的列，另外记实际平均横向代价
\[
\bar c_H=\frac{1}{M_HH}\sum_{\mathrm{selected\ columns}}w\alpha\le1.
\]
沿原文清除分母的计算可保留
\[
\frac{\log|\Delta_H|}{M_HH}
\ge-\bar c_H+\bar b_H-E_{\rm ar},
\]
而不必立即用 1 替代 \(\bar c_H\)。不需额外的选列定理，所有所选列都满足一个严格小于 1 的通用平均代价上界，见第 12 节。若想进一步改善该通用上界，还需利用子式非零带来的结构或证明更强的选列性质。

### 11.5 一个明确但尚未证明的端点比较目标

设每条行展开中选择一项，用 \(\rho\) 标记行，\(a_\rho,\beta_\rho,d_\rho\) 为该项的指标；令 \(\lambda_i=\log(1/\kappa_i)\ge0\)。把两种解析收益及小系数全部保留，可定义
\[
S_H=
\frac{
2\sum_\rho w(a_\rho-\beta_\rho)
+\sum_{\rho,i}\lambda_i(a_{\rho i}-\beta_{\rho i}-d_{\rho i})
+c\sum_aN_a^2
}{M_HH}.
\]
原文同样的有限项展开、绝对值估计和组内抵消给出
\[
\frac{\log|\Delta_H|}{M_HH}
\le E_{\rm an}(2)+o(1)-\inf_{\text{nonzero term choices}}S_H.
\]
所以一个足以导出矛盾的目标是
\[
\liminf_{H\to\infty}
\left(\inf S_H+\bar b_H-\bar c_H\right)
>E_{\rm ar}+E_{\rm an}(2).
\]
参数、分母和中心须先固定，才能取此 \(H\) 极限。如何在仅假设 \(\kappa_i\to0\) 时选择数据并证明这个严格比较，目前没有完成；它可能需要新的插值控制、列选择定理或全局分组估计。这里只给出从原文真实公式提炼出的研究靶子，不把它算作已获得的 π 下界。


## 12. 临界指数 2 的列代价改善与一维障碍（2026-10-28）

本节继续第 11.5 节的比较。所有由论文提供的插值和行列式估计，都按原文的量词顺序使用：先固定维数、参数、权重和中心，再令 \(H\to\infty\)。本节给出的改善不证明 π 的坏逼近。

### 12.1 保留真实系数的收益至多是次幂级

沿用 \(\kappa_i=q_i^2|\pi-a_i/q_i|\le1\)、\(\lambda_i=\log(1/\kappa_i)\) 和 \(w_i=\lceil\log q_i\rceil\)。令
\[
\tau=\max_i \lambda_i/w_i.
\]
对第 11.5 节中每个非零行展开项，\(0\le d_{\rho i}\le a_{\rho i}-\beta_{\rho i}\)，且 \(w a_\rho\le H\)。因此对任意项组合，额外系数收益满足严格的一致界
\[
0\le
\frac{1}{M_HH}
 \sum_{\rho,i}\lambda_i(a_{\rho i}-\beta_{\rho i}-d_{\rho i})
\le\tau.
\]
若使用 OpenAI 主结论，对每个固定 \(\varepsilon>0\)，足够大的分母满足 \(\kappa_i\ge q_i^{-\varepsilon}\)，所以 \(\tau\le\varepsilon\)。因此当所有分母趋于无穷时，这个归一化收益为 \(o(1)\)；对选取有限或变化数量的分母，只要它们的最小值趋于无穷，上界仍成立。

具体含义：如果其他主项已经留下与分母无关的固定负差额，仅让 \(\kappa_i\to0\) 不足以用这个归一化收益跨过它。若能把主项配平或把剩余损失降到同一次幂尺度，\(o(1)\) 收益仍可能有用；这里不判断可否做到，也不交换 \(m\)、权重和 \(H\) 的极限。

### 12.2 不依赖选列定理的通用列代价上界

设列集合为 \(w_0h+w\alpha\le H\)，行数为 \(M_H\)，总列数为 \(P_H\)。第 11.4 节的 \(\bar c_H\) 是任意选取的 \(M_H\) 个不同列的平均横向代价。令
\[
\rho=\lim_{H\to\infty}\frac{M_H}{P_H}
=\frac{K w_0\theta^m}{v_0}\in(0,1).
\]
原文的插值体积条件确保 \(\rho<1\)。

在列单纯形中，横向代价 \(t=w\alpha/H\in[0,1]\) 的极限密度是
\[
f_m(t)=m(m+1)t^{m-1}(1-t).
\]
令 \(t_\rho\in(0,1)\) 唯一满足
\[
\rho=\int_{t_\rho}^1f_m(t)\,dt
=1-(m+1)t_\rho^m+mt_\rho^{m+1}.
\]
在所有可选的 \(M_H\) 个列中，最高代价部分的平均值最大；据此，**任意**这样的列选择均满足
\[
\limsup_{H\to\infty}\bar c_H\le
U_m(\rho):=
\frac{m}{(m+2)\rho}
\left[1-(m+2)t_\rho^{m+1}+(m+1)t_\rho^{m+2}\right]
<1.
\]
证明只用固定维数、固定权重的单纯形 Riemann 和与按代价排序，甚至不需要所选子式非零。这修正第 11 节中“需要先证明选列性质才能获得任何严格改善”的过强表述；如果要比 \(U_m(\rho)\) 更强，仍需进一步利用秩或选择结构。

当 \(m=1\)、\(\rho=1/2\) 时，\(f_1(t)=2(1-t)\)，直接积分得到
\[
U_1(1/2)=1-\frac{\sqrt2}{3}\approx0.52859548,
\qquad
\bar b_H\longrightarrow\theta/3.
\]

### 12.3 一维高指标分支出现正主项，但原碰撞分支仍不闭合

作为可核对的临界指数例子，取
\[
m=1,\quad K=2,\quad
\theta=2/5,\quad A=3/5,\quad B=31/50,\quad
C=39/25,\quad\eta=1/10,
\]
以及任意正有理 \(w_0\)，令 \(v_0=2K\theta w_0\)。这些数据满足插值的两个体积条件 \(K\theta<1\)、\(Kw_0\theta/v_0=1/2<1\)，并满足 \(0<\theta<A<B<1\)、\(B<1/C\)、\(B<C\theta<1\)。横向权重仍需另行选得足够大；本例没有声称选择完成了论文的全部条件。

在指数 \(\nu=2\) 的高指标分支，利用第 11.4、12.2 节，而只看未扣误差的主项，可以保证的差额为
\[
2A(1-\eta)-U_1(1/2)-\theta/3
=\frac{\sqrt2}{3}-\frac4{75}>0.418.
\]
例如原文直接用 \(1\) 和 \(\theta\) 时，同一组参数的主项差额为 \(-0.32\)。因此上一轮发现的平均值与列代价改善，确实消除了**这一分支主项**的符号障碍；这不意味着该分支扣除全部误差后已经成立。

另一分支仍不行。原文碰撞常数为 \(c=(\log2)/4\)，且在本节 \(\rho=1/2\) 下
\[
L=\frac{\eta^2}{4w_0A}.
\]
上述具体参数给出
\[
cL=\frac{\log2}{960w_0}
\qquad\text{但}\qquad
E_{\rm hol}\ge\frac{100K}{w_0}=\frac{200}{w_0}.
\]
所以若低指标分支仍沿原文做法只保留 \(cL\) 而舍弃其他收益，即使取很大的 \(w_0\) 也不可能用这组参数满足比较条件。

这一障碍在一维、半数列的通用上界方法下并非只针对这组数字：若高指标分支的**通用主项下界**为正，则 \(2A(1-\eta)>U_1(1/2)+\theta/3\)，所以 \(A>U_1(1/2)/2\)。于是对任何 \(\eta<1\)、\(K\ge1\)，
\[
cL=\frac{c\eta^2}{4w_0A}
<\frac{c}{2U_1(1/2)w_0}
<\frac{100K}{w_0}
\le E_{\rm hol}.
\]
因此，**在 \(m=1,\rho=1/2\)，使用这个通用列上界、原文两分支单独估计和原文误差上界时，不可能同时完成两分支的严格比较。** 这不排除利用实际非零子式的更小代价、联合保留 \(\kappa_i\) 收益、改进碰撞或解析误差估计，也不排除更高维方案。

### 12.4 有限检查与状态

对 \(m=1,K=2,\theta=2/5,w_0=w_1=1,v_0=8/5\)，直接精确枚举行指标并排序列代价：\(H=100,500,2000,10000\) 时，所选列占比分别约为 \(0.504756,0.500990,0.500249,0.500050\)；最大平均列代价依次约为 \(0.529215,0.528712,0.528624,0.528601\)，趋向 \(1-\sqrt2/3\)。这只检查格点计数，不验证这些有限 \(H\) 的插值或 π 逼近。一般 \(m=1,2,3\) 的排序与积分公式也作过有限检查。新推导尚未经 Lean 核验。


## 13. 半数列、原两分支估计的全维数障碍（2026-10-28）

这一节将第 12.3 节的一维障碍推广到每个固定维数。结论仅涉及以下**具体估计方案**：保留准确行均值和第 12.2 节对任意所选列的通用上界；保留论文对 Taylor 碰撞与解析误差的原始常数；仍按论文的两个分支分别舍弃另一种收益；不使用第 11 节保留的 \(\kappa_i\) 额外因子。它不排除更精细的非零子式选择、联合估计、其他列比例、改进误差界或其他证明路线。

### 13.1 比较条件

固定整数 \(m\ge1\)、\(K\ge1\)，\(0<\theta,A,\eta<1\)，并按原文使用
\[
v_0=2K\theta^m w_0,\qquad w_0>0.
\]
因此行数与总列数的极限比例为 \(1/2\)。记
\[
v=\frac{m}{m+2},\qquad b=v\theta,\qquad U=U_m(1/2),
\qquad c=\frac{\log2}{4}.
\]
第 11.4 节给出 \(\bar b_H\to b\)，第 12.2 节给出 \(\limsup\bar c_H\le U\)。若以这些通用界分别完成指数 2 的高指标、低指标分支，足够的严格比较条件为
\[
2A(1-\eta)-U-b>E,\qquad
cL>U-b+E,
\tag{13.1}
\]
其中
\[
L=\frac{\eta^2}{2(m+1)w_0A^m},\qquad
E=E_{\rm ar}+E_{\rm an}(2).
\]
原文误差公式给出
\[
E\ge E_{\rm hol}\ge\frac{100K}{w_0}+\frac{\log2}{v_0}.
\tag{13.2}
\]
这里使用 \(U-b\) 是为了以对任意子式都成立的列代价上界保证比较；并不声称具体非零子式的实际代价等于 \(U\)。

**结论：无论固定维数 \(m\) 多大，(13.1) 与 (13.2) 都不能同时成立。**

### 13.2 证明

令 \(T\) 为 \(\mathrm{Beta}(m,2)\) 分布，其密度是第 12.2 节的 \(m(m+1)t^{m-1}(1-t)\)。令 \(t_{1/2}\) 为其中位数，则
\[
U=\mathbb E[T\mid T\ge t_{1/2}],\qquad v=\mathbb E T.
\]
上半与下半各占一半，故
\[
U-v=\mathbb E|T-t_{1/2}|
\ge\mathbb E(T-t_{1/2})^2
\ge\operatorname{Var}(T)
=\frac{2m}{(m+2)^2(m+3)}.
\tag{13.3}
\]
第一项比较使用 \(0\le|T-t_{1/2}|\le1\)，第二项使用均值最小化平方误差。这也表明 \(U>v\) 和 \(U-b\ge U-v>0\)。

由 (13.1) 第一项，\(2A(1-\eta)>U+v\theta\)，故
\[
A>\frac{U+v\theta}{2}.
\]
因为 \(U\ge v\) 且 \(0<\theta<1\)，有
\[
\frac{\theta}{A^2}
<\frac{4\theta}{(U+v\theta)^2}
\le\frac4{(U+v)^2}
\le\frac1{v^2}
=\left(1+\frac2m\right)^2.
\]
中间一步可令 \(s=\sqrt\theta\)，观察
\(U+v\theta-s(U+v)=(1-s)(U-vs)\ge0\)。于是
\[
(\theta/A^2)^m<e^4.
\tag{13.4}
\]

由 (13.1) 第二项分别比较 \(U-b\) 与 \(100K/w_0\)，得到
\[
w_0<\frac{c\eta^2}{2(m+1)(U-b)A^m},
\qquad
K<\frac{c\eta^2}{200(m+1)A^m}.
\]
将其代入 \(v_0=2K\theta^m w_0\)，再用 (13.3)、(13.4) 和 \(\eta<1\)，有
\[
\begin{aligned}
v_0
&<\frac{c^2}{200(m+1)^2(U-b)}
      \left(\frac{\theta}{A^2}\right)^m\\
&\le
\frac{c^2e^4(m+2)^2(m+3)}
     {400m(m+1)^2}
\le\frac{9c^2e^4}{400}
<0.037.
\end{aligned}
\]
最后一个有理因子不超过 9，因为
\(((m+2)/(m+1))^2\le(3/2)^2\) 且 \((m+3)/m\le4\)。

因此 \(E\ge\log2/v_0>18\)。但 (13.1) 第一项的左侧小于 2，矛盾。

### 13.3 这对下一步的含义

第 12 节的准确平均值和通用列代价上界确实改善了高指标主项；以上证明表明，仅把它们代入论文原有的**分开估计与原误差常数**，在半数列设置下仍不能闭合，且单纯增加固定维数没有用。

真正需要检查的是至少一处不同的信息：

1. 证明所选非零子式的**实际**平均列代价明显低于对所有列集合的最坏上界；
2. 如第 11.5 节，联合保留碰撞、横向指标以及实际 \(\kappa_i\) 的收益，而不是分别丢弃；
3. 证明适用于该具体子式与中心的更强解析误差或碰撞界；
4. 改动原文选定的 \(M_H/P_H\to1/2\) 比例，并重新推导全套插值与比较条件。

以上四项都是研究方向，不是已经得到的定理。数值计算检查了 \(m=1,2,3,10,100,1000,10000\) 的 (13.3) 并与有限格点枚举相符；证明本身是上述解析推导，尚未经 Lean 核验。


## 14. 非零子式强制的列代价：一个有限阶下界（2026-10-28）

这一节进一步使用原矩阵的实际秩条件。它给的是**列代价的必要下界**，不是第 12 节为改善算术估计而需要的上界，也不证明 π 的逼近性质。

固定原文行指标 \(\rho=(j,s,\beta)\)，其中 \(0\le j<K\)；固定选出的 \(M=M_H\) 个不同列 \(Y^hX^\alpha\)。以原文的截断对数多项式 \(G_i(t)\) 构造矩阵，但将各个中心 \(jr_i\) 同时替换为 \(jzr_i\)，得到关于辅助变量 \(z\) 的方阵 \(M(z)\)。假设其在 \(z=1\) 的子式非零。记
\[
D(z)=\det M(z),\qquad
D_{\rm col}=\sum_{\rm columns}|\alpha|-\sum_{\rm rows}|\beta|.
\]

**结论。**
\[
D_{\rm col}\ge\frac{M(K-1)}2,\qquad
\bar c_H-\bar b_H\ge\frac{w_*(K-1)}{2H}.
\tag{14.1}
\]

**证明。** 按固定的 \((s,\beta)\) 将行分组。每组恰有 \(K\) 行，分别对应 \(j=0,\ldots,K-1\)。由于中心只以 \(jzr_i\) 进入，固定一组的第 \(j\) 行可以写成有限和
\[
\sum_{\ell\ge0}(jz)^\ell V_{s,\beta,\ell},
\]
其中行向量 \(V_{s,\beta,\ell}\) 与 \(j\) 无关。对行列式逐行展开，若同组两行选择相同的 \(\ell\)，对应向量成比例，该项为零。故每个非零项在每组必须使用互异的 \(K\) 个非负次数，其和至少是 \(K(K-1)/2\)。一共有 \(M/K\) 组，因此 \(D(z)\) 在 \(z=0\) 的消失阶至少为 \(M(K-1)/2\)。

另一方面，矩阵元素在 \(\alpha\not\ge\beta\) 时为零，在 \(\alpha\ge\beta\) 时是关于 \(z\) 的至多 \(|\alpha|-|\beta|\) 次多项式。每个行列式项都恰好用每一行和每一列一次，所以
\[
\deg D(z)\le
\sum_{\rm columns}|\alpha|-\sum_{\rm rows}|\beta|=D_{\rm col}.
\]
因为 \(D(1)\ne0\)，它是非零多项式，次数至少是其在零点的消失阶，得到第一个不等式。

对任意非零行列式展开项，其匹配均满足逐坐标 \(\alpha\ge\beta\)。又 \(w_i\ge w_*=\min_iw_i>0\)，故
\[
\sum_{\rm columns}w\alpha-\sum_{\rm rows}w\beta
\ge w_*D_{\rm col}.
\]
除以 \(MH\) 即得第二个不等式。

**尺度限制。** 论文先固定 \(K,w_*\)，再令 \(H\to\infty\)，所以 (14.1) 的归一化下界是 \(O(1/H)\)，不能单独保证一个固定的正代价间隔。若试图让 \(K\) 随 \(H\) 增大，原插值条件 \(K\theta^m<1\) 已经要求其他参数随之变化；原文的固定参数渐近式并未直接提供这种联合极限的统一控制。这指出了一个可研究的尺度问题，但没有解决它。

**有限检查。** 对 \(m=1,K=2,\theta=2/5,w_0=w_1=1,v_0=8/5\)，中心取 \(jz\)，截断多项式取 \(G(t)=t-t^2/2\)，用精确有理数枚举方形子式：\(H=3\) 时 54 个非零子式的最小 \(D_{\rm col}\) 为 5，公式要求至少 3；\(H=4\) 时 3570 个非零子式的最小值为 6，公式要求至少 4。这里仅检查给定有限矩阵，不检查论文的插值权重阈值。上述一般结论来自多项式消失阶证明，尚未经 Lean 核验。
