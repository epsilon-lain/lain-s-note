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
而不必立即用 1 替代 \(\bar c_H\)。原文的满射只保证非零子式存在，不直接保证可选取低代价的列；若要进一步利用这个改善，需要证明带代价的列选择性质。

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
