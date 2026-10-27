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
