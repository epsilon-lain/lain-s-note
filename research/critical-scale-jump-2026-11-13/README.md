# 临界 \(q^{-2}\) 尺度的分母跳跃：从有限见证到构造约束

记录日期：2026-11-13。性质：一个一般性的定量行列式引理、它对连分数的精确特化，以及对已有 runB 有限见证的有理数复核。**不是 P2848 的证明，也不是新的 π 逼近纪录。**

## 1. 本轮选择的工作

仓库当前主目标是寻找既约分数 \(p_n/q_n\)，使

\[
\varepsilon_n:=q_n^2\left|\pi-\frac{p_n}{q_n}\right|\longrightarrow0
\]

沿某个无界分母序列成立。此前的候选构造已经反复表明：只报告误差下降而不同时报告真实约分分母，无法判断是否接近目标。

本轮把“误差下降要付出多少分母增长”写成一个可直接用于筛选构造器的必要条件，并用 2026-11-13 已独立认证的 runB 有限见证作数值实例。

## 2. 定量分母跳跃引理

设 \(\alpha\in\mathbb R\)。设

\[
\frac{p_n}{q_n},\quad \frac{p_{n+1}}{q_{n+1}}
\]

是两个不同的既约有理数，且 \(0<q_n<q_{n+1}\)。定义

\[
\varepsilon_n=q_n^2\left|\alpha-\frac{p_n}{q_n}\right|,
\qquad
R_n=\frac{q_{n+1}}{q_n}>1.
\]

则

\[
\boxed{1\le \varepsilon_nR_n+\frac{\varepsilon_{n+1}}{R_n}.}
\tag{1}
\]

**证明。** 因为两个既约分数不同，

\[
D_n=p_{n+1}q_n-p_nq_{n+1}
\]

是非零整数，所以 \(|D_n|\ge1\)。另一方面，

\[
\begin{aligned}
|D_n|
&=q_nq_{n+1}
 \left|\frac{p_{n+1}}{q_{n+1}}-\frac{p_n}{q_n}\right|\\
&\le q_nq_{n+1}
 \left(
 \left|\alpha-\frac{p_n}{q_n}\right|
 +\left|\alpha-\frac{p_{n+1}}{q_{n+1}}\right|
 \right)\\
&=\varepsilon_nR_n+\frac{\varepsilon_{n+1}}{R_n}.
\end{aligned}
\]

证毕。

由于 \(R_n>1\)，当 \(\varepsilon_{n+1}<1\) 时，(1) 进一步给出

\[
\boxed{R_n\ge\frac{1-\varepsilon_{n+1}}{\varepsilon_n}.}
\tag{2}
\]

因此，若一列既约有理数满足 \(\varepsilon_n\to0\)，则

\[
R_n\to\infty,
\qquad
\varepsilon_nR_n\ge1-\varepsilon_{n+1}\to1.
\tag{3}
\]

这不是 π 特有的性质，对任意实数 \(\alpha\) 都成立。它说明：**临界 \(q^{-2}\) 归一化误差若沿一列相邻候选趋零，分母比就不可能保持有界；其增长至少约为误差的倒数。**

## 3. 连分数把“误差—跳跃”关系写成精确恒等式

令 \(\alpha=[a_0;a_1,a_2,\ldots]\)，其收敛分数为 \(p_n/q_n\)。定义完整商

\[
\alpha_{n+1}=[a_{n+1};a_{n+2},\ldots],
\qquad
\beta_n=\frac{q_{n-1}}{q_n}\in(0,1).
\]

标准收敛分数误差公式给出

\[
\varepsilon_n
=q_n^2\left|\alpha-\frac{p_n}{q_n}\right|
=\frac{1}{\alpha_{n+1}+\beta_n}.
\tag{4}
\]

分母递推则给出

\[
R_n=\frac{q_{n+1}}{q_n}
=a_{n+1}+\beta_n.
\tag{5}
\]

再写

\[
\alpha_{n+1}=a_{n+1}+\theta_n,
\qquad
\theta_n=\frac1{\alpha_{n+2}}\in(0,1).
\]

于是

\[
\boxed{
\frac1{\varepsilon_n}=R_n+\theta_n,
\qquad
\varepsilon_nR_n
=\frac{R_n}{R_n+\theta_n}.
}
\tag{6}
\]

特别地，

\[
\boxed{
1-\frac1{a_{n+1}}
<\varepsilon_nR_n<1,
\qquad
\frac1{a_{n+1}+2}<\varepsilon_n<\frac1{a_{n+1}}.
}
\tag{7}
\]

所以对连分数收敛分数而言，\(\varepsilon_n\) 小与下一部分商 \(a_{n+1}\) 大是同一件事的两种精确表达；同时 \(q_{n+1}/q_n\) 也随之变大。

## 4. 与 \(\mu(\pi)=2\) 的临界窗口

OpenAI 2026-09-24 预印本声称证明 \(\mu(\pi)=2\)。其主结论是：对每个固定 \(\delta>0\)，充分大的 \(q\) 满足

\[
\left|\pi-\frac pq\right|\ge q^{-2-\delta}.
\]

这不是 \(c/q^2\) 型统一下界。相关说明见：

- [OpenAI 预印本](https://github.com/openai/math/blob/main/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/main.tex)
- [官方 Lean 形式化范围说明](https://github.com/openai/math/blob/main/lean/docs/017.md)
- [P2848：π 的连分数部分商无界](https://www.theoremdb.org/statements/P2848/)

P2848 要求部分商无界，即存在 \(n_j\to\infty\) 使 \(a_{n_j+1}\to\infty\)。结合 (4)、(7)，它等价于收敛分数上存在 \(\varepsilon_{n_j}\to0\) 的子序列。

而标准公式

\[
\mu(\pi)
=2+\limsup_{n\to\infty}
\frac{\log a_{n+1}}{\log q_n}
\]

表明，若 \(\mu(\pi)=2\)，则

\[
\frac{\log a_{n+1}}{\log q_n}\longrightarrow0.
\]

因此两个目标若同时成立，所需的不是多项式级别的大部分商，而是

\[
\boxed{
a_{n_j+1}\longrightarrow\infty,
\qquad
a_{n_j+1}=q_{n_j}^{\,o(1)}.
}
\tag{8}
\]

等价地，下一分母比沿该子序列必须趋于无穷，但仍是 \(q_{n_j}^{o(1)}\) 级别。比如 \(a_{n_j+1}\) 像 \(\log q_{n_j}\) 那样增长，在尺度上与 \(\mu=2\) 相容；这只是一个允许的目标尺度，不是 π 已有的渐近结论。

这给构造工作一个明确约束：

1. 不能把 \(q^2|\pi-p/q|\le q^{-c}\)（固定 \(c>0\)）作为无穷成功目标；若有无穷多个这样的既约分数，就会推出 \(\mu(\pi)\ge2+c\)，与预印本主张冲突。
2. 若目标是 P2848，归一化误差必须趋零，但只能处在次幂尺度；同时分母比必须发散。
3. 对一般构造序列，验收时必须同时保存相邻的 \(q_n,q_{n+1}\) 和两项 \(\varepsilon_n,\varepsilon_{n+1}\)，并检查 (1)。只报告一个很小的误差不够。

## 5. runB 有限见证：小误差确实伴随巨大分母跳跃

[runB 独立整数认证](../yohaku-pi-cf-audit-2026-11-12/RUNB-INDEPENDENT.md) 给出：

\[
\begin{aligned}
0.0000127178573727574566035844693935801274
&<\varepsilon_{28420}\\
&<0.0000127178573727574566035844693935801275,
\end{aligned}
\]

且 \(a_{28421}=78629\)。这里 \(p_{28420}\) 已认证为 \(q_{28420}\pi\) 的最近整数，因此该 \(\varepsilon_{28420}\) 与审计文件中的 \(q_{28420}\lVert q_{28420}\pi\rVert\) 完全相同。

由 (7) 可先核对

\[
\frac1{78631}<\varepsilon_{28420}<\frac1{78629}.
\]

再由 (6)，记 \(R_{28420}=q_{28421}/q_{28420}\)，有

\[
78629<R_{28420}<\frac1{\varepsilon_{28420}}
<\frac1{0.0000127178573727574566035844693935801274}
<78629.598.
\]

所以这个有限见证的下一分母比被严格夹在

\[
\boxed{78629<R_{28420}<78629.598.}
\]

其倒数归一化误差约为 \(1/\varepsilon_{28420}\approx78629.597006\)，而分母比与它只相差 \(\theta_{28420}\in(0,1)\)。这正是 (6) 的数值体现。

这是一条很小的有限误差与一次约 \(7.86\times10^4\) 分母增长的配对。它说明跳跃成本不是抽象担忧；但它只涉及一个有限指标，**不能推出部分商无界，也不能推出尾部下极限为零。**

## 6. 精确复核

运行：

~~~sh
python3 check_critical_scale_jump.py
~~~

脚本只用 Python 标准库，从相邻目录 yohaku-pi-cf-audit-2026-11-12/runB-independent.json 读取已经独立认证的有限见证区间，并用 fractions.Fraction 做整数交叉乘法。它核对：

- 见证区间严格包含在 \((1/78631,1/78629)\) 内；
- \(1/\varepsilon_{28420}\) 的严格上界小于 78629.598；
- 由连分数恒等式推出的分母比区间 \((78629,78629.598)\)。

脚本输出保存为 results.json。脚本不重新计算 π，不重新认证 28,422 个连分数项，也不把有限见证升级成无限结论。

## 7. 未解决问题与下一步

仍未解决的核心问题是：能否构造或证明存在无穷多个指标 \(n_j\)，使

\[
a_{n_j+1}\to\infty
\quad\text{同时}\quad
\frac{\log a_{n_j+1}}{\log q_{n_j}}\to0?
\]

本轮没有产生新的 π 逼近分数，没有证明 P2848，也没有证明 π 是坏逼近数。下一轮应继续按仓库的“构造优先”路线设计实际既约分数，并把 (1)、(6)、(8) 作为必要的误差—分母联合验收条件；不能把一个有限的大部分商或单次质量纪录当作无限证明。
