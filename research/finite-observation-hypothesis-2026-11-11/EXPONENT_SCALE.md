# 观测 gcd 的增长指数与临界平方根尺度

2026-11-11。[主文](README.md)的无条件推广。这里不使用或判断任何关于 π 的最新无理指数声明。

固定无理数 α，设 N_B=round(Bα)、g_B=gcd(N_B,B)，B 为趋于无穷的正整数。普通无理指数定义为

\[
\mu(\alpha)=\sup\{\nu:\ |\alpha-p/q|<q^{-\nu}
\text{ 对无穷多个既约分数成立}\}.
\]

由 Dirichlet，μ≥2。

## 1. 无条件增长指数公式

\[
\boxed{\limsup_{B\to\infty}\frac{\log g_B}{\log B}
=1-\frac1{\mu(\alpha)}.}
\tag{1}
\]

μ=∞ 时右端约定为 1。

**上界。** 若 μ<∞，取 ν>μ。主文证明了约分分母 Q_B=B/g_B→∞。按无理指数定义，充分大 B 有

\[
Q_B^{-\nu}\le|\alpha-N_B/B|<1/(2B).
\]

故 g_B<2^{-1/ν}B^{1−1/ν}。先取上极限，再令 ν↓μ。μ=∞ 时用 g_B≤B 即可。

**下界。** 取 1<ν<μ，选取 q→∞ 的既约分数使 δ=|α−p/q|<q^−ν。令 Z=1/(2qδ)、k=floor(Z)、B=kq。Z→∞，所以 k/Z→1；又因 α 无理，Z 不是整数。于是 |Bα−kp|=kqδ<1/2，故 N_B=kp、g_B=k，并有

\[
B=(1+o(1))/(2\delta),\qquad
\frac{\log g_B}{\log B}
=1-\frac{\log q}{\log B}
\ge1-\frac1\nu-o(1).
\]

令 ν↑μ，或在 μ=∞ 时令 ν 任意增大，得到 (1)。

## 2. 指数等于一半不等于统一平方根界

若 μ=2，(1) 意味着：任给 ε>0，充分大 B 有 g_B≤B^{1/2+ε}，而 g_B≥B^{1/2−ε} 会无限次发生。不能写成对每个 B 都有 g_B=B^{1/2+o(1)}，因为 g_B 可以很小。

坏逼近需要更强的 g_B/√B 统一有界。主文的精确恒等式

\[
\limsup_B g_B^2/B=1/(2\lambda(\alpha))
\]

控制的是临界指数内的常数和无界慢增长，而 (1) 只控制幂指数。

## 3. 一个显式对照

取 α=[0;1,2,3,4,…]，第 n 项部分商为 n。收敛分数分母满足 q_n≥n!，且

\[
\frac1{n+3}<\kappa_n
=q_n^2|\alpha-p_n/q_n|<\frac1{n+1}.
\]

前者由 q_n=nq_{n−1}+q_{n−2} 归纳得到。后者由

\[
\kappa_n=\frac1{\alpha_{n+1}+q_{n-1}/q_n},
\quad n+1<\alpha_{n+1}<n+2
\]

直接推出。于是 κ_n→0，α 非坏逼近。但每个 ν>2 下，q_n^{ν−2}>n+3 最终成立，所以收敛分数误差最终大于 q_n^−ν；非收敛的既约分数由 Legendre 判据有误差≥1/(2q²)，同样最终大于 q^−ν。结合 Dirichlet 得 μ=2。

取 k_n=floor(q_n/(2κ_n))、B_n=k_nq_n，则

\[
\frac{n+1}2-\frac1{q_n}
<\frac{g_{B_n}^2}{B_n}<\frac{n+3}2.
\]

所以 g_{B_n}/√B_n 约为 sqrt(n/2)，趋于无穷，而幂指数上极限仍只有 1/2。这是一般无理数的对照，不是 π 的构造或性质结论。

连分数质量恒等式与 Legendre 判据的基础来源见 [Bugeaud, CIMPA 2010](https://irma.math.unistra.fr/~bugeaud/travaux/Bamako2.pdf)，Theorems 1.9、1.11。以上推广在正文给出证明，未作形式化认证，不宣称原创。
