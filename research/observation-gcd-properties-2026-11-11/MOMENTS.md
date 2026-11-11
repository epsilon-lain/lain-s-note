# 观测函数的矩平均与无理指数门槛

2026-11-11。[性质研究主文](README.md)的理论附录。固定无理 α、G(B)=gcd(round(Bα),B)²/B，普通无理指数 μ(α) 沿用[此前定义](../finite-observation-hypothesis-2026-11-11/EXPONENT_SCALE.md)。

## 1. 一阶平均

有以下无条件蕴含：

\[
\boxed{\begin{aligned}
\mu(\alpha)<3&\Longrightarrow
\frac1N\sum_{B\le N}G(B)\to0,\\
\mu(\alpha)>3&\Longrightarrow
\limsup_{N\to\infty}\frac1N\sum_{B\le N}G(B)=\infty.
\end{aligned}}
\]

临界 μ=3 不由本论证决定。这是一般性质，没有使用或判断 π 的最新无理指数声明。

## 2. 任意正阶矩

更一般地，对固定 s>0，临界位置为 2+1/s：

\[
\boxed{\begin{aligned}
\mu<2+1/s&\Longrightarrow N^{-1}\sum_{B\le N}G(B)^s\to0,\\
\mu>2+1/s&\Longrightarrow
\limsup_N N^{-1}\sum_{B\le N}G(B)^s=\infty.
\end{aligned}}
\tag{1}
\]

### 上界证明

主文已证 G→0 依自然密度。因此 0≤G≤1 部分的 s 阶平均趋于零：对任意 ε>0，G≤ε 的贡献至多 εˢ，其余贡献至多异常集合的比例。

G>1 时，真实既约中心必为普通连分数收敛分数 p_n/q_n。记 κ_n=q_n²|α−p_n/q_n|，L_n=floor(q_n/(2κ_n))。该中心在 B=kq_n、1≤k≤L_n 出现，且 G=k/q_n。

对给定 N，令 K=min(L_n,floor(N/q_n))。当 K≥1，整个中心的 s 阶平均贡献至多

\[
F_{n,s}(N)
=\frac1{q_n^sN}\sum_{k=1}^{K}k^s
\le\frac{(K+1)^s}{q_n^{s+1}}
\le C_s\left(\frac1{q_n\kappa_n^s}
+\frac1{q_n^{s+1}}\right).
\tag{2}
\]

第一步用 Σkˢ≤K(K+1)ˢ 和 N≥Kq_n；第二步用 K≤q_n/(2κ_n)。K=0 时贡献为零。每个固定中心 L_n 有限，因此 F_{n,s}(N)→0。

若 μ<2+1/s，选 μ<ν<2+1/s。由指数定义，充分大 n 有 κ_n≥q_n^{2−ν}，所以 (2) 的主要支配项至多 q_n^{s(ν−2)−1}。指数为负，而普通收敛分数分母至少按 Fibonacci 尺度增长，故这些支配项可求和。先截取有限个中心、再控制可求和尾部，得到所有高于 1 的部分平均也趋于零。

### 下界证明

若 μ>2+1/s，选 2+1/s<ν<μ。有 q→∞ 的既约逼近满足 δ=|α−p/q|<q^−ν，其 κ=q²δ→0。充分大时由 Legendre 判据它们是收敛分数。

令 L=floor(q/(2κ))、N=Lq。该中心单独贡献

\[
\frac1N\sum_{k=1}^{L}G(kq)^s
=\frac1{q^{s+1}L}\sum_{k=1}^{L}k^s
\ge\frac{L^s}{(s+1)q^{s+1}}.
\]

由于 L∼q/(2κ)，右端至少为一个仅依赖 s 的正数乘 1/(qκˢ)。又因 κ<q^{2−ν}，该量大于 q^{s(ν−2)−1}，趋于无穷。这证明 (1) 的第二行。

## 3. 与原问题的关系

μ=2 时，每个固定正阶矩平均都趋于零；仍可能有 g_B/√B 无界的稀疏子列。因此即使得到很强的平均衰减，也需另行控制最大值或主文的平方根异常计数。此前[增长指数附录](../finite-observation-hypothesis-2026-11-11/EXPONENT_SCALE.md)提供了 μ=2、非坏逼近的明确对照。

本附录为正文论证，不把有限实验当成无限证明，不宣称文献新颖性。连分数基础见 [Bugeaud, CIMPA 2010](https://irma.math.unistra.fr/~bugeaud/travaux/Bamako2.pdf)。
