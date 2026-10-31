# Global sensitivity coordinate and endpoint asymptotics

Date: 2026-10-31

This note supplements [MIGRATION_GEOMETRY.md](MIGRATION_GEOMETRY.md). Let
\[
F(t)=\sum_{n\ge1}d_nt^{-n},\qquad t>1,
\]
where d_1=3, d_2=0, d_3=1, 0≤d_n≤3, and the word has infinitely many nonzero digits. Define
\[
P_t(n)=\frac{d_nt^{-n}}{F(t)},\quad
\mu(t)=\sum nP_t(n),\quad
V(t)=\sum(n-\mu(t))^2P_t(n).
\]

## Theorem 1: the mean digit position is a global coordinate

For every t>1,
\[
\frac{d\mu}{d\log t}=-V(t)<0.
\]
The variance is strictly positive because both positions 1 and 3 have positive mass.

Moreover,
\[
\boxed{\lim_{t\downarrow1}\mu(t)=+\infty,\qquad
\lim_{t\to\infty}\mu(t)=1.}
\]
Therefore μ is a continuous strictly decreasing bijection from (1,∞) onto (1,∞). In particular, the mean active digit position uniquely determines the base parameter.

### Proof at t→∞

Put x=1/t. Since d_1=3,d_2=0,d_3=1 and d_n≤3,
\[
F(t)=3x+x^3+O(x^4),
\]
\[
\sum n d_nt^{-n}=3x+3x^3+O(x^4).
\]
Dividing gives μ(t)→1. More precisely,
\[
\boxed{\mu(t)=1+\frac{2}{3t^2}+O(t^{-3}).}
\tag{1}
\]
Similarly, the second raw moment is
\[
\frac{\sum n^2d_nt^{-n}}{F(t)}
=1+\frac{8}{3t^2}+O(t^{-3}),
\]
so
\[
\boxed{V(t)=\frac{4}{3t^2}+O(t^{-3}).}
\tag{2}
\]

### Proof at t↓1

Since d_n are nonnegative and infinitely many are nonzero,
\[
F(t)\to\sum_{n\ge1}d_n=+\infty.
\]
For every fixed M,
\[
\sum_{n\le M}P_t(n)
=\frac{\sum_{n\le M}d_nt^{-n}}{F(t)}\longrightarrow0.
\]
Suppose μ(t_j) stayed bounded by C along some sequence t_j↓1. Markov's inequality would give
\[
\sum_{n>M}P_{t_j}(n)\le\frac{C}{M+1}.
\]
Choose M>2C. Then the probability of n≤M would be at least 1/2, contradicting its convergence to zero. Hence μ(t)→∞.

Continuity follows from locally uniform convergence of the defining moment series; strict monotonicity follows from V(t)>0. The endpoint limits prove surjectivity onto (1,∞).

## Corollary: a second coordinate system for base space

The value coordinate F(t) is a strictly decreasing bijection (1,∞)→(0,∞). The sensitivity coordinate μ(t) is another strictly decreasing bijection (1,∞)→(1,∞). Thus one may use either the transported value or the mean digit position as a global coordinate on the base axis.

The relation
\[
\frac{d\log F}{d\log t}=-\mu(t)
\]
can now be read as a first-order differential equation in the global sensitivity coordinate. The variance V controls how quickly that coordinate changes.

This does not classify F(10) as transcendental and gives no new P2848 bound. It does provide a rigorous global parameterization of migration speed, not merely a local sensitivity statistic.
