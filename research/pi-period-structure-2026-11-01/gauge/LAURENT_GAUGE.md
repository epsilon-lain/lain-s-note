# A determinant-preserving Laurent gauge and a joint column norm

2026-11-01. This note reviews and proves the proposed common Laurent gauge. The finite row identity and norm estimate are unconditional under their displayed assumptions. Existence of the required high-dimensional nonzero minors still uses the manuscript's interpolation theorem. Nothing here proves a statement about the bad approximability of π.

## 1. The full-row determinant is unchanged

Fix a complete source packet

\[
\mathcal R_H=\{(j,s,\beta):0\le j<K,\quad v_0s+w\beta/\theta<H\}.
\]

All weights are positive, all indices are nonnegative integers, and the columns are fixed polynomials from the original weighted domain

\[
w_0h+w\alpha\le H.
\]

For any integer $d\ge0$, replace every column $P(Y,X)$ by the Laurent polynomial

\[
P^{[d]}(Y,X)=Y^{-d}P(Y,X).
\]

These Laurent columns are used for estimation. No enlarged polynomial interpolation theorem is being asserted. At every source center the $Y$ coordinate is $1+t$, so the resulting row matrix obeys

\[
A^{[d]}_{j,s,\beta}
=\sum_{\ell=0}^s\binom{-d}{s-\ell}A_{j,\ell,\beta}.
\tag{1}
\]

Every row on the right belongs to the same packet, because only the time index decreases. The transformation is lower triangular in time order and has diagonal entries one. Consequently

\[
\boxed{A^{[d]}=U_dA,\qquad\det U_d=1.}
\tag{2}
\]

Every full-row square minor on the same columns is therefore unchanged. The identity applies before or after fixed column scaling and mixing. It requires closure under lowering time order at each fixed $(j,\beta)$; an arbitrary sparse source subset need not have it.

## 2. Arithmetic integrality and period translation survive

Keep the minimal logarithmic truncations and arithmetic factors from the preceding tail-cancellation note:

\[
T_i=\left\lceil\frac{w_i}{\theta v_0}\right\rceil,
\qquad L_i=\operatorname{lcm}(1,\ldots,T_i-1),
\qquad Q_\alpha=\prod_i(q_iL_i)^{\alpha_i}.
\]

The empty least common multiple is one. With $r_i=2ip_i/q_i$ and $G_i$ the truncated logarithm, a column-and-row scaled Laurent entry is

\[
\binom\alpha\beta[t^s](1+t)^{h-d}
\prod_i[q_iL_i(jr_i+G_i(t))]^{\alpha_i-\beta_i}.
\tag{3}
\]

This is a Gaussian integer. Each factor in brackets has Gaussian-integer coefficients, and $\binom{h-d}{k}$ is an integer for every integer $h-d$ and every $k\ge0$. Only finitely many series coefficients are extracted. Thus the same arithmetic lower bound, including its favorable row factors, holds without an extra denominator cost.

The entire functions now become

\[
f^{[d]}_{a,P}(z)
=\binom\alpha a e^{(h-d)z}z^{|\alpha|-|a|},
\qquad a\le\alpha.
\tag{4}
\]

They are entire even when $h-d<0$. At the true period $\omega=2\pi i$, the integer condition gives $e^{(h-d)j\omega}=1$. Therefore the exact constant-translation identity, after tail cancellation, remains valid with the same scalar factors and with these new functions. Types still satisfy $a\le\alpha$ and $wa\le H$.

The gauge may depend on $H$. Equations (1)–(4) hold separately for every integer $d_H$. A fixed Cauchy radius can still be used, because the norm estimate below explicitly includes the $H$-dependence of $d_H$.

## 3. Optimize the joint column budget

Fix a real Cauchy radius $R>7K$. Define

\[
a_R=\frac R{w_0},\qquad
b_R=\max_i\frac{\ln(q_iL_i)+\ln(2R)}{w_i}.
\]

For the already arithmetically scaled monomial $Q_\alpha Y^{h-d}X^\alpha$, its type-$a$ entire function satisfies, on $|z|\le R$,

\[
\ln|Q_\alpha f^{[d]}_{a,P}(z)|
\le R|h-d|+\sum_i\alpha_i\bigl(\ln(q_iL_i)+\ln(2R)\bigr)
\le R|h-d|+(H-w_0h)b_R.
\tag{5}
\]

The binomial coefficient was bounded by $2^{|\alpha|}$ and the nonnegative residual power of $z$ by $R^{|\alpha|}$. For $0\le d\le H/w_0$, the right side of (5) is a convex function of $h\in[0,H/w_0]$. Its maximum is at an endpoint, giving

\[
\max\{Rd+Hb_R,\ RH/w_0-Rd\}.
\tag{6}
\]

Choose

\[
d_H=\left\lfloor\frac{H(a_R-b_R)_+}{2R}\right\rfloor,
\qquad
G_R=\max\left\{b_R,\frac{a_R+b_R}{2}\right\}.
\]

Then $0\le d_H\le H/(2w_0)$, and (6) is at most

\[
\boxed{HG_R+R.}\tag{7}
\]

Indeed, if $a_R\le b_R$ take $d_H=0$. Otherwise the unrounded value balances the two endpoints at $HG_R$; rounding down decreases the first and increases the second by less than $R$.

The real-valued optimization of (6) is exact. Formula (7) is an upper bound for the actual column functions, not a lower bound on their unavoidable size. It optimizes this particular uniform envelope over common gauges.

For bounded integer mixing of all $N$ scaled columns, the norm is at most

\[
2N\exp(HG_R+R).
\]

The original separate estimate charged $\kappa_{\max}$ plus a holomorphic maximum, where $\kappa_i=\ln(q_iL_i)/w_i$. Here the arithmetic column scaling and analytic growth have been combined within the same weighted column budget. They must not be charged a second time.

## 4. The resulting determinant comparison and rank scope

Let $M$ be the number of source rows and

\[
b_{\kappa,H}=\frac1{MH}\sum_{\rho,i}\beta_{\rho,i}\ln(q_iL_i).
\]

For a nonzero determinant of bounded mixed scaled columns the arithmetic lower bound is $b_{\kappa,H}$. Applying the preceding sharp collision estimate to (4), with

\[
c_R=\frac12\ln\frac R{7K},
\]

shows that the analytic upper bound contains the joint column cost $G_R$, coefficient cost $\ln2/v_0$, and the same scalar-and-collision savings. The extra gauge and mixing remainder is

\[
\frac{R+\ln(2N)}H\longrightarrow0
\]

for fixed parameters and radius. Thus the main comparison threshold is

\[
G_R-b_{\kappa,H}+\frac{\ln2}{v_0}+E_{\rm tr},
\tag{8}
\]

where $E_{\rm tr}$ is the appropriate constant-translation scalar error. This includes the column normalization cost once.

One must distinguish determinant preservation from sparse test-row rank. Equation (2) does **not** say that the determinant of an arbitrary sparse test selection is unchanged. The earlier shell rank certificate can nevertheless be repeated: its selected $(s,a)$ rows lie inside a complete auxiliary one-center packet. Multiplying by $(1+t)^{-d_H}$ acts invertibly on that entire packet, which therefore remains independent; every subset remains independent. The formal-period polynomial argument then certifies rank at the actual period, using the same known transcendence input. A new common bounded mixing may be chosen for the gauged source and gauged test matrices. Its coefficients may depend on $H$.

Thus the gauge is compatible with that conditional construction, but it is not a general preservation theorem for arbitrary sparse test minors on previously fixed columns.

## 5. The same shell family still cannot close this uniform-envelope criterion

The strengthened envelope requires a new radius check; the older bound that separately charged $R/w_0$ cannot simply be reused. Put

\[
\kappa=\max_i\kappa_i>0,
\qquad A=\frac{7K}{w_0},\qquad d_0=\frac1{2v_0}.
\]

Since $b_R\ge\kappa$,

\[
G_R-\kappa
\ge\frac12\max\{0,R/w_0-\kappa\}.
\]

Writing $u=R/(7K)>1$, the extra collision gain above the baseline column cost is at most

\[
\frac{c_R}{v_0}-(G_R-\kappa)
\le d_0\ln u-\frac12\max\{0,Au-\kappa\}.
\tag{9}
\]

Suppose $\kappa>A$ and $\kappa>2d_0=1/v_0$. On $1<u\le\kappa/A$, the right side increases as $d_0\ln u$. Above that point its derivative is $d_0/u-A/2<0$. Consequently

\[
\boxed{
\sup_{R>7K}\left[\frac{c_R}{v_0}-(G_R-\kappa)\right]
\le\frac1{2v_0}\ln\frac{\kappa w_0}{7K}.}
\tag{10}
\]

Along the preceding explicit family $\theta=1/2$, $B=3/5$, $C=13/10$, $K=\lfloor C^m\rfloor$, $w_0=B^{-m}$ and $v_0=2K\theta^mw_0$, one has $A\to0$ and $v_0\to\infty$. If $w_*\to\infty$, then the uniform arithmetic bounds give $\kappa\to1$. The two hypotheses for (10) therefore hold eventually, and its right side is $O(m/v_0)\to0$.

For the mean-quality shell matching, linear cost is at most $\bar\nu J+o(1)$ and the occupancy contribution is at most $c_R/v_0+o(1)$. For $\bar\nu\to2$, $J\to3/20$ and the arithmetic baseline $\kappa-b_{\kappa,H}$ tends to $1/2$ in the prescribed order of limits. Equation (10) shows that jointly optimizing this Laurent gauge and choosing a separate radius for each dimension still leaves the main score below the required threshold. The linear score tends to $0.3$, while the baseline threshold tends to $0.5$.

This concerns the displayed uniform envelope and shell construction. It gives no lower bound on any actual test determinant, does not exclude a sharper grouped norm or cancellation estimate, and does not cover an unrestricted radius depending on $H$. Dimension, weights and radius are fixed before each $H$ limit; only the explicitly controlled integer gauge $d_H$ grows with $H$.

## 6. Exact verification

Run `python3 verify_laurent_gauge.py`. The standard-library verifier checks complete source matrices at several positive gauges, including gauges larger than all original $Y$ degrees; verifies negative-binomial arithmetic and Gaussian-integral scaled entries; and checks the endpoint norm optimization with rational inputs. It also verifies the derivative signs in the plateau-cost radius calculation. These finite checks accompany the general proofs and do not certify the source interpolation thresholds or any new π approximation.
