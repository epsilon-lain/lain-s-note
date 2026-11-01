# Quantitative stability at the actual period

2026-11-01. This note continues the exact tail-cancellation checkpoint. It supplies a computable stability certificate and a quantitative obstruction to assuming that small center errors imply small relative matrix errors. Neither result decides the bad approximability of π. The general finite matrix statements are unconditional; the application to the manuscript's critical high-dimensional family retains its interpolation hypotheses.

## 1. The quantity that must be controlled

Let $A(r_1,\ldots,r_m)$ be a square full-row matrix, on one fixed column family, obtained with the complete logarithmic series. Its entries are polynomials in the center increments $r_i$. Let $\omega=2\pi i$, and put

\[
B=A(\omega,\ldots,\omega),\qquad
A_{\rm rat}=A(2ip_1/q_1,\ldots,2ip_m/q_m).
\]

Assume $B$ is nonsingular. Define the relative perturbation

\[
E=B^{-1}(A_{\rm rat}-B).
\tag{1}
\]

For any submultiplicative matrix norm dominating the spectral radius, if $\eta=\|E\|<1$, then

\[
(1-\eta)^M
\le\frac{|\det A_{\rm rat}|}{|\det B|}
\le(1+\eta)^M,
\tag{2}
\]

where $M$ is the matrix size. Indeed $A_{\rm rat}=B(I+E)$ and all eigenvalues of $E$ have modulus at most $\eta$. This argument does not require $E$ to be normal or diagonalizable.

At the manuscript's normalization $MH$, equation (2) shows that a sufficient condition for asymptotically negligible determinant change is

\[
\eta_H<1,\qquad -\ln(1-\eta_H)=o(H).
\tag{3}
\]

A fixed upper bound strictly below one is more than sufficient; it is not necessary that $\eta_H\to0$. Conversely, small coordinate errors $|r_i-\omega|$ alone do not give (3), because the inverse of the period matrix appears in (1).

The exact tail identity is useful here: the full-row determinant with the rational polynomial truncations equals $\det A_{\rm rat}$. It may therefore supply an arithmetic lower bound while the stability comparison uses the full logarithmic matrices. This does not claim equality of their individual entries.

A coefficient-height estimate alone addresses only part of this question. If every entry of $B$ has modulus at most $L>0$, the adjugate formula and the ordinary determinant expansion give

\[
\|B^{-1}\|_\infty\le\frac{M!L^{M-1}}{|\det B|},
\qquad
\eta\le\frac{M!L^{M-1}\|A_{\rm rat}-B\|_\infty}{|\det B|}.
\tag{3a}
\]

Degree and coefficient-height bounds can control $L$ and the numerator perturbation. Qualitative transcendence gives $\det B\ne0$ for an appropriate nonzero period polynomial, but does not by itself supply the asymptotic lower bound for its modulus needed in (3a). This is not a claim that no quantitative transcendence bounds for π exist; none is assumed or established here. The rational preconditioner method below bypasses this issue for a specified finite matrix, while §3 shows that favorable uniform conditioning can actually fail in the critical regime.

## 2. A finite rational certificate containing the actual period

This construction avoids computing an inverse at a transcendental input. Choose a Gaussian-rational center vector $\sigma$, an invertible Gaussian-rational preconditioner $P$, and write the exact finite Taylor expansion

\[
P A(\sigma+z)-I=F_0+\sum_{\alpha\ne0}F_\alpha z^\alpha.
\]

For a Gaussian rational number define $|a+ib|_1=|a|+|b|$. For a matrix define the corresponding row-sum majorant

\[
\|F\|_{\Box}=\max_j\sum_k |F_{jk}|_1.
\]

It is submultiplicative and dominates the usual complex induced infinity norm. If the desired points lie in the polydisc $|z_i|\le\rho_i$, calculate

\[
\eta_{\Box}=\|F_0\|_{\Box}
+\sum_{\alpha\ne0}\|F_\alpha\|_{\Box}\rho^\alpha.
\tag{4}
\]

All terms in (4) are rational when the radii are rational. If $\eta_{\Box}<1$, every matrix in that polydisc is nonsingular and

\[
|\det P|^{-1}(1-\eta_{\Box})^M
\le |\det A(\sigma+z)|
\le |\det P|^{-1}(1+\eta_{\Box})^M.
\tag{5}
\]

Taking $P=A(\sigma)^{-1}$ makes $F_0=0$. A rational enclosure of π gives a rigorous polydisc containing $(\omega,\ldots,\omega)$; it can also contain a specified rational-center tuple. Formula (4), rather than a floating-point condition number, is the certificate.

The companion script applies this to the complete four-row packet with $K=2$, $\theta=1/4$, $w_0=w_1=w$, $v_0=2w$, and $H=4w$. The rows have $s=0,1$ and $\beta=0$ at both centers. Columns are $1,X,X^3,YX^2$, and

\[
A(r)=
\begin{pmatrix}
1&0&0&0\\
1&r&r^3&r^2\\
0&1&0&0\\
0&1&3r^2&r^2+2r
\end{pmatrix},\qquad
\det A(r)=r^4(1-r).
\tag{6}
\]

Exact Gaussian-rational inversion and a rational Machin-series enclosure of π give:

| Center $\sigma=2ip/q$ | Weight $w=\lceil\ln q\rceil$ | Certified bound in (4) |
|---|---:|---:|
| $2i\cdot22/7$ | 2 | $\eta_{\Box}<1/100$ |
| $2i\cdot355/113$ | 5 | $\eta_{\Box}<1/100000$ |

Both disks contain the true period. The script also verifies the row and column budgets, $w=\lceil\ln q\rceil$, and the familiar error bounds $|\pi-p/q|<q^{-2}$. These are finite source-format examples; the unknown separated-weight thresholds are not asserted. No new π approximation is claimed.

For completeness, the enclosure uses $\pi=16\arctan(1/5)-4\arctan(1/239)$. Twice applying the rational tangent doubling formula and then subtraction gives tangent one for $4\arctan(1/5)-\arctan(1/239)$. Alternating-series bounds place that angle between zero and $4/5<1<\pi/2$; the last inequality follows from $\pi=4\int_0^1(1+x^2)^{-1}dx>2$. Thus the angle is $\pi/4$. Alternating remainder bounds then give the rational interval used in the certificate.

## 3. The arithmetic and analytic bounds force relative instability

The same determinant comparison yields a lower bound on (1). Suppose, for fixed parameters and $H\to\infty$, the common column construction provides

\[
\frac{\ln|\det A_{\rm rat}|}{MH}\ge b_{\kappa,H},
\]

\[
\frac{\ln|\det B|}{MH}
\le\kappa_{\max}-S_{\rm per}+E_{\rm per}+o(1).
\tag{7}
\]

Here the columns have first been scaled by $\prod_i(q_iL_i)^{\alpha_i}$ and then mixed by a bounded integer matrix. The notation is

\[
\kappa_i=\frac{\ln(q_iL_i)}{w_i},\quad
\kappa_{\max}=\max_i\kappa_i,\quad
b_{\kappa,H}=\frac1{MH}\sum_{\rm rows}\sum_i\beta_i\ln(q_iL_i),
\]

\[
R_{\rm ar}=\kappa_{\max}-\frac\theta{m+2}\sum_i\kappa_i.
\]

The first line of (7) is the Gaussian-integral lower bound for the truncated rational-center matrix, transferred by exact tail invariance. The second is any proved analytic bound for the full period matrix; $S_{\rm per}$ denotes its saving and $E_{\rm per}$ the remaining growth cost. The mixed-column scaling is paid exactly once on this analytic side. Define

\[
\Gamma=S_{\rm per}-R_{\rm ar}-E_{\rm per}.
\]

Then

\[
\boxed{
\rho\bigl(B^{-1}(A_{\rm rat}-B)\bigr)
\ge\exp\{H(\Gamma-o(1))\}-1,
}
\tag{8}
\]

where $\rho$ is the spectral radius. In particular (8) holds with any dominating matrix norm on the left.

**Proof.** Equation (7) and the row-simplex mean imply

\[
\left|\det(B^{-1}A_{\rm rat})\right|
\ge\exp\{MH(\Gamma-o(1))\}.
\]

At least one eigenvalue of $B^{-1}A_{\rm rat}=I+E$ therefore has modulus at least $\exp\{H(\Gamma-o(1))\}$. Its corresponding eigenvalue of $E$ has modulus at least that value minus one. This proves (8). ∎

This is invariant under any further common invertible column-basis change: $E$ is then replaced by a similar matrix. It is not merely an artifact of selecting a convenient matrix norm.

Both matrices must be nonsingular on the **same** columns. The manuscript's interpolation theorem supplies full row rank for the rational centers and the exact-period centers; the previously proved bounded mixed-column lemma then supplies such a common square column family. This is a conditional existence statement in the asymptotic regime, not a claim about every raw monomial minor.

## 4. An explicit high-dimensional instability margin

Apply the improved collision estimate with radius $14K$. For the full period packet, each original type $\beta$ appears

\[
n_\beta=K\left\lceil\frac{H-w\beta/\theta}{v_0}\right\rceil
\quad(w\beta<\theta H)
\]

times. Its simplex second moment gives

\[
\frac{\sum_\beta n_\beta^2}{MH}
\longrightarrow\frac{2K}{v_0(m+2)}.
\]

Thus (7) holds with

\[
c=\frac{\ln2}{2},\quad
S_{\rm per}=\frac{2cK}{v_0(m+2)},\quad
E_{\rm per}=\max\left\{\frac{14K}{w_0},
\frac{\ln(28K)}{w_*}\right\}+\frac{\ln2}{v_0}.
\tag{9}
\]

Along the already constructed family

\[
\theta=\tfrac12,\ B_0=\tfrac35,\ C=\tfrac{13}{10},\quad
K=\lfloor C^m\rfloor,\quad
w_0=B_0^{-m},\quad v_0=2K\theta^m w_0,
\]

one has

\[
S_{\rm per}=\frac{c}{m+2}(B_0/\theta)^m
=\frac{\ln2}{2(m+2)}(6/5)^m.
\tag{10}
\]

Choose the transverse weights sufficiently large to meet the interpolation thresholds and make $\ln(28K)/w_*$ small. The weight-independent remaining terms satisfy $K/w_0\to0$ and $v_0\to\infty$. Moreover the minimal-truncation bounds imply

\[
1-1/w_*\le\kappa_i\le1+\frac{4\ln2}{\theta v_0},
\]

so $R_{\rm ar}\to1/2$ when $w_*\to\infty$. Consequently $\Gamma_m\to+\infty$ in this family.

For any fixed dimension with $\Gamma_m>0$, (8) rules out a relative perturbation bound $\|E_H\|=\exp(o(H))$ for these common full-period columns. In particular, no fixed stability radius estimate of the form (3) can be inferred merely from the smallness of the approximation errors. All approximation data, weights and dimension remain fixed before the $H$ limit. There is no contradiction: high-order interpolation can magnify a fixed tiny center displacement exponentially in $H$.

A stronger full-period analytic estimate may replace (9) in the abstract statement (7)–(8); the relative-instability conclusion then strengthens automatically. It is not necessary to assume any unproved quantitative transcendence measure for π.

### Two exact positive margins

The positive gap can be certified without floating-point logarithms. Impose $w_*\ge10^{12}$ in addition to all the original successive interpolation thresholds. Use

\[
\frac{28}{81}<\frac{\ln2}{2}<\frac{25}{72},\qquad
4\ln2<\frac{25}{9}.
\]

These bounds follow by writing $\ln2=2\sum_{k\ge0}3^{-(2k+1)}/(2k+1)$: retain the first two terms for the lower bound, and bound all denominators after the first by 3 for the upper bound. The crude bound $\ln(28K)<28K$ suffices to verify

\[
\frac{\ln(28K)}{w_*}<\frac{14K}{w_0}
\]

in each case below. The arithmetic threshold obeys

\[
R_{\rm ar}\le1-\frac{m\theta}{m+2}
+\frac{4\ln2}{\theta v_0}+\frac\theta{w_*}.
\]

Exact fraction calculations then give:

| Dimension $m$ | $K$ | Certified lower margin | Consequence for all sufficiently large admissible $H$ |
|---:|---:|---:|---|
| 32 | 4427 | $\Gamma>2$ | $\rho(E_H)>e^{2H}-1$ |
| 40 | 36118 | $\Gamma>11$ | $\rho(E_H)>e^{11H}-1$ |

The exact rational lower bounds are about 2.6996 and 11.4452, respectively; the displayed integer margins are the certified conclusions, not rounded tests. These are conditional high-dimensional statements. The bound $10^{12}$ is not a replacement for the unknown interpolation thresholds, and no explicit finite admissible $H$ or mixed-column matrix is supplied by this calculation.

## 5. Fixed genuine centers can be arbitrarily unstable as the packet grows

This phenomenon also has an explicit one-dimensional example, independent of Theorem 2.1. Use the previously established staircase column construction with

\[
K=2,\quad \theta=1/4,\quad v_0=1,\quad w_0=1/4,
\quad w_1=5,\quad H=20n.
\]

At time order $s$, use the complete derivatives $0\le\beta<b_s$, where $b_s=\lceil(H-s)/20\rceil$, and the corresponding staircase columns. Their determinant is the previously proved exact expression

\[
D_n(r)=r^{W_n},\quad
M_n=20n(n+1),\quad
W_n=\frac{10}{3}n(n+1)(2n+1),\quad
\frac{W_n}{M_n}=\frac{2n+1}{6}.
\tag{11}
\]

Take the **fixed genuine rational approximation** $p/q=355/113>\pi$, and let $r=2ip/q$. The weight $w_1=\lceil\ln113\rceil=5$ agrees with that denominator. The centers and weights are fixed as $n$ increases. Both $D_n(r)$ and $D_n(\omega)$ are nonzero, and

\[
\rho(E_n)\ge
\left(\frac{355/113}{\pi}\right)^{(2n+1)/6}-1
\longrightarrow\infty.
\tag{12}
\]

This is a lower bound on a generalized eigenvalue change, even though $|r-\omega|$ remains the same tiny number throughout the sequence.

For a completely finite certificate, take $n=200000000$ and $H=4000000000$. A rational Machin enclosure proves

\[
\frac{355/113}{\pi}-1>\frac1{12000000}.
\]

The sum of the first 41 terms of the positive binomial expansion proves

\[
\left(1+\frac1{12000000}\right)^{2n+1}>101^6.
\]

Hence $\rho(E_n)>100$. The enormous matrix is not materialized; its determinant is certified by (11), and the final inequality uses only finite integer and rational arithmetic. This example is not in the critical high-dimensional regime. Its purpose is to disprove the general inference that a fixed very good center approximation must give a uniformly small relative perturbation as the interpolation packet grows.

## 6. What this advances and what remains open

The two outcomes are complementary. Formula (4) can certify stability for a specified finite matrix at the actual period; the four-row examples do so exactly. Formula (8) shows that, in a regime where the period determinant is analytically much smaller than the arithmetic rational-center determinant, that same type of uniform stability estimate is impossible. The missing ingredient has therefore become an explicit relative spectral-growth question, rather than an appeal to qualitative nonvanishing.

These arguments do not lower-bound the period determinant uniformly in $H$, do not force cancellation of the remaining constant-shift sum, and do not prove π badly or non-badly approximable. They also do not reinterpret a vanishing factor at coalescing centers as a factor of the small approximation error.

Run `python3 verify_stability.py`. The companion JSON records two exact actual-period preconditioner certificates, two rational high-dimensional gap certificates, exact staircase counts and instability comparisons, the finite $355/113$ large-packet witness, and three direct identity-type count checks. Infinite statements rely on the proofs above, not extrapolation from these checks.
