# The actual logarithmic derivative has factorial scale

2026-11-02. This note concerns the same sparse one-dimensional, two-frequency raw-column family as `HANKEL_PERIOD_MAGNITUDE.md` and `LOCAL_PERIOD_DISC.md`. It proves the leading scale of the **actual formal-polynomial derivative** at the period, including its eventual alternating direction. As a consequence, the largest disk on which the polynomial stays below twice its natural period scale has logarithmic radius `-2 n log n + O(n)`. This conclusion is stronger than sharpness of a matrix-norm estimate. No location or uniqueness theorem for nearby roots, and no assertion about the bad approximability of pi, follows here.

## 1. Objects and inherited estimates

For integers n>=2, let

\[
P_n(Z)=\det\bigl([z^s]e^{hz}(jZ+z)^a\bigr),
\qquad j,h\in\{0,1\},\quad 0\le s,a<n.
\]

Rows are ordered by j then s, columns by h then a. Write

\[
t=2\pi,\quad \omega=it,\quad A_n=\prod_{a=0}^{n-1}a!,\quad
\Delta_n=\frac{A_n^4}{A_{2n}},\quad
U_n=t^{n^2}\frac{A_n^2}{A_{2n}}.
\]

The actual exponential determinant is

\[
F_n(Z)=\frac{Z^{n^2}}{A_n^2}\det H_n(Z),\qquad
H_n(Z)_{ij}=\int_0^1u^{i+j}e^{Zu}\,du.
\]

Retain the independent endpoint symbol Y and define

\[
\mu_0(Z,Y)=\frac{Y-1}{Z},\qquad
\mu_k(Z,Y)=\frac YZ-\frac kZ\mu_{k-1}(Z,Y),
\]

\[
G_n(Z,Y)=\frac{Z^{n^2}}{A_n^2}
\det[\mu_{i+j}(Z,Y)].
\]

The exact identities are `F_n(Z)=G_n(Z,e^Z)` and `P_n(Z)=G_n(Z,1)`.
The endpoint derivative moments are

\[
\nu_k(Z)=\partial_Y\mu_k(Z,Y)
=\sum_{\ell=0}^k\frac{(-1)^\ell k!}{(k-\ell)!Z^{\ell+1}},
\qquad V_n(Z)=[\nu_{i+j}(Z)].
\]

Let the probability law on [0,1]^k have density

\[
\frac{\prod_{a<b}(u_b-u_a)^2}{k!\Delta_k},
\qquad X_k=\sum_{j=1}^k u_j-k/2,
\qquad \chi_k(s)=\mathbb E_k\cos(sX_k).
\]

Reflection about 1/2 makes X_k symmetric. The preceding Hankel note proves

\[
\det H_k(is)=\Delta_k e^{iks/2}\chi_k(s),\qquad
\frac16<\chi_k(t)\le1\quad(k\ge2),
\tag{1}
\]

\[
\mathbb E_kX_k^2=\frac{k^2}{4(4k^2-1)}\le\frac1{15}
\quad(k\ge2),\qquad U_n/6<|P_n(\omega)|\le U_n.
\tag{2}
\]

All conclusions below use these already proved uniform estimates. In particular,
`|chi_k'(t)/chi_k(t)| < 6/sqrt(15) < 2` for k>=2.

## 2. One endpoint insertion is an exact Gamma integral

For any polynomial f, the endpoint functional has the exact integral representation

\[
\mathcal V_Z(f)=\frac1Z\int_0^\infty e^{-v}f(1-v/Z)\,dv,
\qquad \mathcal V_Z(u^k)=\nu_k(Z).
\tag{3}
\]

Expand the polynomial and use `integral exp(-v) v^ell dv = ell!` to prove (3).
There is no asymptotic continuation or positivity assumption for this functional.

Put `m=n-1` and

\[
D_n(Y)=\det[\mu_{i+j}(\omega,Y)],\qquad
T_n=\frac{D_n'(1)}{D_n(1)}
=\operatorname{Tr}(H_n(\omega)^{-1}V_n(\omega)).
\]

Multilinearity of the Vandermonde-square integral, with one of its n identical
functionals differentiated, gives

\[
D_n'(1)=\frac1{m!}\mathcal V_\omega^{(x)}
\int_{[0,1]^m} e^{it\sum u_j}
\prod_{a<b}(u_b-u_a)^2\prod_{j=1}^m(x-u_j)^2\,du.
\tag{4}
\]

The factor is `1/m!`, because differentiation changes `n/n!` into `1/m!`.
Consequently

\[
D_n'(1)=\frac{\Delta_m}{it}\int_0^\infty e^{-v}
\mathbb E_m\left[e^{it\sum u_j}
\prod_{j=1}^m(1+iv/t-u_j)^2\right]dv.
\tag{5}
\]

Let V, independent of the u variables, have the Gamma density
`exp(-v) v^(2m)/(2m)!` on v>0. Define

\[
C_m(v,u)=\prod_{j=1}^m\left(1-\frac{it(1-u_j)}v\right)^2.
\]

The two factors `i^(2m)=(-1)^m` and `exp(itm/2)=(-1)^m` cancel. Equation (5)
therefore becomes the exact identity

\[
\boxed{D_n'(1)=-i\,\Delta_m\frac{(2m)!}{t^{2m+1}}
\mathbb E\bigl[e^{itX_m}C_m(V,u)\bigr].}
\tag{6}
\]

The factor C_m is not bounded when V approaches zero. Pointwise convergence or
concentration alone would thus be insufficient in the next step.

## 3. The endpoint insertion has a stable phase in L1

For m>=2t, consider first the event V>=m. Since `t(1-u_j)/V <= 1/2`, define
`L_m=2 sum_j Log(1-it(1-u_j)/V)` using the power-series branch separately for
each factor. Then `C_m=exp(L_m)`, and that power series gives

\[
L_m=-it\frac mV+2it\frac{X_m}V+R,
\qquad |R|\le\frac{2mt^2}{V^2}\le\frac{2t^2}m.
\tag{7}
\]

Here `|log(1+w)-w| <= |w|^2` for `|w|<=1/2`, applied twice per factor, proves the
bound. The exponent outside R is purely imaginary. As `exp(-it/2)=-1`,

\[
|C_m+1|\le e^{2t^2/m}\frac{2t^2}m
+t\left|\frac mV-\frac12\right|+\frac{2t|X_m|}V.
\]

The Gamma variable has mean 2m+1 and variance 2m+1. On V>=m,

\[
\left|\frac mV-\frac12\right|
\le\frac{|V-2m|}{2m},\qquad
\mathbb E|V-2m|\le\sqrt{2m+2}.
\]

Together with (2), this proves

\[
\mathbb E[|C_m+1|;V\ge m]
\le \frac{2t^2e^{2t^2/m}}m
+\frac{t\sqrt{2m+2}}{2m}+\frac{2t}{\sqrt{15}\,m}.
\tag{8}
\]

For the omitted tail, the coarse inequality `|C_m(v,u)| <= (1+t/v)^(2m)` is
useful because it cancels the singular power in the Gamma density:

\[
\begin{aligned}
\mathbb E[|C_m|;V<m]
&\le\frac1{(2m)!}\int_0^m e^{-v}(v+t)^{2m}\,dv\\
&\le e^t\Pr(V\le m+t).
\end{aligned}
\]

Markov's inequality applied to `exp(-V)` gives

\[
\Pr(V\le L)\le e^L\mathbb E e^{-V}=e^L2^{-(2m+1)}.
\]

Applying it once with L=m+t and once with L=m yields

\[
\mathbb E[|C_m+1|;V<m]
\le\frac{e^{2t}+1}{2}\left(\frac e4\right)^m.
\tag{9}
\]

Since e<3, this decays exponentially. Equations (8)--(9) establish the explicit
bound

\[
\boxed{\mathbb E|C_m+1|\le\varepsilon_m:={2t^2e^{2t^2/m}\over m}
+{t\sqrt{2m+2}\over2m}+{2t\over\sqrt{15}\,m}
+{e^{2t}+1\over2}(e/4)^m=O(m^{-1/2}).}
\tag{10}
\]

Thus the remaining oscillatory factor is harmless: its modulus is exactly one,
and symmetry of X_m gives

\[
\left|\mathbb E[e^{itX_m}C_m]+\chi_m(t)\right|\le\varepsilon_m.
\tag{11}
\]

No estimate was multiplied by an uncontrolled oscillatory or unbounded factor.
For one explicit, deliberately loose threshold, `epsilon_m < 1/12` for all
m>=10000. Indeed t<44/7, t^2<40, sqrt(3)<7/4 and sqrt(15)>3 give respectively
the majorants

\[
\frac{80\cdot125}{124m},\qquad\frac{11}{2\sqrt m},
\qquad\frac{13}{3m},\qquad
\frac{3^{13}+1}{2}(3/4)^m.
\]

At m>=10000 their sum is less than 1/12: bound the first three at m=10000,
and the last by its value at exponent 100, which is less than 1/1000. All
inequalities in this last finite certificate are rational.

## 4. Exact factorial scale and the polynomial derivative

Since `D_n(1)=(-1)^n Delta_n chi_n(t)`, equations (6) and (11) prove, as n tends
to infinity,

\[
\boxed{T_n=i(-1)^nK_n\frac{\chi_{n-1}(t)}{\chi_n(t)}
\bigl(1+O(n^{-1/2})\bigr),}
\tag{12}
\]

where the error may initially be complex, and the positive scale is exactly

\[
\boxed{K_n=\frac{[(2n-2)!]^2(2n-1)!}
{[(n-1)!]^4t^{2n-1}}.}
\tag{13}
\]

Indeed, `Delta_m/Delta_(m+1)=(2m)!(2m+1)!/(m!)^4`. The relative error in (12)
is at most `6 epsilon_m`, by the uniform lower bound for chi_m; this use of the
lower bound requires m>=2, or n>=3, as is automatic in the asymptotic regime.

Differentiating the true exponential determinant gives

\[
\frac{F_n'(\omega)}{F_n(\omega)}
=\frac n2-i\left(\frac{n^2}t+\frac{\chi_n'(t)}{\chi_n(t)}\right)=O(n^2).
\tag{14}
\]

The chain rule, with the endpoint variable still independent before
differentiating, is essential:

\[
\frac{P_n'(\omega)}{P_n(\omega)}
=\frac{F_n'(\omega)}{F_n(\omega)}-T_n.
\tag{15}
\]

The scale K_n grows faster than every power of n, by (13) and Stirling's formula.
Thus (12)--(15) give the main derivative theorem

\[
\boxed{\frac{P_n'(\omega)}{P_n(\omega)}
=-i(-1)^nK_n\frac{\chi_{n-1}(t)}{\chi_n(t)}
\bigl(1+O(n^{-1/2})\bigr).}
\tag{16}
\]

In particular,

\[
\boxed{\log\left|\frac{P_n'(\omega)}{P_n(\omega)}\right|
=2n\log n+\bigl[\log(16/\pi^2)-2\bigr]n+O(\log n).}
\tag{17}
\]

To check the linear coefficient, set m=n-1. Stirling applied to (13) gives
`log K_n=2m log m+(6 log 2-2-2 log t)m+O(log m)`. The ratio of the two chi
factors lies between 1/6 and 6, so contributes only O(1) to the logarithm.

The formal moment functional at Y=1 satisfies
`L_(-Z)(f(u))=L_Z(f(1-u))`, as follows directly from its integration-by-parts
recurrence. The change from powers of u to powers of 1-u is triangular with
diagonal entries (-1)^j; applying it on both sides of the moment matrix does
not change its determinant. Multiplication by the prefactor Z^(n^2) therefore
gives `P_n(-Z)=(-1)^nP_n(Z)`. Its coefficients are real.
Hence `P_n'(it)/P_n(it)` is purely imaginary when
the denominator is nonzero. Formula (16) therefore proves that

\[
\operatorname{sgn}\left(\left.\frac{d}{ds}\log|P_n(is)|\right|_{s=t}\right)
=(-1)^n
\tag{18}
\]

for all sufficiently large n. An explicit sufficient threshold is n>=10001:
(10)--(11) imply that the imaginary part of T_n has sign (-1)^n and magnitude
greater than K_n/12. Meanwhile (14) has imaginary part of magnitude at most
`n^2/6+2`. For m>=10000,

\[
K_n\ge\frac{(2m)!}{7^{2m+1}}
\ge\frac{(m/49)^m}{7}
>\frac{200^m}{7}>2(m+1)^2+24,
\]

which proves dominance. The last elementary inequality holds for m>=2 by
induction. This threshold is a proof convenience, not a claim that the sign
first begins there.

## 5. The actual maximal small-value disk has the same leading scale

Define, for this fixed family and this fixed threshold,

\[
\mathcal R_n=\sup\{R>0:\ |P_n(Z)|\le2U_n\text{ for every }|Z-\omega|\le R\}.
\]

The prior disk theorem supplies

\[
\mathcal R_n\ge r_n=
\frac1{16\cdot3^8n^4\,2^{8n}(2n-2)!},
\qquad \log r_n=-2n\log n+O(n).
\tag{19}
\]

Conversely, on any qualifying disk the Cauchy derivative estimate at its center
gives `|P_n'(omega)| <= 2 U_n/R`. By (2),

\[
R\le\frac{2U_n}{|P_n'(\omega)|}
<\frac{12}{|P_n'(\omega)/P_n(\omega)|}.
\tag{20}
\]

Apply (17) and then take the supremum over qualifying radii. Combining (19) and
(20) proves

\[
\boxed{\log\mathcal R_n=-2n\log n+O(n).}
\tag{21}
\]

This is a property of the actual polynomial, independent of the method used to
bound its matrix norm. It says that replacing a coarse matrix-norm estimate
alone cannot enlarge this uniformly small disk beyond its leading scale. It
does not rule out estimates on a chosen direction, a noncircular region,
different columns or a different determinant, nor a method that does not
require uniform smallness on a whole complex disk.

A large first derivative does not by itself locate a zero or control higher
derivative remainders. No root-distance, root-uniqueness or rational-approximation
exponent theorem is asserted here.

## Verification

The companion `verify_period_log_derivative.py` uses exact rational arithmetic
and the Python standard library. It checks the endpoint-insertion factor by
an independent Vandermonde-square expansion, the Gamma normalization by exact
polynomial integration, the factorial scale, and the rational sufficient
constants in the explicit L1/sign threshold. The limiting theorem is proved
by (7)--(17), not extrapolated from numerical samples.


## Companion root-distance corollary

The [nearest-root synthesis](../nearest-root-scale/NEAREST_ROOT_SCALE.md) adds the polynomial degree bound and the existing disk bound to the present derivative theorem. It proves the leading distance scale of actual complex Z roots and, separately, endpoint Y roots. Thus the no-root-inference statements above describe what the derivative alone establishes, not a limitation of the combined package. They do not assert a root on the imaginary axis.

Prior inputs: [period magnitude and Jacobi moments](../../pi-period-mechanism-2026-11-02/raw-root-mechanism/HANKEL_PERIOD_MAGNITUDE.md), and [uniform local disk](../../pi-local-gram-2026-11-02/next-local-period-sensitivity/LOCAL_PERIOD_DISC.md).
