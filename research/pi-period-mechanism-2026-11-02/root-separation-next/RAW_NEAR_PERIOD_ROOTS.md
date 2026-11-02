# Genuine raw minors can have a root extremely close to the actual period

2026-11-02. This note advances the earlier small root-geometry examples to a legal raw-column family with certified roots near \(\omega=2\pi i\). The strongest finite certificate uses twenty complete source rows and places one simple imaginary root within \(6\cdot10^{-12}\) of \(\omega\). It also certifies strong cancellation at the actual period. No infinite root-convergence theorem or rational approximation improvement is asserted.

## 1. A complete packet and two raw frequencies

For an integer \(n\ge4\), choose

\[
K=2,\quad m=1,\quad H=2n,\quad
w_0=1,\quad w_1=2,\quad v_0=2,\quad v_1=2n,
\quad\theta=1/n.
\]

The complete source inequality \(2s+2n\beta<2n\) gives exactly

\[
L_n=\{(j,s,0):j=0,1,\ 0\le s<n\}.
\]

Select the \(2n\) distinct raw columns

\[
C_n=\{Y^hX^a:h=0,1,\ 0\le a<n\}.
\]

They satisfy \(h+2a\le2n-1<H\). Thus no source rows have been removed and no mixed polynomial columns are used. Put

\[
P_n(Z)=\det\left[[z^s]e^{hz}(jZ+z)^a\right]_{(j,s),(h,a)}.
\tag{1}
\]

The order is center then time for rows, frequency then transverse degree for columns. The formal parameter \(Z\) is not assumed to be an exponential period. Only \(Z=\omega\) is the actual period used below.

This family has changing \(\theta\), fixed transverse dimension, and sparse column density. It is not the critical high-dimensional half-density family. In particular, its source count does not satisfy the critical relation \(v_0=2K\theta^m w_0\) as \(n\) grows.

## 2. An exact finite-dimensional operator description

On polynomials of degree below \(n\), in the coefficient basis \(1,x,\ldots,x^{n-1}\), let

\[
E_{sa}=\frac{\mathbf1_{s\ge a}}{(s-a)!},\qquad
T(Z)_{sa}=\binom as Z^{a-s}\mathbf1_{a\ge s}.
\]

Here \(E\) multiplies by \(e^x\) and truncates modulo \(x^n\), while \(T(Z)\) translates a polynomial by \(Z\). Matrix (1) is

\[
\begin{pmatrix}I&E\\T&ET\end{pmatrix},
\qquad
\boxed{P_n(Z)=\det(ET-TE).}
\tag{2}
\]

Thus these period polynomials arise from the failure of truncated exponential multiplication and polynomial translation to commute. If \(D\) is polynomial differentiation, \(u\) is the last coefficient vector, and \(v_a=1/(n-a)!\), the exact boundary identity is

\[
\boxed{DE-E(D+I)=-nuv^\mathsf T.}
\tag{3}
\]

Differentiating a truncated product loses precisely the derivative of its degree-\(n\) coefficient, which proves (3). This identifies a concrete boundary mechanism; it does not by itself bound the determinant or prove a limiting root theorem.

For example,

\[
P_4(Z)=\frac{Z^4(Z^2-12)(Z^2+36)}{144},
\]

\[
P_5(Z)=-\frac{Z^5(Z^2+24)(Z^4+36Z^2-144)}{864}.
\]

The first has the exact root \(6i\); the second has a root \(i\sqrt{18+6\sqrt{13}}\). Later polynomials are stored as exact rational coefficients in the verification JSON.

## 3. Certified finite root proximity

For each \(4\le n\le10\), the verifier recovers the entire polynomial by exact rational interpolation using the degree bound \(n(n-1)\), and checks two additional evaluation points. It then forms the real polynomial

\[
q_n(x)=i^{-n}x^{-n}P_n(ix).
\]

Subtracting each center-zero row from its corresponding center-one row proves that \(Z^n\) divides \(P_n\). Frequency reflection interchanges the two blocks of \(n\) columns, and the general reflection identity therefore gives \(P_n(-Z)=(-1)^nP_n(Z)\). These facts explain why the displayed \(q_n\) is a real rational polynomial; the verifier also checks both coefficient properties directly.

Every displayed root lies in a rational bracket with opposite signs of \(q_n\) at its endpoints. Exact interval evaluation of \(q_n'\) excludes zero throughout the bracket, proving that the root there is simple and unique. The exact root \(6\) for \(n=4\) is checked directly. An independent rational enclosure of \(\pi\) is obtained from Machin's identity and alternating arctangent series. No floating-point root is used as a certificate.

Let \(\rho_n=i r_n\) be the root certified in this way. The following decimals are only readable approximations to the exact rational enclosures:

| n | Complete source rows | \(r_n-2\pi\), approximately |
|---:|---:|---:|
| 4 | 8 | \(-2.8318531\cdot10^{-1}\) |
| 5 | 10 | \(+1.2313692\cdot10^{-2}\) |
| 6 | 12 | \(-3.1338268\cdot10^{-4}\) |
| 7 | 14 | \(+5.6104260\cdot10^{-6}\) |
| 8 | 16 | \(-7.3820113\cdot10^{-8}\) |
| 9 | 18 | \(+7.4258633\cdot10^{-10}\) |
| 10 | 20 | \(-5.8895947\cdot10^{-12}\) |

In particular the exact checks prove

\[
\boxed{0<|\rho_{10}-\omega|<6\cdot10^{-12},
\qquad \frac{|\rho_{10}|}{|\omega-\rho_{10}|}>10^{12}.}
\tag{4}
\]

This supplies a genuine example with a large root-separation statistic from the preceding note. It is not a proof that these distances tend to zero for all \(n\), nor a counterexample to any specified eventual polynomial lower bound.

## 4. Consequence for the exact common-radial aggregate

Every selected column has a single transverse degree, so the preceding reciprocal identity applies on these same columns. With \(D_n=n(n-1)\), define the exact normalized common-radial aggregate

\[
\mathcal F_n(\zeta)
=(1-\zeta)^{D_n}P_n\!\left(\frac\omega{1-\zeta}\right).
\tag{5}
\]

At

\[
\zeta_n^*=1-\frac\omega{\rho_n}
=\frac{r_n-2\pi}{r_n},
\]

equation (5) vanishes exactly. For the certified twenty-row example,

\[
\boxed{0<|\zeta_{10}^*|<10^{-12},\qquad
\mathcal F_{10}(\zeta_{10}^*)=0.}
\tag{6}
\]

In contrast \(\mathcal F_{10}(0)=P_{10}(\omega)\ne0\), since \(P_{10}\) is a nonzero rational polynomial and \(\omega\) is transcendental. Thus a universal zero-free perturbation disk of radius \(10^{-12}\), over all legal raw minors of this general type, is already false.

The displaced center \(\rho_{10}\) is algebraic; it is not asserted to be of the form \(2ip/q\) with integers \(p,q\). Accordingly (6) is not a new rational approximation to \(\pi\). It identifies sensitivity of the aggregate polynomial, not an admissible Diophantine counterexample.

## 5. Actual-period cancellation and normalization

The same exact \(\pi\) enclosure certifies

\[
\boxed{0<|P_{10}(\omega)|<1.3\cdot10^{-16}.}
\tag{7}
\]

Write \(P_{10}(Z)=\sum_k c_k Z^k\). Then

\[
\boxed{\frac{|P_{10}(\omega)|}{\sum_k|c_k|\,|\omega|^k}
<1.8\cdot10^{-32}.}
\tag{8}
\]

Both inequalities concern the specified coefficient-normalized determinant. They demonstrate cancellation of its rational polynomial terms at the actual period. They are not bounds on an unnormalized integer certificate.

Multiplying every source row of time order \(s\) by \(s!\) clears its coefficient factorials. The resulting determinant is multiplied by

\[
B_n=\left(\prod_{s=0}^{n-1}s!\right)^2.
\]

At a rational imaginary center \(Z=2ip/q\), powers of \(q\) must also be cleared; the general safe power is \(q^{D_n}\). Any arithmetic comparison must retain these costs. Small values in (7), or a small normalized root distance, cannot be read as a contradiction to integrality. No improvement in rational approximation quality is established here.

## Verification

`python verify_raw_near_period.py` uses only the Python standard library. It verifies seven complete packets, seven simple-root certificates, 35 exact commutator evaluations, seven rank-one boundary identities, the \(\pi\) enclosure, and the period-value/cancellation bounds. All resulting polynomials and rational enclosures are in `raw_near_period_verification.json`.
