# Nearest complex roots have the factorial sensitivity scale

2026-11-02. This note combines the newly proved logarithmic derivative with the earlier uniform disk estimate. The additional input is the finite degree of the determinant polynomial. Together these give an actual nearby complex zero, rather than only a statement about rapid variation. The leading logarithmic distance to the nearest zero is determined both in the period variable Z and, separately, in the independent endpoint variable Y.

## 1. Inputs and notation

Let n>=2, t=2*pi, omega=it, and

\[
P_n(Z)=\det\bigl([z^s]e^{hz}(jZ+z)^a\bigr),
\quad j,h\in\{0,1\},\quad 0\le s,a<n.
\]

Rows are ordered by j then s; columns by h then a. Put

\[
A_n=\prod_{a=0}^{n-1}a!,\qquad
U_n=t^{n^2}\frac{A_n^2}{A_{2n}},\qquad
L_n=2\cdot3^8n^3 2^{8n}(2n-2)!,\qquad
r_n=\frac1{8nL_n}.
\]

The earlier period-value and local-disk theorems give

\[
\frac{U_n}{6}<|P_n(\omega)|\le U_n,
\qquad |P_n(Z)|<2U_n\quad (|Z-\omega|\le r_n).
\tag{1}
\]

The derivative theorem in the companion `PERIOD_LOG_DERIVATIVE.md` gives

\[
\log Q_n=2n\log n+c n+O(\log n),\qquad
Q_n:=\left|\frac{P_n'(\omega)}{P_n(\omega)}\right|,
\qquad c=\log(16/\pi^2)-2.
\tag{2}
\]

In particular Q_n is nonzero for all sufficiently large n. All root counts and polynomial factorizations below include algebraic multiplicity.

## 2. Finite degree converts the derivative into a nearby zero

Every raw entry in column (h,a) has degree at most a in Z. Each determinant term uses each column once, so

\[
\deg P_n\le 2\sum_{a=0}^{n-1}a=n(n-1).
\tag{3}
\]

At Z=0 the two row blocks coincide, hence P_n(0)=0. By (1), P_n is not identically zero; it is therefore nonconstant. Write D_n=deg P_n, and let its complex zeros be zeta_1,...,zeta_(D_n). The factorization over the complex numbers gives the exact logarithmic derivative

\[
\frac{P_n'(\omega)}{P_n(\omega)}
=\sum_{j=1}^{D_n}\frac1{\omega-\zeta_j}.
\tag{4}
\]

Define the nearest complex-root distance

\[
\delta_n^{Z}:=\min_{P_n(\zeta)=0}|\zeta-\omega|>0.
\]

The triangle inequality in (4) gives Q_n<=D_n/delta_n^Z. Thus, whenever Q_n>0,

\[
\boxed{\delta_n^{Z}\le\frac{D_n}{Q_n}
\le\frac{n(n-1)}{Q_n}.}
\tag{5}
\]

No assertion about the signs or phases of the summands is required. Possible cancellation only weakens the inequality; it cannot invalidate it.

## 3. The disk bound excludes even closer roots

Here is the elementary analytic fact used for the lower bound. If f is analytic on the closed unit disk, |f|<=1 there, and f(a)=0 for |a|<1, then

\[
|f(0)|\le |a|.
\tag{6}
\]

For a!=0, the function g(w)=f(w)(1-conjugate(a)w)/(w-a) has a removable singularity at a and has modulus at most one on the unit circle, because |1-conjugate(a)w|=|w-a| there. The maximum-modulus principle gives |g(0)|=|f(0)|/|a|<=1. The case a=0 is immediate. This is the one-zero case of the Schwarz lemma after a disk automorphism; it uses the full disk bound, not only the first derivative.

Apply (6) to

\[
f(w)=\frac{P_n(\omega+r_n w)}{2U_n}.
\]

By (1), |f|<1 on the closed disk and |f(0)|>1/12. Any root zeta with |zeta-omega|<r_n therefore satisfies |zeta-omega|/r_n>1/12. Roots outside or on that disk satisfy the same lower bound trivially. Consequently, for every n>=2,

\[
\boxed{\delta_n^{Z}>\frac{r_n}{12}.}
\tag{7}
\]

Combining (2), (5), (7), and Stirling's formula for r_n proves

\[
\boxed{\log\delta_n^{Z}=-2n\log n+O(n).}
\tag{8}
\]

More precisely the upper estimate furnished by the derivative is

\[
\log\delta_n^{Z}\le -2n\log n-c n+O(\log n).
\tag{9}
\]

The lower bound has a larger unspecified linear cost, so (9) must not be read as a two-term asymptotic formula for the root distance.

## 4. The nearest independent endpoint root has the same leading scale

For completeness define the endpoint moment matrices

\[
U(Z)_{ij}=(-1)^{i+j}(i+j)!Z^{-i-j-1},\qquad
V(Z)_{ij}=\sum_{\ell=0}^{i+j}
\frac{(-1)^\ell(i+j)!}{(i+j-\ell)!Z^{\ell+1}},
\]

and the endpoint pencil

\[
G_n(Z,Y)=\frac{Z^{n^2}}{A_n^2}\det(YV(Z)-U(Z)).
\]

At omega, H=V(omega)-U(omega) is invertible by (1). Let M=H^(-1)V(omega), and define

\[
\delta_n^Y:=\min_{G_n(\omega,Y)=0}|Y-1|.
\]

The companion endpoint-pencil note proves the exact identities

\[
\frac{G_n(\omega,1+\eta)}{P_n(\omega)}=\det(I+\eta M),
\qquad
\delta_n^Y=\frac1{\rho(M)},
\tag{10}
\]

where rho is the spectral radius. It also proves rho(M)<6L_n, by similarity to its Hilbert-Gram normalization. On the other hand,

\[
|\operatorname{tr}M|\le n\rho(M),
\]

because the trace is the sum of the eigenvalues with algebraic multiplicity, without requiring diagonalizability. Hence, whenever tr M!=0,

\[
\boxed{\frac1{6L_n}<\delta_n^Y
\le\frac n{|\operatorname{tr}M|}.}
\tag{11}
\]

The same endpoint-insertion calculation used in the derivative theorem proves

\[
\operatorname{tr}M
=i(-1)^n K_n\frac{\chi_{n-1}(t)}{\chi_n(t)}
\bigl(1+O(n^{-1/2})\bigr),
\]

where

\[
K_n=\frac{[(2n-2)!]^2(2n-1)!}
{[(n-1)!]^4t^{2n-1}},\qquad
\frac16<\chi_j(t)\le1\quad(j\ge2).
\]

Thus log|tr M|=2n log n+c n+O(log n). Applying this to (11) yields

\[
\boxed{\log\delta_n^Y=-2n\log n+O(n).}
\tag{12}
\]

This improves the earlier standalone endpoint-root estimate, which used the determinant but not the newly established large trace.

## 5. Scope and what the result does not identify

The new conclusion is that **actual complex zeros of P_n approach 2*pi*i at the displayed scale**. It supplements the derivative with two global facts: finite degree for the upper distance bound, and a bounded analytic disk for the lower distance bound. A large derivative by itself would not suffice.

The zeros in (8) are not proved to lie on the imaginary axis. Parity and real coefficients make the zero set symmetric about that axis, but an off-axis conjugate pair is compatible with those symmetries. The theorem does not prove simplicity, uniqueness, a specified side of approach, or convergence of a Newton iteration.

The endpoint zeros in (12) are a different set: Z is fixed and Y varies independently. Their common leading distance scale does not identify them with the Z zeros or impose Y=exp(Z).

No rational approximation p/q has been constructed here. Root distance in a polynomial family does not control the height, degree, or denominator costs needed for pi's irrationality measure. Neither bad approximability nor its negation follows.

## Verification

`verify_nearest_root_scale.py` uses exact rational and Gaussian-rational arithmetic. It constructs small P_n from the compressed polynomial determinant, compares their evaluations with the original raw determinant, checks the degree bound and the zero at the origin, and audits the logarithmic-derivative and trace root-distance inequalities on exact factored examples. The analytic disk lemma and asymptotic conclusions are proved above; finite tests do not establish an asymptotic or a root on the imaginary axis.


## Proof dependencies

- [Actual logarithmic derivative and endpoint trace](../period-derivative-scale/PERIOD_LOG_DERIVATIVE.md).
- [Endpoint eigenvalue/root correspondence](../endpoint-pencil-structure/ENDPOINT_PENCIL_STRUCTURE.md).
- [Prior period magnitude and nonvanishing](../../pi-period-mechanism-2026-11-02/raw-root-mechanism/HANKEL_PERIOD_MAGNITUDE.md).
- [Prior uniform disk and Gram bounds](../../pi-local-gram-2026-11-02/next-local-period-sensitivity/LOCAL_PERIOD_DISC.md).
