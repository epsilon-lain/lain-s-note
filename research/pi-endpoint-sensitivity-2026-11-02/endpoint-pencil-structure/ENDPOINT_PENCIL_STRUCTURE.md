# Endpoint reciprocity and many nearby endpoint roots

2026-11-02. The endpoint perturbation has genuine spectral growth: its large normalized norm is not solely a nonnormal-matrix effect. This note proves exact reciprocal symmetry, fixes the real part of the logarithmic endpoint derivative, and gives nearby roots in the independent endpoint variable. These roots do not yet supply zeros in the period variable or rational approximations to pi.

## 1. Definitions and prior inputs

For n >= 2 put A_n = product_{a=0}^{n-1} a!, omega = 2*pi*i, and t = 2*pi. Use the n by n matrices

\[
E_{sa}=\begin{cases}1/(s-a)!,&s\ge a,\\0,&s<a,\end{cases}
\qquad
T(Z)_{sa}=\begin{cases}\binom as Z^{a-s},&a\ge s,\\0,&a<s.\end{cases}
\]

The two-parameter polynomial and its formal-period specialization are

\[
G_n(Z,Y)=\det(YET(Z)-T(Z)E),\qquad P_n(Z)=G_n(Z,1).
\tag{1}
\]

For Z != 0 define the moment matrices

\[
U(Z)_{ij}=(-1)^{i+j}(i+j)!Z^{-i-j-1},
\]
\[
V(Z)_{ij}=\nu_{i+j}(Z),\qquad
\nu_k(Z)=\sum_{\ell=0}^k
\frac{(-1)^\ell k!}{(k-\ell)!Z^{\ell+1}}.
\]

The earlier formal Hankel identity gives

\[
G_n(Z,Y)=\frac{Z^{n^2}}{A_n^2}\det(YV(Z)-U(Z)).
\tag{2}
\]

At the actual period omega, H := V(omega)-U(omega) equals the matrix with entries integral_0^1 u^(i+j) exp(omega*u) du. The previously proved [two-sided period estimate](../../pi-period-mechanism-2026-11-02/raw-root-mechanism/HANKEL_PERIOD_MAGNITUDE.md) is

\[
\frac16\mathcal U_n<|P_n(\omega)|\le\mathcal U_n,
\qquad
\mathcal U_n=t^{n^2}\frac{A_n^2}{A_{2n}},
\tag{3}
\]
\[
\log |P_n(\omega)|=-n^2\log n+
\left[\frac32+\log(\pi/2)\right]n^2+O(n\log n).
\tag{4}
\]

For J = [1/(i+j+1)] and its real positive inverse square root S = J^(-1/2), the [Gram estimates](../../pi-local-gram-2026-11-02/next-local-period-sensitivity/LOCAL_PERIOD_DISC.md) give

\[
\|SHS\|_2\le1,\qquad
\|SV(\omega)S\|_2\le L_n,
\quad L_n=2\cdot3^8n^3 2^{8n}(2n-2)!.
\tag{5}
\]

Every asymptotic below is for n tending to infinity; constants do not depend on n.

## 2. A determinant-one commutator explains reciprocal symmetry

Both E and T have determinant one, and T(-Z)=T(Z)^(-1). Set

\[
C(Z)=T(Z)^{-1}E^{-1}T(Z)E.
\]

Factoring ET in (1) proves

\[
G_n(Z,Y)=\det(YI-C(Z)),\qquad \det C(Z)=1.
\tag{6}
\]

Thus G is monic of degree n in Y, with constant coefficient (-1)^n. Its Y roots are all nonzero. Direct multiplication gives

\[
C(-Z)=T(Z)C(Z)^{-1}T(Z)^{-1}.
\]

Taking characteristic polynomials proves the exact identity

\[
\boxed{G_n(-Z,Y)=(-1)^nY^nG_n(Z,Y^{-1}).}
\tag{7}
\]

The right side is interpreted as a polynomial, including at Y=0. Because the coefficients in Z are real, at any purely imaginary Z this becomes the self-inversive identity

\[
\boxed{G_n(Z,Y)=(-1)^nY^n
\overline{G_n(Z,1/\overline Y)}.}
\tag{8}
\]

Consequently the endpoint roots pair as lambda <-> 1/conjugate(lambda), with algebraic multiplicity preserved. When n is odd, at least one endpoint root lies on the unit circle. It need not be Y=1: direct expansion gives P_3(Z)=2Z^3, nonzero for every Z != 0.

If Z is purely imaginary and P_n(Z) != 0, differentiate (8) at Y=1. More explicitly, for f(Y)=G_n(Z,Y), the polynomial identity is f(Y)=(-1)^n Y^n conjugate(f)(1/Y), where conjugate(f) means coefficientwise conjugation. Hence

\[
\frac{f'(1)}{f(1)}=n-\overline{\frac{f'(1)}{f(1)}}.
\]

Therefore

\[
\boxed{\operatorname{Re}\frac{\partial_YG_n(Z,1)}{P_n(Z)}=\frac n2.}
\tag{9}
\]

This identity imposes no bound on the imaginary part. It must not be used as an upper bound on the absolute logarithmic derivative.

## 3. Translation preserves the endpoint moment determinant

For a polynomial f introduce the linear functional

\[
\mathscr U_Z[f]=\frac1Z\int_0^\infty e^{-v}f(-v/Z)\,dv.
\]

The integral is a finite sum of convergent factorial moments, valid for every complex Z != 0. Its monomial moments give U, while the moments of f(1+u) give V. In particular, for the unit lower triangular matrix B_{ik}=binom(i,k),

\[
V(Z)=B U(Z) B^T.
\tag{10}
\]

The factorial Hankel determinant is det[(i+j)!] = A_n^2. An elementary proof factors i! from row i and j! from column j; the remaining matrix binom(i+j,i) equals B B^T by Vandermonde's identity. Factoring powers and signs from U now yields

\[
\boxed{\det U(Z)=\det V(Z)=\frac{A_n^2}{Z^{n^2}}.}
\tag{11}
\]

At omega the matrix H is invertible by (2)--(3). Define the relative endpoint operator

\[
M_n=H^{-1}V(\omega).
\]

Equations (2) and (11) imply

\[
\boxed{\frac{G_n(\omega,1+\eta)}{P_n(\omega)}
=\det(I+\eta M_n),\qquad
\det M_n=\det(M_n-I)=\frac1{P_n(\omega)}.}
\tag{12}
\]

Thus M_n has neither 0 nor 1 as an eigenvalue. For an eigenvalue mu of M_n, its corresponding endpoint root is

\[
Y=1-\frac1\mu,\qquad \eta=-\frac1\mu.
\tag{13}
\]

No diagonalizability assumption is needed: det(I+eta*M_n) is the product of 1+eta*mu over eigenvalues counted with algebraic multiplicity. Equations (8), (9), and (13) give

\[
\boxed{\mu\longleftrightarrow1-\overline\mu,
\qquad \operatorname{Re}\operatorname{tr}M_n=\frac n2.}
\tag{14}
\]

Large eigenvalues can therefore occur in pairs whose real parts sum to one. The real trace identity does not exclude a large spectral radius.

## 4. Genuine spectral growth and a nearby endpoint zero

The Hilbert determinant is det J = A_n^4/A_{2n}. Therefore (2)--(3) give

\[
|\det(SHS)|=\frac{|P_n(\omega)|}{\mathcal U_n}>\frac16.
\]

All singular values of SHS are at most one by (5). Their product exceeds 1/6, so each exceeds 1/6; in particular ||(SHS)^(-1)||_2 < 6. The matrix

\[
\widetilde M_n=(SHS)^{-1}(SVS)=S^{-1}M_nS
\]

is similar to M_n. Its norm is less than 6L_n, hence

\[
\boxed{|P_n(\omega)|^{-1/n}\le\rho(M_n)<6L_n.}
\tag{15}
\]

The lower bound uses |det M_n| <= rho(M_n)^n and (12). In logarithmic terms,

\[
n\log n+O(n)\le\log\rho(M_n)\le2n\log n+O(n).
\tag{16}
\]

This proves growth of actual eigenvalue moduli; a large norm alone would not do so.

Let delta_n be the distance from Y=1 to the nearest root of G_n(omega,Y). Equation (13) gives delta_n = 1/rho(M_n), so

\[
\boxed{\frac1{6L_n}\le\delta_n
\le |P_n(\omega)|^{1/n}\le\mathcal U_n^{1/n}.}
\tag{17}
\]

Thus a zero of the independent endpoint pencil approaches Y=1 at least as fast as exp(-n log n+O(n)). This is not yet a matching sharp endpoint-root scale.

## 5. A positive fraction of the endpoint roots approach one

Fix 0 < alpha < 1. Let k_n(alpha) count eigenvalues of M_n with

\[
|\mu|>T_n,\qquad T_n=\exp(\alpha n\log n),
\]

always counting algebraic multiplicity. By (15),

\[
|\det M_n|\le(6L_n)^{k_n}T_n^{n-k_n}.
\]

For all sufficiently large n the denominator in the next expression is positive. Taking logarithms and using (4),

\[
k_n\ge\frac{-\log|P_n(\omega)|-n\log T_n}
{\log(6L_n)-\log T_n}
=\left(\frac{1-\alpha}{2-\alpha}+O_\alpha(1/\log n)\right)n.
\]

Consequently

\[
\boxed{\liminf_{n\to\infty}\frac{k_n(\alpha)}n
\ge\frac{1-\alpha}{2-\alpha}.}
\tag{18}
\]

Equivalently, that fraction of the endpoint roots satisfies

\[
|Y-1|<\exp(-\alpha n\log n).
\]

For alpha=1/2 the guaranteed fraction is 1/3 in the limit inferior. The count does not assert that these roots are distinct. A companion exterior-power estimate can improve this count by controlling products of the largest eigenvalue moduli; (18) is the standalone consequence requiring only the determinant and the previously established single norm bound.

There is also a finite sufficient certificate involving only rational arithmetic. For even n set

\[
T_n=n^{n/2},\qquad
\mathcal U_n^{+}=(44/7)^{n^2}\frac{A_n^2}{A_{2n}}.
\]

If, for an integer 1 <= k <= n,

\[
\boxed{\mathcal U_n^{+}(6L_n)^{k-1}T_n^{n-k+1}<1,}
\tag{19}
\]

then at least k endpoint roots satisfy |Y-1| < n^(-n/2), with algebraic multiplicity. Otherwise at most k-1 eigenvalues would exceed T_n; their product would be at most (6L_n)^(k-1) T_n^(n-k+1), contradicting
\(|\det M_n|=1/|P_n(\omega)|>1/\mathcal U_n^{+}\).
The companion verifies (19) for (n,k) = (64,1), (96,3), (128,5). These deliberately conservative finite counts do not estimate the actual counts or replace the asymptotic argument.

## 6. What this resolves and what remains open here

The endpoint matrix has an exact multiplicative determinant constraint, a reciprocal-conjugate symmetry, and many spectrally large directions. Its large scale cannot be dismissed solely as an artifact of nonnormality, nor does its fixed real trace show that only a few directions matter.

All the nearby zeros above are in **Y at fixed Z=omega**. They are not zeros of P_n(Z)=G_n(Z,1). They also need not lie on the actual endpoint relation Y=exp(Z). Converting this endpoint information into control along either of those curves requires an additional coupling argument. No bound of order q^(-2), no new rational approximation record, and no conclusion that pi is badly approximable follows from this note.

## Verification

`verify_endpoint_pencil.py` uses only Python's standard library and exact Gaussian rational arithmetic. It checks the commutator inverse-conjugacy, the two determinant representations, moment translation congruence, reciprocal and self-inversive identities, both determinant identities for M_n, and its real trace at finite imaginary rational parameters. It also produces sufficient finite eigenvalue-count certificates using t < 44/7 and the explicit L_n. These finite checks audit the algebra; the all-size conclusions are proved above.


## Subsequent refinement in this package

The [exterior profile](../exterior-endpoint-cost/EXTERIOR_SINGULAR_PROFILE.md) strengthens the asymptotic endpoint-root count to `(1-alpha)n-O_alpha(n/log n)`. The [nearest-root synthesis](../nearest-root-scale/NEAREST_ROOT_SCALE.md) also uses the independently proved large trace to sharpen the nearest endpoint-root distance to `log delta_n=-2n log n+O(n)`. Its separate argument locates nearby complex Z roots; the endpoint-root correspondence alone does not do this.
