# The entire endpoint singular spectrum and its exterior cost

2026-11-02. This note strengthens the preceding normalized endpoint norm theorem to all singular values. Exact Pascal factorizations show that a positive fraction of directions has factorial size. Exterior powers improve the bound obtained by taking powers of the largest singular value. Together with the exact determinant, they also force a stronger count of actual endpoint-pencil zeros near the period parameter. None of these statements proves that pi is badly approximable.

## 1. Definitions and conclusions

Indices run from 0 to n-1; n>=2. Put t=2*pi, omega=i*t, and

\[
J=[1/(r+s+1)],\quad S=J^{-1/2},\quad
u_k(Z)=(-1)^k k!Z^{-k-1},\quad U(Z)=[u_{r+s}(Z)],
\]
\[
\nu_k(Z)=\sum_{\ell=0}^k\frac{(-1)^\ell k!}{(k-\ell)!Z^{\ell+1}},
\qquad V(Z)=[\nu_{r+s}(Z)].
\]

The actual moment matrix satisfies

\[
H(\omega)=\left[\int_0^1u^{r+s}e^{\omega u}\,du\right]=V(\omega)-U(\omega).
\]

Write

\[
\mathsf A=SH(\omega)S,\qquad \mathsf B=SV(\omega)S,
\qquad \mathsf C=\mathsf A^{-1}\mathsf B.
\]

Singular values are in decreasing order. Our first result is

\[
\boxed{\log\sigma_j(\mathsf B)=2(n-j)\log n+O(n),
\quad 1\le j\le n,}
\tag{1}
\]

with an error uniform in j. The same holds for C. Consequently, uniformly for 1<=k<=n,

\[
\boxed{\log\|\wedge^k\mathsf B\|_2
=k(2n-k-1)\log n+O(nk),}
\tag{2}
\]

again with the same result for C. These are statements about singular values and operator norms, not eigenvalue asymptotics or noncancellation of individual determinant coefficients. Many directions are large on an absolute scale; the theorem does not rule out a low-rank approximation whose error is measured relative to the largest singular value, and does not assert proportional stable rank.

## 2. Two exact Pascal factorizations

Let P be the lower triangular Pascal matrix and D_f the factorial diagonal:

\[
P_{rs}=\binom rs\quad(s\le r),\qquad
(P^{-1})_{rs}=(-1)^{r-s}\binom rs,
\qquad D_f=\operatorname{diag}(r!).
\]

Finite binomial convolution gives

\[
\boxed{V(Z)=P U(Z)P^T.}
\tag{3}
\]

Indeed the (r,s) entry of the right side is
sum_l binom(r+s,l)(-1)^l l! Z^(-l-1), exactly nu_(r+s).
Likewise Vandermonde's binomial identity gives

\[
F=[(r+s)!]=D_f P P^T D_f.
\tag{4}
\]

To verify (4), divide row r and column s by r! and s!: the required identity is sum_l binom(r,l)binom(s,l)=binom(r+s,r).

Let Q=diag(i^r), D_t=diag(t^(-r)), and D_q=diag(q_r), where q_r=r!/t^r. Then

\[
\boxed{U(it)=\frac1{it}Q K Q,
\qquad K=D_t F D_t=(D_qP)(D_qP)^T.}
\tag{5}
\]

Here Q is unitary and K is real positive definite. It follows that the singular values of U(it) are t^(-1) times the squared singular values of D_qP. All these identities hold for arbitrary real t>0, not only t=2*pi.

## 3. Explicit bounds for every singular value

The squared Frobenius norms of P and P^(-1) are equal, and

\[
\|P\|_2^2,\ \|P^{-1}\|_2^2
\le\sum_{r=0}^{n-1}\binom{2r}{r}
\le\sum_{r=0}^{n-1}4^r<4^n.
\tag{6}
\]

The Hilbert estimates from the preceding Gram proof are

\[
\|J\|_2\le\operatorname{tr}J\le n,
\qquad \|J^{-1}\|_2\le2n^2 2^{8n}.
\tag{7}
\]

For completeness the inverse diagonal equals
(2r+1) binom(n+r,n-r-1)^2 binom(2r,r)^2, bounded by 2n*2^(8n); positivity bounds the norm by the trace.

Let a_1>=...>=a_n be the decreasing rearrangement of the positive q_r. For invertible X,Y, singular-value min-max gives

\[
\frac{\sigma_j(M)}{\|X^{-1}\|_2\|Y^{-1}\|_2}
\le\sigma_j(XMY)\le\|X\|_2\|Y\|_2\sigma_j(M).
\tag{8}
\]

Apply this first to D_qP, then use (5), then the Pascal congruence (3), and finally Gram normalization. The result is the fully explicit estimate

\[
\boxed{
\frac{a_j^2}{t n 2^{4n}}
\le\sigma_j(\mathsf B)
\le\frac{2n^2 2^{12n}}t a_j^2.
}
\tag{9}
\]

The two factors of 2^(2n) before Gram normalization come respectively from squaring the D_qP comparison and from the additional congruence by P. They must both be retained.

For fixed t>1 and 0<=r<n, an entirely rationalizable bound is

\[
(3t)^{-n}n^r\le q_r\le n^r.
\tag{10}
\]

The single positive exponential-series term n^r/r! is at most e^n<3^n. Hence r!>n^r/3^n; also t^r<=t^n. The upper bound follows from r!<=n^r and t>1. Sorting preserves coordinatewise inequalities between vectors, so (10) gives

\[
(3t)^{-n}n^{n-j}\le a_j\le n^{n-j}.
\tag{11}
\]

It is not necessary, and for small indices not correct, to assume that q_r is increasing. Equations (9)--(11) prove (1) with an absolute error bounded by a constant times n, uniformly in j. Taking the product of the k largest singular values proves (2), since
2*sum_(j=1)^k(n-j)=k(2n-k-1).

## 4. The actual period controls inversion

The earlier actual-period moment estimate proves

\[
\|\mathsf A\|_2\le1,\qquad
\frac16<|\det\mathsf A|\le1.
\tag{12}
\]

The first inequality is the Gram compression bound for multiplication by e^(itu). The second is precisely the normalized two-sided period theorem in `../../pi-period-mechanism-2026-11-02/raw-root-mechanism/HANKEL_PERIOD_MAGNITUDE.md`. Because every singular value of A is at most 1, its smallest singular value is at least |det A|. Thus A is invertible and ||A^(-1)||_2<6. Equation (8) now gives

\[
\boxed{\sigma_j(\mathsf B)\le\sigma_j(\mathsf C)<6\sigma_j(\mathsf B).}
\tag{13}
\]

This is a comparison in the fixed Gram geometry. The matrix C is similar to H(omega)^(-1)V(omega), but similarity by itself would not preserve singular values.

## 5. Determinants, coefficients, and a useful envelope

Set a factorial-product symbol distinct from the matrices:

\[
\mathcal F_m=\prod_{r=0}^{m-1}r!,\qquad
\Delta_n=\det J=\frac{\mathcal F_n^4}{\mathcal F_{2n}},\qquad
\mathcal U_n=t^{n^2}\frac{\mathcal F_n^2}{\mathcal F_{2n}}.
\]

The Pascal determinant is 1. From (3)--(5), or direct factorial determinants,

\[
|\det V(it)|=|\det U(it)|=\frac{\mathcal F_n^2}{t^{n^2}},
\qquad
\boxed{|\det\mathsf B|=\mathcal U_n^{-1}.}
\tag{14}
\]

Consequently

\[
\mathcal U_n^{-1}\le|\det\mathsf C|<6\mathcal U_n^{-1},
\qquad
\log|\det\mathsf C|=n^2\log n+O(n^2).
\tag{15}
\]

The endpoint perturbation factorizes as

\[
\det(\mathsf A+\eta\mathsf B)
=\det\mathsf A\det(I+\eta\mathsf C),
\quad
\det(I+\eta\mathsf C)
=\sum_{k=0}^n\eta^k\operatorname{tr}(\wedge^k\mathsf C).
\tag{16}
\]

In particular

\[
|\operatorname{tr}(\wedge^k\mathsf C)|
\le\binom nk\|\wedge^k\mathsf C\|_2.
\tag{17}
\]

The trace-norm version gives the stronger complete product bound

\[
\boxed{|\det(I+\eta\mathsf C)|
\le\prod_{j=1}^n(1+|\eta|\sigma_j(\mathsf C))
\le\prod_{r=0}^{n-1}(1+|\eta|K_n n^{2r}),
\quad K_n=\frac{12n^2 2^{12n}}t.}
\tag{18}
\]

Indeed the singular values of the k-th exterior power are all products of k singular values; its trace norm is their sum. Summing these trace-norm bounds in (16) proves the first product inequality. Equations (9), (11), and (13) prove the second.

The sharper single-norm baseline is log(||C||^k)=2k(n-1) log n+O(nk). Against that baseline, the exterior estimate saves k(k-1) log n in the leading exponent. This materially improves high-order upper bounds. At k=1 there is no saving. Lower bounds for exterior norms do not imply lower bounds for traces: cancellation remains possible.

## 6. A stronger count of actual endpoint-pencil zeros

Let mu_1,...,mu_n be the eigenvalues of C, ordered by decreasing absolute value and counted with algebraic multiplicity. They are nonzero by (15). Every product of k eigenvalues is an eigenvalue of the k-th exterior power. Hence the spectral-radius inequality and (2) imply

\[
\log\prod_{j=1}^k|\mu_j|
\le k(2n-k-1)\log n+O(nk).
\tag{19}
\]

Fix 0<alpha<1, and let K be the number of eigenvalues with |mu_j|>=exp(alpha*n*log n). The remaining eigenvalues have smaller absolute value. Combining (15) and (19), including the trivial k=0 case, gives

\[
n^2\log n+O(n^2)
\le K(2n-K-1)\log n+(n-K)\alpha n\log n+O(nK).
\]

Writing y=(n-K)/n and absorbing lower-order terms gives
y(y-alpha)<=O(1/log n). If y>alpha, then y>=alpha, so y-alpha<=O_alpha(1/log n). Therefore

\[
\boxed{K\ge(1-\alpha)n-O_\alpha(n/\log n).}
\tag{20}
\]

This is stronger than the estimate obtained from the determinant and largest singular value alone. It does not assert that the individual eigenvalues have the full singular-value profile (1).

The unnormalized endpoint pencil is

\[
\mathcal H_n(\omega,Y)=\det(YV(\omega)-U(\omega)).
\]

The [companion pencil note](../endpoint-pencil-structure/ENDPOINT_PENCIL_STRUCTURE.md) gives
\(G_n(Z,Y)=Z^{n^2}\mathcal H_n(Z,Y)/\mathcal F_n^2\), so their roots in Y agree for Z!=0. Since H(omega)=V(omega)-U(omega), equations (16) and the similarity of C with H^(-1)V give the exact root correspondence

\[
Y_j=1-\mu_j^{-1}.
\tag{21}
\]

Thus, for fixed 0<alpha<1, at least (1-alpha)n-O_alpha(n/log n) roots of this degree-n polynomial, counted with multiplicity, satisfy

\[
\boxed{|Y_j-1|\le e^{-\alpha n\log n}.}
\tag{22}
\]

This conclusion concerns the independent complex endpoint variable Y with Z fixed at omega. It does not locate roots of the one-variable formal-period polynomial P_n(Z), and it does not directly produce a rational approximant to pi. The relation between Y and Z must still be tracked in any such application.

## Verification and scope

The standard-library companion `verify_exterior_profile.py` checks exact Pascal inverses and congruences, the factorial-Hankel factorization, the phase factorization at rational imaginary centers, Hilbert inverse and determinant identities, exact endpoint determinants, and exterior-power algebra in finite dimensions. It uses rational and Gaussian-rational arithmetic, without approximating pi or numerically extracting singular values. The uniform estimates, rearrangement step, exterior bounds, and asymptotic root count follow from the proof above; finite tests are only checks against algebraic and indexing errors.

These results refine the mechanism within the existing two-frequency sparse family. They do not extend it to the critical high-dimensional parameter range, prove a matching lower bound for a determinant coefficient, or settle pi's irrationality measure or bad approximability.
