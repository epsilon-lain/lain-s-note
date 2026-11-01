# Escape-count Gram identities and a spectral form of the period-minor envelope

2026-11-01. The statements below are finite linear algebra, followed by explicitly qualified fixed-parameter asymptotics. They organize all escape counts, including a vanishing escape fraction. They do not estimate the required period spectrum, and do not prove a Diophantine property of pi.

## 1. Normalized matrices and the invertibility hypothesis

After removing the inside unit-triangular transport block, write the complete matrix as

\[
A=B+RC,
\]

where B is the M-by-M source-period matrix on one fixed family of raw columns, C is the n-by-M matrix of allowed outside period rows, and R is the M-by-n boundary transport matrix. All matrices are complex. Every matrix here uses the same ordinary-Taylor row convention; arbitrary factorial row rescaling must be carried into both factors.

For the source-relative formulas assume B is invertible, and define

\[
X=CB^{-1},\qquad A=(I_M+RX)B.
\tag{1}
\]

This hypothesis is a genuine restriction. An arbitrary column family can have a singular B, even if A is nonsingular. No interpolation theorem is asserted to supply B's invertibility in the desired high-dimensional family. When B is singular the first, unnormalized Gram identity below remains valid, but the spectral normalization by X is unavailable.

Index source rows by L={1,...,M} and outside rows by O={1,...,n}. A full target set J replaces some source subset I of size k by an outside subset S of size k. Right multiplication by B^{-1} turns the full period-row matrix into [I_M; X]. Therefore

\[
\frac{\det B_J}{\det B}=\pm\det X_{S,I}.
\tag{2}
\]

The sign is the row-order permutation sign. Here B_J denotes the chosen rows from [B;C], not a submatrix of B alone. The target's escape count is k=|J\setminus L|.

## 2. Exact generating functions for every escape count

For u>0, Cauchy--Binet applied to [B;sqrt(u) C] gives

\[
\boxed{
\sum_{|J|=M}u^{|J\setminus L|}|\det B_J|^2
=\det(B^*B+uC^*C)
=|\det B|^2\det(I_M+uX^*X).}
\tag{3}
\]

The first equality does not require B invertible. The second uses (1). The sum is over every target set in this finite ambient row set; it does not impose nonzero transport coefficients or bucket-count constraints.

Let sigma_1>=...>=sigma_M>=0 be the singular values of X, padded by zeros if needed. Let e_k denote the kth elementary symmetric polynomial, with e_0=1. Comparing powers of u in (3) gives the exact identity

\[
\boxed{
\sum_{|J\setminus L|=k}|\det B_J|^2
=|\det B|^2 e_k(\sigma_1^2,\ldots,\sigma_M^2).}
\tag{4}
\]

Likewise, the maximal minors of [I_M,R] satisfy

\[
\sum_{|J\setminus L|=k}|\det[I_M,R]_J|^2
=e_k(\lambda_1,\ldots,\lambda_M),
\tag{5}
\]

where lambda_i are the eigenvalues of RR^*. This includes zero eigenvalues. The total escape-count generating function is det(I_M+u RR^*).

Introduce a formal marker eta in the actual determinant:

\[
\det(I_M+\eta RX)=\sum_{k=0}^{\min(M,n)}d_k\eta^k.
\]

The same expansion, or Cauchy--Binet on the kth exterior power, yields

\[
d_k=\sum_{|I|=|S|=k}\det R_{I,S}\det X_{S,I}.
\tag{6}
\]

Consequently a separate Cauchy--Schwarz inequality at each count gives

\[
\boxed{|d_k|\le
\sqrt{e_k(RR^*)e_k(X^*X)},\qquad
|\det A|\le|\det B|
\sum_k\sqrt{e_k(RR^*)e_k(X^*X)}.}
\tag{7}
\]

The notation e_k(G) means the elementary symmetric polynomial of G's eigenvalues, equivalently the sum of its k-by-k principal minors. No separate factorial count or rowwise Hadamard estimate has been introduced.

For comparison, a single weighted Cauchy--Schwarz step gives, for every u>0,

\[
\boxed{
|\det A|^2\le|\det B|^2
\det(I_M+u^{-1}RR^*)\det(I_M+uX^*X).}
\tag{8}
\]

Applying scalar Cauchy--Schwarz to the sum in (7) proves that (7) is at least as strong as every bound (8), including its infimum over u. Neither bound retains the relative phases and common support of the individual matching minors in (6).

### A version that allows a singular source matrix

Assume only that the full stacked matrix [B;C] has column rank M, so G=B^*B+C^*C is positive definite. Set

\[
Q=G^{-1/2}B^*BG^{-1/2},\qquad 0\le Q\le I,
\]

and let lambda_1,...,lambda_M be its eigenvalues. Factoring G^(1/2) on both sides gives

\[
\boxed{
\sum_{k=0}^M p_k u^k
=\det G\prod_{i=1}^M[\lambda_i+u(1-\lambda_i)],
\qquad
p_k=\sum_{|J\setminus L|=k}|\det B_J|^2.}
\tag{8a}
\]

Thus p_k/det G is exactly the mass at k of a sum of independent Bernoulli variables with success probabilities 1-lambda_i. This is an algebraic description of the coefficients, not a probabilistic hypothesis about pi. Eigenvalues lambda_i=0 force an outside selection, while lambda_i=1 contribute no outside selection. Their multiplicities are M-rank(B) and M-rank(C), respectively, so p_k=0 outside M-rank(B)<=k<=rank(C). If the full stacked matrix has rank less than M, all p_k vanish and every period minor is zero.

For N_k=binom(M,k)binom(n,k), direct maximum-versus-sum bounds give

\[
\sqrt{p_k/N_k}\le\max_{|J\setminus L|=k}|\det B_J|
\le\sqrt{p_k}
\tag{8b}
\]

whenever such sets exist. Their logarithmic gap is O(M), uniformly k in the fixed-weight packet. This supplies a singular-source-safe form of the unrestricted envelope. It still requires estimates on the actual period Gram pencil; whitening has not supplied those estimates. When B is invertible, (3)--(4) recover the source-relative form without whitening.

## 3. The unrestricted period-minor envelope is a cumulative singular-value sum

For 0<=k<=min(M,n), define

\[
F_k=\max_{|J\setminus L|=k}\frac{|\det B_J|}{|\det B|},
\quad N_k=\binom Mk\binom nk,
\quad S_k=\prod_{i=1}^k\sigma_i,
\]

with empty products equal to one. If k exceeds rank(X), F_k=S_k=0. For every other k,

\[
\boxed{\frac{S_k}{\sqrt{N_k}}\le F_k\le S_k.}
\tag{9}
\]

For the upper bound, every k-by-k minor is a matrix entry of the kth exterior power of X in orthonormal coordinate bases. Its operator norm is the product of X's k largest singular values: this follows immediately by applying the singular-value decomposition to the exterior power. A matrix entry cannot exceed that norm. For the lower bound, (4) and the spectral identity give

\[
N_kF_k^2\ge e_k(\sigma^2)\ge S_k^2.
\]

These are finite bounds, uniform over all counts. In the previous fixed-weight packet, M is asymptotic to a positive constant times H^(m+1), and M+n=O(H^(m+1)). Thus log N_k<= (M+n) log 2=O(M), uniformly k. It follows that

\[
\boxed{
\frac{\log\max_{|J\setminus L|=k}|\det B_J|}{MH}
=\frac{\log|\det B|}{MH}
+\frac1{MH}\sum_{i=1}^k\log\sigma_i(X)+O(H^{-1})}
\tag{10}
\]

uniformly for k<=rank(X). The error is one-sided before being written as O(H^{-1}): the left side is no larger than the displayed spectral expression.

Equation (10) concerns all period target sets. For the smaller class with nonzero transport coefficients, it gives only an upper bound. The maximizing coordinate minors may violate center/bucket conservation or may have zero Pascal transport determinant. No lower bound for the transport-compatible maximum follows from (10).

Combining this upper bound with the positive-fraction transport estimate in the preceding research packet gives a concrete sufficient upper target for the complete sum:

\[
\frac{\log|\det A|}{MH}
\le\frac{\log|\det B|}{MH}+\frac{\log2}{w_*}
+\max_{0\le k\le\operatorname{rank}X}
\left\{\frac1{MH}\sum_{i=1}^k\log\sigma_i(X)
-\frac{\tau}{(m+1)\mathcal A}(k/M)^2\right\}+o(1).
\tag{11}
\]

Here w_*, tau and mathcal A have exactly the fixed-parameter meanings from [POSITIVE_FRACTION_LAW.md](../pi-positive-fraction-2026-11-01/POSITIVE_FRACTION_LAW.md). This is for the normalized determinant after the inside factor is removed; any earlier diagonal arithmetic factors must be restored separately. Its identity term can be split off to avoid charging log(2)/w_* to k=0.

## 4. A one-parameter convex transform records the same spectrum

For positive sigma_i put a_i=log(sigma_i)/H. For a zero singular value set a_i=-infinity. At u=exp(-2tH), the normalized relative Gram logarithm is

\[
\frac{\log\det(I_M+uX^*X)}{2MH}
=\frac1M\sum_i(a_i-t)_++\epsilon_H(t),
\qquad 0\le\epsilon_H(t)\le\frac{\log2}{2H}.
\tag{12}
\]

This uniform bound follows term by term from
0<=log(1+exp(y))-max(y,0)<=log 2. The zero singular values contribute zero. Also, for k<=rank(X),

\[
\sum_{i=1}^k a_i
=\inf_{t\in\mathbb R}\left(kt+\sum_i(a_i-t)_+\right).
\tag{13}
\]

Choose t between the kth and (k+1)st sorted values to obtain equality, and compare terms on either side for the opposite inequality. The endpoint cases k=0 and k=rank(X) follow by placing t beyond the finite values. Equations (12)--(13) provide a scalar convex route to the same unrestricted envelope, but do not estimate the spectrum itself.

## 5. A precise compatibility loss that scalar weighting cannot repair

The following is a generic linear-algebra obstruction, not a constructed period packet. Let n=2M, r,s>0, B=I_M and

\[
R=[rI_M\ \ 0],\qquad
X_{\rm aligned}=\begin{bmatrix}sI_M\\0\end{bmatrix},
\qquad
X_{\rm disjoint}=\begin{bmatrix}0\\sI_M\end{bmatrix}.
\]

Both choices have exactly the same RR^*=r^2 I and X^*X=s^2 I. They therefore have the same escape-count Gram data in (4)--(5), the same spectral envelope, and the same bounds (7)--(8). But

\[
\det(I+RX_{\rm aligned})=(1+rs)^M,\qquad
\det(I+RX_{\rm disjoint})=1.
\tag{14}
\]

The fixed-count bound (7) equals (1+rs)^M; the optimized scalar-weight bound (8) has the same value, attained at u=r/s. Thus even the full collection of escape-count squared norms cannot distinguish an aligned example attaining the bound from a disjoint-support example whose entire positive-degree part is zero.

For r=2^(-H), s=2^(2H), the lost logarithm divided by MH is log(1+2^H)/H, tending to log 2. The loss occurs on the main MH scale despite a small transport norm. A method using only these two Gram spectra cannot improve its universal upper bound below the aligned example. An improvement for actual period rows must use their additional structure or the common matching-minor support in (6).

## 6. A vanishing escape fraction can still matter without spectral control

Here is a separate generic warning, again not a period-packet counterexample. For positive integers M>=2,H, let

\[
B=\operatorname{diag}(2^{-MH},1,\ldots,1),
\quad C=e_1^T,\quad R=2^{-H}e_1.
\]

Every B and C entry has modulus at most one, and R has norm 2^(-H). Yet X=2^(MH)e_1^T has a single large singular value. Directly,

\[
\det(B+RC)=2^{-MH}+2^{-H},\qquad
\frac{\det(B+RC)}{\det B}=1+2^{(M-1)H}.
\tag{15}
\]

Only one outside row is available, so this increase comes from escape fraction 1/M. As M,H grow, its normalized logarithmic increase tends to log 2. Bounded entries, small R, or a vanishing escaped-row fraction do not suffice to control a source-relative determinant.

A sufficient way to exclude this phenomenon is a bound log^+(sigma_1(X))<=C H independent of H: then the cumulative positive gain for k=o(M) is at most C k/M and vanishes. More generally one needs suitable uniform integrability of the positive normalized logarithmic singular values. Such a condition has not been established for the required actual-period family here. This source-relative warning does not contradict an absolute upper bound based on the Hardy evaluation Gram matrix: a chosen B can be far smaller than that absolute bound.

## 7. A support-aware diagonal refinement

The generic disjoint-support example is an obstruction to scalar escape weighting or the separate fixed-count squared norms, not to every weighted Gram method. Let E and D be positive diagonal matrices of sizes M and n. Replace R and X by

\[
\widetilde R=E^{1/2}RD^{-1/2},\qquad
\widetilde X=D^{1/2}XE^{-1/2}.
\]

Their product is similar to RX. Applying (8) at u=1 gives

\[
\boxed{
|\det(I+RX)|^2\le
\det(I+E^{1/2}RD^{-1}R^*E^{1/2})
\det(I+E^{-1/2}X^*DXE^{-1/2}).}
\tag{16}
\]

The logarithm of the right side is convex in the logarithms of the diagonal entries of E and D. To prove this without any matrix-convexity assertion, expand each determinant into its principal minors and then use Cauchy--Binet. Every term is a nonnegative squared minor times a monomial in those diagonal entries. Taking logarithms yields a log-sum-exp of affine functions, which is convex; the sum of the two such expressions remains convex. The identity term ensures neither determinant vanishes.

This refinement retains which source and outside coordinates participate in each minor through their separate weights. In the disjoint example (14), take E=I, D=a I on the first M outside coordinates and D=b I on the last M. Let a tend to infinity and b to zero. Both determinants in (16) tend to one, exactly recovering the actual value. The analogous scalar weighting could not do this.

Equation (16) is therefore a possible support-aware variational target. It still discards complex phases, still depends on X=CB^{-1}, and has not been estimated in the intended arithmetic regime. No claim that its infimum equals the actual determinant for general matrices is made.

## 8. What is new relative to the earlier Hardy Gram bound

The earlier grouped Hardy factorization retained correlations among evaluation rows in an infinite coefficient space. The present finite identities instead grade all actual compatible-column period minors by their number of outside rows. They identify the exact object behind the unresolved maximum: cumulative singular values of C B^{-1}, when B is invertible. They also separate the fixed-count Cauchy--Schwarz estimate from the weaker single-weight estimate.

This is a reduction of the missing estimate, not a resolution. The same inverse-period conditioning cost remains visible in X. The generic examples show precisely which kinds of information are discarded if one uses only separate transport and period Gram spectra. They say nothing by themselves about whether those obstructions occur in the arithmetic-period structure.

Run `python verify_escape_gram.py`. The standard-library verifier checks exact rational generating functions, all matching minors, marker coefficients, weighted inequalities, a rational matrix with known singular values, the arbitrary alignment loss, and the single-exit spectral spike. It does not infer an infinite-period estimate from finite computations.
