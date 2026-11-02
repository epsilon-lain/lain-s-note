# Exterior column energy sharpens the genuine thermal kernel bound

2026-11-02. This note refines the full weighted-column kernel introduced in the preceding checkpoint. It gives an exact determinant-level energy bound, its weighted-simplex asymptotics, and a genuine complete-jet example in which the old thermal cutoff loses a leading MH exponent while the refinement removes that loss. It does not settle the near-total-escape estimate for the critical pi construction.

## 1. Keep the sum of the selected column degrees

Let W=(w_0,...,w_m) be positive weights. A raw column is indexed by an exponent vector nu=(h,alpha) in the nonnegative integer lattice, with cost c(nu)=W dot nu. Let v_nu be that column's vector on any fixed M-row target selection J. Define

\[
G_H[J]=\sum_{c(\nu)\le H}v_\nu v_\nu^*,\qquad
G_\infty(q)[J]=\sum_\nu q^{c(\nu)}v_\nu v_\nu^*,
\quad 0<q<1,
\]

where the second sum and the required finite jets converge. In the actual-period polynomial setting this is exactly the product kernel from `DIAGONAL_PROFILES_AND_COLUMN_KERNEL.md`, with its cumulative 1/(1-q) factor removed after integralizing the weights. Rational weights can instead be used directly in q^(c(nu)).

For a finite set C of M distinct columns, put E(C)=sum_{nu in C}c(nu). Cauchy--Binet gives

\[
\det G_\infty(q)[J]
=\sum_{|C|=M}q^{E(C)}|\det B_{J,C}|^2.
\tag{1}
\]

This follows first for finite truncations and then by monotone convergence of the positive sums and convergence of the finite Gram matrices. Consequently

\[
\boxed{|\det B_{J,C}|^2\le q^{-E(C)}\det G_\infty(q)[J].}
\tag{2}
\]

Define E_max(H,M) to be the sum of the M largest costs among all columns with c(nu)<=H. When there are fewer than M columns all M-row minors vanish. Otherwise the stronger full-Gram inequality is

\[
\boxed{\det G_H[J]\le
q^{-E_{\max}(H,M)}\det G_\infty(q)[J].}
\tag{3}
\]

Indeed, each admissible C has E(C)<=E_max and hence q^(E(C))>=q^(E_max). Restricting the sum (1) to admissible C proves (3). This is a determinant inequality; it does not assert a corresponding matrix Loewner inequality with q^(-E_max/M). The earlier matrix comparison q^(-H)G_infinity gives the looser determinant cost q^(-MH).

Equations (2)--(3) use the actual raw-column energy once. If columns have been arithmetically scaled or mixed, those operations must be included in v_nu and in the same Cauchy--Binet expansion; a mixed column does not have a fictitious single raw degree.

### Anisotropic version

For a positive vector lambda=(lambda_0,...,lambda_m), replace q^(c(nu)) by exp(-lambda dot nu), in its convergence domain. The analogue of (2) has prefactor exp(sum_{nu in C}lambda dot nu). The uniform prefactor in (3) is the maximum of this sum over allowed M-column sets. Its logarithm is a finite maximum of linear functions of lambda. Also log det G_infinity(lambda) is a log-sum-exp of affine functions by (1). Thus their sum is convex in lambda whenever the Gram determinant is positive. This gives a column-weight optimization that retains genuine column compatibility; it does not estimate the optimum in the near-total-escape regime.

## 2. Maximum exterior energy in a weighted simplex

Write d=m+1 for the number of raw exponent coordinates, and let N_H count the columns of cost at most H. At fixed positive rational weights,

\[
N_H\sim\frac{H^d}{d!\prod_{i=0}^m w_i}.
\]

Suppose M_H/N_H tends to r in (0,1]. Then

\[
\boxed{\frac{E_{\max}(H,M_H)}{M_HH}
\longrightarrow
e_d(r):=\frac{d}{d+1}
\frac{1-(1-r)^{(d+1)/d}}{r}.}
\tag{4}
\]

The normalized cost x=c(nu)/H has limiting distribution with density d*x^(d-1) on [0,1]. The most expensive fraction r consists of x above a=(1-r)^(1/d). Integrating x over that tail gives (4). Weighted-lattice cubes give an O(H^(d-1)) count error and an O(H^d) energy error, so ties at the cutoff do not affect the limit. The formula also holds at r=1 by continuity or by summing the whole simplex.

In the prior critical source-count choice v_0=2K theta^m w_0, the number of source rows satisfies M_H/N_H->1/2. Hence

\[
e_d(1/2)=\frac d{d+1}(2-2^{-1/d})
=1-\frac{1-\ln2}{d}+O(d^{-2}).
\tag{5}
\]

The exterior-energy saving over MH is positive at each fixed dimension. In this subsequent large-d comparison it is only of order 1/d. No claim that this alone supplies the missing endpoint margin is made.

## 3. A complete genuine-column model with an exact thermal gap

Take d=m+1 variables (Y,X_1,...,X_m), one center j=0, all raw weights and source weights equal to one, and H=n a positive integer. The complete source packet consists of all ordinary time/transverse jets with s+|beta|<n. The raw column space consists of all monomials Y^h X^alpha with h+|alpha|<=n. Thus

\[
M=\binom{n+d-1}{d}.
\]

The actual local substitution is Y=e^z, X_i=z+u_i. Change local variables to T=Y-1 and X_i. On all jets of total order less than n this is an invertible triangular coordinate change of determinant one. One can see this by first applying the linear shear u_i->u_i+z and then the one-variable series e^z-1=z+O(z^2); each induced truncated-jet transformation has diagonal one in a compatible ordering. Therefore all full-row Gram determinants are unchanged by this local coordinate change.

In the direct coordinates (T,X), a row (s,beta) on a column Y^h X^alpha is delta_(alpha,beta)*binom(h,s). Put

\[
n_\beta=n-|\beta|\quad(|\beta|<n).
\]

Each type beta therefore contributes an n_beta-by-(n_beta+1) Pascal block to the finite column matrix. Its maximal minors are binom(n_beta,k), 0<=k<=n_beta, up to permutation signs. To verify this, adjoin its last row to obtain the square unit-triangular Pascal matrix; its inverse has entries (-1)^(h-s)binom(h,s). The complementary-minor formula gives the displayed minors. Their squared sum is binom(2n_beta,n_beta). Thus

\[
\boxed{\det G_n=
\prod_{|\beta|<n}\binom{2n_\beta}{n_\beta},\qquad
1\le\det G_n\le4^M.}
\tag{6}
\]

In particular log(det G_n)/(2Mn) tends to zero at fixed d.

For the infinite geometric column kernel, the direct-coordinate block at beta is q^(|beta|) times the Gram of binom(h,s), h>=0, with weights q^h. The confluent Cauchy identity gives

\[
\det\left[\sum_{h\ge0}q^h\binom hs\binom ht\right]_{0\le s,t<l}
=\frac{q^{l(l-1)/2}}{(1-q)^{l^2}}.
\]

For a direct proof, the generating kernel is [1-q(1+x)(1+y)]^(-1). Apply the ordinary Cauchy determinant formula to l distinct nearby evaluation points, divide out their Vandermonde factors, and let them coalesce at one. The resulting normalized derivative determinant has the stated numerator and denominator, with no derivative factorials left over.

Multiplying the blocks yields

\[
\boxed{\det G_\infty(q)=q^{E_0}(1-q)^{-T},\quad
E_0=\frac{d(n-1)}{d+1}M,\quad
T=\frac{2n+d-1}{d+1}M.}
\tag{7}
\]

Here E_0 is the total degree of all monomials of degree less than n. Also T=sum_beta n_beta^2. The identities follow by the standard composition count: E_0=d*binom(n+d-1,d+1) and T=M+2*binom(n+d-1,d+1).

The old thermal bound q^(-Mn)det G_infinity(q) is minimized at

\[
q_* =\frac{n+d}{3n+2d-1}.
\]

Its leading normalized logarithm is

\[
\boxed{\lim_{n\to\infty}\frac1{2Mn}
\ln\inf_{0<q<1}\left[q^{-Mn}\det G_\infty(q)\right]
=\frac{3\ln3-2\ln2}{2(d+1)}>0.}
\tag{8}
\]

Thus even after optimizing the scalar thermal parameter, the old cutoff loses a factor exp(c_d MH) against this exact genuine-column Gram determinant. This is a complete fixed-weight polynomial/jet example, rather than an arbitrary matrix. However theta=1 and v_i=v_0 here, so there is no positive boundary transport cost. It is not a counterexample to the critical separated-weight pi construction or to a refined kernel estimate there.

## 4. Exterior energy removes the leading loss in this model

There are binom(n+d-1,d-1) new columns of exact degree n. Starting from all columns of degree below n, their total added degree is exactly dM. Since excluded column costs are nonnegative,

\[
E_{\max}(n,M)\le E_0+dM.
\]

Apply (3) with q=1/(n+1). Since T=sum n_beta^2<=nM and ln(1+1/n)<=1/n,

\[
\boxed{\frac1{2Mn}\ln\left[
q^{-E_{\max}}\det G_\infty(q)\right]
\le\frac{d\ln(n+1)+1}{2n}\longrightarrow0.}
\tag{9}
\]

The corrected thermal estimate therefore agrees with the exact leading exponent in (6). The decisive change is accounting for M distinct columns' actual total energy, rather than assigning the maximum possible degree to every column independently. It is not a new cancellation theorem.

## 5. What this contributes to the near-total-escape problem

The genuine product kernel remains a useful exact object, but its cutoff penalty should be applied at exterior level as in (2)--(3). The repair can matter on the MH scale and can remove a whole leading loss in a complete analytic model. In the earlier critical count ratio r=1/2, the uniform distinct-column improvement is quantified by (5); it is not by itself an endpoint proof.

The remaining determinant in (3), after using these sharper column energies, still needs a uniform estimate for the actual near-total-escape rows. Its relation to error powers and the same columns' arithmetic degree differences must be kept. All statements about pi's bad approximability remain open in this work.

Run `python verify_exterior_energy.py`. The standard-library verifier checks the full local jet coordinate changes, exact finite and infinite Gram determinants, the optimized old kernel formula, the repaired exterior bound, and finite weighted-simplex energy selections. Infinite limits follow from the written formulas, not extrapolation from the checks.
