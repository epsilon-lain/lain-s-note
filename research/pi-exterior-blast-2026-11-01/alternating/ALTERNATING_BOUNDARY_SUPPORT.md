# Weighted alternating boundary sums: exact support and a sharp exponential limitation

2026-11-01. This note analyzes the scalar alternating sum in the exact boundary transport. It does not estimate the compatible mixed determinants and does not prove an approximation property of pi.

Let gamma=(gamma_1,...,gamma_m) be a nonzero nonnegative integer multi-index, let c_i>0 for every active coordinate gamma_i>0, and put

\[
N=|\gamma|,\qquad L=c\cdot\gamma,\qquad
S_\gamma(d)=\sum_{\substack{0\le\kappa\le\gamma\\c\cdot\kappa<d}}
(-1)^{|\kappa|}\binom\gamma\kappa.
\]

The strict inequality matches the source row packet. The boundary coefficient is

\[
W_{\rho,\gamma}=\binom{\beta+\gamma}{\beta}(-\zeta)^\gamma S_\gamma(d_\rho).
\]

Zero-cost active coordinates were already shown to cancel identically. Here all active costs are positive.

## 1. Exact localization to one thin slab

Choose any active coordinate i. Splitting off one of its gamma_i copies gives

\[
\boxed{
S_\gamma(d)=
\sum_{0\le\kappa\le\gamma-e_i}
(-1)^{|\kappa|}\binom{\gamma-e_i}{\kappa}
\mathbf 1_{\{c\cdot\kappa<d\le c\cdot\kappa+c_i\}}.}
\tag{1}
\]

This follows by pairing subsets that exclude or include that copy. Their indicators differ by
\(\mathbf1_{x<d}-\mathbf1_{x+c_i<d}=\mathbf1_{x<d\le x+c_i}\).
Thus the entire interior alternating sum is supported on a slab of width c_i at the threshold. In particular it vanishes if d lies outside the union of these half-open intervals. This is an exact zero test, not an inequality after taking absolute values.

Choose i of smallest active cost. View the remaining N-1 copies as distinct binary choices, with weights at least c_i. The subsets contributing to (1) have total weight in [d-c_i,d), so they form an antichain: a proper inclusion raises the weight by at least c_i and cannot remain in that interval. The elementary antichain bound gives

\[
\boxed{|S_\gamma(d)|\le
\binom{N-1}{\lfloor(N-1)/2\rfloor}.}             \tag{2}
\]

For completeness, each subset A of an antichain lies on |A|!(n-|A|)! permutation-generated maximal chains, where n=N-1. Each maximal chain meets the antichain at most once. Counting all n! permutations gives
\(\sum_A1/\binom n{|A|}\le1\), and hence the displayed cardinality bound. Repeated equal weights cause no problem because the copies, rather than their costs, are distinct.

Equal costs attain (2) at the middle threshold, by the partial alternating-binomial identity. Compared with the earlier 2^N bound this is a polynomial improvement of order sqrt(N), not a uniform extra exponential factor. Inserting (2) yields the valid coefficient bound

\[
|W_{\rho,\gamma}|
\le\binom{\beta+\gamma}{\beta}
\binom{N-1}{\lfloor(N-1)/2\rfloor}|\zeta|^\gamma.
\tag{3}
\]

## 2. Exact zero bands under multiplicity-dependent separation

Suppose h is the largest active coordinate, the lower coordinates contain at least one active copy, and

\[
L_<:=\sum_{i<h}c_i\gamma_i<c_h.
\tag{4}
\]

Let k=ceil(d/c_h)-1 and e=d-kc_h, so 0<e<=c_h. All slices with kappa_h<k contain the entire lower cube and cancel; those with kappa_h>k are empty. Therefore

\[
\boxed{S_\gamma(d)=(-1)^k\binom{\gamma_h}{k}
S_{\gamma_<}(e)}                                \tag{5}
\]

when 0<=k<=gamma_h, and it is zero otherwise. In particular its support is contained in

\[
\bigcup_{k=0}^{\gamma_h}(kc_h,kc_h+L_<].          \tag{6}
\]

This can be iterated if each next largest active cost dominates the total lower **multiplicity-weighted** cost. It supplies large exact gaps in finite problems with bounded multiplicities.

The distinction in (4) matters. Successively separated fixed weights c_h>>c_{h-1} do not imply (4) when gamma grows with H. At fixed weights the lower multiplicities are unbounded as H tends to infinity. Thus (5) cannot be applied uniformly to all boundary indices in the manuscript's order of limits merely from weight separation.

## 3. Generating polynomial and persistent exponential coefficients

Assume now c_i are positive integers. For integer d=1,...,L,

\[
\boxed{S_\gamma(d)=[z^{d-1}]P_\gamma(z),\qquad
P_\gamma(z)=\frac{\prod_i(1-z^{c_i})^{\gamma_i}}{1-z}.}
\tag{7}
\]

The numerator is divisible by 1-z, so P_gamma is a polynomial of degree L-1. Formula (7) is simply the cumulative-sum generating function for the coefficients of the numerator. It also proves the reflection identity

\[
\boxed{S_\gamma(L+1-d)=(-1)^{N-1}S_\gamma(d).}    \tag{8}
\]

If all active c_i are odd, evaluation at z=-1 gives

\[
\sum_{d=1}^L(-1)^{d-1}S_\gamma(d)
=P_\gamma(-1)=2^{N-1}.
\]

Consequently

\[
\boxed{\max_{1\le d\le L}|S_\gamma(d)|
\ge\frac{2^{N-1}}{L}.}                          \tag{9}
\]

By (8) a maximizing threshold can be chosen with d>=(L+1)/2. For fixed positive odd costs and gamma=n gamma^(0), N grows linearly in n and L grows linearly in n. Thus the maximal alternating coefficient has exponential rate exactly log(2) per copy:

\[
\lim_{n\to\infty}\frac1{n|\gamma^{(0)}|}
\log\max_d|S_{n\gamma^{(0)}}(d)|=\log2.           \tag{10}
\]

The upper bound follows already from 2^N, and the lower bound is (9). This holds for arbitrarily widely separated positive odd costs, even with every coordinate active. Hence weight separation alone cannot produce a uniform bound of the form exp((log2-epsilon)N) for all thresholds and multiplicities, for any fixed epsilon>0.

For rational costs the same statements hold whenever c_i=a C_i for a positive common unit a and positive odd integers C_i, with d=a D. Here L in (7)--(9) is the integer sum C dot gamma, and the physical support is scaled by a.

## 4. These large scalar coefficients are not outside the source budget

An explicit admissible family takes theta=1/2, v0=1 and any positive integer weights w_i. Then c_i=2w_i-1 are odd and may be chosen successively as large as desired. Pick gamma, put N=|gamma| and L=c dot gamma, and choose a maximizing upper-half integer d from (9). Set

\[
\beta=0,\quad s=N,\quad H=N+d.
\]

The source row has cost s=H-d<H; its transported target (0,gamma) has cost
\(2w\gamma=N+L\ge H\). The target is within the available polynomial type support because

\[
w\gamma=(N+L)/2\le N+d=H.
\]

Also |gamma|=s. Thus the witnesses obey the time-order, boundary, and column-support restrictions used in the scalar coefficient bound. The corresponding monomial column X^gamma has a nonzero target coefficient. This does not assert that a full high-dimensional interpolation minor with these entries has any specific determinant size, or that these parameters certify the manuscript's rank threshold.

The example generalizes at the scalar-budget level to any rational v0=A/B>0 with theta=1/2: write A=2^t A_odd and choose integer w_i divisible by 2^t. Then
\(c_i=(2^t/B)(2B(w_i/2^t)-A_odd)\), so sufficiently large chosen weights have positive odd costs in a common unit. The choices can meet any prescribed successive lower thresholds. This does not show that the denominator-dependent weights ceil(log q_i) of an actual pi approximation family fall into those arithmetic progressions; that extra claim is not made.

## 5. Implication for the proposed explosion

There are exact local zero bands and a sharp antichain improvement. However, there is no uniform exponential cancellation from the scalar alternating sum alone under mere separation of weights. The large scalar coefficients in (9) can still cancel when multiplied by compatible period rows and assembled into the determinant. That determinant-level compatibility is where an additional decisive mechanism must occur; (9) does not rule it out.

`verify_alternating_support.py` checks the thin-slab identity, antichain bound, dominance reduction, generating polynomial, reflection, odd-cost lower bound and feasible upper-half source witnesses with exact integer/rational arithmetic. Infinite limits are established by the displayed inequalities, not inferred from the finite checks.

## 6. In one transverse dimension the entire boundary block is Cauchy

This section gives a determinant-level formula for the transport coefficients themselves. It does not replace the compatible period row determinants multiplying them.

Let m=1 and c=w/theta-v0>0. Fix the preserved total degree n=s+beta. The source indices on this layer are beta=0,...,r, where

\[
r=\min\left(n,\left\lceil\frac{H-v_0n}{c}\right\rceil-1\right).
\]

Consider a nonempty layer with 0<=r<n; its outside indices are a=r+1,...,n. The one-dimensional boundary formula simplifies to

\[
\boxed{R_{\beta,a}
=(-1)^{a+r}\zeta^{a-\beta}
\frac{(r+1)\binom r\beta\binom a{r+1}}{a-\beta}.}
\tag{11}
\]

To see this, the source row's first-exit cutoff is r-beta, so its boundary entry is
\(\binom a\beta(-\zeta)^{a-\beta}(-1)^{r-\beta}
\binom{a-\beta-1}{r-\beta}\).
Cancelling factorials gives (11). Thus R is a diagonal matrix times the Cauchy matrix 1/(a-beta) times another diagonal matrix.

For increasing row indices B=(beta_1,...,beta_k) in {0,...,r} and increasing target indices A=(a_1,...,a_k) in {r+1,...,n}, write V(B) and V(A) for their positive Vandermonde products. Then

\[
\boxed{
\begin{aligned}
\det R_{B,A}
={}&(-1)^{\sum_j a_j+kr+k(k-1)/2}
\zeta^{\sum_j a_j-\sum_i\beta_i}(r+1)^k\\
&\quad\times
\left(\prod_i\binom r{\beta_i}\right)
\left(\prod_j\binom{a_j}{r+1}\right)
\frac{V(B)V(A)}{\prod_{i,j}(a_j-\beta_i)}.
\end{aligned}}
\tag{12}
\]

All displayed factors other than the sign and power of zeta are positive. For zeta nonzero, every such minor is nonzero. The sign (-1)^{k(k-1)/2} comes from the Cauchy determinant when both beta and a are listed increasingly; reversing the row order removes it. A sign-free total-positivity assertion without specifying this ordering would be incorrect.

For all r+1 source rows, the factors simplify further. Put
\(E=\sum_{j=0}^r a_j-r(r+1)/2\). Then

\[
\boxed{\det R_{\{0,\ldots,r\},A}
=(-\zeta)^E\frac{V(A)}{\prod_{k=0}^r k!}.}         \tag{13}
\]

One can also prove (13) without the Cauchy formula: R=U^{-1}T_{L,O}, det U=1, and
\(T_{\beta,a}=\binom a\beta(-\zeta)^{a-\beta}\).
Factoring the zeta powers leaves the determinant of the degree-beta polynomials binom(a,beta), whose leading coefficients are 1/beta!. Their evaluation determinant is V(A)/prod beta!.

If n>=2r+1, choose the consecutive outside indices a_j=r+1+j. Their Vandermonde is prod_{k=0}^r k!, and (13) becomes

\[
\boxed{\det R_{\{0,\ldots,r\},\{r+1,\ldots,2r+1\}}
=(-\zeta)^{(r+1)^2}.}                            \tag{14}
\]

This is a pure approximation-error power with unit integer coefficient. Grouping all boundary assignments into this transport minor does not annihilate it or produce an additional small scalar factor. It may still cancel against other compatible period-minor contributions in the original determinant. Equations (12)--(14) provide an exact object for evaluating that remaining possibility.
