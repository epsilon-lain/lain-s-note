# Active transverse spans, normal profiles, and necessary weight separation

2026-11-02. This note extends the earlier single diagonal-quotient capacity test. With equal transverse weights it decomposes the entire source matrix into normal-profile blocks and sums their capacity deficits. It also gives a two-dimensional quotient test for the original displaced centers with unequal weights. These are rank obstructions and necessary conditions, not sufficient interpolation thresholds or a conclusion about pi.

## 1. The matrix and its active transverse subspace

Use ordinary coefficient rows of

\[
A_{j,s,\beta}(P)=[z^s u^\beta]P(e^z,c_j+z\mathbf1+u),
\quad 0\le j<K,
\]

where the vector of ones is denoted by \(\mathbf1\). The column space and complete source packet are

\[
V_H=\operatorname{span}\{Y^hX^\alpha:w_0h+w\cdot\alpha\le H\},
\qquad
L_H=\{(j,s,\beta):v_0s+w\cdot\beta/\theta<H\}.
\tag{1}
\]

All weights are positive and fixed before H tends to infinity; \(0<\theta\le1\). The coefficients are ordinary normalized Taylor coefficients. Complete logarithmic time packets can be changed to these packets by determinant-one row operations for each fixed transverse type.

The common-period matrix has \(c_j=jZ\mathbf1\). Its entries are formal polynomials in Z after the center exponential has been removed. They are the actual period rows at \(Z=2\pi i\); arbitrary rational Z used in a verifier is not asserted to be an exponential period.

The displaced matrix has \(c_j=jr\), with \(r=(r_1,\ldots,r_m)\). Define the active subspace

\[
E=\operatorname{span}\{\mathbf1,c_0,\ldots,c_{K-1}\}
\subseteq\mathbb C^m,\qquad d=\dim E.
\tag{2}
\]

At a common period, d=1. For the displaced centers, d<=2; if K>=2 and the r_i are not all equal, d=2. When K=1 the only center is zero and d=1 regardless of r.

## 2. Exact full-profile decomposition when all transverse weights agree

Assume \(w_1=\cdots=w_m=w\), put \(v=w/\theta\), and let n=m-d. Choose linear transverse coordinates

\[
(T_1,\ldots,T_d,U_1,\ldots,U_n)
\]

such that the U coordinates vanish on E. Both the column budget and each fixed-time source derivative space are invariant under this coordinate change: the change preserves transverse total degree, and the packet contains every transverse derivative of each permitted total degree. This statement has no condition on the relative sizes of v_0 and v.

Write a column as \(Y^hT^aU^\mu\), and a row as \((j,s,\ell,\nu)\). At every center and along the entire local time direction the normal coordinates vanish. Consequently

\[
\boxed{
A^{\rm ad}_{j,s,\ell,\nu;h,a,\mu}
 =\mathbf1_{\nu=\mu}
[z^s t^\ell]e^{hz}\bigl(T(c_j)+zT(\mathbf1)+t\bigr)^a.}
\tag{3}
\]

The simultaneous row and column changes are invertible. Thus the full matrix is a direct sum of active blocks \(B_\nu\), plus columns belonging to profiles with no source rows. A profile of normal total degree k has

\[
r_k=K\#\{(s,\ell)\in\mathbb Z_{\ge0}\times\mathbb Z_{\ge0}^d:
v_0s+v|\ell|<H-vk\},
\]

\[
c_k=\#\{(h,a)\in\mathbb Z_{\ge0}\times\mathbb Z_{\ge0}^d:
w_0h+w|a|\le H-wk\}.
\tag{4}
\]

For n>=1 its multiplicity is \(g_n(k)=\binom{k+n-1}{n-1}\). The exact rank identity and resulting capacity bound are

\[
\boxed{\operatorname{rank}A_H=\sum_\nu\operatorname{rank}B_\nu
\le\sum_{k:\,vk<H}g_n(k)\min(r_k,c_k).}
\tag{5}
\]

If n=0 there is only the zero profile and the same statement is interpreted accordingly. The minimum in (5) is an upper bound, not a formula for each block's rank. The appendix gives an explicit 15-by-15 active block of rank 14, so this distinction is necessary.

For the common-period case one can choose \(T=X_1\), \(U_i=X_{i+1}-X_1\). Then the coordinate changes on the finite column and coefficient-row spaces are integer and determinant one. The active block is the genuine bivariate matrix

\[
B_{j,s,\ell;h,a}(Z)=\binom a\ell[z^s]e^{hz}(jZ+z)^{a-\ell},
\tag{6}
\]

with the source and column budgets reduced as in (4). This extends the earlier single-type capacities and zero-normal-profile quotient to every normal profile at once.

### Determinants and the choice of column metric

In the adapted column basis, the full Gram matrix is exactly block diagonal, so

\[
\det(A^{\rm ad}(A^{\rm ad})^*)=\prod_\nu\det(B_\nu B_\nu^*).
\tag{7}
\]

A square minor using adapted columns is zero unless it takes exactly r_k columns from each individual source profile; if those counts match, its determinant factors into the corresponding block minors, up to ordering signs.

These statements do not identify the adapted Gram with the raw-column Gram. If \(A^{\rm ad}=S A C\), then exactly

\[
AA^*=S^{-1}A^{\rm ad}Q(A^{\rm ad})^*S^{-*},
\qquad Q=C^{-1}C^{-*}.
\tag{8}
\]

The positive metric Q generally couples the profiles. A determinant-one column change is not a unitary change. Equations (3)--(5) are valid rank statements for the raw matrix regardless of this metric distinction.

## 3. A positive fraction of missing rank, then an almost total loss

Suppose n=m-d>=1. Put

\[
\lambda=\lim_{H\to\infty}\frac{M_H}{N_H}
=\frac{K\theta^m w_0}{v_0},\qquad
A_d=\frac{K\theta^d w_0}{v_0}.
\]

The fixed-parameter leading volumes in a profile with y=vk/H are proportional to \((1-y)^{d+1}\) for rows and \(A_d^{-1}(1-\theta y)^{d+1}\) for columns, after dividing both by the row coefficient. Summing (5), including the profile multiplicity, gives

\[
\boxed{
\limsup_{H\to\infty}\frac{\operatorname{rank}A_H}{M_H}
\le \mathcal C_{m,d}:=
\frac{\displaystyle\int_0^1y^{n-1}
\min\{(1-y)^{d+1},A_d^{-1}(1-\theta y)^{d+1}\}\,dy}
{\displaystyle B(n,d+2)}.}
\tag{9}
\]

Here B is the beta integral. Fixed rational weights give the required lattice estimates; counts on a boundary layer are lower order. In particular, if A_d>1 and theta<1, every sufficiently small y is overloaded, and the integral in (9) is strictly less than one. The crossover is

\[
y_* =\frac{A_d^{1/(d+1)}-1}{A_d^{1/(d+1)}-\theta}.
\]

Thus the earlier single quotient obstruction can represent a positive fraction of all source rows, rather than only an O(H^(d+1)) exceptional family.

A simpler bound is especially informative in high dimension. Every column with normal total degree k>=H/v is invisible to all source rows. The fraction of all columns with \(wk/H<\theta\) tends to

\[
I_\theta(n,d+2)=
\frac{\int_0^\theta x^{n-1}(1-x)^{d+1}\,dx}{B(n,d+2)}.
\]

It follows directly that

\[
\boxed{\limsup_{H\to\infty}\frac{\operatorname{rank}A_H}{M_H}
\le\min\{1,\lambda^{-1}I_\theta(m-d,d+2)\}.}
\tag{10}
\]

This bound counts all visible columns, even when a visible profile is itself underfilled or has further analytic dependencies, so it is generally weaker than (9).

On the critical count path \(v_0=2K\theta^m w_0\), one has lambda=1/2. For the common-period matrix (d=1), (10) becomes

\[
\boxed{\limsup_H\frac{\operatorname{rank}B_L}{M_H}
\le \min\{1,\theta^{m-1}[
(1-\theta)^2m^2+(1-\theta^2)m+2\theta]\}.}
\tag{11}
\]

For the displaced centers with d=2 and m>=3 it becomes

\[
\boxed{\limsup_H\frac{\operatorname{rank}A^r_L}{M_H}
\le\min\{1,2\theta^{m-2}[
\binom{m+1}{3}(1-\theta)^3+
\binom{m+1}{2}\theta(1-\theta)^2+
(m+1)\theta^2(1-\theta)+\theta^3]\}.}
\tag{12}
\]

For every fixed theta<1, the bounds in (11) and (12) tend to zero as m tends to infinity, after the separate H limit at each m. This does not require a bound on K across that subsequent dimension comparison, since K has canceled on the critical count path.

For theta=1/2, the untruncated bounds are respectively

\[
\frac{m^2+3m+4}{2^{m+1}},\qquad
\frac{m^3+3m^2+8m+12}{6\,2^m}.
\]

At m=10 these are 67/1024 and 29/128. Hence equal transverse weights in this regime cannot provide a full-rank source matrix even after the centers leave the common diagonal. The obstruction at arbitrary displaced centers is consequential: one cannot take the equal-weight model's rank loss at the period and presume that a nonzero arithmetic source minor is available nearby.

## 4. Unequal weights: a quotient test with two greedy bases

The simultaneous full-profile decomposition above need not preserve unequal weighted filtrations. The following quotient test does not assume such preservation.

For the active space E in (2), let \(\xi_i=X_i|_E\in E^*\). Process the coordinates in increasing order of w_i, retaining a restriction only when it increases the span. Let the d selected weights be \(a_1,\ldots,a_d\). Process them in decreasing order instead, giving weights \(b_1,\ldots,b_d\). Ties can be broken consistently in any order. Both selections are bases of E*.

**Column image.** In the cheapest basis, every restriction xi_i is a linear combination of retained restrictions with weights at most w_i. Substitution therefore preserves the column budget. Conversely the chosen basis restrictions come from the corresponding actual coordinate monomials. Thus restriction to E has image exactly

\[
\operatorname{span}\{Y^h T^\alpha:w_0h+\sum_{q=1}^d a_q\alpha_q\le H\},
\]

where alpha is the exponent vector; equivalently its dimension is

\[
D_H=\#\{(h,\alpha)\in\mathbb Z_{\ge0}^{d+1}:
w_0h+\sum_q a_q\alpha_q\le H\}.
\tag{13}
\]

**Source directions.** Let e_q in E be dual to the most-expensive basis, retaining its descending order. Every coordinate of weight greater than b_q lies in the span of restrictions selected before q, so it annihilates e_q. Thus e_q is supported only on coordinates whose weights are at most b_q. A product of normalized directional derivatives of multi-order ell therefore uses only source derivatives beta with

\[
w\cdot\beta\le\sum_q b_q\ell_q.
\]

For each center, all tests with \(v_0s+\sum_q b_q\ell_q/\theta<H\) belong to the original source row span. Their coefficient vectors in the original row-label space are independent: the direction vectors e_q are independent, their symmetric monomials give independent differential symbols, and distinct centers or time orders use disjoint row-label groups.

All these tests factor through restriction to E, since both the center and time direction lie in E. Define

\[
L_H^E=\#\{(s,\ell)\in\mathbb Z_{\ge0}^{d+1}:
v_0s+\sum_q b_q\ell_q/\theta<H\}.
\]

We obtain the finite and asymptotic necessary conditions

\[
\boxed{K L_H^E>D_H\quad\Longrightarrow\quad
\text{the full source matrix is row-rank deficient},}
\tag{14}
\]

\[
\boxed{\text{full source rank for unbounded H requires}\qquad
K\theta^d w_0\prod_q a_q\le v_0\prod_q b_q.}
\tag{15}
\]

The reverse weak inequality in (15) is only necessary. Neither it nor the finite count inequality in the opposite direction proves independence.

### Common-period and displaced-center specializations

At the common period, every coordinate restriction to E=span(1) is the same one-dimensional functional. Thus a_1=w_min and b_1=w_max. Formula (15) recovers exactly the earlier diagonal quotient condition; this special case is not claimed as new.

For m>=2, K>=2 and pairwise-distinct r_i, the active plane is span(1,r), and coordinate i restricts to the linear form \(x+r_i y\). Any two such forms are independent. Ordering weights increasingly therefore gives

\[
(a_1,a_2)=(w_1,w_2),\qquad
\{b_1,b_2\}=\{w_{m-1},w_m\}.
\]

Consequently full source rank at the displaced centers requires

\[
\boxed{\frac{w_{m-1}w_m}{w_1w_2}
\ge \frac{K\theta^2w_0}{v_0}.}
\tag{16}
\]

On the critical count path this becomes

\[
\boxed{\frac{w_{m-1}w_m}{w_1w_2}
\ge\frac1{2\theta^{m-2}}.}
\tag{17}
\]

When some r_i coincide, the greedy bases must be used rather than automatically taking the two smallest or largest weights. For example, the second cheapest pivot is the cheapest coordinate whose r_i differs from the first pivot's value. These tests constrain the original displaced source matrix, not only a period specialization. Strongly separated weights can satisfy them; no contradiction to an interpolation theorem with such hypotheses is asserted.

## 5. Scope and verification

The new equal-weight conclusions are the exact all-profile rank decomposition and the aggregate rank-fraction bounds, including the displaced active plane. The new unequal-weight conclusion is the active-plane quotient test with its two greedy bases. The earlier diagonal quotient is recovered as a special case. The appendix identifies a regime with exact active-block ranks and separately shows why minimum row/column counts cannot be substituted for ranks in general.

The companion verifier checks coordinate changes and exact profile identities, finite rank sums and adapted Gram products, finite profile counts, the rational beta-integral expressions, and the greedy restriction/directional-support construction. These finite computations accompany the general proofs; neither asymptotic full rank nor any statement about pi is inferred from numerical ranks.
