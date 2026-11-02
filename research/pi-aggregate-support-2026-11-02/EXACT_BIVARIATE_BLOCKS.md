# Exact bivariate block ranks in a shear-invariant packet

2026-11-02. This appendix gives an exact block-rank formula under a stated time-weight condition and a counterexample outside that condition. It complements the normal-profile decomposition; it does not assert that every normal block has the smaller of its row and column counts as its rank.

Fix K>=1 and Z!=0, positive weights w_0,w,v_0,v, and cutoffs B_s,B_c>0. The source rows and raw columns of one bivariate block are

\[
\mathcal L=\{(j,s,l):0\le j<K,\ v_0s+vl<B_s\},\qquad
\mathcal C=\{(h,a):w_0h+wa\le B_c\}.
\]

Its normalized ordinary coefficient matrix is

\[
B_{j,s,l;h,a}(Z)
=\binom al[z^s]e^{hz}(jZ+z)^{a-l},
\tag{1}
\]

with the entry defined as zero if a<l. At the actual period Z=2 pi i these are the genuine local rows of Y^h T^a. At a different Z, (1) is a formal polynomial specialization after removing the center exponential, not an assertion that e^Z=1.

## 1. Exact rank when v_0>=v

For every nonnegative integer h define

\[
n_h=\max\left\{0,\left\lceil\frac{B_s-v_0h}{v}\right\rceil\right\},
\qquad
d_h=\max\left\{0,1+\left\lfloor\frac{B_c-w_0h}{w}\right\rfloor\right\}.
\]

Only finitely many n_h and d_h are positive.

**Proposition.** If v_0>=v, then for every nonzero Z,

\[
\boxed{\operatorname{rank}B(Z)
=\sum_{h\ge0}\min\{K n_h,d_h\}.}
\tag{2}
\]

Consequently this entire bivariate source packet has full row rank if and only if

\[
d_h\ge K n_h\quad\hbox{for every }h\hbox{ with }n_h>0.
\tag{3}
\]

This is an exact finite condition, including the strict source cutoff and the non-strict column cutoff.

**Proof.** At center j introduce direct local coordinates x=Y-1 and t=T-jZ. The ordinary local substitution in (1) is

\[
x=e^z-1,\qquad t=z+u.
\]

The coefficient of z^s u^l in the direct monomial x^r t^a is

\[
\binom al[z^{s-a+l}](e^z-1)^r.
\tag{4}
\]

A nonzero coefficient requires a>=l and r<=s-a+l. Therefore

\[
v_0r+va
\le v_0s+vl-(v_0-v)(a-l)
\le v_0s+vl<B_s.
\]

Thus every direct coefficient involved in an allowed old row is itself an allowed direct row. Order both finite row lists by increasing time order. In (4), r<=s, and r=s requires a=l; its coefficient then equals one. The resulting square row-change matrix is lower triangular with diagonal one. It is invertible and preserves rank. This argument proves the filtration property directly rather than assuming that a shear preserves an anisotropic jet packet.

Replace the raw columns Y^h T^a by (Y-1)^h T^a. The binomial expansion only lowers h, so it preserves the full weighted column space and is an invertible triangular column change. On the direct local rows, the new entries are

\[
\mathbf1_{r=h}\binom al(jZ)^{a-l}.
\]

The matrix is therefore a direct sum over h. Its h-block consists of n_h consecutive normalized Taylor tests at each of the K distinct points 0,Z,...,(K-1)Z, applied to all univariate polynomials of degree below d_h. A polynomial vanishing to order n_h at all K points is divisible by their product of powers, whose degree is K n_h. Conversely Hermite interpolation at these distinct points realizes arbitrary jet data with degree below K n_h. These two facts give rank min(K n_h,d_h) and prove (2). ∎

The row changes above have determinant one. The full-row Gram determinant therefore also factors into the corresponding Hermite Gram determinants **when the columns are the shifted family (Y-1)^h T^a**. The raw-column Gram matrix cannot be identified with that shifted-column Gram matrix merely because the column change is unimodular; it need not be unitary. The rank identity (2) is independent of this norm issue.

In the normal-profile setting with equal transverse weights, a profile of total normal degree k uses B_s=H-vk and B_c=H-wk. Formula (2) applies to that profile only when v_0>=v. The intended critical time anisotropy generally has v_0<v, so this formula must not be used there without a new argument.

## 2. Row and column counts alone do not give the remaining block ranks

Take a single center j=0 and

\[
K=1,\quad w_0=w=2,\quad v_0=1,\quad v=3,
\quad B_s=B_c=8.
\]

There are 15 source rows s+3l<8 and 15 raw columns h+a<=4. Nevertheless the matrix has rank exactly 14. The nonzero polynomial

\[
P(Y,T)=(T-Y+1)^4
\]

belongs to the column space, and its local substitution is

\[
P(e^z,z+u)=(u+z+1-e^z)^4.
\]

The series z+1-e^z starts in degree two. A term with u-degree l therefore starts in z-degree at least 2(4-l). Its row weight is at least 2(4-l)+3l=8+l, so every tested coefficient with s+3l<8 vanishes. This gives a nonzero kernel vector without any numerical rank estimate.

For an exact lower-rank certificate, order the rows (s,l) lexicographically in s then l, and the columns (h,a) lexicographically in h then a. Delete the row (0,1) and the column (0,0). The resulting 14-by-14 determinant is

\[
\boxed{1/640.}
\]

The companion rational verifier checks this determinant, the full kernel vector, and the rank. No formal-period issue occurs in this example, since its only center is zero.

More generally, at B_s=B_c=2n with these same weights, (T-Y+1)^n is always an invisible column combination: each tested term has weight at least 2n. This observation supplies a kernel vector, not a claim about the complete nullity for arbitrary n.

## Verification

Run `python verify_bivariate_blocks.py`. The standard-library verifier checks both exact coordinate changes, the direct Hermite blocks, their rank formula at several nonzero rational formal centers, and the explicit 15-by-15 rank counterexample. It also checks the displayed invisible-polynomial family at finite n. The universal rank formula follows from the filtration and Hermite arguments above, not from extrapolating the finite examples.
