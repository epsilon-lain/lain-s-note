# A nonzero packed transport coefficient can have deficient period rank

2026-11-01. This is a finite obstruction to inferring period rank from the local packing theorem in [the earlier local packing theorem](../pi-exterior-blast-2026-11-01/packing/BUCKET_BASE_PACKING.md). Every source row is assigned a target, all low-order rows are retained, every selected target row is individually nonzero, and every bucket Pascal minor is nonzero. Nevertheless the selected target rows have rank 8 instead of 10. It does not rule out a different successful selection, or an asymptotic common-base theorem under stronger hypotheses.

## 1. Full weighted column space and formal-period rows

Take one transverse variable, two centers, and

\[
H=5,\quad w_0=1,\quad w_1=2,\quad
v_0=1,\quad v_1=6=w_1/\theta,\quad\theta=1/3.
\]

Use the entire twelve-dimensional column space

\[
V=\operatorname{span}_{\mathbf Q}\{Y^hX^\alpha:h+2\alpha\le5\}.
\]

For a formal variable \(Z\), define factorial-scaled time rows

\[
R_{j,s,b}(Z)P
=\left.\frac{d^s}{dt^s}\right|_{t=0}
[u^b]P(e^t,jZ+t+u).
\tag{1}
\]

The ordinary coefficient-normalized jet is \(L_{j,s,b}=R_{j,s,b}/s!\). Row ranks are unchanged by this scaling. These are the formal-period matrix entries after the period identity has already been used. No assertion that \(e^Z=1\) for an arbitrary formal or numerical \(Z\) is needed. At the actual period \(Z=2\pi i\), (1) gives the appropriate actual-period rows.

The complete source packet is

\[
\mathcal S=\{(j,s,0):j=0,1;\ 0\le s\le4\},
\]

since \(s+6b<5\) forces \(b=0\).

## 2. A selection satisfying the local packing constraints

Here \(c=v_1-v_0=5\), \(t=2\), and the sufficient bulk threshold in the packing note is

\[
N_0=\frac{t(H+c)}{c+tv_0}=\frac{20}{7}.
\]

Every source bucket has \(B_N=\{0\}\), for \(0\le N\le4\). Keep identity transport for \(N=0,1,2\). At each of \(N=3,4\), take target type \(a=1\) at center 0 and \(a=2\) at center 1. Both belong to

\[
A_N=\{a\ge0:5a<2(5-N+5),\ a\le N\}.
\]

They are disjoint across centers, satisfy \(2a\le H\), and each bucket Pascal minor is the nonzero singleton \(\binom a0=1\). In the packing note's coefficient normalization, the complete transport coefficient is

\[
(-\zeta)^{(1+2)+(1+2)}=\zeta^6\ne0
\quad(\zeta\ne0).
\]

The full selected target packet, in the order used below, is

| Index | Center \(j\) | Time order \(s\) | Type \(b\) | Conserved order \(N=s+b\) |
|---|---:|---:|---:|---:|
| 0 | 0 | 0 | 0 | 0 |
| 1 | 1 | 0 | 0 | 0 |
| 2 | 0 | 1 | 0 | 1 |
| 3 | 1 | 1 | 0 | 1 |
| 4 | 0 | 2 | 0 | 2 |
| 5 | 1 | 2 | 0 | 2 |
| 6 | 0 | 2 | 1 | 3 |
| 7 | 1 | 1 | 2 | 3 |
| 8 | 0 | 3 | 1 | 4 |
| 9 | 1 | 2 | 2 | 4 |

## 3. Exact source rank and exact target deficiency

Order all raw columns by

\[
(h,\alpha)=(0,0),(0,1),(0,2),(1,0),(1,1),(1,2),
(2,0),(2,1),(3,0),(3,1),(4,0),(5,0).
\]

On the first ten columns, the source rows ordered first by \(j\) and then by \(s\) have determinant

\[
\boxed{\det R_{\mathcal S,\,0:9}(Z)=-9216Z^5.}
\tag{2}
\]

Thus the complete source packet has rank 10 for every \(Z\ne0\), including the actual period. This conclusion is verified directly; it does not invoke an interpolation theorem.

Write \(R_k\) for target row \(k\) in the table. On the entire twelve-dimensional space \(V\), there are two exact polynomial row relations:

\[
\boxed{R_4-R_5+ZR_6+Z^2R_7=0,\qquad R_7-R_9=0.}
\tag{3}
\]

For the second relation, the only allowed monomials with transverse degree 2 are \(X^2\) and \(YX^2\). Extracting \([u^2]\) leaves \(1\) and \(e^t\). Their first and second time derivatives at zero agree. Consequently the two factorial-scaled rows are equal. In the ordinary normalization the statement is \(L_{1,1,2}=2L_{1,2,2}\), so those rows are proportional, not equal.

For the first relation, regard \(D_t^2P(e^t,t+x)|_{t=0}\) as a polynomial in \(x\) of degree at most 2. Its exact Taylor identity between \(x=0\) and \(x=Z\) gives

\[
R_{1,2,0}=R_{0,2,0}+ZR_{0,2,1}+Z^2R_{0,2,2}.
\]

The type-2 rows are independent of the center, and the preceding observation identifies \(R_{0,2,2}=R_{1,1,2}\). This proves (3) without calculation of a determinant.

All ten target rows are individually nonzero for every Z: for a row of type b, the permitted raw monomial Y X^b has row value one. Moreover, target rows \(0,1,2,3,4,5,6,8\), restricted to the first eight columns, have determinant

\[
\boxed{16Z^4.}
\tag{4}
\]

Hence the target rank is exactly 8 for every \(Z\ne0\), including \(2\pi i\). No alternative choice of columns or common column mixing can raise this fixed selection's rank.

## 4. Scope and remaining lemma

This example satisfies the local packing theorem's lower-set, dilation, conserved-order, disjointness, Pascal-minor, and transverse column-support conditions. It even satisfies the elementary scalar inequalities \(K\theta=2/3<1\) and \(Kw_0\theta/v_0=2/3<1\). It is not asserted to satisfy the manuscript's successive large-weight hypotheses. Therefore it is not a counterexample to that manuscript's interpolation theorem or to an existence theorem for better target selections.

It proves that these local packing conditions alone do not imply global period independence, even if all low-order rows are retained and zero target rows are excluded. The obstruction depends on the time directions available in the shared weighted column space: a type can be supported while having too few independent time profiles.

The outstanding asymptotic claim must select a low-cost common base of the bucketwise Pascal matroid and the actual-period row matroid on the same full \(P_W(H)\), under the intended parameter hypotheses. A bound on the number of unseparated low-order rows is insufficient without a rank-preserving repair argument. The finite successful selection in `BULK_PERIOD_CERTIFICATE.md` is compatible with this obstruction: some selections work, while others do not.

## 5. Reproducible exact verification

Run `python verify_period_obstruction.py`.

The standard-library verifier checks all local bounds, every bucket coefficient, nonzero individual target rows, and the coefficientwise identities (3). It proves (2) and (4) by exact rational determinant evaluations at more points than the respective polynomial degree bounds. These checks establish polynomial identities, not numerical approximations or samples of an unbounded claim.
