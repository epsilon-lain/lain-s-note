# The active vectors of the Euclidean lattice flow are convergent denominators

Date: 2026-11-20  
Status: **proved theorem** for every irrational target \(\alpha\), plus a finite pi experiment. This identifies the arithmetic source of the lower-envelope itinerary.

## 1. Setup

For irrational \(\alpha\), let
\[
E_q=\|q\alpha\|=\min_{p\in\mathbb Z}|q\alpha-p|,
\qquad q\ge1.
\]
The shortest candidate for denominator \(q\) is the nearest-integer numerator. After setting \(x=e^{4t}\), minimizing the squared lattice-vector length is equivalent to minimizing
\[
\ell_q(x)=q^2+E_q^2x,\qquad x\ge1.
\]
The vertical vector contributes \(\ell_0(x)=x\).

Call \(q\) active if \(\ell_q(x)\) is the strict minimum among all candidates for at least one \(x\ge1\).

## 2. Active implies a strict record of the second-kind error

**Proposition.** If \(q\) is active, then
\[
E_q<E_r\qquad\text{for every }1\le r<q.
\tag{1}
\]

**Proof.** If some \(r<q\) satisfied \(E_r\le E_q\), then for every \(x\ge1\),
\[
\ell_r(x)=r^2+E_r^2x
<q^2+E_q^2x=\ell_q(x),
\]
contradicting activity of \(q\). \(\square\)

Thus every active denominator is a strict record denominator for \(\|q\alpha\|\), i.e. a best approximation of the second kind.

## 3. Strict record denominators are convergent denominators

Let
\[
\alpha=[a_0;a_1,a_2,\ldots]
\]
be the regular continued fraction, with convergents \(p_n/q_n\). The classical best-approximation theorem of the second kind says:
\[
0<q<q_{n+1}
\quad\Longrightarrow\quad
\|q\alpha\|\ge\|q_n\alpha\|,
\]
with equality cases handled by the usual initial-index conventions. In particular, if \(q\) is a strict record for \(\|q\alpha\|\), then \(q=q_n\) for some convergent denominator.

Combining this with Proposition 1 gives:

\[
\boxed{\text{Every active shortest-vector denominator is a continued-fraction convergent denominator.}}
\tag{2}
\]

This is a theorem, not a numerical observation. The converse is false: a convergent can be skipped by the Euclidean lower envelope because its point \((q_n^2,E_{q_n}^2)\) lies above the lower convex hull.

Intermediate convergents are therefore not needed as active vectors for this particular Euclidean norm. They may still be relevant for other norms or other deformation paths.

## 4. Exact switching formula between active convergents

Suppose active denominators \(q_n<q_m\) are consecutive vertices of the lower envelope, with errors \(E_n,E_m\). Their switching time is
\[
\boxed{
x_{n,m}=\frac{q_m^2-q_n^2}{E_n^2-E_m^2},
\qquad
t_{n,m}=\frac14\log x_{n,m}.
}
\tag{3}
\]
Since an active transition has \(q_m>q_n\) and \(E_m<E_n\), the denominator and numerator are positive.

The vector \(q_n\) is active on the interval between its entry and exit intersections. This interval is nonempty precisely when its two adjacent intersection abscissae are ordered. Equivalently, the points \((q_n^2,E_n^2)\) are successive vertices of the lower convex hull with the appropriate restriction \(x\ge1\).

The individual balance time
\[
t_n^*=\frac12\log(q_n/E_n)
\]
is not the switching time (3); it minimizes one line's corresponding two-coordinate norm, whereas (3) equates two different vectors' lengths.

## 5. Pi experiment, interpreted arithmetically

For \(\alpha=\pi\), the finite search \(q\le200000\) produced the active sequence
\[
1,\ 7,\ 113,\ 33215,\ 66317,\ 99532.
\]
These are all regular continued-fraction convergent denominators. The convergent denominators \(106\) and \(33102\) are skipped by the Euclidean lower envelope. The next convergent denominator \(364913\) appears when the cutoff is extended to one million.

The exact theorem (2) explains why no semiconvergent denominator can appear in this active itinerary. It does not determine which convergents survive the convex-hull test; that depends on the slopes between consecutive points \((q_n^2,E_n^2)\).

## 6. A useful continued-fraction criterion for a skipped convergent

For a convergent denominator \(q_n\), define
\[
P_n=(q_n^2,E_n^2).
\]
It is active if and only if there exists \(x\ge1\) such that
\[
q_n^2+xE_n^2
<
q_k^2+xE_k^2
\quad\text{for every }k\ne n.
\tag{4}
\]
Only convergents need be checked in (4), by (2). Thus the infinite rational search has been reduced to the continued-fraction sequence, followed by a one-dimensional convex-hull test.

A convergent may be discarded if it is dominated by a convex combination of two other convergent points in the \((q^2,E^2)\)-plane. This is the precise source of the low-redundancy effect.

## 7. What remains open for pi

The active-denominator theorem is universal and classical in its continued-fraction input. It does not give a lower bound on the minimum lattice length for \(\pi\). The desired statement remains
\[
\inf_{t\ge0}\lambda_1(\Lambda_\pi(t))>0,
\]
equivalently \(\inf_{q\ge1}q\|q\pi\|>0\).

The next useful task is to express the convex-hull survival condition (4) directly in terms of the partial quotients \(a_n\), and to determine whether the survival pattern for \(\pi\) carries any information beyond the already-known continued-fraction data. No new lower bound for the partial quotients of \(\pi\) is claimed here.

## 8. Attribution

The best-approximation property of continued-fraction convergents is classical; see Cassels, *An Introduction to Diophantine Approximation*, Chapter I. The contribution here is the short deduction that the active vectors of this particular diagonal lattice flow must be among those convergents, together with the explicit lower-envelope criterion (4).
