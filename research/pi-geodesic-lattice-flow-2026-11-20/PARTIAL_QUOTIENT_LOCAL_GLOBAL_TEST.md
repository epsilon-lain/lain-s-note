# Partial-quotient simplification of the lattice-flow convexity test

Date: 2026-11-20  
Status: **proved local criterion and exact global weighted-average criterion**. The result explains precisely why partial quotient 1 is the only possible source of local deletion.

## 1. Setup and indexing

Let \(\alpha\) be irrational, with regular continued fraction
\[
\alpha=[a_0;a_1,a_2,\ldots],
\]
convergents \(p_n/q_n\), and errors
\[
E_n=|q_n\alpha-p_n|.
\]
For \(n\ge1\), define
\[
r_n=\frac{q_{n-1}}{q_n}\in(0,1),
\qquad
z_n=\frac{E_{n+1}}{E_n}\in(0,1),
\qquad a=a_{n+1}.
\]
The recurrences are
\[
q_{n+1}=a q_n+q_{n-1},\qquad
E_{n-1}=a E_n+E_{n+1}.
\]
The lattice-flow line associated with \(q_n\) is
\[
\ell_n(x)=q_n^2+E_n^2x,\qquad x=e^{4t}\ge1.
\]

## 2. Adjacent slopes and exact simplification

Define adjacent positive drop ratios
\[
D_{n-1,n}=\frac{E_{n-1}^2-E_n^2}{q_n^2-q_{n-1}^2},
\qquad
D_{n,n+1}=\frac{E_n^2-E_{n+1}^2}{q_{n+1}^2-q_n^2}.
\]
Substitution of the recurrences gives
\[
D_{n-1,n}
=\left(\frac{E_n}{q_n}\right)^2
\frac{(a+z_n)^2-1}{1-r_n^2},
\]
\[
D_{n,n+1}
=\left(\frac{E_n}{q_n}\right)^2
\frac{1-z_n^2}{(a+r_n)^2-1}.
\]
A necessary local condition for \(P_n=(q_n^2,E_n^2)\) to be a strict lower-hull vertex is
\[
D_{n-1,n}>D_{n,n+1}.
\tag{1}
\]

### Case A: \(a_{n+1}\ge2\)

Since \(0<r_n,z_n<1\),
\[
\frac{(a+z_n)^2-1}{1-r_n^2}>a^2-1\ge3,
\]
whereas
\[
\frac{1-z_n^2}{(a+r_n)^2-1}
<\frac1{a^2-1}\le\frac13.
\]
Therefore (1) is automatic.

**Conclusion:** a convergent followed by a partial quotient at least 2 is always locally convex. It can still fail the global support test because a distant chord may lie below it.

### Case B: \(a_{n+1}=1\)

Condition (1) becomes
\[
\frac{(1+z_n)^2-1}{1-r_n^2}>
\frac{1-z_n^2}{(1+r_n)^2-1}.
\]
All denominators are positive. Clearing them and simplifying yields the exact equivalent condition
\[
\boxed{
(r_n+z_n)^2+2r_nz_n(1+r_n+z_n)>1.
}
\tag{2}
\]
Thus the only convergents that can fail even the local convexity test are those followed by a partial quotient 1. For these, (2) is necessary and sufficient for local convexity.

A useful immediate sufficient condition is \(r_n+z_n\ge1\), which forces (2). If \(r_n+z_n<1\), the additional positive term \(2r_nz_n(1+r_n+z_n)\) decides the outcome.

## 3. The global criterion as comparison of weighted averages

For \(i<n\), define
\[
D_{i,n}=\frac{E_i^2-E_n^2}{q_n^2-q_i^2}.
\]
For \(j>n\), define
\[
D_{n,j}=\frac{E_n^2-E_j^2}{q_j^2-q_n^2}.
\]
Let
\[
D_k=D_{k,k+1}.
\]
Then telescoping gives
\[
D_{i,n}
=\frac{\sum_{k=i}^{n-1}(q_{k+1}^2-q_k^2)D_k}
{\sum_{k=i}^{n-1}(q_{k+1}^2-q_k^2)},
\tag{3}
\]
and
\[
D_{n,j}
=\frac{\sum_{k=n}^{j-1}(q_{k+1}^2-q_k^2)D_k}
{\sum_{k=n}^{j-1}(q_{k+1}^2-q_k^2)}.
\tag{4}
\]
So every nonadjacent secant drop is a weighted average of adjacent drops, with positive weights
\[
w_k=q_{k+1}^2-q_k^2.
\]

The exact global support interval from the preceding note can therefore be rewritten as follows. Put
\[
d^-_n=\inf_{i<n}D_{i,n},
\qquad
d^+_n=\sup_{j>n}D_{n,j}.
\]
Then
\[
L_n=\max(1,1/d^-_n),\qquad U_n=1/d^+_n.
\]
Hence
\[
\boxed{
q_n\text{ is strictly active for some }t\ge0
\iff
d^+_n<\min(1,d^-_n).
}
\tag{5}
\]
Equations (3)–(5) show why checking only adjacent curvature is insufficient: a long weighted average from the right can exceed a short adjacent drop, or a long average from the left can be unusually small.

## 4. Pi diagnostic

For the initial convergents of \(\pi\), the local criterion gives:

| Denominator \(q_n\) | Next partial quotient \(a_{n+1}\) | Local test |
|---:|---:|---|
| 7 | 15 | automatic |
| 106 | 1 | fails |
| 113 | 292 | automatic |
| 33102 | 1 | fails |
| 33215 | 1 | passes |
| 66317 | 1 | passes |

For \(q=106\),
\[
(r_n+z_n)^2+2r_nz_n(1+r_n+z_n)
\approx0.0053066721<1.
\]
For \(q=33102\), the same expression is approximately \(0.3417179371<1\). For \(q=33215\) and \(q=66317\), it is approximately \(7.0174422\) and \(1.4021788\), respectively.

These are finite numerical diagnostics. The universal implication for \(a_{n+1}\ge2\) and the equivalence (2) are exact proofs.

## 5. What has and has not been reduced

The continuous lower-envelope problem now has three exact layers:

1. An active denominator must be a convergent denominator.
2. A convergent followed by \(a_{n+1}\ge2\) automatically passes local convexity; one followed by \(a_{n+1}=1\) is decided locally by (2).
3. Actual global activity is decided by the weighted-average comparison (5), involving all earlier and later convergents.

This gives a fast local filter, but not a finite-memory rule: (5) still involves arbitrarily long weighted averages. The \(\pi\) bad-approximability question remains equivalent to a uniform lower bound for the shortest-vector orbit, and is not solved here.
