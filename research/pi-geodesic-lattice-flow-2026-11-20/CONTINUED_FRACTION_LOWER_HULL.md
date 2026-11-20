# Continued-fraction-only lower-envelope experiment

Date: 2026-11-20  
Status: finite numerical verification of the theorem in ACTIVE_VECTORS_AND_CONTINUED_FRACTIONS.md.

## 1. Why this is a substantial reduction

The theorem says every active denominator is a convergent denominator. Thus the lower envelope does not need all \(q=1,2,\ldots,Q\). It is enough to generate the continued-fraction convergents and test their points
\[
(q_n^2,\|q_n\pi\|^2).
\]
This changes the candidate count from \(Q\) to \(O(\log Q)\) for typical continued-fraction growth, though no universal logarithmic bound is asserted for every irrational.

## 2. Exact error formula from complete quotients

Write
\[
\pi=[a_0;a_1,a_2,\ldots]
\]
and let \(p_n/q_n\) be its convergents. If \(\alpha_{n+1}=[a_{n+1};a_{n+2},\ldots]\) is the complete quotient, then
\[
\pi=\frac{p_n\alpha_{n+1}+p_{n-1}}
{q_n\alpha_{n+1}+q_{n-1}}.
\]
Using \(p_nq_{n-1}-p_{n-1}q_n=(-1)^{n-1}\),
\[
\boxed{
|q_n\pi-p_n|=
\frac1{q_n\alpha_{n+1}+q_{n-1}}.
}
\tag{1}
\]
Therefore the point used in the hull is
\[
\boxed{
P_n=\left(q_n^2,
\frac1{(q_n\alpha_{n+1}+q_{n-1})^2}\right).
}
\tag{2}
\]
The switching abscissa between two surviving convergents \(n<m\) is still the line-intersection formula (3) of the active-vector note.

## 3. Finite pi result

Generating convergents only, then taking their lower convex hull, reproduces the active denominators found by the full scan \(q\le200000\):
\[
1,\ 7,\ 113,\ 33215,\ 66317,\ 99532.
\]
The next active convergent is \(364913\), which appears when the denominator cutoff is raised to one million.

Skipped convergents in the first range include \(106\) and \(33102\). They are genuine strict records of \(\|q\pi\|\), but their points lie above the lower convex hull. Thus “best approximation of the second kind” is necessary for activity, not sufficient.

## 4. Reproducibility and scope

The accompanying script computes the convergents at 100 decimal digits, evaluates their errors, constructs the lower hull, and compares it with a direct denominator scan for a configurable finite cutoff. Numerical agreement is a regression check; the infinite statement that active denominators are convergents comes from the theorem, not the scan.
