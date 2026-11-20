# A sharp local filter: only partial quotient 1 can cause a local hull failure

Date: 2026-11-20  
Status: **proved exact simplification** of the adjacent-curvature test, with a finite pi diagnostic.

## 1. Setup

Let \(\alpha\) be irrational, with convergents \(p_n/q_n\), errors
\[
E_n=|q_n\alpha-p_n|,
\]
and lower-hull points
\[
P_n=(q_n^2,E_n^2).
\]
Set
\[
r_n=\frac{q_{n-1}}{q_n}\in(0,1),
\qquad
z_n=\frac{E_{n+1}}{E_n}\in(0,1).
\]
The adjacent drop-slope test from the preceding note is
\[
D_{n-1,n}>D_{n,n+1}.
\tag{1}
\]

## 2. Partial quotient at least two: automatic local convexity

Write \(a=a_{n+1}\). The exact slope formulas are
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
If \(a\ge2\), then \(0<r_n,z_n<1\) gives
\[
\frac{(a+z_n)^2-1}{1-r_n^2}>a^2-1\ge3,
\qquad
\frac{1-z_n^2}{(a+r_n)^2-1}<\frac1{a^2-1}\le\frac13.
\]
Thus (1) is automatic.

## 3. Partial quotient one: a one-line threshold

If \(a_{n+1}=1\), clearing positive denominators in (1) and factoring gives
\[
(1+2r_n)z_n^2+2r_n(2+r_n)z_n+r_n^2-1>0,
\]
and the left side factors exactly as
\[
(1+2r_n)
\left(z_n-\frac{1-r_n}{1+2r_n}\right)(z_n+1+r_n).
\]
Since \(z_n+1+r_n>0\), the local convexity condition is equivalent to
\[
\boxed{
z_n>\frac{1-r_n}{1+2r_n}.
}
\tag{2}
\]

Now use the continued-fraction identities
\[
z_n=\frac{E_{n+1}}{E_n}=\frac1{\alpha_{n+2}},
\qquad
r_n=\frac{q_{n-1}}{q_n}=\frac1{\beta_n},
\]
where
\[
\alpha_{n+2}=[a_{n+2};a_{n+3},\ldots]
\]
is the forward complete quotient and
\[
\beta_n=\frac{q_n}{q_{n-1}}=[a_n;a_{n-1},\ldots,a_1]
\]
is the finite reversed quotient.

Substituting into (2) yields the especially simple exact criterion
\[
\boxed{
a_{n+1}=1
\quad\Longrightarrow\quad
\left[
P_n\text{ is locally strictly convex}
\iff
\alpha_{n+2}<\frac{\beta_n+2}{\beta_n-1}
\right].
}
\tag{3}
\]
Thus the local fate of a convergent followed by a partial quotient 1 is decided by a comparison between a forward complete quotient and a finite reversed continued fraction.

Equivalently, the only possible local failures occur at indices with \(a_{n+1}=1\), and there the exact threshold is (3). This does not settle global activity: a locally convex point may still lie above a chord joining distant convergents.

## 4. Pi diagnostic

For the initial convergents of \(\pi\), the threshold test gives:

| \(q_n\) | \(a_{n+1}\) | \(\beta_n=q_n/q_{n-1}\) | \(\alpha_{n+2}\) | Local result |
|---:|---:|---:|---:|---|
| 106 | 1 | 15.142857… | 292.634591… | fail |
| 33102 | 1 | 292.938053… | 1.736659… | fail |
| 33215 | 1 | 1.00341369… | 1.357479… | pass |
| 66317 | 1 | 1.996601… | 2.797366… | pass |
| 265381 | 1 | 2.666286… | 3.935009… | fail |

For \(q_n=265381\), the local test fails, consistent with this convergent being popped from the finite hull at cutoff \(10^6\). The first two failures are also the previously observed skipped convergents.

The decimals are finite numerical diagnostics. Formula (3) is exact for every irrational \(\alpha\).

## 5. Consequence for the streaming algorithm

The stack algorithm can use (3) as a cheap prefilter:
- if \(a_{n+1}\ge2\), the new point always passes the adjacent local test;
- if \(a_{n+1}=1\), evaluate the single threshold in (3);
- points passing the local test still enter the global stack, where distant weighted averages may remove them.

This does not make the algorithm fixed-memory, but it explains why partial quotient 1 is the only source of local curvature failures.

## 6. Scope

The result is a local convexity theorem for the Euclidean lattice-flow hull. It does not prove a uniform lower bound for \(q\|q\pi\|\), does not classify all globally active convergents of \(\pi\), and does not assert that local convexity is sufficient for global activity.
