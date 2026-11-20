# Global active-convergent criterion and partial-quotient curvature

Date: 2026-11-20  
Status: **proved exact criterion**. It separates a local necessary test from the genuinely global convex-hull condition.

## 1. Geometry of the lower envelope

Let \(\alpha\) be irrational, and let \(p_n/q_n\) be its regular continued-fraction convergents. Put
\[
E_n=|q_n\alpha-p_n|=\|q_n\alpha\|,
\qquad X_n=q_n^2,\quad Y_n=E_n^2.
\]
The lattice-flow candidates are the lines
\[
\ell_n(x)=X_n+xY_n,\qquad x=e^{4t}\ge1.
\]
The active convergent is one for which \(\ell_n(x)\) is strictly smaller than every other candidate for at least one \(x\ge1\).

By the active-vector theorem, no nonconvergent denominator needs to be considered.

## 2. Exact global criterion

For indices \(i<n<j\), define the positive secant-drop ratios
\[
D_{i,n}=\frac{Y_i-Y_n}{X_n-X_i},
\qquad
D_{n,j}=\frac{Y_n-Y_j}{X_j-X_n}.
\]
They are positive because convergent denominators increase and their errors strictly decrease.

The inequality \(\ell_n(x)<\ell_i(x)\), for \(i<n\), is equivalent to
\[
x>\frac{X_n-X_i}{Y_i-Y_n}=\frac1{D_{i,n}}.
\]
For \(j>n\), \(\ell_n(x)<\ell_j(x)\) is equivalent to
\[
x<\frac{X_j-X_n}{Y_n-Y_j}=\frac1{D_{n,j}}.
\]
Consequently define
\[
L_n=\max\left(1,\ \sup_{i<n}\frac1{D_{i,n}}\right),
\qquad
U_n=\inf_{j>n}\frac1{D_{n,j}}.
\]
(For the first/last index of a finite truncation, the missing side is interpreted as \(0\) or \(+\infty\), respectively.)

**Theorem.** The convergent denominator \(q_n\) is strictly active for some \(t\ge0\) if and only if
\[
\boxed{L_n<U_n.}
\tag{1}
\]
If (1) holds, its exact active interval in the \(x\)-parameter is
\[
\boxed{L_n<x<U_n,}
\tag{2}
\]
and the corresponding time interval is
\[
\boxed{\frac14\log L_n<t<\frac14\log U_n.}
\tag{3}
\]
If \(L_n=U_n\), the line can only touch the lower envelope at a tie and is not strictly active.

This is the global condition. It checks all earlier and later convergents, not merely the two adjacent ones.

## 3. Partial-quotient formula for the adjacent curvature test

Write the convergent recurrence as
\[
q_{n+1}=a_{n+1}q_n+q_{n-1},
\qquad
E_{n-1}=a_{n+1}E_n+E_{n+1}.
\]
Set
\[
r_n=\frac{q_{n-1}}{q_n}\in(0,1),
\qquad
z_n=\frac{E_{n+1}}{E_n}\in(0,1),
\qquad a=a_{n+1}.
\]
Then
\[
q_{n+1}/q_n=a+r_n,\qquad E_{n-1}/E_n=a+z_n.
\]
The adjacent drop ratios are therefore
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
Hence a necessary local convexity condition for \(P_n=(X_n,Y_n)\) to be a strict lower-hull vertex is
\[
\boxed{
\frac{(a+z_n)^2-1}{1-r_n^2}
>
\frac{1-z_n^2}{(a+r_n)^2-1}.
}
\tag{4}
\]
The common factor \((E_n/q_n)^2\) cancels.

The two variables are themselves continued-fraction quantities:
\[
r_n=[0;a_n,a_{n-1},\ldots,a_1],
\qquad
z_n=\frac1{\alpha_{n+2}},
\]
where \(\alpha_{n+2}=[a_{n+2};a_{n+3},\ldots]\) is the complete quotient.

Condition (4) is only a local necessary test. It is not sufficient for global activity: a point can be locally convex relative to its immediate neighbors and still lie above a chord joining two more distant points. The exact sufficient-and-necessary criterion remains (1).

## 4. Why the global test is naturally an envelope test

For any \(i<n<j\), the point \(P_n\) lies strictly below the chord joining \(P_i,P_j\) precisely when
\[
D_{i,n}>D_{n,j}.
\]
Thus \(P_n\) is a strict lower-hull vertex exactly when it admits a supporting line of slope \(-1/x\) for some \(x\ge1\). Equations (1)–(3) record the same condition without ambiguity about the restricted range \(x\ge1\).

This formulation also handles skipped convergents correctly. It avoids the incorrect shortcut “the adjacent slopes turn the right way, therefore the point is active.”

## 5. Pi diagnostic

For the initial convergents of \(\pi\), the local test (4) rejects the skipped denominators \(106\) and \(33102\):
- at \(q=106\), the left adjacent drop ratio is smaller than the right one;
- at \(q=33102\), the same failure occurs.

The local test accepts the first active denominators \(1,7,113,33215,66317,99532\) in the computed range. Acceptance alone is not a global proof; the full envelope criterion (1) is still required.

The finite experiment supplies the numerical illustration. No infinite statement about the active sequence of \(\pi\) follows from this diagnostic.

## 6. What this accomplishes

The original one-parameter lattice flow has now been reduced in two exact stages:

1. all integer denominators \(\to\) continued-fraction convergent denominators (best approximation of the second kind);
2. all convergents \(\to\) those satisfying the global support interval (1).

The first reduction is classical. The second is an exact convex-geometric criterion, now expressed in terms of the continued-fraction errors and, locally, the partial quotients.

The remaining arithmetic difficulty is unchanged:
\[
\inf_{q\ge1}q\|q\pi\|>0
\]
is still unknown. The criterion identifies which convergents control the continuous lattice orbit, but does not bound their values \(q_nE_n\) away from zero.

## 7. Attribution

The continued-fraction recurrences and best-approximation theorem are classical. The global interval (1) is the direct supporting-line criterion for the particular lattice-flow lower envelope; formula (4) is its adjacent partial-quotient expansion. No claim is made that the classical bad-approximation problem for \(\pi\) is resolved.
