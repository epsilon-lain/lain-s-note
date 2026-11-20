# Finite but nonuniform lookahead: exact future-veto stopping theorem

Date: 2026-11-20  
Status: **proved theorem**. For each fixed convergent, its infinite future-veto infimum is attained after finitely many future convergents, with an explicit stopping certificate. The required lookahead is data-dependent and has no uniform bound established here.

## 1. Two-sided activity thresholds

Let \(\alpha\) be irrational, with convergents \(p_n/q_n\), errors \(E_n=|q_n\alpha-p_n|\), and points
\[
P_n=(X_n,Y_n)=(q_n^2,E_n^2).
\]
The candidate line is \(\ell_n(x)=X_n+xY_n\), \(x=e^{4t}\ge1\).

For a fixed index \(n\), define the past-admission threshold
\[
P_n^*=\max\left(1,\max_{i<n}
\frac{q_n^2-q_i^2}{E_i^2-E_n^2}\right),
\]
where the maximum over an empty set is omitted, and the future-veto threshold
\[
F_n=\inf_{j>n}
\frac{q_j^2-q_n^2}{E_n^2-E_j^2}.
\]
Then \(n\) is strictly active iff \(P_n^*<F_n\).

## 2. The future-veto infimum is attained

For every \(j>n\),
\[
0<E_j<E_n,
\]
so
\[
0<E_n^2-E_j^2<E_n^2.
\]
Consequently
\[
\frac{q_j^2-q_n^2}{E_n^2-E_j^2}
>
\frac{q_j^2-q_n^2}{E_n^2}.
\]
As \(j\to\infty\), \(q_j\to\infty\), hence the right-hand side tends to \(+\infty\). Therefore the sequence of future crossing scales tends to \(+\infty\), and its infimum is a minimum attained at some finite index:
\[
\boxed{
F_n=\min_{j>n}
\frac{q_j^2-q_n^2}{E_n^2-E_j^2}.
}
\tag{1}
\]
This already gives a precise answer to one concern: although the future is formally infinite, any fixed candidate is vetoed (or survives) by a finite future witness.

## 3. Explicit stopping certificate

Suppose we have examined \(j=n+1,\ldots,m\), and let
\[
U_m=\min_{n<j\le m}
\frac{q_j^2-q_n^2}{E_n^2-E_j^2}.
\]
For every unexamined \(j\ge m+1\),
\[
\frac{q_j^2-q_n^2}{E_n^2-E_j^2}
>
\frac{q_j^2-q_n^2}{E_n^2}
\ge
\frac{q_{m+1}^2-q_n^2}{E_n^2}.
\]
Thus, if
\[
\boxed{
\frac{q_{m+1}^2-q_n^2}{E_n^2}\ge U_m,
}
\tag{2}
\]
then no future convergent can improve the current minimum. We may stop and certify
\[
F_n=U_m.
\]

Condition (2) is a finite, checkable stopping certificate. Since \(q_{m+1}\to\infty\), it must eventually hold. If the minimum is tied, the non-strict inequality still certifies its value; strict activity is then decided by comparing \(P_n^*\) with \(F_n\).

The important distinction is:
- **Finite lookahead for each fixed \(n\): proved.**
- **A universal bound on how many future convergents are needed, independent of \(n\) and \(\alpha\): not proved and generally not expected.**

## 4. A posteriori decision algorithm

For a fixed \(n\):

1. Compute \(P_n^*\) from the finitely many earlier convergents.
2. Initialize \(U\) with the crossing scale from \(j=n+1\).
3. Examine successive future convergents, updating \(U\) whenever a smaller crossing is found.
4. After each update, test the lower bound in (2), using \(q_{m+1}\).
5. Stop when (2) holds. Then \(F_n=U\), and the exact two-sided criterion decides activity.

This is not a fixed finite-memory online rule: it may need arbitrarily many future terms. But it is a terminating per-candidate certification algorithm.

For rigorous computation with a specified transcendental \(\alpha\), use interval bounds for \(E_n\) and \(E_j\). The stopping test must be verified with outward-rounded intervals; floating-point approximations alone are not a proof.

## 5. Pi: a concrete future veto in the finite experiment

For \(\alpha=\pi\), the convergent \(q=265381\) (numerator \(833719\)) has, in the high-precision finite experiment:
- past admission threshold imposed by \(q=99532\):
  \[
  P_n^*\approx1.9749138\times10^{22};
  \]
- future veto threshold imposed by \(q=364913\) (numerator \(1146408\)):
  \[
  F_n\approx1.2536577\times10^{22}.
  \]
Thus the later convergent crosses below it before it can enter the lower envelope:
\[
F_n<P_n^*.
\]
This is a numerical illustration of genuine future veto. A formal certificate for these decimal values would require interval arithmetic; the theorem above is exact independently of this example.

By contrast, \(q=99532\) has
\[
P_n^*\approx9.5922976\times10^{19},
\qquad
F_n\approx1.5275923\times10^{22},
\]
so it survives in the computed hull.

## 6. Interpretation: future changes the record, not the past

The candidate at index \(n\) is first admitted by the past at scale \(P_n^*\). It remains viable only until the future reaches scale \(F_n\). If the future threshold lies earlier, the candidate is retrospectively removed from the global record.

The streaming stack is the online implementation: it provisionally retains a candidate, then a later point may pop it. The stopping theorem gives the complementary offline fact: for each fixed candidate, a finite future prefix eventually certifies that no still-later point can change its verdict.

This is a mathematically precise two-boundary selection mechanism. It resembles a future-and-past constrained variational problem, but it is ordinary global convex optimization, not physical retrocausality or a Feynman path integral.

## 7. What remains open

1. Find a useful explicit bound for the stopping index in terms of \(q_n\), \(E_n\), and partial quotients.
2. Determine whether the lookahead can be bounded for \(\pi\) by a function of the local partial quotients.
3. Obtain certified interval arithmetic for the displayed pi future-veto example.
4. Investigate whether the finite future-veto certificate can be combined with the BBP arithmetic information without assuming a global bound on all partial quotients.

No lower bound for \(q\|q\pi\|\), and hence no proof that \(\pi\) is badly approximable, is claimed.
