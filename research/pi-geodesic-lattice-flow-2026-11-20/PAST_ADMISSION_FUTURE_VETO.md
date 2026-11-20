# Past admission and future veto: a two-sided principle for the continuous lattice flow

Date: 2026-11-20  
Status: **exact theorem and explicit counterexample to bounded-lookahead rules**. This is a mathematical “future changes the past” mechanism in global optimization, not physical retrocausality.

## 1. The lower envelope is a global boundary-value object

For an irrational target \(\alpha\), let the candidate convergent points be
\[
P_n=(X_n,Y_n)=(q_n^2,E_n^2),\qquad E_n=|q_n\alpha-p_n|.
\]
At scale \(x=e^{4t}\ge1\), the candidate cost is
\[
\ell_n(x)=X_n+xY_n.
\]
The global cost is the lower envelope
\[
\boxed{M(x)=\inf_{n\ge0}\{X_n+xY_n\}.}
\tag{1}
\]
As an infimum of affine functions, \(M\) is concave and nondecreasing. Its active slope is one of the values \(Y_n\), and the active index can only move forward as \(x\) increases.

The decision whether index \(n\) is active is not determined by its past alone: it must admit a scale \(x\) that beats every earlier and every later candidate.

## 2. Past admission threshold

For \(i<n\), the condition \(\ell_n(x)<\ell_i(x)\) is equivalent to
\[
x>\frac{X_n-X_i}{Y_i-Y_n}.
\]
Define the threshold imposed by the entire past:
\[
\boxed{
P_n=\sup_{i<n}\frac{X_n-X_i}{Y_i-Y_n},
}
\tag{2}
\]
with \(P_0=0\). The point \(n\) cannot become active before \(x>P_n\). Include the physical domain \(x\ge1\) by setting
\[
P_n^\ast=\max(1,P_n).
\]
Interpretation: the past determines the earliest scale at which this candidate can enter the lower envelope.

## 3. Future veto threshold

For \(j>n\), the condition \(\ell_n(x)<\ell_j(x)\) is equivalent to
\[
x<\frac{X_j-X_n}{Y_n-Y_j}.
\]
Define
\[
\boxed{
F_n=\inf_{j>n}\frac{X_j-X_n}{Y_n-Y_j}.
}
\tag{3}
\]
Interpretation: the future determines the latest scale at which the current candidate can remain active.

Combining all inequalities gives the exact two-sided criterion:
\[
\boxed{
n\text{ is strictly active}
\iff
P_n^\ast<F_n.
}
\tag{4}
\]
When active, its exact interval is
\[
\boxed{P_n^\ast<x<F_n.}
\tag{5}
\]
Thus the past supplies an admission boundary and the future supplies a veto boundary. A future point \(j\) retroactively removes \(n\) precisely when its crossing scale
\[
x_{n,j}=\frac{X_j-X_n}{Y_n-Y_j}
\]
falls at or below \(P_n^\ast\).

This is the precise sense in which the future “decides the past” in our construction: the complete path is selected by a two-sided global minimization. No signal travels backward in physical time.

## 4. Explicit example: same past, opposite future verdict

Consider the common prefix
\[
A=(0,10),\qquad B=(1,9).
\]
For \(x\ge1\),
\[
\ell_A(x)=10x,\qquad \ell_B(x)=1+9x,
\]
so \(B\) is strictly active for every \(x>1\) in this two-point prefix.

Now consider two possible futures.

### Future I: preserve B
Append
\[
C_+=(2,8.5).
\]
Then
\[
\ell_{C_+}(x)=2+8.5x.
\]
The crossings are \(A/B\) at \(x=1\), and \(B/C_+\) at \(x=2\). Thus \(B\) remains active on
\[
1<x<2.
\]

### Future II: veto B
Instead append
\[
C_-=(2,7).
\]
Then
\[
\ell_{C_-}(x)=2+7x.
\]
The crossing of \(A\) and \(C_-\) is at
\[
x=\frac{2-0}{10-7}=\frac23.
\]
For every \(x\ge1\), \(C_-\) is already below \(B\), so \(B\) is never active.

The two sequences have exactly the same past through \(B\), but the future changes B's status from active to inactive.

The example can be delayed arbitrarily far: insert any prescribed finite number of intermediate points on a shallow decreasing-slope segment after \(B\), then append a sufficiently steep final point. Hence no rule that inspects only a fixed number of future candidates can decide activity for arbitrary monotone candidate sequences.

This last statement is about arbitrary monotone point sequences. It does **not** prove that no finite-window theorem exists for the much more rigid continued-fraction sequence of \(\pi\).

## 5. Connection to the streaming stack and weighted PAVA

The streaming lower-hull algorithm is a causal way to maintain the *current finite-prefix answer*. When a new point arrives, it may pop a suffix of previously retained points. This is exactly the algorithmic form of the future veto in (3).

The two-sided criterion (4) is the offline, infinite-horizon characterization. The stack is its online finite-prefix implementation. The weighted PAVA identity explains how each veto merges adjacent slope blocks.

There is no contradiction between the two:
- a point may be active for every prefix seen so far;
- a later point may lower its future threshold below its past admission threshold;
- once removed, it can never return, because the chord witnessing its removal remains present.

## 6. A Feynman-like analogy, with a strict boundary

The mathematical structure resembles a two-boundary variational problem: a candidate is accepted only if it is compatible with constraints from both earlier and later data. The action-like quantity is \(\ell_n(x)\), and the selected trajectory is the global lower envelope.

This is only an analogy to global variational selection. It is not Feynman's path integral, Wheeler–Feynman electrodynamics, or a physical advanced signal. Our model is deterministic convex optimization over integer approximation candidates.

## 7. What to investigate next

1. For continued fractions, determine whether the future veto \(F_n\) can be bounded using a short block of future partial quotients plus a uniform tail interval.
2. Determine whether the past admission \(P_n^\ast\) has a similar finite-state description from reversed partial quotients.
3. Compare the two bounds to obtain certified active intervals without constructing the full hull.
4. For \(\alpha=\pi\), test the bounds against the convergents that survive or are removed by the streaming stack.

The key unresolved question is now precise: can the two infinite boundary functionals \(P_n^\ast,F_n\) be compressed to finite continued-fraction data for this special target? No such compression is claimed yet.
