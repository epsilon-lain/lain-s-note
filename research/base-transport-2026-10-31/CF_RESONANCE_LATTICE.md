# The continued-fraction resonance lattice of the migrating π word

Date: 2026-10-31

This note adds a geometric interpretation to [CONTINUED_FRACTION_MIGRATION.md](CONTINUED_FRACTION_MIGRATION.md). It distinguishes two different continued-fraction processes:

1. **Read the continued fraction of the moving value** Ψ(t) as t changes. Its partial quotients are piecewise constant and jump when Ψ(t) crosses a rational.
2. **Transport the rational convergents of π back to the base axis** by the inverse map σ=Ψ^{-1}. These are resonance points t_n=σ(p_n/q_n).

The second process gives a direct geometric encoding of the partial quotients of π.

## 1. Setup and derivative at the native point

Let
\[
\Phi(t)=\sum_{j\ge1}d_jt^{-j},\qquad
d_1=3,\ d_2=0,\ d_3=1,\quad
\Phi(\pi)=1,
\]
and define
\[
\Psi(t)=t\Phi(t),\qquad \Psi(\pi)=\pi.
\]
Put
\[
m=\sum_{j\ge1}j d_j\pi^{-j}.
\]
Termwise differentiation gives
\[
\boxed{\Psi'(\pi)=1-m<0.}
\tag{1}
\]
The map Ψ is strictly decreasing on (1,∞), with range (3,∞), so it has a real-analytic inverse σ:(3,∞)→(1,∞). In particular,
\[
\boxed{\sigma'(\pi)=-\frac1{m-1}.}
\tag{2}
\]

## 2. Convergents of π become resonance points on the base axis

Let
\[
\frac{p_n}{q_n}=[a_0;a_1,\ldots,a_n]
\]
be the ordinary continued-fraction convergents of π, with q_n>0. For all sufficiently large n, p_n/q_n>3, so define
\[
\boxed{t_n=\sigma(p_n/q_n),\qquad
\Psi(t_n)=p_n/q_n.}
\tag{3}
\]
These are the base parameters at which the migrating curve hits the rational convergents of π.

Since p_n/q_n→π,
\[
t_n\to\pi.
\]
The inverse-function Taylor expansion and (2) yield
\[
\boxed{
t_n-\pi=-\frac{p_n/q_n-\pi}{m-1}
+O\left((p_n/q_n-\pi)^2\right).
}
\tag{4}
\]
Because convergents satisfy \(|\pi-p_n/q_n|<1/q_n^2\), multiplying by q_n² gives
\[
\boxed{
q_n^2|t_n-\pi|
=\frac{q_n^2|\pi-p_n/q_n|}{m-1}
+O(q_n^{-2}).
}
\tag{5}
\]
Thus the normalized approximation quality of the convergents is transported to the base axis with the fixed scale factor 1/(m−1).

In particular,
\[
\boxed{
\liminf_n q_n^2|t_n-\pi|=0
\iff
\liminf_n q_n^2|\pi-p_n/q_n|=0.
}
\tag{6}
\]
This is a geometric reformulation of the bad-approximation question along the resonance sequence. It is not a proof that either liminf is zero.

## 3. Adjacent resonance gaps encode the next partial quotient

Consecutive convergents satisfy
\[
|p_{n+1}q_n-p_nq_{n+1}|=1.
\]
Therefore
\[
\left|\frac{p_{n+1}}{q_{n+1}}-\frac{p_n}{q_n}\right|
=\frac1{q_nq_{n+1}}.
\]
By the mean value theorem, for some ξ_n between t_n and t_{n+1},
\[
\boxed{
|t_{n+1}-t_n|
=\frac1{q_nq_{n+1}|\Psi'(\xi_n)|}.
}
\tag{7}
\]
Since both endpoints tend to π, ξ_n→π. Hence
\[
\boxed{
q_nq_{n+1}|t_{n+1}-t_n|
\longrightarrow\frac1{m-1}.
}
\tag{8}
\]

The denominator recurrence is
\[
q_{n+1}=a_{n+1}q_n+q_{n-1}.
\]
Consequently
\[
q_n^2|t_{n+1}-t_n|
=
\frac{1+o(1)}
{(m-1)\left(a_{n+1}+q_{n-1}/q_n\right)}.
\tag{9}
\]
A large next partial quotient creates an unusually small gap between adjacent resonance points. Thus the partial quotient is visible as a compression factor in the base-parameter geometry.

## 4. Semiconvergents split a parent resonance interval into a chain

Fix n and let a=a_{n+1}. Define the intermediate fractions
\[
\frac{p_{n,k}}{q_{n,k}}
=\frac{k p_n+p_{n-1}}{k q_n+q_{n-1}},
\qquad k=0,1,\ldots,a.
\]
The endpoints k=0 and k=a are p_{n-1}/q_{n-1} and p_{n+1}/q_{n+1}. For sufficiently large n all these fractions exceed 3, so their resonance points
\[
t_{n,k}=\sigma(p_{n,k}/q_{n,k})
\]
are defined.

Adjacent semiconvergents have determinant one:
\[
|p_{n,k+1}q_{n,k}-p_{n,k}q_{n,k+1}|=1.
\]
Therefore each neighboring base gap obeys the exact mean-value formula
\[
\boxed{
|t_{n,k+1}-t_{n,k}|
=\frac1{q_{n,k}q_{n,k+1}|\Psi'(\xi_{n,k})|}.
}
\tag{10}
\]
There are exactly a_{n+1} such subintervals. Their rational lengths telescope:
\[
\sum_{k=0}^{a-1}
\left|\frac{p_{n,k+1}}{q_{n,k+1}}-
\frac{p_{n,k}}{q_{n,k}}\right|
=
\sum_{k=0}^{a-1}\frac1{q_{n,k}q_{n,k+1}}
=
\frac{a}{q_{n-1}q_{n+1}}.
\tag{11}
\]
The right side is exactly the width of the parent Farey interval between p_{n-1}/q_{n-1} and p_{n+1}/q_{n+1}, because their determinant has absolute value a.

Thus the partial quotient a_{n+1} counts the number of elementary resonance gaps into which the parent interval is split. Under σ, the individual gaps are distorted by the local factor 1/|Ψ'|, but their order and adjacency are preserved (with orientation reversed).

This is the most concrete geometric interpretation obtained here:
- convergents are resonance points;
- semiconvergents are intermediate resonance points;
- the next partial quotient counts the number of subgaps;
- large partial quotients create many closely spaced subgaps.

## 5. Two different CF questions must not be conflated

The integer-base theorem in [CONTINUED_FRACTION_MIGRATION.md](CONTINUED_FRACTION_MIGRATION.md) describes the continued fraction of Ψ(b) for integer b→∞. For example, its first partial quotients have explicit polynomial and residue-class patterns.

The resonance construction above instead fixes the continued fraction of π and transports its rational convergents to the base axis. The first question studies the itinerary of Ψ(t) at a selected t; the second follows a fixed sequence of rational targets through the inverse map.

They are related by the same rational-hit equation Ψ(t)=p/q, but they involve different limits:
- integer-base partial quotient laws: fixed CF depth, base b→∞;
- π resonance laws: fixed target π, CF depth n→∞.

The fixed-depth eventual quasi-polynomial theorem does not imply any bound on the partial quotients of π. Conversely, the resonance-gap formula does not prove a law for Ψ(b) as b grows.

## 6. Consequence for P2848 and remaining obstacle

For convergents of π, (5) shows that the normalized base resonance distance and the usual normalized rational approximation error have the same liminf up to a positive constant. Since for irrational x the condition
\[
\liminf_{q\to\infty}q\|qx\|=0
\]
is equivalent to unbounded partial quotients in its ordinary continued fraction, P2848 can be read geometrically as asking whether the resonance sequence has arbitrarily deep compression.

This does not simplify away the hard arithmetic: the sequence of partial quotients is exactly the information being tested. What the migration model adds is a geometric observable for it—the gaps between rational-hit base parameters.

No proof or disproof of P2848 is claimed. No transcendence claim for Φ(10) is made.

## 7. Reproducibility and literature

The identities (4)–(11) follow from the inverse function theorem, the determinant identity for convergents, and the mean value theorem. A numerical certificate may verify selected resonance intervals, but the asymptotic equivalences are analytic statements.

The use of formal Laurent-series continued fractions and specialization is established; see:
- A. J. van der Poorten and J. Shallit, [A Specialised Continued Fraction](https://doi.org/10.4153/CJM-1993-058-5), *Canadian Journal of Mathematics* 45 (1993), 1067–1079.
- A. J. van der Poorten, [Specialisation and reduction of continued fractions of formal power series](https://doi.org/10.1007/s11139-005-0827-3), *Ramanujan Journal* 9 (2005), 83–91.
- Qian Hu, [Certain Hidden Stabilities of Continued Fraction Arithmetic](https://scholars.unh.edu/dissertation/2988), University of New Hampshire dissertation (2026), whose abstract studies quasi-polynomial regularity of partial quotients for sequences governed by rational asymptotic expansions.

The present note does not claim those general methods as new. Its contribution is the resonance-point and semiconvergent-gap interpretation for the specific π-native transport Ψ.
