# Crossing networks: Descartes bounds and certified finite-prefix comparisons

Date: 2026-10-31

This note extends the spectral-crossing theory. It focuses on a question not answered by isolatedness alone: how many times can two transport curves cross, and how can an infinite π-generated curve be compared rigorously with finite objects?

Let a,b be finitely supported nonnegative sequences and
\[
Z_a(t)=\sum_{n\ge1}a_nt^{-n},\qquad
Z_b(t)=\sum_{n\ge1}b_nt^{-n}.
\]
Put u=1/t and
\[
H(u)=\sum_{n\ge1}(a_n-b_n)u^n.
\]

## 1. Descartes crossing bound

Let V(a-b) be the number of sign changes in the nonzero coefficient list
\[
a_1-b_1,a_2-b_2,\ldots
\]
after zero coefficients are deleted.

**Theorem.** The number of crossings Z_a(t)=Z_b(t) for t>1, counted with multiplicity, is at most V(a-b).

**Proof.** H is a real polynomial. Descartes' rule of signs bounds the number of its positive roots, counted with multiplicity, by the sign changes in its coefficient sequence. Every crossing t>1 corresponds bijectively to a positive root u=1/t in (0,1). Restricting to (0,1) cannot increase the count. ∎

Consequences:
- If the coefficient difference has no sign change, distinct curves never cross.
- One sign change permits at most one crossing.
- A root of multiplicity k consumes at least k units of the sign-change budget.
- The order-k construction H_k(u)=u(1-2u)^k has exactly k sign changes and one root u=1/2 of multiplicity k, so the bound is sharp for every k.

This turns the crossing graph of any finite family of finite objects into a controlled arrangement: each pairwise edge has a crossing count bounded by the sign variation of that pair's coefficient difference.

## 2. A finite-prefix certificate for an infinite object

Now let a be any bounded nonnegative sequence with 0≤a_n≤A, and let b be a finite or infinite bounded nonnegative sequence with 0≤b_n≤B. Define
\[
H(u)=\sum_{n\ge1}(a_n-b_n)u^n,
\qquad
H_N(u)=\sum_{n=1}^N(a_n-b_n)u^n.
\]
For 0<u<1, the tail satisfies
\[
\boxed{
|H(u)-H_N(u)|
\le \max(A,B)\frac{u^{N+1}}{1-u}.
}
\tag{1}
\]
Indeed |a_n-b_n|≤max(A,B), and the geometric tail sums exactly.

Thus, at any specified u_0∈(0,1), if
\[
|H_N(u_0)|>
\max(A,B)\frac{u_0^{N+1}}{1-u_0},
\tag{2}
\]
then H(u_0) has the same sign as H_N(u_0). With rational interval bounds for u_0 and the finite polynomial, (2) becomes a fully certifiable sign test.

On a compact interval u∈[u_-,u_+]⊂(0,1), the uniform tail bound is
\[
|H-H_N|\le \max(A,B)\frac{u_+^{N+1}}{1-u_+}.
\tag{3}
\]
If the finite polynomial H_N has a certified sign margin larger than this bound throughout the interval, the full infinite curves cannot cross there.

For the π-native digit object d_n∈{0,1,2,3} compared with another digit object b_n∈{0,...,B}, one may take A=3. This gives an exact finite-prefix route to exclude crossings on any compact base interval bounded away from 1.

## 3. Cumulant contact and sign variation fit together

At a crossing t_0, the multiplicity of the positive root u_0=1/t_0 of H(u) is exactly the contact order of Z_a and Z_b as functions of log t. It is also the first order j at which their cumulants differ.

Therefore for finite objects:
\[
\boxed{
\text{first differing cumulant order}
=\text{root multiplicity}
\le V(a-b).
}
\tag{4}
\]
This is a direct bridge between three descriptions:
- spectral geometry (curve contact);
- probability statistics (cumulants);
- discrete arithmetic (sign changes of integer coefficients).

In particular, if two finite nonnegative integer objects have coefficient difference with only r sign changes, their spectra cannot have contact order exceeding r at any crossing.

## 4. Examples

### Example A: no crossing

If a_n≥b_n for all n and at least one inequality is strict, then H(u)>0 for every u>0. The curves never cross. This is the zero-sign-change case.

### Example B: one crossing

Let a=(1,0,0,...) and b=(0,2,0,...). Then
\[
H(u)=u-2u^2=u(1-2u).
\]
There is exactly one crossing at u=1/2, i.e. t=2.

### Example C: sharp quadratic contact

Let a=(1,0,4,0,...) and b=(0,4,0,...). Then
\[
H(u)=u-4u^2+4u^3=u(1-2u)^2.
\]
There are two sign changes and one double root at u=1/2. The curves have exactly quadratic contact.

## 5. What is still missing for Φ(10)

For the π-generated word d, comparison with a finite candidate b can be reduced to a polynomial prefix plus the explicit tail bound (1). This can certify that no crossing occurs in a chosen interval, or that a crossing must occur between two certified sign points.

It does not by itself prove transcendence of Φ(10): a transcendence proof must exclude every algebraic relation, not merely equality with one finite transport curve. Nor does it prove or disprove P2848.

The immediate computational project is to build a certified crossing-network explorer:
1. choose a finite library of source words (π, √2, √3, integer-base words, and designed controls);
2. compute exact finite-prefix difference polynomials;
3. use Descartes' rule to bound possible crossings;
4. use rational interval arithmetic and (1) to isolate actual crossings;
5. record contact multiplicity and the first differing cumulant.

This gives a finite, auditable map of the spectral geometry without confusing numerical intersections with proven ones.
