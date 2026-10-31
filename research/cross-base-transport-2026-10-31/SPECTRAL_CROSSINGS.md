# Spectral crossings of digit objects: projective invariants and exact contact order

Date: 2026-10-31

Let a=(a_n), b=(b_n) be nonzero bounded nonnegative sequences. Define
\[
Z_a(t)=\sum_{n\ge1}a_nt^{-n},\quad
P_{a,t}(n)=\frac{a_nt^{-n}}{Z_a(t)},\quad t>1.
\]
Write x=log t and let κ_{a,j}(t) denote the j-th cumulant of P_{a,t}.

## 1. Crossings are isolated unless the objects coincide

The difference H(u)=Σ(a_n-b_n)u^n is analytic for |u|<1. If Z_a(t)=Z_b(t) on a set with an accumulation point in (1,∞), then H has zeros accumulating at an interior point u=1/t. The identity theorem gives H≡0 and a_n=b_n for every n.

Thus distinct objects have isolated crossings in t>1 and only finitely many crossings on every compact interval [1+ε,T]. Infinite crossings can accumulate only at t↓1 or escape to infinity.

## 2. Contact order is the first differing cumulant

Set
\[
R(x)=\log Z_a(e^x)-\log Z_b(e^x).
\]
Then for every j≥1,
\[
R^{(j)}(x)=(-1)^j(κ_{a,j}(e^x)-κ_{b,j}(e^x)).
\tag{1}
\]
At a crossing x_0, R(x_0)=0. If j is the smallest order at which the cumulants differ, R has a zero of exact order j. Since
\[
Z_a(e^x)-Z_b(e^x)=Z_b(e^x)(e^{R(x)}-1),
\]
the transport curves have exactly the same contact order.

If all cumulants agree, the moment generating functions of the tilted position laws agree near zero, so the laws agree. At a crossing the total values also agree; hence a_n=b_n for all n. Distinct objects cannot have infinite-order contact.

## 3. Every finite contact order occurs

For each integer k≥1, define
\[
H_k(u)=u(1-2u)^k
=\sum_{j=0}^k\binom{k}{j}(-2)^j u^{j+1}.
\]
Split the positive and negative coefficients:
\[
a_{j+1}=\max\{\binom{k}{j}(-2)^j,0\},\quad
b_{j+1}=\max\{-\binom{k}{j}(-2)^j,0\}.
\]
Then a,b are distinct finitely supported nonnegative integer objects and
\[
Z_a(t)-Z_b(t)=H_k(1/t).
\]
At t=2 the difference has a zero of exact order k, and both values are positive.

Examples:

- **Transverse crossing:** a=(1,0,0,...), b=(0,2,0,...). Then Z_a=t^{-1}, Z_b=2t^{-2}; they cross at t=2 with means 1 and 2.
- **Quadratic tangency:** a=(1,0,4,0,...), b=(0,4,0,...). Then Z_a-Z_b=u(1-2u)^2. At t=2 both values are 1. The first tilted law has mass 1/2 at positions 1 and 3, hence mean 2 and variance 1. The second is concentrated at position 2, hence mean 2 and variance 0. The curves have exact quadratic contact.

Thus arbitrary finite contact order is possible, while infinite-order contact forces equality of the objects.

## 4. A cross-object projective invariant

When a_n,b_n,a_m,b_m are positive,
\[
\frac{P_{a,t}(n)}{P_{b,t}(n)}
=\frac{a_n}{b_n}\frac{Z_b(t)}{Z_a(t)}.
\]
Taking the ratio at positions n and m cancels the partition functions:
\[
\boxed{
\frac{P_{a,t}(n)P_{b,t}(m)}
{P_{b,t}(n)P_{a,t}(m)}
=\frac{a_nb_m}{b_na_m}.
}
\tag{2}
\]
This is independent of t. The relative log-odds profile of two objects is transported rigidly.

For one object, the three-position invariant is
\[
\frac{P_{a,t}(j)^{k-i}}
{P_{a,t}(i)^{k-j}P_{a,t}(k)^{j-i}}
=\frac{a_j^{k-i}}{a_i^{k-j}a_k^{j-i}},
\quad i<j<k,
\]
also independent of t.

## 5. The full spectrum fingerprints the object

At a known base t_0, the moment generating function of the tilted law is
\[
M_{a,t_0}(s)=\frac{Z_a(t_0e^{-s})}{Z_a(t_0)}.
\]
Its logarithm has Taylor coefficients κ_{a,j}(t_0). The complete cumulant sequence determines M near zero, hence the tilted probability law. The total value then recovers
\[
a_n=Z_a(t_0)P_{a,t_0}(n)t_0^n.
\]
Therefore the total value plus the complete cumulant spectrum uniquely identifies the object. A finite truncation gives only a finite jet and cannot generally recover the infinite word.

## 6. Relevance to the π-generated word and remaining gap

For d=(3,0,1,1,0,2,...), the spectrum (Φ(t),κ_1(t),κ_2(t),...) is a complete fingerprint at a known t>1. The crossing and contact theorems apply when comparing d with any other bounded nonnegative digit object.

These results classify the geometry of curve intersections, but do not prove that Φ(10) is transcendental. Nonperiodicity gives irrationality only. A transcendence proof still needs explicit arithmetic information about the π-generated coefficients, such as effective block complexity, recurrence obstructions, or modular constraints.

No proof or disproof of P2848 is claimed.
