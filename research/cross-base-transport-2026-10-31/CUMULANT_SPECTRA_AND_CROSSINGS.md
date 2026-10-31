# Cumulant spectra, object identification, and crossings of base-transport curves

Date: 2026-10-31

This note studies a general fixed digit object
\[
a=(a_n)_{n\ge1},\qquad 0\le a_n\le A,
\]
through its base-transport function
\[
Z_a(t)=\sum_{n\ge1}a_nt^{-n},\qquad t>1,
\]
assuming a is not identically zero. Put x=log t and define the tilted position law
\[
P_{a,x}(n)=\frac{a_ne^{-nx}}{Z_a(e^x)}.
\]
Its cumulants are denoted κ_{a,j}(x), with κ_1 the mean and κ_2 the variance.

## 1. The full cumulant spectrum identifies the digit object

**Theorem.** Fix a known base t>1. The data
\[
\mathcal M_a(t)=\bigl(Z_a(t),\kappa_{a,1}(t),\kappa_{a,2}(t),\ldots\bigr)
\]
uniquely determine the entire coefficient sequence a_n.

**Proof.** The moment generating function of P_{a,t} is
\[
M_{a,t}(s)=\sum_{n\ge1}P_{a,t}(n)e^{sn}
=\frac{Z_a(te^{-s})}{Z_a(t)}.
\]
It is analytic for |s|<log t (indeed for Re(s)<log t). Its logarithm has Taylor expansion
\[
\log M_{a,t}(s)=\sum_{j\ge1}\kappa_{a,j}(t)\frac{s^j}{j!}
\]
near s=0. Hence the complete cumulant sequence determines M_{a,t} in a neighborhood of zero. Uniqueness of analytic continuation, or uniqueness of the moment generating function, determines the probability law P_{a,t}(n). Finally
\[
\boxed{a_n=Z_a(t)P_{a,t}(n)t^n.}
\tag{1}
\]
Thus the spectrum plus the total value reconstructs every digit exactly.

In fact, because cumulants determine the normalized law, Z_a(t) is the one additional scalar needed to recover the absolute coefficient scale. If the coefficients are required to lie in a fixed digit alphabet with gcd 1, this scale can sometimes be inferred from the normalized law and integrality, but no such inference is needed for (1).

**Important qualification.** A finite list of cumulants does not generally determine an infinite digit object. The theorem requires the full infinite sequence and a known base t. Nor does the scalar value Z_a(t) alone identify the object.

## 2. Spectral representation and invariants

The characteristic function
\[
\chi_{a,t}(\theta)=\sum_{n\ge1}P_{a,t}(n)e^{in\theta}
\]
is the Fourier spectrum of the tilted position law. Its Fourier coefficients are precisely P_{a,t}(n). Thus the full cumulant spectrum and the full characteristic spectrum encode the same normalized positional information.

For any two active positions m,n,
\[
\frac{P_{a,t}(n)}{P_{a,t}(m)}
=\frac{a_n}{a_m}t^{m-n}.
\]
Therefore
\[
\frac{d}{d\log t}\log\frac{P_{a,t}(n)}{P_{a,t}(m)}=m-n.
\]
For any three active positions i<j<k, the projective invariant
\[
\boxed{
I_{i,j,k}(a):=
\frac{P_j^{k-i}}{P_i^{k-j}P_k^{j-i}}
=\frac{a_j^{k-i}}{a_i^{k-j}a_k^{j-i}}
}
\tag{2}
\]
is independent of t. It is also unchanged if all coefficients are multiplied by the same positive scalar. These invariants are part of the projective digit profile.

The support S(a)={n:a_n>0} is invariant under transport. A qualification is necessary: on a common support, the complete family of triple invariants determines positive coefficients only up to b_n=C·λ^n a_n, with C,λ>0. Indeed, the logarithms of b_n/a_n must be affine in n, and this affine freedom also preserves every triple invariant. For example, the finite objects (1,1,1) and (1,2,4) have identical triple invariants but different coefficient ratios. Fixing two reference coefficients removes both freedoms. If the common support is infinite and both objects have bounded positive integer coefficients there, then λ≠1 would force one coefficient ratio to tend to zero or infinity, so λ=1; only in that additional setting do the invariants determine all coefficient ratios without two anchors. Singleton support is trivial.

## 3. Crossings of two transport curves are governed by cumulant differences

Let a,b be two nonzero bounded nonnegative digit objects and define
\[
R(x)=\log Z_a(e^x)-\log Z_b(e^x).
\]
Then
\[
\boxed{
R^{(j)}(x)=(-1)^j\bigl(\kappa_{a,j}(x)-\kappa_{b,j}(x)\bigr),
\qquad j\ge1.
}
\tag{3}
\]
Suppose the curves cross at t_0=e^{x_0}, so Z_a(t_0)=Z_b(t_0). Their logarithmic contact is classified by the first cumulant at which they differ:
- if their means differ, the crossing is transverse;
- if means agree but variances differ, the curves have quadratic contact;
- if the first difference is in the j-th cumulant, their log-ratio has a zero of exact order j at x_0.

If all cumulants agree at x_0, the two tilted probability laws agree by the moment-generating-function uniqueness theorem. Since the total values also agree, (1) gives a_n=b_n for every n. Therefore **two distinct digit objects cannot have infinite-order contact at a positive base**.

This is a strong rigidity result: a crossing can be arbitrarily flat to any prescribed finite order in principle, but exact equality of every cumulant forces equality of the underlying object.

The slope of the log-ratio is
\[
R'(x)=-(\mu_a(x)-\mu_b(x)).
\]
Thus a crossing is transverse exactly when the two objects have different mean active positions at that base.

## 4. Finite-jet matching and what it cannot prove

The first J cumulants describe the first J derivatives of log Z at a chosen base. Matching these gives contact of order at least J+1 when the values also agree. But finite-jet matching does not imply equality of the digit sequences.

This separates two tasks:
- exact reconstruction requires the full infinite cumulant spectrum;
- numerical or finite symbolic experiments can only test finite-order contact.

For the π-generated word d, the spectrum is well-defined at every t>1. Its first cumulant is the global coordinate μ(t) studied in [GLOBAL_SENSITIVITY_COORDINATE.md](GLOBAL_SENSITIVITY_COORDINATE.md); the second is the variance controlling its speed; higher cumulants describe increasingly fine shape information.

## 5. Arithmetic type and the explicit value problem

The transport curve for d is real analytic and strictly decreasing. Its algebraic-valued parameter set is countable dense, and its transcendental-valued parameter set has full measure and is comeager. These generic statements do not decide the arithmetic type of a specified value such as Φ(10).

The full cumulant spectrum is a complete encoding of the digit object, but it does not by itself supply an arithmetic theorem about Φ(10). To prove transcendence at t=10, one needs an additional explicit arithmetic property of the coefficient word d or of the value series. Non-eventual periodicity proves irrationality only.

Potential routes:
1. derive effective block-complexity bounds for d from its π-remainder recursion;
2. use integer recurrences, Hankel minors, or modular obstructions to constrain a hypothetical algebraic relation for Φ(10);
3. compare the π word with other source-base words and identify invariants that distinguish their algebraic relations;
4. study whether finite cumulant jets at t=10 admit certified lower bounds incompatible with algebraicity. No such lower bound is currently proved.

No proof of transcendence of Φ(10), and no proof or disproof of P2848, is claimed.

## 6. Reproducible finite experiments

A finite verifier should:
- generate a certified prefix of each source word;
- compute Z and tilted moments with rational interval bounds at a selected algebraic t;
- convert moments to cumulants using the standard triangular moment-cumulant relations;
- compare crossing values and the first differing cumulant;
- report a finite-jet match only as such, never as equality of infinite objects.

The exact theorems above require no numerical experiment. Finite calculations can illustrate crossing order and estimate which cumulants distinguish two objects, but cannot prove transcendence or infinite spectral equality.
