# Algebraic and transcendental values along the cross-base transport curve

Date: 2026-10-31

This note continues the fixed-word transport
\[
F(t)=\sum_{n\ge1}d_nt^{-n},\qquad t>1,
\]
where d_n are the greedy π-base digits of 1:
\[
d_n\in\{0,1,2,3\},\qquad
1=\sum_{n\ge1}d_n\pi^{-n}.
\]
No decimal expansion and no continued fractions are used.

## Theorem: the transport curve sweeps all positive real values

The function F is real analytic on (1,∞), and
\[
F'(t)=-\sum_{n\ge1}n d_nt^{-n-1}<0,
\]
because d_1=3. Thus F is strictly decreasing.

As t→∞, dominated convergence (or the power-series estimate) gives F(t)→0.

The digit sequence has infinitely many nonzero terms. Otherwise (1) would be a finite polynomial relation for π. Since d_n are nonnegative integers, monotone convergence gives
\[
\lim_{t\downarrow1}F(t)=\sum_{n\ge1}d_n=+\infty.
\]
Consequently F is a continuous strictly decreasing bijection
\[
\boxed{F:(1,\infty)\longrightarrow(0,\infty).}
\tag{1}
\]
In particular, it is a homeomorphism.

## Theorem: algebraic-valued bases form a countable dense set

Let
\[
\mathcal A_F=\{t>1:F(t)\text{ is algebraic over }\mathbb Q\}.
\]
The positive algebraic numbers form a countable dense subset of (0,∞). Since F is a bijection, each positive algebraic number has exactly one preimage. Therefore
\[
\boxed{\mathcal A_F=F^{-1}(\overline{\mathbb Q}\cap(0,\infty))
\text{ is countable and dense in }(1,\infty).}
\tag{2}
\]
Its complement
\[
\mathcal T_F=\{t>1:F(t)\text{ is transcendental}\}
\]
is uncountable, has full Lebesgue measure, and is comeager (a dense G_delta set). Indeed, the complement of a countable set has full measure, and the complement of a countable set is comeager because singletons are nowhere dense.

Thus, along this one analytic cross-base curve, transcendental values occur for almost every base parameter, even though algebraic-valued parameters are dense. The two kinds of parameter values are interwoven: there is no interval on which all values are transcendental, because algebraic values occur densely; nor an interval on which all values are algebraic.

This is a rigorous sense in which arithmetic type changes along the base parameter. It is not a continuous transition from algebraic to transcendental: algebraicity is a pointwise, discontinuous arithmetic property.

## What this does not identify

The theorem does not tell us whether F(10) is transcendental. It proves only that F(10) is irrational, using the non-eventual periodicity of d and the rational-expansion criterion in integer bases. The unique parameter t_y for a specified algebraic value y is well-defined, but its own arithmetic type is not classified.

Nor does this prove that the base parameter t is algebraic or transcendental whenever F(t) is transcendental. There is no such implication from the theorem.

## A useful inverse coordinate

Because F is a decreasing homeomorphism, define the inverse base coordinate
\[
\tau(y)=F^{-1}(y),\qquad y>0.
\]
Then
\[
\tau(1)=\pi.
\]
For every positive algebraic y, τ(y) is an algebraic-valued point of the base axis in the sense F(τ(y))=y; this does not say τ(y) is algebraic. The set of all such τ(y) is countable and dense.

The inverse derivative exists everywhere because F'(t)<0:
\[
\tau'(y)=\frac1{F'(\tau(y))}<0.
\]
This gives a second vector field, now using the transported value y as coordinate:
\[
\frac{d\tau}{dy}
=-\left(\sum_{n\ge1}n d_n\tau(y)^{-n-1}\right)^{-1}.
\]
It may be useful when studying how a chosen arithmetic target value determines its unique base.

## Relation to the Laplace moment flow

With t=πe^s and w_n=d_nπ^{-n}, the previous note gives
\[
G(s)=F(\pi e^s)=\sum_{n\ge1}w_ne^{-ns}.
\]
The present theorem adds global range information to the local moment-flow identities: as s runs over (-log π,∞), G decreases continuously from +∞ to 0. The probability-measure representation is valid for s>−log π; the limit at the left endpoint is an infinite-mass boundary after exponential tilting, not a probability distribution at that endpoint.

## Consequence and next question

The fixed-word transport construction is not merely an example where one integer value becomes irrational at a larger integer base. It gives a real-analytic bijection whose algebraic-value parameter set is countable dense and whose transcendental-value parameter set has full measure and is comeager.

The next arithmetic question is now more precise: can one identify even one explicit parameter t (for example t=10) for which F(t) is provably transcendental? The measure-category theorem cannot answer this explicit-value question. A route through digit complexity would require an effective complexity theorem for the specific π-generated word d; no such bound is proved here.

No proof or disproof of P2848 is claimed.
