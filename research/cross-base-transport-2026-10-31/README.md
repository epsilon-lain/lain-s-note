# Cross-base transport of a native π digit word: a Laplace flow and arithmetic migration

Date: 2026-10-31

## Status and scope

This note fixes the digit word generated intrinsically by the greedy π-expansion of 1, then transports that *same symbolic object* to every real base t>1. It does not recompute the digits in each base and does not use decimal digits or continued fractions.

The main rigorous conclusions are:
1. the transported value is a strictly decreasing, real-analytic, completely monotone function in logarithmic base-time;
2. its logarithmic derivatives are moments of one fixed probability measure on the positive integers;
3. the value is exactly 1 at t=π, but is irrational at every integer base b≥4, in particular at 10 and every larger integer base.

No transcendence claim is proved for the transported values at b≥4. In particular, irrationality at base 10 does not imply transcendence.

## 1. The fixed symbolic object

Let
\[
r_0=1,\qquad d_n=\lfloor \pi r_{n-1}\rfloor,\qquad
r_n=\pi r_{n-1}-d_n.
\]
Then d_n is an integer in {0,1,2,3}, 0<r_n<1, and
\[
1=\sum_{n\ge1}d_n\pi^{-n}.
\tag{1}
\]
The word d=d_1d_2... is not eventually periodic. Indeed, if it were eventually periodic, the right side of (1) would be a rational function of π with rational coefficients; clearing its denominator would give a nonzero polynomial over Q vanishing at π, contradicting the transcendence of π.

The word d is now held fixed. For every real t>1 define its cross-base evaluation
\[
\boxed{F(t)=\sum_{n\ge1}d_nt^{-n}.}
\tag{2}
\]
This is a transport map from a symbolic word to a value in base t. It is not the greedy t-expansion of 1 in general.

## 2. Exact migration from an integer to irrational values

Since d_n≥0 and at least one d_n is positive,
\[
F'(t)=-\sum_{n\ge1}n d_n t^{-n-1}<0.
\]
Thus F is strictly decreasing on (1,∞). Equation (1) gives
\[
F(\pi)=1.
\]
Consequently,
\[
1<t<\pi\implies F(t)>1,\qquad
t>\pi\implies 0<F(t)<1.
\tag{3}
\]
This is an exact crossing: the same symbolic object has integer value 1 in its native base and a value strictly between 0 and 1 in every larger base.

For an integer b≥4, the digits d_n belong to the standard base-b alphabet {0,...,b−1}. Hence (2) is a base-b expansion of F(b). A rational number has an eventually periodic base-b expansion. Since d is not eventually periodic, we obtain
\[
\boxed{F(b)\notin\mathbb Q\quad\text{for every integer }b\ge4.}
\tag{4}
\]
In particular F(10) is irrational, despite F(π)=1.

The theorem is only irrationality. The nonperiodicity of a digit word alone does not prove that its value is transcendental. Establishing transcendence of F(10), or of F(b) for any integer b≥4, is an open question in this project.

## 3. Logarithmic base-time as a Laplace transform

Write the base as
\[
t=\pi e^s,\qquad s> -\log\pi.
\]
Define weights
\[
w_n=d_n\pi^{-n}.
\]
By (1), w_n≥0 and \sum_n w_n=1. Thus
\[
\mu=\sum_{n\ge1}w_n\,\delta_n
\]
is a probability measure on the positive integers. The transported value becomes
\[
\boxed{G(s):=F(\pi e^s)=\sum_{n\ge1}w_ne^{-ns}
=\int e^{-sx}\,d\mu(x).}
\tag{5}
\]
So cross-base transport is exactly the Laplace transform of one fixed discrete probability measure.

For every integer k≥0 and s>−log π, termwise differentiation is justified by exponential convergence and gives
\[
(-1)^kG^{(k)}(s)=\sum_{n\ge1}n^kw_ne^{-ns}>0.
\tag{6}
\]
Therefore G is completely monotone. Its derivatives encode the moments of the fixed measure after exponential tilting.

At the native base s=0,
\[
G(0)=1,\qquad
(-1)^kG^{(k)}(0)=M_k:=\sum_{n\ge1}n^kd_n\pi^{-n}.
\tag{7}
\]
All M_k are finite. In particular, the constant m from the previous intrinsic-arithmetic note is exactly M_1. The local expansion is
\[
G(s)=1-M_1s+\frac{M_2}{2}s^2-\frac{M_3}{6}s^3+\cdots.
\tag{8}
\]
The first derivative is not an arbitrary sensitivity parameter: it is the mean position of the digit-mass distribution. Higher derivatives are higher position moments.

The moments are invariants of the symbolic object under this base-transport construction. The base parameter changes only the exponential tilt e^{-ns}; it does not change the weights w_n or the underlying support positions n.

## 4. The associated infinite-dimensional vector field

Define the tilted moments
\[
M_k(s)=\sum_{n\ge1}n^kw_ne^{-ns}.
\]
Then
\[
\boxed{\frac{d}{ds}M_k(s)=-M_{k+1}(s),\qquad k=0,1,2,\ldots.}
\tag{9}
\]
This is an infinite moment-flow vector field: the derivative of each coordinate is the negative of the next coordinate. The transported value is the zeroth coordinate, G(s)=M_0(s).

Normalize the tilted weights by
\[
w_n(s)=\frac{w_ne^{-ns}}{G(s)}.
\]
They form a probability distribution. Its mean is
\[
\mathbb E_s[N]=\frac{M_1(s)}{M_0(s)}
=-\frac{d}{ds}\log G(s),
\]
and differentiating once more gives
\[
\frac{d^2}{ds^2}\log G(s)
=\operatorname{Var}_s(N)\ge0.
\tag{10}
\]
Thus log G is convex. Equality in (10) would require the tilted distribution to be supported at one position; here at least d_1=3 and d_3=1, so the variance is strictly positive for every admissible s.

This gives a concrete invariant package:
- the fixed support {n:d_n>0};
- the native weights w_n=d_nπ^{-n};
- the full moment sequence (M_k) at s=0;
- the differential hierarchy (9);
- strict log-convexity (10).

The measure μ uniquely determines G, and G on any open interval uniquely determines μ by uniqueness of the Laplace transform. Thus no information is lost if the entire transport curve is retained; information is lost when only one value F(t) is retained.

## 5. What this says—and does not say—about arithmetic migration

We now have a rigorous example of a fixed symbolic object whose value migrates as follows:
\[
\begin{array}{c|c}
\text{base }t&F(t)\\ \hline
\pi&1\text{ (integer)}\\
b\in\mathbb Z,\ b\ge4&\text{irrational in }(0,1)\\
10&\text{irrational decimal value in }(0,1)
\end{array}
\]
The path t↦F(t) is analytic and strictly decreasing, while the arithmetic type of its value can change at special parameter values. The crossing at t=π is exact; irrationality at integer b≥4 follows from nonperiodicity.

For noninteger t, rationality is not classified here. For algebraic t, transcendence of F(t) is not established. For integer b≥4, transcendence of F(b) is also not established. A possible future route is to study the block-complexity function p_d(n), the number of distinct length-n blocks in d. Results of Adamczewski–Bugeaud show strong complexity restrictions on expansions of algebraic irrational numbers in integer bases; to use such a criterion here, one would first need a rigorous upper or lower bound for p_d(n). The present note proves no such bound.

This distinction is important: the word d is non-eventually periodic, which suffices for irrationality in integer bases, but it does not by itself settle transcendence.

## 6. Research tasks opened by the transport model

1. Determine effective bounds for the moments M_k and for the cumulants of μ. In particular, compare M_1 with the constant m already appearing in the return-count asymptotic.
2. Study the block complexity p_d(n) from the π remainder dynamics. Can one prove a complexity estimate strong enough to force transcendence of F(10)?
3. For algebraic or rational noninteger t>1, investigate whether F(t) can be rational, algebraic irrational, or transcendental. The fixed word d is admissible as a digit sequence only when its digits fit the chosen alphabet; the series (2) remains defined regardless.
4. Compare two transports: (i) keep the word d fixed and vary t, as here; (ii) at each t regenerate the greedy expansion of 1 in base t. The second operation changes the word itself. The difference between these two maps may be measured by symbolic distance and the associated change in the digit-mass measure.
5. Use the moment-flow invariants as a bridge back to P2848 only if a quantitative theorem connects them to the integer return counts or to q||qπ||. Complete monotonicity alone does not imply a Diophantine approximation bound.

## 7. Relation to established theory

The use of a fixed digit word and its shift/β-expansion coding belongs to the established theory of β-expansions and symbolic dynamics. Rényi's β-expansions and Parry's symbolic characterization are background, not claimed as new here. See:
- [β-Expansions and symbolic dynamics](https://doi.org/10.1016/0304-3975(89)90038-8), *Theoretical Computer Science* 65 (1989), 131–141.
- [Quasi-periodic β-expansions and cut languages](https://doi.org/10.1016/j.tcs.2018.02.028), for the relation between integer-base rationality and eventual periodicity.
- [Adamczewski–Bugeaud, Dynamics for beta-shifts and Diophantine approximation](https://adamczewski.perso.math.cnrs.fr/AB_ETDS.pdf), for complexity constraints on algebraic numbers' expansions.

The Laplace-transform and moment-flow identities (5)–(10) are direct consequences of the fixed probability measure μ; they do not invoke IUT or any unverified analogy with it.

## 8. Reproducibility note

All displayed identities in this note follow from (1), positivity, and termwise differentiation of exponentially convergent series. The irrationality assertion (4) uses the standard theorem that rational numbers have eventually periodic expansions in an integer base. No numerical experiment is needed for these proofs. No claim of P2848 being proved or disproved is made.

## 9. Extension to the whole word space and carry jets

The companion [continuous-base geometry note](../base-transport-2026-10-31/README.md) transports the entire closed π-greedy word space. It proves the dimension formula logπ/logt, identifies c_n as the number of depth-n gaps, and connects their total length to D·C=1. Further results include the critical gap-depth distribution, generic algebraic independence of t and F(t), and the correspondence between powers of the carry kernel and Taylor contact order. Exact finite certificates and a figure accompany the proofs. None of these statements settles transcendence of F(10) or P2848.
