# Rational-base displacement of the π-unit word and continued-fraction quality

Date: 2026-10-31

This note studies the user's proposed observable: take the fixed π-unit digit word
\[
1=\sum_{n\ge1}d_n\pi^{-n},\qquad d=(3,0,1,1,0,2,\ldots),
\]
and evaluate the *same word* in a rational base r:
\[
F(r)=\sum_{n\ge1}d_nr^{-n}.
\]
The central quantity is its displacement from 1:
\[
\Delta(r)=F(r)-1.
\]
This is not the ordinary expansion of the real number 1 in base r; it is the transport of a chosen symbolic representation.

## 1. Exact factorization of the displacement

Put
\[
w_n=d_n\pi^{-n},\qquad \sum_{n\ge1}w_n=1,\qquad m=\sum_{n\ge1}n w_n.
\]
For any r>1, let x=\pi/r. Since
\[
x^n-1=(x-1)(1+x+\cdots+x^{n-1}),
\]
we obtain the exact identity
\[
\boxed{
F(r)-1=\frac{\pi-r}{r}
\sum_{n\ge1}w_n\sum_{j=0}^{n-1}\left(\frac{\pi}{r}\right)^j.
}
\tag{1}
\]
The series converges for every r>1: its terms are bounded by 3n·max(π⁻¹,r⁻¹)^n, since d_n≤3. This bound also covers r>π, where a bound by a fixed constant times n r⁻ⁿ need not hold.

Define the migration secant factor
\[
\boxed{
\mathcal R(r)=\frac{F(r)-1}{\pi-r}
=\frac1r\sum_{n\ge1}w_n\sum_{j=0}^{n-1}(\pi/r)^j.
}
\tag{2}
\]
It is positive for all r>1, including on both sides of π. At the native base, dominated/local uniform convergence gives
\[
\boxed{\mathcal R(\pi)=\frac m\pi.}
\tag{3}
\]
Thus
\[
F(r)-1=\left(\frac m\pi+o(1)\right)(\pi-r)
\qquad(r\to\pi).
\tag{4}
\]
This is stronger than a first-order Taylor slogan: (1) gives an explicit positive weighted-geometric representation of the secant factor.

For r≥π, 1≤Σ_{j=0}^{n-1}(π/r)^j≤n, hence
\[
\boxed{\frac1r\le\mathcal R(r)\le\frac mr.}
\tag{5}
\]
On any compact neighborhood of π contained in (1,∞), continuity and positivity give constants 0<c<C<∞ with c≤R(r)≤C.

## 2. The key Diophantine consequence

Let p_k/q_k be any sequence of rational approximants to π with q_k→∞ and p_k/q_k→π. Then (4) implies
\[
\boxed{
q_k^2\left|F(p_k/q_k)-1\right|
=
\left(\frac m\pi+o(1)\right)
q_k^2\left|\pi-\frac{p_k}{q_k}\right|.
}
\tag{6}
\]
Therefore
\[
\boxed{
\liminf_k q_k^2|F(p_k/q_k)-1|=0
\iff
\liminf_k q_k^2|\pi-p_k/q_k|=0.
}
\tag{7}
\]
For the full denominator sequence, the standard continued-fraction theorem says the latter liminf is zero exactly when the partial quotients of π are unbounded. Equivalently, π is badly approximable exactly when its partial quotients are bounded. See the cited standard references below.

This is the central answer to the proposed idea:
- the rational-base displacement of the transported π-unit word is a faithful local gauge of π's approximation quality;
- it preserves the exponent and the zero/nonzero liminf classification;
- it does **not** create an independent Diophantine invariant or make the hard question disappear. Locally it is a smooth positive rescaling of |r−π|.

For convergents, write the complete quotient after p_k/q_k as α_{k+1}=[a_{k+1};a_{k+2},...]. The exact error formula gives
\[
q_k^2\left|\pi-\frac{p_k}{q_k}\right|
=\frac1{\alpha_{k+1}+q_{k-1}/q_k}.
\]
Consequently
\[
\boxed{
q_k^2|F(p_k/q_k)-1|
=
\frac{m/\pi+o(1)}
{\alpha_{k+1}+q_{k-1}/q_k}.
}
\tag{8}
\]
Large next partial quotients create small transported displacements. The displacement is a weighted, smooth encoding of the same continued-fraction gene.

## 3. The concrete case r=22/7

Since 22/7>π, equation (1) gives F(22/7)<1. Using the certified π-digit prefix beginning
\[
30110211100202211300010200021022221221202012120200
\]
and the universal tail bound d_n≤3, the first 50 digits give
\[
\sum_{n=1}^{50}d_n(7/22)^n
\le F(22/7)
\le
\sum_{n=1}^{50}d_n(7/22)^n+
\frac{3(7/22)^{51}}{1-7/22}.
\]
The tail width is less than \(1.906\times10^{-25}\). Numerically this yields
\[
F(22/7)\approx0.9995539515839704267337230261,
\]
and
\[
\boxed{
1-F(22/7)\approx0.0004460484160295732662769739.
}
\]
The normalized comparison is
\[
\frac{1-F(22/7)}{22/7-\pi}
\approx0.35274986316372215,
\]
while the limiting secant factor at π is
\[
\frac m\pi\approx0.35291755888141614.
\]
The difference between these two slopes is expected: 22/7 is close to, but not equal to, π. The finite-prefix plus geometric-tail formula above is the certification method; the displayed decimals are rounded.

Since 22/7 is the first nontrivial convergent [3;7] of π, its transported displacement is a concrete first sample of the same approximation quality. Multiplying by \(7^2\) gives
\[
49\,|F(22/7)-1|\approx0.02185637238544909.
\]
This number alone does not diagnose the infinite liminf; the sequence of convergents is essential.

## 4. Why the rational base matters—and what it does not imply

For r=p/q, the value can be written as
\[
F(p/q)=\sum_{n\ge1}d_n(q/p)^n.
\]
The rationality of the base makes each finite truncation rational and gives exact integer numerator recurrences after clearing powers of p. The infinite limit, however, need not be rational merely because the base is rational. Eventual periodicity is a rationality criterion for *integer* bases with the standard digit alphabet; it should not be applied without proof to a noninteger rational base p/q.

The displacement identity is valid for every real r>1, so rationality is not needed for the analytic comparison. Rational bases are special because they provide a countable, arithmetically structured sampling set.

A possible further arithmetic question is whether the sequence of rational truncations at p/q has modular or gcd structure that can be combined with the factor (1). This would be genuinely additional information; the secant identity alone is only a smooth rescaling.

## 5. Generalization to any irrational native base β

Nothing in (1) uses special properties of π beyond the existence of a bounded nonnegative unit digit word
\[
1=\sum_{n\ge1}d_n^{(\beta)}\beta^{-n}.
\]
For any β>1 with such a word, define
\[
F_\beta(r)=\sum_{n\ge1}d_n^{(\beta)}r^{-n},\qquad
m_\beta=\sum_{n\ge1}n d_n^{(\beta)}\beta^{-n}.
\]
Then
\[
\boxed{
\frac{F_\beta(r)-1}{\beta-r}
\longrightarrow\frac{m_\beta}{\beta}
\quad(r\to\beta).
}
\tag{9}
\]
Thus every native base carries its own positive local conversion factor between rational-base displacement and ordinary rational approximation error.

Comparing these constants \(m_\beta/\beta\) across β may reveal differences in the positional geometry of the native words. But the zero-liminf property for rational approximants to β is preserved, not altered, by this local transport.

## 6. Research status

Established:
1. exact weighted-geometric factorization (1);
2. positive secant factor and native slope m/π;
3. equivalence of normalized displacement liminf and the ordinary Diophantine liminf;
4. exact relation to the next partial quotient for convergents;
5. a finite-prefix numerical enclosure for the 22/7 example.

Not established:
- whether π's partial quotients are unbounded;
- whether P2848 is true or false;
- whether F(22/7) or F(10) is transcendental;
- any new independent lower or upper bound for π's irrationality measure.

The main conceptual result is a transport dictionary: a rational approximation p/q to π is sent to a rational-base displacement F(p/q)-1, with asymptotically constant scale m/π. This is a useful new observable, but not yet a stronger approximation theorem.

## 7. References

- [Continued fraction convergent error bounds](https://crypto.stanford.edu/pbc/notes/contfrac/converge.html): exact formula and bounds involving consecutive convergent denominators.
- [Badly approximable numbers and bounded partial quotients](https://androma.org/theorems/5518): standard equivalence between bounded partial quotients and a uniform \(q^{-2}\) lower bound.
- [MIT continued-fractions notes](https://math.mit.edu/~rmd/IAP/continuedfractions.pdf): standard continued-fraction approximation background.

No proof or disproof of P2848 is claimed.

## 8. Integer carry continuation

[Rational carry escape](RATIONAL_CARRY_ESCAPE.md) identifies the secant factor with the native-remainder series, proves exact denominator and two-sided modular identities, and translates approximation quality into excess integer survival depth. It gives a finite integer-strip certificate equivalent to P2848 in the unbounded-depth limit, together with a non-eventually-periodic √10 control whose excess depth is uniformly bounded. No unboundedness theorem for π is claimed.
