# π approximation scope, Flint Hills mechanisms, and claim audit

Research snapshot: 2026-11-01. This note separates a complete **metadata scope screen** from selective proof reading. It does not certify every paper in the repository or a Lean build.

## Scope that was actually checked

The pinned repository is `openai/math@4eaee87548250b2682affda7d05d947c05554bae`. The local catalogue contains 721 manuscripts in 721 distinct directories; its inventory has 751 directories. The 30 entries in `uncatalogued.json` form exactly the complement of those 721 directories. All 721 titles, abstracts, and family descriptions were screened, and all 30 extra directory names were screened. A deliberately broad search produced 32 candidate metadata entries. Manual subject classification leaves one direct π approximation manuscript:

- Family 017, [The irrationality exponent of pi is 2](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/paper.pdf), directory `The-irrationality-exponent-of-pi-is-2-September-24-2026`.

No additional revision of this manuscript is identified in the pinned inventory. The scope is **rational approximation to the real constant π, irrationality exponents, and the associated Flint Hills series**. Papers with π only in geometric/physical constants, billiard-angle normalization, product-gate notation `ΣΠΣΠΣ`, or homotopy/fundamental-group notation are outside that subject. The other 31 broad hits and all 30 extra directory names are individually recorded in `scope.json`.

The 30 extra directories were not full-text screened. Thus the defensible claim is that the pinned metadata identifies one direct π approximation manuscript, not that 751 full texts have been independently searched or read. Supporting external papers in the π manuscript’s bibliography are separate from the 721-paper corpus.

## The Lean distinction: selected statement, actual solution, and verification

The official [017 scope document](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/lean/docs/017.md) selects the eventual lower bound for every exponent ν>2 and the exact supremum characterization of μ(π)=2. It expressly excludes the Flint Hills consequence from that **selected comparator statement**.

The linked [comparator target](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/lean/ComparatorChallenges/PiExponent.lean) ends with `by sorry`. That is the challenge statement, and must not be mistaken for the actual solution. Its adjacent [configuration](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/lean/ComparatorChallenges/PiExponent.json) names `OAI.NumberTheory.PiExponent.Main` as the solution module. The actual [Main.lean](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/lean/OAI/NumberTheory/PiExponent/Main.lean) contains proof terms for the selected statement and also a theorem named `flint_hills_summable`.

Main derives the π lower bound from `AdmissibleMatrixInterpolation.globalInterpolation`, and then derives the exponent statement and summability. Its direct imports connect geometry/interpolation to the analytic aggregate. The scope document, target, configuration, solution Main, and two direct imports were read. The transitive import closure was not audited, and no Lean compiler, kernel export, or comparator command was run. Therefore: **solution source exists, but independent executable verification has not been performed here**. The catalogue’s `formalization_source_listed=false` means only that this exact preprint PDF is not directly listed by the snapshot parser in `formalization.yaml`; it does not imply that no solution source exists.

## Why individual-term estimates lose the important information

Write

\[
 F=\sum_{n\ge1}\frac1{n^3\sin^2 n},\qquad
 \alpha=\frac1\pi,\qquad \|x\|=\min_{m\in\mathbb Z}|x-m|.
\]

For every real y,

\[
 2\|y\|\le |\sin(\pi y)|\le\pi\|y\|.
\]

Consequently F converges exactly when

\[
 D(\alpha)=\sum_{n\ge1}\frac1{n^3\|n\alpha\|^2}
\]

converges. The singularity is a distance to the nearest integer; trigonometry only changes fixed factors.

An irrationality bound μ(π)=2 would give, for every ε>0 and all sufficiently large n,

\[
 |\sin n|\gg_\varepsilon n^{-1-\varepsilon},\qquad
 \frac1{n^3\sin^2 n}\ll_\varepsilon n^{-1+2\varepsilon}.
\]

This upper bound is not summable. It treats every integer as though it were an exceptionally good approximation. That is why termwise comparison alone fails to extract Flint Hills convergence even from μ(π)=2. A series can have terms tending to zero and still diverge.

The missing ingredient is **spacing**. Two extremely accurate approximations with nearby denominators would subtract to give an even more accurate approximation with the smaller difference denominator. This is the mechanism in Meiburg’s density argument. The continued-fraction formulation below makes the same obstruction more concrete.

## A reconstruction of the spacing argument

Let p_k/q_k be the convergents of α and put β_k=||q_k α||. After discarding finitely many initial indices, the standard continued-fraction facts are

\[
 \frac1{q_{k+1}+q_k}<\beta_k<\frac1{q_{k+1}},\qquad
 0<m<q_{k+1}\Longrightarrow\|m\alpha\|\ge\beta_k.
\]

Consider a block q_k≤n<q_{k+1}. Any two distinct indices differ by a positive integer smaller than q_{k+1}; therefore the corresponding points nα modulo 1 are β_k-separated on the circle. If their distances to zero are ordered as x_1≤x_2≤…, then x_1≥β_k and x_j≥(j−1)β_k/2 for j≥2. Hence

\[
 \sum_{q_k\le n<q_{k+1}}\|n\alpha\|^{-2}
 \le (1+4\zeta(2))\beta_k^{-2}.
\]

The single index n=q_k supplies the opposite bound. Multiplying by the weight n^-3 gives

\[
 \beta_k^{-2}q_k^{-3}
 \le \sum_{q_k\le n<q_{k+1}}
       \frac1{n^3\|n\alpha\|^2}
 \le (1+4\zeta(2))\beta_k^{-2}q_k^{-3}.
\]

Thus the largest spike controls the entire block up to an absolute factor. This is a statement about both height and frequency, not just height. Using

\[
 q_{k+1}=a_{k+1}q_k+q_{k-1},
\]

we obtain the exact convergence criterion

\[
 \boxed{\quad
 F<\infty
 \iff\sum_k\frac{q_{k+1}^2}{q_k^3}<\infty
 \iff\sum_k\frac{a_{k+1}^2}{q_k}<\infty.
 \quad}
\]

This reconstruction checks the core block argument in López Zapata v5 and the comparable reduction in Dan v1. The constants are unimportant for convergence; the structural fact is that the cost of block k is approximately a_(k+1)^2/q_k.

## The 5/2 threshold and why equality is different

The continued-fraction identity

\[
 \mu(\alpha)=1+\limsup_k\frac{\log q_{k+1}}{\log q_k}
 =2+\limsup_k\frac{\log a_{k+1}}{\log q_k}
\]

and μ(1/π)=μ(π) yield:

- If μ(π)<5/2, choose δ>0 so that eventually q_(k+1)^2/q_k^3≤q_k^-δ. Denominators grow at least geometrically, so the sum over blocks converges.
- If μ(π)>5/2, infinitely many blocks have q_(k+1)^2/q_k^3 growing without bound. The summands at n=q_k fail to tend to zero.
- At μ=5/2, the exponent records only the power scale. It does not record logarithmic or index-dependent corrections, which can decide summability.

For an explicit critical family, set q_0=1, a_1=2, and recursively

\[
 a_{k+1}=\left\lceil\frac{\sqrt{q_k}}{(k+1)^{r/2}}\right\rceil,
 \qquad q_{k+1}=a_{k+1}q_k+q_{k-1}.
\]

The quantity inside the ceiling tends to infinity because q_k grows at least exponentially. Therefore q_(k+1)^2/q_k^3∼(k+1)^-r, while log q_(k+1)/log q_k→3/2. Every resulting irrational has μ=5/2. For r>1 its generalized Flint Hills series converges. For 0<r≤1 the series diverges even though its individual summands tend to zero. For r=0 the convergent-denominator summands do not tend to zero. This is Dan’s Proposition 5.1 specialized to (s,t)=(3,2), and its proof follows directly from the displayed recursion.

Nothing here proves that the actual continued fraction of 1/π has any particular tail. Every finite prefix can be extended into either kind of critical behavior. A long numerical computation, or an almost-everywhere theorem, does not settle membership of this one prescribed number.

## Primary-paper reading and status

| Paper | Claim-relevant proof reading | What it contributes |
|---|---|---|
| [Alekseyev 2011, v1](https://arxiv.org/html/1104.5100v1) | Full 5-page paper; Theorem 2 pp2–3, Corollary 4 pp3–4, Theorem 5 p4 | Convergence forces μ(π)≤5/2. Its termwise sufficient condition is weaker and does not solve the classical series. |
| [Meiburg 2022, v1](https://arxiv.org/html/2208.13356v1) | Full 12-page paper; density proof pp7–9, refined partition argument pp9–10, critical divergent construction pp10–11 | Strict sufficient condition μ(π)<5/2, obtained by controlling how often large terms occur. Equality is not covered. |
| [Dan 2026, v1](https://arxiv.org/html/2609.12338v1) | Sections2–5 and beginning6; reduction pp3–5, dimension p6, critical construction pp7–8. Later gauge proofs not fully audited. | General weighted denominator criterion, and critical examples of both behaviors. Its own text leaves the classical π point undecided. |
| [López Zapata 2026, v5](https://arxiv.org/html/2603.09719v5) | Core Sections1–3, pp1–4; later result statements checked | Exact partial-quotient criterion. Introductory attribution of strict sufficiency to Alekseyev is inaccurate; that improvement is Meiburg’s. |
| [Carella 2022, v10](https://arxiv.org/html/1902.08817v10) | Proof routes in Sections3–8, pp4–14; not a complete audit of all later material | Claims μ(π)=2, but the displayed proof routes have explicit obstructions below. It is not used as an established upper bound. |
| [Mantzakouras–López Zapata 2025, withdrawn v2](https://arxiv.org/abs/2502.03474v2) | Official withdrawal notice; original v1 selectively read at pp2,13–16,24–27,38–42 of its 70-page PDF | Withdrawn on 2026-09-14. The notice disclaims the convergence proof and μ(π)≤2.5. It is historical claim material, not a bound. |

The older sources’ phrases “best known” describe their dates, not the current accepted state. In particular, neither Carella’s assertion, the withdrawn claim, nor the new 017 manuscript becomes an accepted bound merely by appearing in a bibliography. Current-bound assessment belongs to the separate primary-source audit of the rational-integral papers.

## Carella: concrete sign and cancellation obstructions

In the first proof, PDF p5 Eq(8) inserts a lower bound of order 1/q_n for the sine error, although the preceding facts provide an upper bound of that order and comparability between the sine and the approximation error. Two upper bounds do not produce that lower bound.

The bounded-partial-quotient proof on PDF p13 is more directly checkable. It defines

\[
 e=p_{n+1}-\pi q_{n+1},\qquad C=e-1/q_n,\qquad D=1/q_n.
\]

For convergents, |e|<1/q_(n+1)<1/q_n. Therefore C<0. For all sufficiently large n, C is small and negative, so **sin C<0**. Eq(55) instead gives sin C a positive lower bound asymptotic to 1/(2q_n). The sign is wrong. Moreover C+D=e, so

\[
 \cos C\sin D+\cos D\sin C=\sin e.
\]

The two pieces can almost cancel. The reverse triangle inequality does not justify replacing both magnitudes by separate lower bounds inside their difference. In general A≥a and B≥b do not imply |A−B|≥|a−b|. This explains precisely where the purported lower bound of order 1/q_n disappears.

A different route, Lemma5.2 on PDF p10, treats near-maximal endpoint values of a continuous oscillatory function as a lower bound throughout its interval. Continuity supplies no such conclusion. These are failures of specific displayed arguments; the present note does not claim a complete audit of every later statement in the 35-page PDF.

## Withdrawn paper: why a change of representation does not remove the arithmetic

The original v1, PDF p41, claims that floor(x) is Hölder continuous with a positive exponent. Take x=m−ε and y=m: the jump is 1 while |x−y|^β=ε^β tends to zero. No fixed Hölder constant can make the assertion true for any β>0. The Young-integral argument therefore lacks a stated hypothesis. Even a finite-interval existence result would not itself establish a finite improper integral at infinity.

The Bessel pathway, PDF p27, also does not give the claimed decay. On the imaginary axis the half-order Bessel expression retains the oscillation of sin n; the second exponential cannot simply be discarded. Even in the displayed one-exponential expressions, e^(-3in)/(e^(-in))^3=1, so there is no exponential decay. The exact half-order identity turns the relevant ratio back into a constant times n·sin(3n)/sin(n)^3. The triple-angle identity then recovers the same cosecant-square singularity. Thus this representation alone does not bound the exceptional spikes. The official v2 withdrawal is consistent with these concrete obstructions.

## Coverage limits

This note provides a complete pinned-metadata scope screen and a proof-focused reading of the cited Flint Hills pathways. It does not claim full reading of the 721-paper corpus, of the 30 uncatalogued directory bodies, of all 70 pages of the withdrawn original, of Dan’s later Hausdorff gauge arguments, or of Lean’s transitive source closure. No existing full-reading ledger was modified.
