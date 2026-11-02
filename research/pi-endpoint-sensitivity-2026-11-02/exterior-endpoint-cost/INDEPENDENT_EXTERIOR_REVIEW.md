# Independent review of the singular profile and exterior cost

Reviewed 2026-11-02 by a separate research agent. Scope: `EXTERIOR_SINGULAR_PROFILE.md`, including its stronger endpoint-root count. Verdict: the uniform singular-value estimate, exterior estimate, and root-count argument are correct, conditional on the cited two-sided actual-period estimate. No unresolved mathematical blocker was found.

## Exact factorization and constants

The two distinct Pascal factorizations both hold:

\[
V=PUP^T,\qquad [(r+s)!]=D_fPP^TD_f.
\]

The first is a moment translation; the second is Vandermonde's binomial identity. At imaginary `Z=it`, the phase factorization is `U=(it)^{-1}QKQ` with `Q=diag(i^r)` and positive definite `K=(D_qP)(D_qP)^T`. A conjugate transpose is not needed in the phase identity: left and right multiplication by unitary matrices both preserve singular values. Positivity of the real `K` then identifies its singular values with squared singular values of `D_qP`.

The Frobenius bounds give `||P||,||P^{-1}||<2^n`. Squaring the comparison for `D_qP` costs `2^(2n)`, and applying the separate congruence `V=PUP^T` costs another `2^(2n)`. Neither is omitted. The Hilbert normalization costs at most `n` on the lower side and `2n^2 2^(8n)` on the upper side. Consequently equation (9) has the correct constants:

\[
\frac{a_j^2}{t n 2^{4n}}
\le \sigma_j(SVS)
\le \frac{2n^2 2^{12n}}{t}a_j^2.
\]

This uses the individual-index min-max inequality, so it holds for every `j`, including the smallest singular value.

## Rearrangement and uniformity

The sequence `q_r=r!/t^r` is not monotone in its initial range. The note correctly avoids replacing its sorted entries directly by `q_(n-j)`. Its elementary bounds

\[
(3t)^{-n}n^r\le q_r\le n^r
\]

hold for all `0<=r<n` at the fixed `t=2*pi>1`: use `n^r/r!<=e^n<3^n`. Order statistics preserve coordinatewise inequalities, which proves the sorted bounds in equation (11). There is therefore no missing small-index or crossover case. Combining these bounds with equation (9) gives a single constant in the `O(n)` error, independent of `j`.

Summing the top `k` singular logarithms gives the exact leading sum `k(2n-k-1) log n` and error `O(nk)`, uniformly for `1<=k<=n`. At `k=n`, the difference between `n(n-1) log n` and the exact determinant's `n^2 log n` is only `n log n`, which is within the asserted `O(n^2)` error. There is no determinant inconsistency at the full exterior power.

The comparison from `B=SVS` to `C=(SHS)^{-1}B` is valid at every singular-value index, since `||SHS||<=1` and its inverse norm is less than six. Singular values of the unnormalized similar matrix are not implicitly substituted for these values.

## Coefficients and actual eigenvalues

The exterior trace identity and trace-norm inequality prove the determinant product envelope in equation (18). They give coefficient upper bounds. They do not give coefficient lower bounds from the corresponding exterior norm lower bounds; the note correctly retains this cancellation warning.

For the root count, the product of the largest `K` eigenvalue moduli is bounded by `||wedge^K C||`, without normality or diagonalizability assumptions. Combining that bound with the exact total determinant and bounding the other `n-K` moduli by `exp(alpha n log n)` gives, for `y=(n-K)/n`,

\[
y(y-\alpha)\le C/\log n.
\]

For fixed positive `alpha`, this implies `y<=alpha+O_alpha(1/log n)` and hence

\[
K\ge(1-\alpha)n-O_\alpha(n/\log n).
\]

The `K=0` case is covered by the empty exterior product. Multiplicities are algebraic. This count does not claim that individual eigenvalues share the singular-value profile. Under `Y=1-1/mu`, it yields the stated count of nearby endpoint roots.

## Corrections incorporated and verification

During review, two presentation corrections were incorporated:

- The unnormalized moment pencil was renamed `mathcal H_n`, and its scaling to the companion's `G_n` was stated, avoiding a conflicting definition of `G_n`.
- The comparison with the sharper single-norm baseline now gives a saving of `k(k-1) log n`; at `k=1` there is no claimed saving.

I reran `python verify_exterior_profile.py`; all its exact Gaussian-rational checks passed. Separately, 170 exact symbolic checks shared with the pencil audit tested sizes `n=2,...,6` at rational and imaginary-rational parameters. Relevant checks here included the factorial-Hankel factorization, endpoint translation congruence, phase factorization, and endpoint determinants. The uniform estimates and root count were reviewed analytically, not inferred from numerical singular values or those finite tests.

The nearby roots remain roots in independent `Y` at fixed `Z=omega`. The note makes no unsupported inference to period-variable roots, irrationality measures, bad approximability of pi, or the unaddressed critical high-dimensional regime.
