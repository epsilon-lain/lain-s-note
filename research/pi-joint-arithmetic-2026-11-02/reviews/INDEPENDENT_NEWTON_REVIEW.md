# Independent audit of the unique near-period root

2026-11-02. Reviewed input:
`imported-next-work/root_geometry/UNIQUE_NEAR_PERIOD_ROOT_LATEX.md`.
The supplied internal review was not used. The input files were not changed.

## Verdict

**The claimed theorem is valid for all sufficiently large n**, conditional
on the previously established period-value and logarithmic-derivative
theorems. The new all-order mixed-endpoint bound supplies the higher-order
control that the earlier first-derivative argument lacked. It proves one
simple nearest root on the imaginary axis, its side of approach, and the
stated relative Newton error `O(n^(-2))`.

There is no substantive correction required. Two harmless clarifications are
to state `0<=k<n` beside the recurrence and to define `K_1=1/t` by the same
factorial formula. The uniform local estimate is used only for `n>=2`.

In fact, retaining one factor from the exact recurrence improves the relative
Newton error to **`O(n^(-3))`**; this strengthening is proved below separately.
No rational-approximation conclusion follows.

## 1. Mixed coefficient normalization and recurrence

For `m=n-k`, expansion of the determinant as a Vandermonde-square integral
gives the mixed coefficient the factor `1/(m! k!)`. After taking absolute
values, the ordinary-variable integral contributes
`m! Delta_m exp(m max(Re Z,0))`. The endpoint-variable integral, after the
shift `y_a=v_a+|Z|`, is bounded by

\[
e^{k|Z|}\int_{[0,\infty)^k}e^{-\sum y_a}
\prod_a y_a^{2m}\prod_{a<b}(y_b-y_a)^2\,dy.
\]

Its value is `k! A_k product_(j=0)^(k-1) (2m+j)!`. The powers of `|Z|`
from cross differences, endpoint differences, and endpoint measures total
`2mk+k(k-1)+k=k(2n-k)`. Thus there is no missing factorial or binomial
coefficient in the claimed majorant

\[
|C_{n,k}(Z)|\le
6e^{m\max(\Re Z,0)+k|Z|}B_{n,k}(|Z|).
\]

The displayed monic Rodrigues polynomials have squared norms
`j!(j+2m)!`: integration by parts `j` times proves both orthogonality and
the norm. Boundary terms vanish also at `m=0`. This verifies the Laguerre
integral, including the pure-endpoint case `k=n`.

Direct factorial cancellation, for `0<=k<n` and hence `m>=1`, gives

\[
\frac{B_{n,k+1}(t)}{B_{n,k}(t)}
=\frac{k!K_{n-k}}{(2n-2k)_k}.
\]

At `k=0` the quotient is `K_n`; at `k=n-1` it is `K_1/n`.
Also `B_(n,n)=A_n^2/(Delta_n t^(n^2))`, the correct final-endpoint scale.
Since `(2m)_k>=k!`, the geometric envelope follows once
`max_(1<=j<n) K_j=K_(n-1)`. This holds already for `n>=4`: the displayed
successive ratio is greater than one for `j>=2`, and
`K_3/K_1=4320/t^4>4320/7^4>1`.

## 2. Uniformity in Z and in the number of insertions

On `|Z-omega|<=1/n`, with `n>=2`, one has
`log(t/|Z|)<=1/(5n)`. The three exponential costs in the coefficient bound
are at most

\[
1+kt+k/n+2k/5\le1+k(t+1).
\]

Consequently the bound `|C_(n,k)(Z)|<=6e c^k B_(n,k)(t)`,
`c=exp(t+1)`, is uniform for **every** `0<=k<=n`; its constant is not
silently growing with `n` or `k`.

The coefficients are analytic on this disk because it excludes zero.
Cauchy's derivative estimate on its inner half gives
`|C_(n,1)(omega+w)-T_n|<=12ec nK_n |w|`. The ordinary integral gives
`|C_(n,0)(omega+w)-1|<=6n|w|exp(n|w|)`. Their denominator is the fixed
nonzero value `det H_n(omega)`, so these arguments do not assume that
`det H_n(Z)` remains nonzero elsewhere.

For `|w|<=L/K_n`, the exact substitution is
`eta=1-exp(omega+w)=1-exp(w)`, with
`|eta|<=2|w|` and `|eta+w|<=|w|^2` eventually. In particular the linear
term has sign `-wT_n`. Put `rho_n=K_(n-1)/K_n=O(n^(-2))`.
The geometric-series ratio is at most `2cL rho_n`, eventually less than
one half. This controls the full tail through `k=n`, giving the stated
`O_L(n^(-2))` remainder. The ordinary-coefficient variation, linear-coefficient
variation, and error in replacing `eta` by `-w` contribute only
`O_L(n/K_n)`. Thus the normal form is valid uniformly on the entire disk.

## 3. Rouché, imaginary location, and Newton displacement

On `|w|=48/K_n`, the bound `|T_n|>=K_n/12` gives
`|1-wT_n|>=3`. The uniform remainder tends to zero, so Rouché gives exactly
one zero, counted with multiplicity. The nonzero prefactor
`Z^(n^2)/A_n^2` does not change this zero count: the disk excludes `Z=0`
eventually. There are no zeros on its boundary.

At the zero, `wT_n=1+O(n^(-2))`, so `|w|<=24/K_n` eventually.
Every other root is outside the radius-`48/K_n` disk. Thus this root is
simple and uniquely nearest. The symmetry `Z -> -conjugate(Z)` preserves
this disk and fixes its center. Uniqueness forces the root onto the imaginary
axis; symmetry by itself would not have sufficed.

The exact chain rule gives `P_n'/P_n=-T_n+O(n^2)`. Since
`n^2/K_n` is smaller than every fixed inverse power of `n`, replacing
`T_n^(-1)` by the exact Newton displacement `-P_n/P_n'` preserves the
relative `O(n^(-2))` error. The known phase of `T_n` then places the root
below `omega` for even `n` and above it for odd `n`. The asserted logarithmic
distance follows from the established Stirling expansion.

At fixed `Z=omega`, the same coefficient estimates give one simple endpoint
root with `|Y_n-1|<=24/K_n`. Its reciprocal-conjugate image has distance
`|Y_n-1|/|Y_n|<48/K_n` eventually. The endpoint symmetry and uniqueness
therefore imply `|Y_n|=1`. This is a separate zero set from the period roots.

## 4. A valid sharpening available from the same proof

The exact recurrence at `k=1` retains an extra factor:

\[
B_{n,2}(t)=\frac{K_nK_{n-1}}{2n-2}.
\]

For `k>=2`, the remaining ratios are still at most `K_(n-1)`. Therefore
the whole higher-order tail satisfies the stronger bound

\[
\sum_{k=2}^n|C_{n,k}(\omega+w)\eta^k|
\le\frac{6ec^2B_{n,2}(t)|\eta|^2}
 {1-cK_{n-1}|\eta|}
=O_L(n^{-3}).
\]

The other errors remain `O_L(n/K_n)=o(n^(-3))`. Repeating the same Rouché
and chain-rule argument consequently proves

\[
\boxed{\zeta_n-\omega
=-\frac{P_n(\omega)}{P_n'(\omega)}(1+O(n^{-3})).}
\]

The corresponding endpoint displacement has relative `O(n^(-3))` error as
well. This is an improvement of the proved bound, not a claim that the
exponent three is optimal.

## Exact checks executed

The supplied verifier's exact routines were imported with bytecode writing
disabled and rerun without executing or relying on its numerical root
examples. All checks passed: 1377 coefficient recurrences, 81 endpoint
identities, 1377 geometric bounds, 49 Laguerre determinants, and 27 synthetic
Gaussian-rational mixed coefficient bounds. Results are recorded in
`reviews/newton_exact_check_results.json`. These finite checks audit the
algebra; the uniform all-size argument above establishes the theorem.

For the separate `O(n^(-3))` corollary, additional exact checks passed all
81 identities `B_(n,2)=K_n K_(n-1)/(2n-2)` and 1296 strengthened geometric
envelopes, for `n=4,...,30` and the rational values `t=6,44/7,7`.

To reproduce both exact suites without the optional numerical examples, run
`python -B reviews/verify_newton_exact.py`. The portable verifier finds
the staged `imported-next-work/root_geometry/verify_unique_root.py`, with a
local-input fallback; `--source PATH` can select a source file or directory.
It writes only its own JSON report. Local-default, explicit-file,
explicit-directory, and staged-default invocation modes were checked, with
the input file contents unchanged and no bytecode cache created there.
