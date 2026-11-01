# Batch 06 reading checkpoint

Pinned repository: `openai/math`, commit `4eaee87548250b2682affda7d05d947c05554bae`.

This is source-body reading, not independent proof verification. Downloaded source is not counted as read. Recursive directory trees were fetched without relevance exclusions, README and main TeX were inspected, and each input dependency was resolved. Files and exact numeric reading ranges are recorded in `records.json`. External papers and Lean scope documentation were not inspected in this batch.

## Position 16 — family 333, smooth isometric immersion

**Full body read.** Main TeX, all nine body inputs, the nested geometry figure, and bibliography were read, with no omitted appendix. Total text closure including README: 3,168 lines.

The claim is a smooth isometric immersion into R^4 for every closed smooth Riemannian surface, without orientability or curvature restrictions. The proof prepares a short spherical immersion with one global negative radial normal, decomposes the metric deficit into finitely many positive primitive additions, realizes each by an oscillatory circle, and preserves strict geometric conditions needed for the remaining additions.

The analytical separation is useful: finite mode solves give arbitrary finite accuracy; a perturbation exponent is fixed independently of that accuracy; smoothing then converts sufficient finite accuracy into an exact smooth correction. The iteration smooths only its inputs and adds the resulting increment to the unsmoothed map. Its exact quadratic error identity keeps higher derivatives linear, and the block schedule waits finitely long before controlling the next derivative order. The large primitive construction and small C^2 correction have distinct roles. The exact Gauss equation recovers an unscaled component that the large oscillation alone does not control.

Candidate interface to the pi project: fix overhead exponents before choosing cancellation order, and carry a rank witness plus explicit occupancy/cost invariants through a finite schedule. The bounded coloring of many phase cells suggests avoiding weights that grow with total family size when only boundedly many supports interact. Exact premises and required arithmetic bridges appear in the JSON record. None supplies symbolic rank over Q(Z), an admissible integer specialization, or the same determinant's lower and upper bounds.

## Position 17 — family 362, relativistic Vlasov–Maxwell

**Full body read.** All ten body sections, three nested figure/caption sources, bibliography, main and README were read. Total text closure: 4,019 lines.

The claim concerns arbitrary smooth one-species data in three spatial and three momentum dimensions: compactly supported nonnegative particle density, bounded smooth finite-energy fields and the two divergence constraints. Its central estimate is on signed momentum increments, with the norm taken after integration. This leads to a divergent sum of lower bounds for momentum-doubling times, then bounded-momentum continuation yields global smoothness.

The proof first establishes baseline occupation and direct estimates under the signed bootstrap, then derives an exact source-trajectory cancellation including the transport term. First use of that cancellation controls direction changes from baseline information; the resulting improved occupation is used only afterward in final selection. This order avoids a circular use of the improved estimate.

Particularly relevant methodological details:

- A single energy budget indexed by source momentum and receiver angle is shared by all source-angle and distance bins. Bounded overlap is proved before summation, and the budget remains inside Cauchy–Schwarz instead of being charged once per bin.
- Near/far estimates are balanced by dyadic geometric sums, including explicit exceptions. Expensive selected bins carry a summable weight.
- If a selected receiver coefficient were too large, exact rational exponent inequalities force the bin into the improved-occupation window, where they become inconsistent. The stronger estimate is required only in that window.
- Source and geometric cutoff errors are transferred to unselected neighbors using the exact partition identity after factoring a common multiplier. Time endpoints are restored using signed increments, preserving the gain absent from absolute-force control.

For the determinant project these are candidate methods for one shared occupancy budget, a cheap/expensive signature split and exact cancellation before taking absolute values. They do not transfer the PDE conservation law or timelike geometry. A project-specific resource budget, exact determinant identity, nonzero symbolic minor and integral normalization are still required.

## Position 18 — family 002, exact BSD from low Selmer corank

**Partial body read.** README, main, and complete sections 01–06 have been read. This covers 4,355 body lines, or 4,458 lines including main and README, out of a 7,523-line text closure. All sources are cached and their Git blob SHA values match the pinned recursive tree. Sections 07–10 and references.tex remain unread by this agent. Next is `build/sections/07-theta.tex`, lines 1–1082. Remaining files and exact ranges are:

| File | Unread original-source lines |
|---|---:|
| build/sections/07-theta.tex | 1–1082 |
| build/sections/08-extraction.tex | 1–1033 |
| build/sections/09-unit.tex | 1–449 |
| build/sections/10-center.tex | 1–271 |
| build/references.tex | 1–230 |

The stated theorem is full BSD for every E/Q whose full q-power Selmer group has corank zero or one at some prime q. The read introduction explicitly separates companion rank/finiteness and two-primary inputs from its new odd-primary comparisons. Read section 02 retains exact periods, full point-lattice indices, torsion and Manin factors before taking valuations.

Read section 03 constructs finite free models retaining cup maps, homotopies, duality and derived base change before inverting p; determinant switches retain exact volume information. Integral character division requires bounded series and a residual nonvanishing hypothesis. A central Gaussian pivot formula need not survive specialization: a suitable surviving minor is used instead.

Read section 04 constructs first-order lift obstructions. The available tests span a matrix space; the text does not assert every matrix in that span is a single Galois evaluation. A linear combination with zero left kernel gives a full-rank first-order differential, hence a nonzero initial term of an actual maximal minor. Alternating parity is separately established where invertibility of an alternating matrix is needed.

Read section 05 constructs integral Siegel classes, calibrates their exact period coordinates, and proves a normalized determinant coordinate integral by line switches and a separate removal of smoothing/Euler divisors. Integrality of the smoothed class alone does not permit division by its smoothing factor. In rank one, a nonzero tame-height direction, a cup identity and an exact spectral-height identity determine the central logarithmic coordinate. Integral Poitou–Tate torsion lengths and local Haar factors cancel omitted Euler factors and yield `X_p(E)=v_p(U(0))>=0` for non-CM E. The rank-one comparison is rationally exact up to a common orientation sign, not merely a unit-valued proportionality.

Read section 06 constructs integral depleted logarithmic primitives on the ordinary modular source, without assuming ordinary reduction of the target E. Bounded CM disk measures provide exact finite-character and central interpolation. The analytic ideals descend, and vertical primitivity is proved separately via tame-stalk independence, formal subtori and monodromy. At high characters, local Kummer lines and horizontal square switches preserve determinant valuations and prove divisibility after inverting p. A separate residual determinant condition and integral Weierstrass division then yield `U_w=B_w/L_w` in the integral ring with no factor p. On removing added tame variables the numerator and denominator have the same nonzero squared factor, canceled after an exact marked determinant comparison; that factor need not be a central unit. Horizontal quotient divisors remain and require the unread sections 07–09.

These are candidates for a symbolic row-rank certificate over Q(Z), with explicit formal variables and a separate controlled integral specialization. They also warn against replacing integral nonvanishing by generic rational nonvanishing or assuming that a span statement produces one admissible configuration. Sections 07–10 contain the theta, extraction, unit and center arguments; the JSON architecture labels their description as the roadmap stated in the read introduction, not inspected proof content.

## Exact-source normalization audit

The initial local write appended one extra terminal LF to each cached text file. This added only a final blank line and shifted no nonblank source line. It was removed from all cached files; every resulting text cache was verified against its pinned Git blob SHA and byte length. `records.json` uses the exact original source line counts and clamps read-range endpoints to those counts. The verified totals are 3,168/3,168 for position 16, 4,019/4,019 for position 17, and 4,458/7,523 for position 18.

## Scope of the speculative transfer

The common project interface remains: construct a specific symbolic minor of full row rank over Q(Z), specialize it integrally without violating signature caps, track the quadratic occupancy contribution to column cost, and establish the small upper bound for that same nonzero integer determinant. No source reading in this batch establishes those bridges or resolves pi's BA/nonBA status.
