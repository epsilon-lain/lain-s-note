# Full reading batch 03 — positions 7–9

Pinned repository: `openai/math`, commit `4eaee87548250b2682affda7d05d947c05554bae`.

All three assignments have status **full_body_read**. The recursive directory trees were fetched with `recursive=1`; none was truncated. Every line of the main source, every repository-local included mathematical body/figure file, the bibliography, and the README was read. There are no unread local body or appendix files. These articles have no separate local appendices. The cache contains 28 text files, totalling 6,216 upstream source lines: 2,248 for position 7, 2,681 for position 8, and 1,287 for position 9. Text caches were checked against their Git blob SHAs after restoring the exact terminal newline bytes.

**This records reading coverage, not independent verification of the proofs or their imported theorems.** The binary PDFs were not separately inspected. Standard package source and cited papers are outside article body closure. In position 9, the two external `glyphtounicode` inputs are PDF glyph maps, not mathematical text. Its spectral companion is an explicit substantive external dependency, and is not silently counted as read here.

Full per-file SHAs, upstream line counts, numeric read ranges, claims, proof architectures, assumptions for construction transfers, and limitations are in `records.json`. Exact pinned source text and recursive trees are in `/workspace/scratch/66870cb2f937/source_context/fullread_03/{7,8,9}`.

## Position 7 — A translational tile with no fully periodic tiling in dimension three

Claim as read: a finite tile T in Z³ admits translational tilings but none with a finite-index period subgroup. Its unit-cube thickening admits no full-rank-lattice-periodic tiling in R³ even with arbitrary real translations. This excludes full periodicity, not every nonzero period, and does not assert connectedness.

The proof preserves two different obligations at every encoding: one common solution exists, and every periodic encoded solution decodes to a forbidden periodic object. Keeping the *same* translation set across all simultaneous equations is central. Separate favorable solutions would be useless.

The starting arithmetic object is a strip array obeying a two-digit p-adic affine word rule. Sparse exceptional cells permit an exception-free 4×4 square. Vanishing second differences along horizontal and two diagonal directions give an affine approximation on that square, then a propagation argument gives it everywhere. Nonconstant columns force a nonzero vertical coefficient. A minimal-period argument normalizes the affine form to Bm, samples rows pm, proves that nonconstancy survives the rescaling, and produces a smaller period. Existence is supplied by the last nonzero p-adic digit.

Symbols are encoded by active input digits of finite-valued functions. A graph tiling equation gives one point above every base input. Other equations remove unwanted dependencies. Cycle tests say that a finite map of output coordinates must be a permutation: if one selected digit is inactive, invert successively around the broken cycle; if all are active, choose shift maps that collide two different input tuples.

The especially useful construction is the activation shift list. Its multiplicity function is
`μ=1+δ(0,0)−δ(1,0)−δ(0,1)+δ(1,1)`.
The joint histogram is nonuniform while each coordinate marginal is uniform. The nonuniformity forces ordinary activity, while the uniform marginals admit inverse orderings that make an explicit common solution possible.

Seed activation uses more than averaging. If a seed were nowhere active, its nonnegative integer output histogram Q has mass D=p⁴ and satisfies Q*ω=0. Zero mixed difference forces Q(α,ξ)=u(α)+v(ξ), with u,v nonnegative integral after choosing a minimizing row. Summation gives
`D=b Σu + a Σv`.
Taking b>D kills Σu, so a divides D, contradicted by a=2 or 3 and p>200 prime. Both seeds share one auxiliary output; the explicit construction checks inverse fibres rather than merely checking marginals.

Stacking uses a fresh prime coordinate partitioned into classes whose self-difference sets are the whole cyclic group. There can be at most one contribution per color in any target fibre; total coverage then forces one of every color. The final marked-point construction is equally exact: a large common difference set forces centers of a hypothetical periodic Euclidean tiling onto a grid, and marked-cell coverage is `1−b(a)+b(a−w)=1`. This gives kernel invariance and permits quotient descent.

Transfer leads for the π determinant/transversal project:

1. Seek a nonuniform joint signature budget with balanced marginal row/column counts. Such a budget could detect simultaneous dependencies missed by marginal capacity estimates. This requires an exact coverage/permutation representation, not only an inequality.
2. Encode forbidden simultaneous functional dependencies by cycle collisions, provided outputs genuinely separate into channels with one-digit witnesses.
3. Look for a nonnegative *integer* histogram decomposition whose total mass and an auxiliary modulus force divisibility. The method fails if integrality, nonnegativity, or the exact mixed-difference identity is missing.
4. At every construction stage explicitly keep one common surviving row transversal and prove reverse decoding. Independent choices at different constraints do not establish a global construction.

These are speculative transfers. No admissibility channel model or relevant integer histogram identity for Q(Z) has been established by this reading.

## Position 8 — Positivity of Serre's Intersection Multiplicity

Claim as read: complementary-dimensional nonzero finite modules over every Noetherian regular local ring, with finite-length tensor product, have strictly positive Tor Euler characteristic, including ramified mixed characteristic.

The proof reduces to domain quotients D,E and two finite Galois closure extensions over regular normalization bases. It constructs length algebras whose parameter quotients have positive normalized multiplicity length and whose positive parameter Koszul homology has normalized length zero. Equal characteristic uses perfection and a trace denominator; mixed characteristic imports specific perfectoid results.

There are two separate norm comparisons. The first shows
`|G_D| a_D a_E χ(D,E) = Σσ cσ`.
It does not assert equality of ordinary and normalized lengths for each conjugate. Rational supported K-theory descent supplies equality on invariant norm sums.

Fixing σ, a parameter quotient U of the first length algebra has positive normalized length. It is approximated by finite B-modules M_n, all killed by one fixed power m_B^a, retaining the relevant normalized Tor data and the existing actions. This fixed annihilator gives
`length(M_n) ≤ β₀(M_n) length(B/m_B^a)`,
hence a strictly positive normalized generator density.

Minimal resolutions have negligible normalized ranks above dim E. Brutal truncations become supported after quotienting by those ranks. The article takes care with uniformly supported lifts, supported roofs, and nonlifting idempotents, constructing nonnegative homology evaluations in the quotient.

A warning worth retaining is the rank-one resolution of T/(z^n) with normalization n: it is negligible in rank but has normalized raw homology length one. Thus small rank alone cannot justify throwing away a term with uncontrolled support/annihilator. Koszul tensoring with one fixed power restores uniform bounds.

The second norm comparison over E′ gives
`|G_E| cσ = Σγ bσγ`.
Nonnegative homology coefficients and parameter Koszul vanishing force this second evaluation into degree zero. Minimality supplies a surjection onto β₀(M_n) copies of a positive-length residue quotient, yielding an explicit positive lower bound. This proves each cσ>0 and then positivity of the first average.

Transfer leads:

1. Positive total mass plus a *uniform* per-generator capacity forces positive type density. In row/signature language, one needs a capacity independent of the scale. A growing capacity only gives a potentially vanishing density.
2. Before discarding low-rank pieces in an integer determinant construction, control their coefficient height or arithmetic mass uniformly. Negligible rank need not mean negligible determinant contribution.
3. When symmetries do not preserve an individual evaluation, compare complete orbit sums after proving a descent/comparison mechanism. Termwise invariance is a stronger claim and cannot be assumed.
4. An averaged signed quantity can become strictly positive after a second evaluation with nonnegative components and a quantitative degree-zero lower bound. The exact comparison identities are indispensable.

These are structural analogies, not an established reduction from the π problem to supported K-theory. The novel categorical interfaces and cited perfectoid/vanishing inputs have been read as stated, not independently reconstructed or verified.

## Position 9 — Random Walks on Critical FK–Ising Maps and Liouville Brownian Motion

Claim as read: under the deterministic metric-measure-limit premise and the companion inputs stated in the paper, stationary half-edge random walks on critical spherical FK–Ising maps converge jointly with the same unit-area √3-LQG sphere to Liouville Brownian motion. Cable paths are accelerated by exactly the number n of original edges. Every fixed finite number of conditionally independent walks is retained with its time parameter. Arbitrary starting vertices are not covered.

The article explicitly says metric-measure convergence alone is insufficient. Its spectral companion supplies a common conformal coordinate, a corner tour, local unit-conductivity energy recovery, and centered inverse convergence. Those proofs are not repeated in this article.

The normalization is checked exactly:
`μ_n(v)=deg(v)/(2n)`, `G_n=n(I−P_n)`, and
`E_n(u,v)=½ Σ_e d_eu d_ev`.
Parallel edges remain separate summands. Loops add degree and attempted steps but zero energy. The companion uses unhalved operators 2G_n and 2G; the paper carries this distinction through functional calculus.

Centered inverse convergence yields positive-time semigroup and resolvent convergence. Uniform tour tracking makes multiplication operators converge, giving stationary finite-dimensional distributions. This is still insufficient for paths: rapid excursions between sampled times must be ruled out.

Buffered, clipped cutoff recoveries are uniformly accurate and have no excess energy. Resolvent smoothing gives bounded generators. A reversible forward–backward martingale estimate controls the supremum of a test-function error along an entire stationary trajectory by its L² norm and Dirichlet energy. Optional sampling then provides a stopping-time increment bound; finite separating cutoff families and vanishing jumps give continuous path tightness.

A metric bridge construction embeds all cable graphs and the limit in one compact space without extending the conformal coordinate to edge interiors. An independent exponential attempt clock has a maximal L² deviation tending to zero; continuous-limit modulus control removes this clock. Loop attempts are counted in the comparison. Finally, conditional product laws on each environment and bounded convergence yield the joint random-environment result.

Transfer leads:

1. If one can construct a reversible stationary chain on admissible Q(Z), small quadratic energy of a signature/determinant statistic could potentially give whole-search-trajectory control. This would need the chain, its stationary capacity law, and all estimates; none is supplied here.
2. Finite collections of favorable constraints do not automatically produce one infinite construction. A compactness/tightness bridge needs uniform control through refinements, analogous to the distinction between finite-dimensional laws and whole paths.
3. Carry multiplicities and normalizations through one common representation. Unweighted proxies can erase the factors that determine the final clock, or in a determinant setting the final integer lower/upper comparison.

These probabilistic ideas are speculative construction principles. The strongest concrete connection in this batch remains the tiling paper's dependence tests and nonnegative integer capacity obstruction.

