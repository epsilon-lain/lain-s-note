# Continuation 09: position 18 — Exact Birch–Swinnerton-Dyer Formula from Low Selmer Corank

Commit: `4eaee87548250b2682affda7d05d947c05554bae` (`openai/math`).

**Assigned parts complete:** sections 07–10 and `references.tex`, every source line (3,065 total). This agent does not claim the other sections were read. Full source reading is not independent proof verification.

## Source receipt

| Source | Original lines read | Git blob SHA |
|---|---:|---|
| `build/sections/07-theta.tex` | 1–1082 | `76a08c1065b59964f60f9395806a9ef937c16bd9` |
| `build/sections/08-extraction.tex` | 1–1033 | `b2858e5f47069310a766cf403b032704ce5af382` |
| `build/sections/09-unit.tex` | 1–449 | `ee25c7a9e7dc38f8f1d298ca0708c80a2d2ccf24` |
| `build/sections/10-center.tex` | 1–271 | `4ed4109705bcb3a623c587d95baaee4b9db71ac6` |
| `build/references.tex` | 1–230 | `a7d51622016221f621129266ef64ed6903ed2f7b` |

All five exact cached byte sequences match both the source SHAs in batch_06 position 18 and the pinned directory tree. Initially each cache had one extra terminal LF; agent `fullread_06` removed that byte and corrected original line counts. A second Git-blob SHA audit passed for every file. Nonblank line numbering did not change.

## Assigned-part result

For the integral pair quotients U_w=B_w/L_w, prove U_w is a unit, compute the central strict/full determinant on the full Tate lattice, obtain the odd-prime pair discrepancy identity, and combine the non-CM one-sided inequality or CM rank-zero input with the stated low-corank converse and two-primary companion formula to finish the exact BSD formula for rational elliptic curves of full Selmer corank zero or one at some prime.

## Proof architecture read

1. 07: Start from a horizontal zero of the undepleted analytic series b_w^0 and suppose both relevant primitive strict field-test Selmer groups vanish. Choose an algebraic interior point preserving their vanishing and local acyclicity, with safe auxiliary spanning/separation denominators. Two disjoint auxiliary-prime sets and independent exponent vectors ensure one safe set.
2. 07: Build theta sections on canonical ordinary (2,1,1) levels whose depth grows with p-conductor. The text identifies these levels with finite etale flags/frames of the multiplicative character group, proves cusp saturation without depth-dependent precision loss, and calculates both normalized positive p-operators to be 1.
3. 07: Prove a uniform bound for arbitrary bounded locally constant cuts of unary position moments: use compatible upper-triangular cusp changes, coefficientwise group-algebra bounds, squares and saturated integral expansions to remove vertical denominators, then parity projectors and polynomial density. This supplies bounded measures and joint weight limits even for arbitrarily deep cuts.
4. 07: Use Kirillov Fourier inversion at p with a Mellin-one test, and identify upper-unipotent translations with integral-frame Serre–Tate translations. A finite Fourier sum cancels p^{-b} exactly against the orthogonality sum, leaving a bounded unit-restricted measure integral. Toric limits reproduce b_w^0 (up to fixed nonzero scalar/group units) while the companion period remains bounded.
5. 07: Apply bilinear binary Siegel–Weil/Rallis and Waldspurger comparisons with all exceptional raw factors retained. Matrix/toric gamma factors at p cancel including conductors; moving fake Euler and raw additive factors are cleared, fixed Euler factors cancel as rational identities before specialization, and archimedean factorials cancel. The squared cleared identity forces all constant Fourier–Jacobi terms uniformly small. Stagewise choose weight precision first, ordinary weight size next, then increase group depth.
6. 08: Obtain a positive coefficient of bounded valuation at that same nonaugmentation point. One-index unfolding separates a fixed nonzero stabilizer period and two unary periods. Split each moving character into two slot characters, retain integral open-pivot tests, and avoid conductor-dependent averaging.
7. 08: Show the first unary period remains a nonzero analytic function of degree at prescribed t_0 by pushing fixed position measures to the true character group and composing a nonzero measure transform with the nonconstant degree-character map. Restore moving unit cuts as nonzero degree factors, use an ultraproduct-valued analytic limit for varying coefficient fields, and choose degrees/slot parameters where both periods are nonzero. Restore controlled tests using safe finite spanning denominators.
8. 08: Uniform cusp lifting on finite etale canonical branches, Hasse powers and ample twists gives cusp forms without depth-dependent loss. A bounded positive coefficient cancels integral operator relations, producing finite operator orders mapping to O_k/(p^{M_i}) with M_i tending to infinity. The two p-operators retain value 1.
9. 08: Transfer local Galois controls at each classical stage: Hodge–Newton contacts give local ranks 2,3,4 and last quotient/first conjugate line. Inertia- or Frobenius-separated binary projections give inertia-fixed binary allocations; bounded allocation filtrations kill operator nilpotents by a fixed power of a nonzero separation factor, so the congruences descend to reduced eigensystem orders.
10. 08: Form bounded-generator ultraproduct matrix algebras, localize/henselize at the limit map, and retain a faithful reduced ambient realization. The limiting blocks are f=V(-2)chi^{±1}, x=epsilon^{-1}, y=epsilon^{-2}. Lift corner idempotents/matrix units, and use integral contractions of marked finite models to identify admissible cocycles with field-test Selmer classes.
11. 08: Eliminate the cyclotomic y-to-x corner modulo paths through f: local controls make the resulting Kummer class strict at all p-places and unramified elsewhere; the CM biquadratic field has unit rank one and no nonzero global unit with all p-adic logarithms zero. Then a highest-weight argument excludes a globally stable last line at every ambient field evaluation, forcing the y-to-f corner to have zeroth Fitting ideal zero and nonzero residue fiber.
12. 08: A nonzero functional on that corner yields a nonzero extension in V chi^{±1}, strict on the required p-half. Separated binary projections and parahoric rank-one inertia show unramifiedness at each fixed bad place. Shapiro/descent from L to K splits the class between V and its quadratic twist; at least one primitive field test jumps, contradicting the chosen point. This proves the theta implication.
13. 09: At every horizontal DVR factor the strict determinant into fixed singular Euler determinants times the primitive differential determinant. A quotient divisor therefore forces b_w^0 to vanish, hence a primitive jump. This makes each quotient factor horizontal with nonzero reduction modulo (p,t); p-primitivity separately excludes p. Integral Weierstrass division upgrades high-character associations to equality (U_w)=(U_barw).
14. 09: Prove that augmented alternating matrices [A|B] with two extra columns have rank-defect codimension at least two, as does simultaneous defect in two independent even alternating blocks. A tangent-slice kernel inclusion transfers this codimension bound to actual maximal-minor ideals by initial forms in a regular local ring.
15. 09: Add finitely many tame directions at every old quotient divisor. Cup-lift obstruction tests and joint class evaluations span the required independent augmented-alternating or paired-alternating tangent slices. Parity forces one-sided jump dimensions even; two one-sided jumps for one representation lie in the both-strict jump locus. Thus neither both sides for one representation nor one side for both representations can jump together in codimension one.
16. 09: Exact integral specialization in the added directions preserves quotient identities up to units. Any enlarged height-one quotient factor specializes to a nonunit old factor and survives in its localized enlarged ring. Association makes it divide both side quotients; the two theta implications then force one of the excluded simultaneous jump patterns. Hence enlarged and original quotients are units; cancel the same nonzero restored Euler factors and remove artificial variables.
17. 10: At the true-line center compute the strict/full determinant by derived specialization of the full Tate lattice. Begin with Kummer volume 2n_P+2tau_g−s_K, subtract twice the local logarithmic coordinate plus p-torsion, then subtract non-p local torsion twice for split places. Compare its negative valuation with the analytic logarithm expression; exact Haar identities convert surviving local terms into every component number. This yields 2(n_P+tau_g)=s_K+2v_p(c_E)+2sum_q v_p(c_q(E)).
18. 10: Substitute the full group index n_P+tau_g into the earlier discrepancy identity to get X(E)+X(E^D)=0. For non-CM curves both terms are nonnegative by the single-curve inequality. For CM rank zero use the stated Burungale–Flach input after a Galois-stable lattice isogeny to maximal-order CM and descent from the CM field; CM rank one pairs with a rank-zero twist. Thus every odd valuation is exact. Finally the companion converse supplies rank and whole-Sha finiteness, and the two-primary companion supplies positive rationality and the 2-valuation; equality at all prime valuations gives the real exact formula.

## Construction candidates

### 1. Use a low-dimensional tangent slice to prove that actual rank defects cannot have a divisorial component.

Source premises:

- build/sections/09-unit.tex:136–168 proves codimension >=2 for rank([A|B])<s, A alternating and B having two independent extra columns, by an incidence variety with projective kernel vector.
- build/sections/09-unit.tex:170–228 proves that a tangent injection whose first-order kernel is contained in the augmented-matrix kernel forces height >=2 for the actual maximal-minor ideal; it also treats two independent even alternating blocks.
- build/sections/09-unit.tex:288–382 supplies these slices from cup-lift obstructions and independent commutator evaluations, with parity and field-extension compatibility made explicit.

**Speculative interface.** For a symbolic row-rank construction, seek an admissible formal parameter family whose first-order obstruction has an augmented alternating slice or another explicitly codimension-two rank-defect slice. A height >=2 minor ideal is a stronger certificate than finding one generically nonzero minor: it can prevent a single common polynomial factor from explaining every failure.

**Required bridge.** Exhibit a regular local parameter ring for the actual construction, an injected slice compatible with signature caps and column costs, and the kernel inclusion over every extension of the residue field. The arithmetic span and parity in this paper do not supply those facts for an arbitrary determinant family.

### 2. Upgrade integral divisibility to a unit by adding variables, controlling simultaneous failures, and specializing exactly.

Source premises:

- build/sections/09-unit.tex:37–127 separates integral quotient regularity from mutual divisibility and uses the nonzero reduction modulo (p,t) needed for bounded integral series division.
- build/sections/09-unit.tex:252–389 constructs finitely many added directions for every old possible divisor and excludes simultaneous jump patterns in codimension one.
- build/sections/09-unit.tex:396–449 specializes an alleged enlarged factor to a nonunit old factor, keeps it alive in S_(f,y), uses side association plus two implications to force an excluded jump, and cancels equal nonzero restored factors even if they are nonunits at the center.

**Speculative interface.** In a determinant comparison, track the numerator and denominator integrally, identify the finite list of possible common factors, and add controlled variables to make the two failure loci meet only in codimension >=2. Exact specialization back to the original variables can then rule out hidden factors that generic rational comparisons or high-parameter tests miss.

**Required bridge.** Prove an actual exact specialization identity up to a unit, ensure an enlarged irreducible factor specializes nontrivially, establish association of the two relevant quotient ideals, and show every factor forces the controlled failure pattern. These are separate obligations; generic full rank alone proves none of them.

### 3. Use a bounded witness coefficient to cancel scalar relations without losing precision with the size of the family.

Source premises:

- build/sections/08-extraction.tex:73–159 unfolds a coefficient into a fixed nonzero stabilizer period times two unary periods and keeps local averaging denominators fixed.
- build/sections/08-extraction.tex:317–363 restores prescribed local controls using cross-multiplied finite spanning denominators and integral sums, retaining bounded positive content.
- build/sections/08-extraction.tex:452–487 applies an integral polynomial relation to a lifted section and evaluates a coefficient of valuation at most C, costing only the fixed C before obtaining an order-to-residue-ring map.

**Speculative interface.** When proving scalar determinant identities from large families of integral vectors or rows, isolate one coefficient or minor with an explicit lower bound. Apply relations to the full integral object and cancel only against that witness, so the loss depends on its bound rather than the number of terms or family size.

**Required bridge.** Produce the witness for the actual determinant and prove a uniform lower bound under all admissible specializations. Over real or Archimedean estimates, long sums can grow and cancellations matter; the paper's ultrametric coefficientwise bounds do not transfer automatically.

### 4. Estimate a Fourier projection after exact orthogonality cancellation, and choose stage precisions in the dependency order.

Source premises:

- build/sections/07-theta.tex:699–744 identifies CM translations with integral connected frames and proves p^{-b} times the finite double Fourier sum equals a bounded measure integral over units.
- build/sections/07-theta.tex:756–829 pays finite-stage inverse-Fourier losses by weight approximation before increasing conductor, while proving the limiting group polynomials have a common bound.
- build/sections/07-theta.tex:1060–1082 orders the choices: finite group/depth, weight precision, ordinary size, then increasing group orders with slower retained precision.

**Speculative interface.** For a structured averaging or interpolation argument, compute the complete kernel sum before applying triangle bounds; an apparent denominator may cancel exactly through orthogonality. If finite-stage condition numbers grow, separate existence of sufficiently accurate input at each stage from a false uniform operator-norm claim.

**Required bridge.** Identify exact finite orthogonality and the norm appropriate to the determinant problem. The real absolute value lacks the nonarchimedean bound on long sums, and arbitrary-stage approximations must still satisfy the project's integer height, signature, and cost constraints.

## Limits

- Only the assigned sections 07–10 and bibliography were read by this agent. This continuation does not independently assert that the whole paper source closure was read; that merge belongs to the parent.
- This is full-source reading and reconstruction of the stated architecture, not an independent verification of exact BSD, the local estimates, or the cited companion theorems.
- The increasing-depth argument depends on the paper's claimed canonical finite-etale geometry, expansion saturation, arbitrary-cut moment bounds, and cusp lifting without depth-dependent loss; these claims were read but not independently rebuilt.
- The extraction uses a substantial companion Galois package, local Langlands/monodromy controls, bounded-generator matrix ultraproducts, henselian corner constructions, and admissible-cochain comparisons. Their hypotheses and asserted transitions were followed; external results and companion proofs were not checked here.
- The new-direction tangent slices rely on the finite-model, evaluation detection, parity, cup-primitive, and finite-precision Chebotarev results in earlier sections, assigned to another reader. This agent did not independently verify those arithmetic constructions.
- The bibliography was read as source metadata. Its links to companion preprints are main-branch links; no external paper or link target was retrieved as part of this continuation.
- The CM rank-zero assembly invokes Burungale–Flach with a lattice-isogeny and restriction-of-scalars argument; the original external theorem scope and the claimed nonmaximal-order extension were not independently verified.
- Candidate construction transfers are methodological only. They establish neither symbolic row rank over Q(Z), signature-cap feasibility, an integer specialization of controlled height, a column-cost bound, nor a contradiction for pi.

The original batch_06 receipt was not edited by this continuation. The parent is responsible for merging the disjoint reading coverage.
