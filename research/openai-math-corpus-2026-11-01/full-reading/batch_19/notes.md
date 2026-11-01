# Queue 27 — Eventual categoricity for abstract elementary classes

Complete first reading at pinned `openai/math` commit `4eaee87548250b2682affda7d05d947c05554bae`. The preflight directory tree is `40a644eae73a7039c4a7da684df4d69ee61051d7`. Every cached source file matches its upstream Git blob SHA and byte size. Reading status is **full_body_read**, with complete recursive source mapping and no unread body or appendix source. This status does not certify the proof independently.

## Coverage

18 files; 9,773 source lines; 475,881 bytes. Actual main includes fourteen body sections, the introduction includes one external TikZ figure, and the bibliography is `build/references.bib`. The return-tree figure in section 09 is inline and was read with that section. No appendix input is present. The bibliography style `alpha-linked` is a formatting dependency absent from this paper directory, not an unread body input.

| Source within paper directory | Full range read | Git blob SHA |
|---|---:|---|
| `README.md` | 1–17 | `2cd5f23ffe746f2abf89666a3a6e2341aaaab52b` |
| `build/figures/proof-route.tex` | 1–39 | `d5ebc250bbe36f23ab4c309f68fd41918401ff97` |
| `build/main.tex` | 1–89 | `e4c4e47a7ed322d0a759c1c7289b0d0b000706c5` |
| `build/references.bib` | 1–207 | `2c42f5393d95b7a0484cd49440844c0e9da0a0a8` |
| `build/sections/01-introduction.tex` | 1–308 | `6c1140977719729fad75324f05f97562f37309ac` |
| `build/sections/02-reduction-orders.tex` | 1–387 | `61f2dbf9c113ff7b5d1b851df3572eace7b20f92` |
| `build/sections/03-reversal.tex` | 1–718 | `63293af935405088cc4e3818376077916103db18` |
| `build/sections/04-diagrams.tex` | 1–614 | `bce9ffbc665f2c2c239b5d45d1dd0c91b4586ae9` |
| `build/sections/05-working-schemes.tex` | 1–782 | `5c9db5bc01ee95dae7558aec9179a21d47ec1141` |
| `build/sections/06-calculus.tex` | 1–893 | `80b8af231796f19cc12c1ec7a4e9a10061d48562` |
| `build/sections/07-isolation.tex` | 1–724 | `d140a33d17dd1ad7ef5a7187e1d9826a1fd41536` |
| `build/sections/08-marked-increments.tex` | 1–898 | `36c00399d1739173f63cc3601632bd15298776e9` |
| `build/sections/09-marked-boundaries.tex` | 1–789 | `a59ec84b3a070e8e35a120127e8650eb668f2471` |
| `build/sections/10-geometry.tex` | 1–607 | `b1fb37fd0e7dbba8653182f30b7582957184d673` |
| `build/sections/11-geometric-systems.tex` | 1–905 | `f5ee8d3053c7e00a5e24609703c018b68b9f46e5` |
| `build/sections/12-geometric-boundaries.tex` | 1–514 | `56a9b975471ac95204c16a4a3844fc7b50e15b4a` |
| `build/sections/13-uniform-closure.tex` | 1–850 | `4dc89c2af007a99a3b397dca6360668bd4f706cc` |
| `build/sections/14-cardinal-transfer.tex` | 1–432 | `7a8f4925a781ec74b8fc05e232d21e7c7c2b04d1` |

All source was freshly fetched through the GitHub connector at the pin. Each connector string had one additional final LF after local saving; that LF was removed only after checking that the resulting Git blob SHA exactly equaled the preflight SHA. The full byte audit is in `byte_audit.json`; the include graph is in `source_closure.json`. Granular successful display ranges are in `source_context/fullread_round02_19/reading_chunks.json`. The one initial truncated segment, section 02 lines 12–76, was explicitly reread within lines 1–95.

## Claim as stated

For every infinite cardinal lambda there is an existential cardinal mu(lambda)>=lambda such that every AEC K with LS(K)<=lambda, if categorical in any one kappa>=mu(lambda), is categorical in every kappa-prime>=mu(lambda). The same threshold appears in premise and conclusion. The paper claims a ZFC proof for arbitrary AECs without amalgamation, joint embedding, tameness, no-maximal-models, arbitrarily-large-models, successor-start or regularity assumptions. It gives no explicit numerical threshold. Its qualitative equivalent is: categoricity in unboundedly many cardinals implies categoricity on a tail, including existence.

## Proof architecture

1. Recover each LS-bounded AEC and its strong relation from a set of small model/pair codes. Once a qualitative tail is known for each unbounded categoricity spectrum, Replacement over that set gives one existential threshold also exceeding every bounded spectrum.

2. Extract a coherent finite-support Ehrenfeucht–Mostowski order presentation from arbitrarily large models by successive finite-arity partition thinnings; labels and finite patterns determine literal structure data and strong order substitutions.

3. Build pure orders H(delta)=Q^(Z^(delta)) with balanced nonprincipal cuts and controlled short-sequence orbits. A bounded palette of support records and orbit-arrow profiles, repeated at a small and a larger categorical size, yields simultaneous alignment of every short tuple in a small subpresentation.

4. Eliminate any nonreversible orbit arrow by theta+ rounds, processing all a^theta small-model tuples in a larger mu reservoir. Persistent exclusions contradict the forbidden orbit in the categorical small union. This proves reversal at all sufficiently large categoricity cardinals without a cardinal-exponentiation equality or cofinality restriction.

5. Stabilize surjective tuple-class projections using their bounded palette; countably cofinal closure under the thresholds gives actual bounded-diagram homogeneity. Uniform testing from ordinal-vector witnesses and finite-support event-budget compact placement yield amalgamation and singleton type counts at selected levels.

6. Package finite systems with explicit structural unions, mixed presentations/downward closure, natural bounded diagrams, saturation, and a strong-submodel exhaustion criterion. Extra basis arrow conditions remain separate obligations. Fixed set-state back-and-forth tests provide a uniform diagram code only after the package has been proved.

7. Construct majority approximation sequences with at most |Z|+aleph0 exceptions per query. Their bounded whole diagrams transport unique stationary free extensions. Unequal approximation lengths prove symmetry of products; local character is proved for finite tuples or tuples shorter than the chain cofinality. Universal strips give pointwise-cofinite approximate retractions.

8. Relative isolation either has a uniform bounded fragment or admits a nonfree joint extension with free marked marginal. A large stationary matrix plus Fodor and local character obtains isolation through one tuple length. A common family of retractions then descends isolation for every bounded length simultaneously, including tuples spanning a countable construction. This gives prime all-dominated extensions and whole types for long product-independent lists.

9. Extract marked-pair templates and their bounded comparison tools before proving new chain closure. Two-part testing and scheduled approximate retractions into every limit base handle spanning tuples and yield the full marked working package without assuming it prematurely.

10. Projection type correspondence permits finite marked iteration. Return trees retain full arrows only along branches and boundary maps in a shared accumulator; orbit bounds force all joint boundary comparisons, including actual cross-face equalities. First obtain completion into an extension, then separately exact isomorphisms onto a prescribed target using bounded boundary-fixing data.

11. A singleton splitting tree and type count yield a minimal nonalgebraic type. Forced containment in strong couplings constructs a finitary pregeometry and change-of-base theory. Bounded anchoring makes all root prescriptions a fixed set before choosing a common marked level. All finite exact marked boundaries glue over finite subsets of a large index set; a hypothetical proper no-point-growth pair would produce a large categorical model with too few points. Alternating LS shrinkage extends no-growth and rank=size to arbitrary larger cardinalities.

12. For recursively admissible finite geometric shapes, auxiliary basis pieces give relative-independence arrows. Tiled/buffered order presentations prove independent decomposition and full relative rank before categoricity. A paired induction establishes system classification by root/nonroot basis maps, then the actual working package, then extension classification by arbitrary prescribed complementary basis maps, allowing small or empty complements.

13. Fresh types are unique and lift freely; dominated hulls have finite complementary rank cost for finite tuples and bounded cost for bounded tuples. Exact geometric boundaries first prove actual intersection consistency by lower-dimensional return trees, then use rationally tiled testing and an unoccupied top-only strip to realize partial matches inside the given target.

14. Uniform finite-depth set recursions, one fixed adequacy formula, a set of candidate descriptors/templates, and a fixed finite register of witnessed bound requests yield one majorant Phi. Nested G_(j+1)=G_j intersect Lim(G_j) absorb finite refinements; the unbounded G_omega supplies all finite depths simultaneously. Marked and geometric uses are separate, with geometry chosen after the no-growth threshold.

15. Starting from one common geometric gamma, induct on maximum basis-piece cardinality simultaneously for all finite shapes and monotone unequal assignments. Subcubes retain smaller vertices, and a new axis records the preceding full isomorphism as boundary data while the maximum size decreases. Coherent unions handle successor, regular-limit and singular cardinals. The zero-axis case gives pointed categoricity; LS plus categoricity at one fixed smaller rho puts the fixed point pattern in every larger model. Order presentations prove existence, and small-code Replacement completes uniformity.

## Construction ideas and possible interfaces

These are transferable construction patterns, including distant ideas. Every interface to the pi project is speculative. The current arithmetic target still requires **one and the same normalized integer determinant** to be both nonzero and sufficiently small; neither source cardinality arguments nor generic rank preservation supplies that contradiction.

### 1. Count only genuinely new support events when turning coherent local realizations into a global placement.

**Required hypotheses.** A fixed finite-support template with natural diagrams; at most sigma new support positions; disjoint block transpositions preserve diagrams; every coherent prescription on a parameter set of size at most chi is locally realizable; chi^sigma=chi and chi dominates the fixed language and arities. The non-well-ordered version additionally needs all bounded local prescriptions and closure under the entire sigma+-stage profile-growth recursion.

**Mechanism.** Over each current parameter set, consider all support profiles and one failed query per no-new-event continuation. A final local realization must incur a new old-position hit or a new separation of two new positions at every stage. There are at most sigma such first events, contradicting sigma+ stages. The conclusion preserves complete joint diagrams, including equalities, rather than just separate finite tests.

**Speculative interface.** Speculative interface to the pi project: maintain a finite ledger of row-signature hits, row/column identifications, and separations while refining a candidate common-basis family for Q(Z). If all failures can be charged to first events and local prescriptions are jointly feasible, an adaptive search might terminate with one globally compatible row selection. This is a combinatorial candidate-selection principle, not a determinant estimate.

**Missing step.** One must prove that the actual rational row-rank/common-basis constraints admit the required coherent local realizations and that no-event moves preserve the same formal minor. The source does not control integer specialization height, common denominator, column cost, quadratic occupancy, or the nonzero and small bounds for one normalized integer determinant.

**Source anchors.** `build/sections/04-diagrams.tex:411` (lev:compact-placement); `build/sections/06-calculus.tex:139` (calc:compact-placement).

### 2. Separate a small object that supplies exact alignment from a larger reservoir that pays for enumerating all its constraints.

**Required hypotheses.** The source has a fixed palette of theta-tuple support records and named strong-map profiles, categoricity at repeated sizes a<mu, a^theta<=mu, large cofinality of the ambient index orders, and an aligned a-sized subpresentation computing every theta-tuple record exactly. Strong unions, coherence, and persistence of inability to reach an orbit are essential.

**Mechanism.** Repeat profiles at two sizes. During theta+ rounds enumerate all a^theta tuples of the small model in the larger model, push every tuple that can reach a forbidden orbit onward to an irreversible orbit, and retain aligned small submodels. Their small union still realizes the forbidden orbit by categoricity, contradicting the persistent elimination. No closure a^theta=a or mu^theta=mu is needed.

**Speculative interface.** Speculative interface: distinguish a compact signature/core object used to certify rank from a larger auxiliary sampling or interpolation reservoir used to process many basis constraints. The larger reservoir can absorb enumeration without demanding an unrealistic self-reproduction property of the final determinant size.

**Missing step.** For pi, the palette would have to include exact formal-rank information, the literal chosen row/column map, and arithmetic normalization. A repeated coarse occupancy profile need not determine a nonzero minor. Reservoir size can itself destroy the height/cost inequality, so this architecture offers no quantitative gain until both resources are budgeted.

**Source anchors.** `build/sections/03-reversal.tex:365` (found:orbit-bound); `build/sections/03-reversal.tex:510` (found:alignment); `build/sections/03-reversal.tex:593` (found:reversal).

### 3. Replace ambient-size-dependent obstruction searches by a set of monotone ordinal-vector witnesses.

**Required hypotheses.** The coordinates form one fixed set C. Valid configurations are uniformly set-theoretically definable using fixed set parameters. Their variable portions are ordinal lengths of well-ordered slots, and coordinatewise domination between valid configurations gives an embedding preserving designated terms and natural comparisons. For rational tilings, freeze every exceptional countable fiber hit or split by bounded candidate supports, the initial order, and all support comparisons first.

**Mechanism.** If no set of witnesses generates the class upward, choose a sufficiently long sequence of vectors each dominating none of the earlier ones. Color each pair by a coordinate that decreases; the partition argument produces an infinite descending ordinal sequence. Bounding the sizes of a set of valid unequal witnesses yields a parameter-test bound independent of the ambient base size.

**Speculative interface.** Speculative interface: encode ordered column blocks or multiplicity capacities by a fixed vector of lengths and seek a finite or explicitly bounded set of minimal rank/compatibility obstructions. The finite-coordinate analogue suggests a Dickson-style organization of obstruction search and could separate support shape from ambient multiplicity.

**Missing step.** The determinant constraints must actually be monotone under the chosen embeddings. Adding columns or changing multiplicities can change cancellations, common bases, and cost, so monotonicity is not automatic. The paper supplies existential set bounds, not computable generators or low-height witnesses; a useful pi application requires explicit quantitative bounds.

**Source anchors.** `build/sections/04-diagrams.tex:324` (lev:ordinal-vectors); `build/sections/04-diagrams.tex:359` (lev:testing); `build/sections/12-geometric-boundaries.tex:217` (exact:tiled-testing).

### 4. Use two different approximation lengths to justify interchanging majority evaluations.

**Required hypotheses.** For each bounded query, approximation families have a uniform exception bound e on the relevant row/column comparisons, independent of the lengths of the families. In the source e is an infinite cardinal, kappa=e+ and lambda=kappa+ lie below the eligible level, and diagrams are preserved when queries move over the full bounded code sequences.

**Mechanism.** For the comparison array, all but at most e rows have one majority S, with at most e exceptional columns per good row. Their union has size at most kappa*e=kappa. Choose a column among lambda>kappa outside this union and the bad-column set; its column majority R and more than e good-row entries force R=S. Equal sequence lengths would not support this argument.

**Speculative interface.** Speculative interface: if several row-rank or specialization tests each fail on a controlled set of parameter values, deliberately choose unequal sampling budgets so that one sample survives their combined exceptions. A finite implementation needs explicit uniform failure counts and can be audited as a common-witness selection argument.

**Missing step.** Neither a cardinal majority argument nor typical full rank controls coefficients or proves that the same integer determinant is small. For pi one still needs actual finite root/failure bounds, an admissible specialization range, denominator clearance, and the column/occupancy cost of the common surviving choice.

**Source anchors.** `build/sections/06-calculus.tex:378` (calc:approximations); `build/sections/06-calculus.tex:572` (calc:products).

### 5. Descend whole bounded configurations using one family of pointwise almost-identity maps, even when the configuration spans a limit construction.

**Required hypotheses.** Same-level working objects A<=U admit arrows r_i:U->A such that each a in A is fixed by all but finitely many indices. Smaller compatible sections A-prime<=U-prime must be closed under all selected retraction images. The query size, majority sequence, diagram tests, and approximation costs must remain strictly below the smaller eligible level.

**Mechanism.** For a bounded query Z, all but at most |Z|+aleph0 maps fix Z pointwise. Hence (r_i(d)) computes the actual type of d over A by majority, and the average-descent theorem identifies it with the free lift of the smaller type over A-prime. This handles every bounded d separately without assuming it appears in one stage of a countable exhaustion.

**Speculative interface.** Speculative interface: construct variable substitutions or support-deletion maps that preserve a rank witness while fixing each relevant coefficient/test outside a controlled exceptional set. Close a smaller coefficient/support model under their images before descending a joint row configuration.

**Missing step.** Arbitrary specialization or deletion maps do not preserve formal rank and are not literal retractions on the required determinant. A finite arithmetic version needs uniform exception counts, injectivity on chosen witnesses, and exact compatibility with a common denominator and normalization. The source gives no analytic smallness estimate.

**Source anchors.** `build/sections/06-calculus.tex:828` (calc:retractions); `build/sections/06-calculus.tex:851` (calc:retraction-descent); `build/sections/07-isolation.tex:514` (iso:all-lengths).

### 6. Retain boundary compatibility as the shared invariant while postponing incompatible full realizations.

**Required hypotheses.** Finite system schemes have a proved working package, exact projection/type correspondence, uniform parameter-orbit bounds, and completion in lower dimensions. Boundary data are occurrence-tagged and must preserve joint diagrams, including actual equality on extra overlaps. Full coherent arrows are retained on each branch, but only compatible boundary maps are retained in the common accumulator.

**Mechanism.** A return tree amplifies an alleged boundary disagreement to more than chi distinguishable leaf images over at most chi shared tests. Temporary full maps are used locally and discarded; full maps of leaves are installed only after split tests are fixed. This avoids forcing contradictory full identifications into the accumulator. A separate argument is still needed for onto lifting to the prescribed target.

**Speculative interface.** Speculative interface: in a divide-and-conquer search for compatible row common bases, share exact boundary signatures and overlap identifications while allowing different full interior bases in different branches. Postpone committing an entire candidate until the boundary constraints have been checked jointly. This could make a search certificate more precise about which row equalities and column occupancies are genuinely shared.

**Missing step.** The finite row/basis problem has no established analogue of lower-dimensional completion or the source orbit bound. Compatible marginal data do not themselves imply a common full basis or preservation of the determinant minor. A usable bridge must prove those extension properties and account for quadratic occupancy across faces.

**Source anchors.** `build/sections/09-marked-boundaries.tex:330` (iter:boundary-completion); `build/sections/09-marked-boundaries.tex:589` (iter:boundary-fixed-parameters); `build/sections/09-marked-boundaries.tex:705` (iter:exact-boundary).

### 7. Add one dimension to remember an old full map while performing an induction on a different size parameter.

**Required hypotheses.** Exact boundary lifting is already available at one common uniform level for every finite admissible shape. Pointed models have a finitary pregeometry detecting every proper strong extension; basis cardinality equals model size. Piece-cardinality assignments are monotone under label inclusion, allowing unequal cardinalities. Compatible subcubes use subsets of one fixed global basis partition.

**Mechanism.** At a growth step below a maximum cardinal kappa, include old systems and requested entries in smaller subcubes while gaining exactly chi fresh coordinates on the upward-closed set of growing labels. Treat the old and new cubes as two layers on one additional axis. The old full isomorphism becomes a prescribed boundary map of this enlarged shape, whose maximum piece cardinal is still below kappa; apply the simultaneous cardinal induction. Coherent unions handle singular as well as regular limits.

**Speculative interface.** Speculative interface: enrich a finite linear-algebra state by one extra block/axis that explicitly records a previous row-basis correspondence, so a subsequent extension lemma is required to preserve it rather than selecting an unrelated fresh basis. The source bounded-rank-cost argument also suggests tracking added coordinate rank independently of ambient model size.

**Missing step.** A matrix analogue needs a proved exact extension theorem for these enriched states. General common-basis systems need not have the required gluing property; preserving rank does not bound column cost or quadratic occupancy. The source finite-rank-cost result is qualitative and does not provide the constants needed for the pi determinant contradiction.

**Source anchors.** `build/sections/11-geometric-systems.tex:750` (sys:rank-cost); `build/sections/11-geometric-systems.tex:840` (sys:bounded-prime); `build/sections/14-cardinal-transfer.tex:56` (tail:support-subcubes); `build/sections/14-cardinal-transfer.tex:187` (tail:all-cardinals).

## Dependencies and limits to retain

- This is a complete first reading of the pinned source, not independent proof validation. No theorem is certified correct merely by the reading status; no external cited proof was checked and no PDF compilation or page-by-page PDF comparison was performed.

- The existential threshold and the many set-theoretic testing/closure bounds are not numerical or effective bounds. They cannot be substituted for explicit height and determinant-cost estimates.

- Selected eligible/common levels are unbounded but need not form a tail. The paper uses a separate all-cardinal induction to obtain the tail; arbitrary unboundedness of an eligible class is explicitly insufficient for the intermediate refinements.

- New system packages require actual structural union, homogeneity, saturation and any extra relative-basis arrow condition. Diagram coding is conditional on already established packages and does not prove those properties by itself.

- Finitary products alone do not automatically identify long tuple types or justify limit continuity. The source supplies witnessed chains and a separate argument for independent lists; local character has the stated finite/cofinality restriction.

- Compatible maps on named subfaces do not automatically preserve actual additional overlaps. Both boundary arguments explicitly prove joint equalities before gluing, and completion into an enlarged target is distinct from exact onto lifting.

- No concrete improvement to the pi project is proved here. The bridges are speculative interfaces requiring exact Q(Z) rank/common-basis lemmas, explicit column/occupancy budgets, and nonzero and small estimates for the same normalized integer determinant.

## Scope

Fresh, nonselective, complete source reading at pinned commit 4eaee87548250b2682affda7d05d947c05554bae: README, actual build/main.tex, every line of all fourteen included body sections, the recursively included proof-route figure, inline return-tree TikZ, and complete references.bib. All 18 cached files (9773 lines, 475881 bytes) match the preflight tree blob SHA and byte size. Extra connector LF was removed only when it restored the expected Git blob SHA. No appendix or additional body input occurs in the complete directory tree. The main names bibliography style alpha-linked, a formatting dependency absent from the paper directory tree; it is not an unread body source. Source closure and reading do not imply independent validity, citation verification, successful compilation, or a quantitative pi result.
