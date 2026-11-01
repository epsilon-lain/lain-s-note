# Full-reading batch 02

Pinned repository: openai/math at 4eaee87548250b2682affda7d05d947c05554bae. All 42 authored/provenance source files fetched and byte-verified against Git blob SHAs. Full-body reading is not independent theorem validation. No manuscripts were excluded by relevance.

| Queue | Status | Source lines read | Source files read |
|---|---|---:|---:|
| 4 | full_body_read | 1727/1727 | 10/10 |
| 5 | full_body_read | 4373/4373 | 14/14 |
| 6 | partial_body_read | 4376/12204 | 9/18 |

Queue 06 remaining 03–06 is assigned to fullread_07; remaining 10–11 is assigned to fullread_05. This file claims only this reader’s own exact ranges.

## Queue 4: The symmetric Mahler conjecture and its equality cases

Status: full_body_read. Pinned authored main source and all eight included body sections read completely. Figure captions read in source; rendered figure PDFs, bibliography entries, third-party TeX implementation packages and separate Lean formalizations are outside this full-body reading scope. Full reading is not independent theorem validation.

Claim as stated: For every origin-symmetric convex body K in R^n, |K||K°|≥4^n/n!, with equality exactly for invertible linear images of Hanner polytopes. The introduction also derives the even functional Mahler and symmetric entropy–transport inequalities from the claimed all-dimensional geometric result.

Proof architecture:

1. Construct a convex conformal lens whose uniform boundary angle gives an independent uniform sign and a uniform imaginary coordinate; establish strict curvature and λ(1−η)∼(2/π)η log(1/η).
2. Use Stokes around an isolated homogeneous holomorphic zero to lower-bound total squared Jacobian-minor mass, then apply it to high powers of inverse-lens coordinate maps. Rescaled radial variables yield Σ_I P_I≥1 over independent row bases.
3. Convert each boundary feasibility probability into a real integral with factor |det B_I|/4^d. Extra ties are null, and unique maximizers make feasible signed simplices disjoint in interior, giving 1≤d!4^(−d)∫|Σ_X|≤d!4^(−d)|A||A°|. Approximate arbitrary symmetric bodies.
4. For equality lift by one segment, approximate the polar product by inner polytopes, and turn vanishing integrated missing volume into nearby feasible point–direction pairs using positive-volume product neighborhoods. Compactness gives one common feasibility witness even for degenerate limiting row representations.
5. Compare all positive representations using that common witness. Let the point approach the lifted apex; the endpoint ratio of λ yields support inequalities and a median for each triple in the polar norm.
6. Convert medians to the three-ball intersection property and invoke Hansen–Lima's finite-dimensional classification. Polarity and norm-sum formulas identify precisely the linear Hanner equality bodies.

Construction interfaces (speculative):

- A global minor-mass budget can locate a nonzero feasible basis before optimizing its cost. For candidate row transversals Q(Z), attach strip constraints to the actual row vectors and view Σ_I P_I as a distribution of independent bases. Test a cost filtration by column degree/denominator size and seek a mass inequality concentrated below the analytic smallness threshold. This could separate proving some minor nonzero from charging its columns, but the existing lemma alone gives no cheap basis and cannot repair deficient rank. Assumptions: The candidate row family spans the intended determinant dimension, over the field relevant to the integer determinant.; A holomorphic constraint construction can be made whose minor factors are the same candidate determinants and whose compact sublevels and isolated zero are proved.; An additional cost-weighted estimate, absent from the manuscript's ΣP_I≥1 statement, controls the total degree/height of a feasible basis.
- A single dual feasibility witness can certify all competing row costs, including singular limits. Represent a target signature or budget direction by costed candidate rows. Find a dual vector Y giving c_i·Y=cost_i on selected rows and h·Y≤cost_h on every candidate. Then every competing representation has cost at least ξ·Y. The manuscript's integrated deficit and degenerate-limit treatment suggest how to retain such a certificate when candidate bases coalesce. Apply first to a finite model of signature capacities/Q(Z), keeping this optimization certificate separate from the needed integer nonzero determinant. Assumptions: The proposed determinant/column costs admit a convex positive-representation or matroid base formulation.; The same chosen witness is valid for every candidate constraint; competitors cannot each choose a different witness.; Any zero-deficit-to-pointwise upgrade uses positive-volume neighborhoods and actual compactness, rather than assuming pointwise convergence.

Limits:

- No independent checking of the claimed Mahler resolution or its external classification inputs was performed.
- ΣP_I≥1 concerns an already spanning finite row family; it does not prove rank for the user's actual determinant.
- The lens includes π in its normalization, so reusing the analytic device in a π rational-approximation argument requires a noncircular definition and explicit integrality/denominator accounting.
- No degree/height estimate suitable for BA or non-BA π is supplied by this geometric proof.

Exact continuation:

None in the authored body closure.

## Queue 5: The Unique Games Theorem

Status: full_body_read. Pinned authored main source, all eight included body sections, complete finite-consequences appendix and all three included figure sources read completely. Bibliography entries, third-party TeX implementation packages and separate Lean formalizations are outside this full-body reading scope. Full reading is not independent theorem validation.

Claim as stated: For each fixed ε,δ∈(0,1/2), a deterministic polynomial-time reduction from 3SAT to explicit unweighted simple bipartite Unique Games over a fixed K=F_2^s has completeness≥1−ε and soundness≤δ; every edge is a translation constraint. This is a manuscript claim, not independently validated.

Proof architecture:

1. Construct latent additive noise using characteristic-two quadratic finite-field blocks. Uniform genericity over every linear extension map and random orientations permit adaptive restrictions; the harmonic potential H_r retains positive rank while nonlinear change decays. Quotient shift kernels and add injection-indexed copies to realize any sufficiently large fixed output alphabet.
2. Preprocess a weighted nearly satisfiable parity gap to three distinct variable IDs per occurrence. Form affine table keys that retain support, actual-bit coordinates and exact intercepts; identify keys through projection and fold free constant translations. The resulting rank-one equality test has near completeness from global assignment evaluation.
3. Fourier comparison turns high latent-test acceptance into ordinary shortcode acceptance. Grassmann expansion yields a bounded row/column slice. Erase all good row fibers and apply concentration uniformly over every slice–target description to force positive good-advice mass; folding forces a nonempty short set of actual answers.
4. Compare the full sparse-projection padded table law with unrestricted sampling by a product chi-squared bound: kβ²→0. Transfer the selected witness using a positive probability difference, then a private Fourier decoder reaches the projected short answer set with fixed positive probability. The decoder never receives the other side’s witness.
5. Choose β=k^(−2/3), so kβ→∞ clean singleton coordinates remain. Condition only on the clean mask, advice intercept and outside data, retaining the independent equation–variable law. Each side reconstructs its original local input; projection-game repetition makes full agreement at most exp(−c k^(1/3)), contradicting the decoded constant lower bound.
6. Largest-remainder rounding, four-edge subdivision and alphabet-independent repetition give explicit unweighted simple bipartite translation games with prescribed errors. Check quantifier order and polynomial encoding. The appendix supplies simple strong-left games and a rational finite ordering composition using local smoothing, low-influence estimates, influence decoding and conditioning on distinct neighbors.

Construction interfaces (speculative):

- Use a harmonic rank potential to compare useful contraction against loss of independence. In a costed-column determinant construction, track H_r=Σ_{j≤r}1/j instead of raw remaining rank. If collisions among candidate Q(Z) row transversals have probability O(rθ), expected H_r loss is O(θ), permitting many smallness-improving steps while retaining positive rank with explicit probability. A finite simulation should measure rank loss conditioned on signature occupancy, then test whether candidate selection obeys a uniform extension-map bound. Assumptions: A concrete family of permitted row/column restrictions has rank drop at most one per step, with conditional loss probability≤c r θ at current rank r.; The rank-loss estimate holds for every extension map that earlier choices can induce; assuming independence from an adaptive construction is insufficient.; The operation improves the analytic determinant size on the same object and preserves its integer structure.
- Prove genericity uniformly over adversarial coefficients using an augmented nonzero minor. The latent proof handles all maps g simultaneously by showing a nonlinear target evaluation vector is outside a matrix column span, then exhibits one augmented minor as a nonzero bounded-degree polynomial. For quadratic occupancy/signature constraints, formulate potential zero terms as target-in-span conditions; prove a formal augmented minor nonzero via monomial/parity support, then specialize integer parameters away from its zero set while costing the specialization. This is a route to certify that structural zero terms do not consume every candidate, rather than counting signatures alone. Assumptions: The evaluation vectors have algebraic parameter dependence and an algebraic obstruction to column-span membership.; A nonzero polynomial minor is established symbolically before specializing parameters.; Finite-field genericity must be replaced or lifted with explicit integer heights and rational-approximation errors for the actual π problem.
- Keep nonlinear stability and linear detection as separately quantified requirements. Treat the analytic smallness of an integral or Padé-type column as the stable nonlinear observable, and determinant rank as detection by a linear family. Search for a perturbation of row choices that changes analytic errors rarely but is detected by any large annihilator/signature space. The manuscript shows that these goals can coexist abstractly; it does not construct such a perturbation for π. Assumptions: The two properties are verified under the exact same noise or perturbation law.; Nonlinear error probability and linear detection threshold are fixed before ambient dimension/number of columns grows.; Any lifting to rational approximation controls denominator growth from products and quotients.
- Use a sparse perturbation regime with many useful coordinates but asymptotically unchanged full law. For candidate row transversals Q(Z), design a mixture that marks a coordinate with probability β. A full-law bound of order kβ² can vanish while kβ diverges, giving many special rows under nearly the same candidate distribution. β=k^(−2/3) balances these requirements. This is directly relevant to quadratic occupancy: a proposed capacity improvement should be checked under the actual joint row law and include the cost of conditioning on distinct/signature-safe selections. It remains speculative until the local likelihood ratio and determinant nonvanishing event are established. Assumptions: Independent block contributions have a bounded chi-squared divergence between the ordinary and special-coordinate laws.; The ambient alphabet and all local dimensions are fixed before the number k of blocks grows.; A nonzero determinant or capacity event is transferred using the joint table distribution; agreement of separate marginals is insufficient.

Limits:

- No independent verification of the claimed Unique Games resolution or of the exact cited external expansion/repetition/invariance hypotheses was performed.
- The finite-field genericity, characteristic-two trace and square-root calculations do not automatically transfer to characteristic zero.
- The enormous but finite recursive gadget and all constants are fixed before input size; the proof supplies no usable integer-height or denominator bounds for rational approximation.
- A distributional or combinatorial analogy supplies neither an integer determinant nonzero certificate nor an analytic estimate below one for that same determinant.
- No theorem about π, irrationality measure, BA or non-BA follows.

Exact continuation:

None in the authored body closure.

## Queue 6: Uniform bounds for planar polynomial limit cycles

Status: partial_body_read. Main source, macros, introduction, entire finite-logarithmic-flags section 02, one-rate regular/stable section 07, two-scale mixed-passage section 08, finite-geometric-subdivision section 09, and proof-roadmap figure source read. Sections03–06 and 10–11 plus remaining three figure sources are fetched but unread by this reader. Root assigned 03–06 to fullread_07 and 10–11 to fullread_05. No relevance exclusion. Primitive proofs in07–08 are read but their earlier abstract packet/calculus hypotheses await the other reading.

Claim as stated: For every degree d≥1 there is a finite B(d) bounding the number of limit cycles of every real planar polynomial vector field of degree≤d throughout the plane, without coefficient/stability restrictions; no effective formula for B(d) is claimed.

Proof architecture:

1. Encode an alleged escaping solution sequence in a finite saturated logarithmic flag Z_1≫…≫Z_m≫1 with exact triangular log relations; preserve existing free coordinates when promoting vanishing ordinary coordinates.
2. Construct complex paths through far anchors independent of ordinary backgrounds. Quantitative angular rates produce nested good contours and buffered sides where a selected exponential scale dominates later ones; finite-jet controls retain these margins under changes of background.
3. The introduction states that two lateral asymptotic packet trees, with full domain/transition/error data, yield a first nonzero term or exact vanishing via weighted Cauchy corrections and maximum principles. Differentiated implicit calculus and elimination would then contradict isolated escaping zeros.
4. The full one-rate section 07 and mixed-clock section 08 construct actual holomorphic passage maps and finite-charge descriptions. They use formal normalization, fixed-one-side tail convergence, Gevrey optimal truncation, exact clusters at colliding residue costs, safe meromorphic columns, compatible endpoint charts and ordered retirement of numerical cutoffs. Earlier abstract packet/calculus hypotheses and the additive section06 remain unread by this reader.
5. The fully read geometric section 09 uses simultaneous subanalytic preparation of rational fields, including every pole modulus and numerator balance radius. Compact pole-free rational families have finite jet-order bounds; the induction potential 2m+e terminates box refinements, and a one-redo argument controls nonzero resonances. Jordan-curve crossing signs bound boundary events, yielding a finite passage alphabet and bounded cyclic words uniformly over coefficient degeneracies.
6. An independent absolute-to-uniform projection-counting principle is stated to bound fiber components, hence hyperbolic cycles. Signed rotation and spatial rescaling are stated to include nonhyperbolic isolated cycles throughout the plane.

Construction interfaces (speculative):

- Protect zero-term analysis by an exact scale flag and distinguish zero germs from cancellation on a sequence. For Δ_N(q,p) built from costed columns near p/q≈π, encode factorial, exponential and powers of q scales in a flag and reduce exact dependencies before ordering determinant terms. Label the zero-term loss by the first genuinely nonzero signature coefficient. Preserve free approximation scale variables when inserting −log|qπ−p|. This could expose where quadratic occupancy arguments count only formally distinct terms that actually cancel. Assumptions: Every determinant contribution/error term can be expressed in a finite ordered collection of actual scales, with exact dependencies and bounded ordinary backgrounds.; After cancellation, an asserted leading coefficient is a nonzero analytic germ, not merely nonzero at selected sample points.; Any claim that all asymptotic terms vanishing implies exact zero requires a proved domain/separation theorem; real flatness alone is inadequate.
- Use finite-request safe intervals to avoid all finitely many bad denominator or pole lattices without asserting equidistribution. When coefficient transformations create finitely many sine, gamma or reciprocal factors, select a common detuned parameter interval by bounding time spent near each bad lattice and union bounding; no independence or incommensurability of sweeps is needed. Then separately prove an allowed integer index exists in that interval and charge the interval width. This could stabilize candidate-row transversals near singular factors but is not itself a Diophantine approximation theorem. Assumptions: Relevant sweeps are eventually increasing with derivatives tending to infinity and |x''/(x')²|≤C/x.; Only a fixed finite list of bad lattices is excluded at a given determinant size.; The application needs continuous parameter choices or a separate integer-point argument; safe real columns do not automatically supply integer N/q.
- Maintain explicit buffers and fixed derivative-loss exponents before strengthening accuracy. For a near-vanishing determinant perturbation or parameter correction, write a ledger of analytic decay versus derivative and denominator losses. Reserve a positive buffer that survives the row/column transformations. The flags proof's hierarchy and high finite-jet sensitivity suggest a disciplined way to prevent a later approximation accuracy requirement from silently increasing the column-degree budget past the small-integer threshold. Assumptions: The analytic upper bound can be continued to complex domains with quantitative buffers, not just real estimates.; Only finitely many derivative/coherence orders are needed per determinant and all cost exponents are fixed before accuracy is increased.; Dependent scales are reevaluated at current background values after implicit changes.
- Combine colliding terms exactly before separating costs or taking absolute values. The mixed-passage proof keeps the analytic residue cluster (exp(−C0 L)−exp(−C1 L))/(C1−C0), whose equal-cost limit is L exp(−C0 L), instead of treating two singular residue terms separately. In the determinant expansion for π, group row transversals with the same signature/exponent before capacity counting. A signature may have many combinatorial occupants but a zero total coefficient. Conversely near collisions may contribute polynomial prefactors rather than independent exponential terms. Derive an exact cluster formula first, then cost its polynomial factor in the analytic smallness ledger. Assumptions: A finite cluster identity is proved at coinciding exponents or signatures, with its removable singularities retained.; Every clustered coefficient is added using exact algebra before calling it nonzero or estimating its size.; The same clustered object remains an integer after denominator clearing; this must be checked separately.

Limits:

- Sections03–06 and 10–11 remain unread by this reader, including the abstract separation/calculus/elimination, additive construction, terminal transfer reduction and final counting argument. Other named readers have those continuations assigned; their coverage is not asserted in this record.
- The packet constructions in07–08 reference earlier hypotheses not yet read by this reader; reading their text is not independent verification of those interfaces.
- Sequence-dependent flags are explicitly not a finite covering of parameter space.
- Exact vanishing requires the later complex-domain separation hypotheses; no real-flatness shortcut is valid.
- No π rational-approximation bound or determinant-rank theorem is established.

Exact continuation:

- build/figures/coupled-bands.tex: [[1, 50]]
- build/figures/matching-arguments.tex: [[1, 46]]
- build/figures/retained-charges.tex: [[1, 32]]
- build/sections/03-separation.tex: [[1, 1428]]
- build/sections/04-calculus.tex: [[1, 1144]]
- build/sections/05-elimination.tex: [[1, 1147]]
- build/sections/06-additive.tex: [[1, 1748]]
- build/sections/10-terminal.tex: [[1, 1201]]
- build/sections/11-counting.tex: [[1, 1032]]
