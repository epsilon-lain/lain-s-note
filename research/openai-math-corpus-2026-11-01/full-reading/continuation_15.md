# Position 46 — FOUR rigidity continuation

Pinned source: `openai/math@4eaee87548250b2682affda7d05d947c05554bae`.

**Scope:** Fresh continuation for position 46 only. Every line of the six assigned build/sections/four/rigidity files was actually read in the listed untruncated outputs (3,831 lines / 171,677 bytes). The cache was not modified. Git blob SHAs and line counts were recomputed and matched both source_closure.json and batch_24/records.json. This receipt does not read the other manuscript modules, regenerate the ledger, change base records, compile/render the PDF, or independently validate any theorem. The full 721-manuscript nonselective queue remains the parent task; this assignment changes no selection rule.

## Actual reading coverage

| Exact source path | Git blob SHA | Lines actually read | Output chunks |
| --- | --- | --- | --- |
| `preprints/Spacetime-Penrose-inequalities-enclosing-area-charge-and-rigidity-October-5-2026/build/sections/four/rigidity/01-enclosures.tex` | `979086b53210bfc6d618c3cdfb8fc441cb4490fc` | 1–543 | 1–543 |
| `preprints/Spacetime-Penrose-inequalities-enclosing-area-charge-and-rigidity-October-5-2026/build/sections/four/rigidity/03-variation.tex` | `82e3bd729a77425b9cb2c7088cff1dacd209720d` | 1–845 | 1–420, 421–845 |
| `preprints/Spacetime-Penrose-inequalities-enclosing-area-charge-and-rigidity-October-5-2026/build/sections/four/rigidity/04-adjoint.tex` | `996ec7cd4c6a31b7f8495295eedfbab30be1db88` | 1–673 | 1–350, 351–673 |
| `preprints/Spacetime-Penrose-inequalities-enclosing-area-charge-and-rigidity-October-5-2026/build/sections/four/rigidity/05-static-base.tex` | `6a97d57b83ac36a0ce9019792004cbf4812a7942` | 1–506 | 1–270, 271–506 |
| `preprints/Spacetime-Penrose-inequalities-enclosing-area-charge-and-rigidity-October-5-2026/build/sections/four/rigidity/06-classification.tex` | `95896f4f2a777db2ef4a53d715f1f01374c6fbd7` | 1–682 | 1–350, 351–682 |
| `preprints/Spacetime-Penrose-inequalities-enclosing-area-charge-and-rigidity-October-5-2026/build/sections/four/rigidity/07-recovery.tex` | `314b8517257edf6daba3af1a84c6fcaf3021e416` | 1–582 | 1–310, 311–582 |

Total: **3,831 lines / 171,677 bytes**; all six hashes and line counts match the base receipt and source closure.

## Claim as stated

Under the manuscript’s four-dimensional horizon and equality hypotheses (one-ended connected exterior with connected outermost outer area-minimizing future MOTS), these six files claim to recover the original (g,K) as a smooth spacelike hypersurface in the regular maximal Schwarzschild–Tangherlini extension with geometric mass m=sqrt(E²−|P|²). Positive boundary lapse gives a smooth future-horizon section; zero boundary lapse gives the bifurcation sphere. The converse is restricted to hypersurfaces satisfying every stated exterior, decay and horizon hypothesis; they have invariant ADM mass M0 and horizon three-volume ω3(2M0)^(3/2), hence equality, including nonzero K or ADM momentum. This records the source claim rather than independent validity.

## Sections actually read

- 01-enclosures.tex 1–543: full-boundary obstacle perimeter; uniform outer truncation; BV existence and C1,1 obstacle regularity; connected exterior; smooth-cut approximation; strict inner collar; compact all-minimizer tangent measures; right derivative as the minimum of all tied area derivatives; common tangent-plane continuity.
- 03-variation.tex 1–845: fixed supporting ADM charge; full-ball DEC derivative convention; strict mixed conformal direction and differentiated tail/flux estimates; independent end prototypes; nonlinear feasible path and integrability; positive-measure separation retaining all enclosures; local Gaussian-time normal tests with Bianchi cancellation through μ=0; localization to the original horizon via weak contact MOTS regularity and outermostness.
- 04-adjoint.tex 1–673: causal moment measures, distributional constraint adjoint and elliptic smoothness; full metric matter term; finite-jet prolongation and boundary-atom exclusion; causal constant asymptotics and five-probe normalization; boundary flux conditions and lapse dichotomy; constructed stationary metric, null stress, Ric(ξ,·)=0 and vacuum on N>0; both timelike collar expansions.
- 05-static-base.tex 1–506: smooth Killing norm/current identities through arbitrary interior null sets; transverse logarithmic energy; twist cutoffs near null loci, infinity and both horizon-lapse cases; local static equations; connectedness of the positive set; infinite base distance to interior zeros; complete regular boundary attachment and even-metric/odd-lapse reflection.
- 06-classification.tex 1–682: reduced static equations, harmonic coordinates and full multipole improvements; positive monopole; C² scalar-flat conformal double including inversion and possible extra ends; compact smoothing and mixed-exhaustion conformal deformation preserving the distinguished-end flux; cited arbitrary-ends positive-mass contradiction and Euclidean volume rigidity; elimination of interior null regions; round boundary and global Schwarzschild–Tangherlini classification.
- 07-recovery.tex 1–582: exact graph recovery of original g and K; positive-lapse logarithmic time and even angular-coordinate horizon regularity; zero-lapse bifurcation extension; proper global embedding; intrinsic-decay derivation of a spacelike affine end via causal jets and proper covering projection; gauge flux cancellation; explicit four-dimensional boost ADM calculation; horizon three-volume and converse equality.

## Proof architecture

- Build compact full-perimeter minimizing enclosures, retaining original-boundary mass. The right derivative is the minimum over every tied minimizer, avoiding differentiable selection.
- Use a strict conformal DEC/trapping direction and independent end charge prototypes to justify nonlinear feasible paths. Separate their image from a positive cone, obtaining a probability measure over minimizers and positive constraint/boundary multipliers.
- Construct compact normal graph tests with spatial Einstein tensor J⊗J/(μ+1/j); Bianchi and DEC active-ray cancellation survive μ=0 uniformly. Common tangent planes convert aggregate trace cancellation to MOTS frontiers; outermostness and contact regularity localize the probability measure to S.
- Turn positive active-ray measure moments into a smooth causal lapse/shift. Derive both adjoint equations with the vector-identification correction, prolong their jets, remove boundary atoms, and normalize asymptotic constants using independent charge fluxes. Boundary tests give X=uν, surface gravity κ, and the positive/zero lapse dichotomy.
- Construct the stationary development with exactly original K and a possible null stress tensor. A degenerate Killing-current identity and transverse log estimate remove every twist cutoff cost, including arbitrary compact interior null loci; local staticity then gives a connected complete base with only S at finite distance.
- Improve static asymptotics by harmonic coordinates and multipoles, glue the two conformal factors, compactify the minus end in C², and retain possible complete extra ends. If Ricci were nonzero, compact smoothing/variation and a mixed conformal solve force negative mass at the distinguished end; the cited positive-mass input forbids it. Ricci-flat volume rigidity makes the double Euclidean, excluding the null set and identifying the entire static base.
- Integrate the closed time one-form on the now simply connected exterior and use exact metric/second-form identities. Handle the original horizon smooth structure separately in each lapse case using Kruskal coordinates. For the converse derive the boosted affine end from intrinsic decay, cancel linear gauge fluxes exactly, and directly integrate the four-dimensional boost charges.

## Conditional construction ideas

### 1. Retain every tied optimizer before extracting a dual certificate.

**Premises.** A fixed compact space of actual feasible candidates; continuity of each first-variation functional; an attained right derivative equal to the minimum over all current minimizers; a strict direction on every active constraint; and tests that localize the resulting positive measure.

**Mechanism.** The enclosure derivative is min_T a′_T(h), with no differentiable choice of T. A separating positive functional becomes a probability measure on the entire minimizing space, a positive active-constraint measure, and a boundary multiplier. Compact normal tests annihilate active derivatives; the common tangent plane makes the signed trace a single function of the aggregate area measure. Outermostness then identifies almost every selected frontier with the original horizon.

**Possible interface.** For a proposed Q(Z) rank/common-basis construction, retain all admissible tied bases or minors in one candidate family before attempting to extract one common obstruction functional. This can expose a silent change of optimizer between a nonvanishing argument and an analytic bound.

**Missing bridge.** The source candidates are real BV sets and the separation is in a real function space. No compact arithmetic basis family, strict integral correction, row rank over Q(Z), or positive measure supported on integer nonzero minors is supplied. A common optimizer must still yield the very same denominator-cleared nonzero integer determinant that is analytically small; probability averaging cannot replace that witness.

**Source anchors:** `build/sections/four/rigidity/01-enclosures.tex:317` (`rig4:enc:area-derivative`); `build/sections/four/rigidity/01-enclosures.tex:464` (`rig4:enc:common-plane`); `build/sections/four/rigidity/03-variation.tex:493` (`rig4:var:multiplier`); `build/sections/four/rigidity/03-variation.tex:741` (`rig4:var:area-localization`).

### 2. Use finite-jet rigidity plus independent normalization tests to preserve the original object.

**Premises.** A smooth homogeneous prolonged first-order system for the field and its first jet; positivity or causality ruling out affine growth; a Hessian bound with exponent q0>1; and a finite family of independent flux probes normalized in the same frame.

**Mechanism.** The adjoint moment measures become smooth by a diagonal Laplace system. Prolongation turns the first jet into a linear ODE along curves, giving one-sided boundary extension and unique determination. Causality removes the linear term in the end expansion. Five explicit end prototypes then fix the asymptotic lapse and shift from their independent energy/momentum derivatives. Recovery uses an exact metric identity, so both original g and original K are retained.

**Possible interface.** Before comparing arithmetic and analytic descriptions of a determinant, specify one exact coordinate/basis identity and enough independent normalization probes to rule out an unnoticed rescaling, different row module, or basis change. A finite list of exact tests can be an audit device if it determines the object.

**Missing bridge.** The finite-jet system and causal inequality have no established counterpart for Q(Z). Five real ADM probes do not establish full row rank over Q(Z), integral basis compatibility, or denominator bounds. Any arithmetic normalization must prove exact equality of the integer witness and explicitly control its column heights and quadratic occupancy.

**Source anchors:** `build/sections/four/rigidity/04-adjoint.tex:239` (`rig4:adj:boundary-regularity`); `build/sections/four/rigidity/04-adjoint.tex:330` (`rig4:adj:causal-jet-asymptotics`); `build/sections/four/rigidity/04-adjoint.tex:388` (`rig4:adj:asymptotics`); `build/sections/four/rigidity/07-recovery.tex:32` (`rig4:rec:interior-graph`).

### 3. Control a coupled singular error before removing a bad locus.

**Premises.** A globally smooth nonnegative N and degenerate divergence identity; a nonnegative tensor C_b with known transverse eigenvalues; uniform compact coefficient bounds; a logarithmic energy estimate; and a current contraction orthogonal to the degenerating direction.

**Mechanism.** Testing the smooth norm equation by ζ²/(N+ε) proves a uniform weighted log-gradient bound. The twist-current contraction cancels its component parallel to X, leaving only d_D log N. Its cutoff cost on ε<N<2ε is bounded by the integral of 1+|d_D log N|² over that moving band and tends to zero by dominated convergence. Separate boundary-lapse expansions and q0>1 end decay remove the remaining cutoff terms. No measure or regularity of the interior null set is assumed.

**Possible interface.** For the determinant strategy, seek an exact coupled error/occupancy identity and an integrable bound on the actual bad-configuration cost before discarding rare or high-occupancy configurations. The controlling quantity must retain cross terms and the same column-height normalization.

**Missing bridge.** No arithmetic analogue of the current cancellation or logarithmic energy is established. Column-height control alone cannot substitute for the quadratic occupancy estimate, and separate marginal estimates need not control the coupled cost. This method does not prove an integer determinant nonzero or small, and it does not classify π as BA or non-BA.

**Source anchors:** `build/sections/four/rigidity/05-static-base.tex:94` (`rig4:sta:log-estimate`); `build/sections/four/rigidity/05-static-base.tex:182` (`rig4:sta:twist-vanishes`); `build/sections/four/rigidity/05-static-base.tex:209` (`rig4:sta:current-contraction`); `build/sections/four/rigidity/05-static-base.tex:247` (`rig4:sta:null-cutoff-cost`).

### 4. Preserve one exact target flux while allowing uncontrolled auxiliary ends.

**Premises.** A complete scalar-flat C² space with a distinguished end, compactly supported perturbation/smoothing, uniform core estimates anchored at that end, and finite exhaustions with Dirichlet data at the chosen outer sphere and homogeneous Neumann data at every other artificial boundary.

**Mechanism.** Anchored radial/Poincaré control yields a coercive conformal mixed problem despite arbitrary other ends. The finite-domain divergence identity has zero flux at the auxiliary Neumann boundaries. Passing to the limit on a fixed sphere and fixed curvature support preserves the distinguished-end flux exactly; no unproved zero-flux assertion at limiting ends is used. Its positive sign would give negative mass, contradicting the cited smooth positive-mass nonnegativity theorem; Ricci flatness and volume comparison then classify the double.

**Possible interface.** In a family of basis changes or auxiliary corrections, isolate the actual determinant functional at the finite stage and force auxiliary contributions to vanish by exact algebraic conditions before taking limits. A limit argument should transfer the same witness, not assume discarded terms vanish.

**Missing bridge.** The conformal problem is a real PDE construction with geometric fluxes. It supplies no integral solution or lattice-preserving correction and no common basis for Q(Z). The required nonzero determinant, analytic smallness of that same integer, and simultaneous column-height plus quadratic-occupancy estimates remain separate unproved obligations.

**Source anchors:** `build/sections/four/rigidity/06-classification.tex:329` (`rig4:cla:zero-mass`); `build/sections/four/rigidity/06-classification.tex:381` (`rig4:cla:single-end-flux`); `build/sections/four/rigidity/06-classification.tex:399` (`rig4:cla:anchored-estimate`); `build/sections/four/rigidity/06-classification.tex:474` (`finite-domain fixed-sphere flux passage`).

## Limits

- Actual complete source reading is distinct from independent theorem verification. Obstacle regularity, Allard graphical compactness, contact MOTS regularity, weighted harmonic-coordinate/multipole arguments, low-regularity doubling and smoothing, arbitrary-ends positive mass, and cited external results were read as used in the manuscript, not independently established here.
- The six-file receipt depends on horizon/setup and numerical theorem hypotheses defined elsewhere; it does not claim these six files are the whole manuscript. The numerical inequality and original-data rigidity have different hypotheses. The converse explicitly retains outer area minimization and all theorem conditions.
- No global vacuum assumption is used at the start: the source constructs a causal stationary field with possible null stress and handles compact interior N=0 loci before the classification excludes them. Positive versus zero boundary lapse require distinct smooth structures and cutoff estimates; conflating them loses original horizon regularity.
- No result about π BA/non-BA is proved. BA would require a uniform positive lower bound on q²|π−p/q| for all integer p,q>0; non-BA requires a sequence with that quantity tending to zero. This remains unresolved.
- Outstanding arithmetic obligations remain simultaneous: a nonzero integer determinant and an analytic smallness estimate for that exact same denominator-cleared integer; full row rank over Q(Z); one compatible common basis/lattice with controlled denominators; and sufficient quantitative column-height and genuinely quadratic occupancy bounds. Real multiplier separation, causal positivity, PDE existence, geometric compactness and flux identities do not supply them.
