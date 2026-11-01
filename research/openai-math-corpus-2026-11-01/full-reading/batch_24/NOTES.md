# Queue 46 — lead reader: source mapping and neutral numerical module

Pinned source: `openai/math@4eaee87548250b2682affda7d05d947c05554bae`, family 260, *Spacetime Penrose inequalities: enclosing area, charge, and rigidity*, October 5, 2026. This is a source-reading receipt, not an independent proof verification.

## Source closure and coverage

The actual PDF source root is `build/paper.tex`. It inputs `preamble`, `sections/00-introduction`, and `modules`, then uses `references.bib`. The modules manifest orders the neutral numerical proof, neutral many-end bridge, neutral rigidity, charge, four-dimensional numerical and rigidity proof, and local AdS proof. Recursive input tracing finds all 49 non-PDF text blobs in this directory and no missing or unreferenced text blobs, including the three included TikZ figures. The README is retained for identity/version metadata. The compiled text closure is 1,025,372 bytes and 23,468 source lines.

Exact bytes are cached under `source_context/fullread_round03_24/46/`. Every cached file's Git blob hash was recomputed as SHA1(`blob <byte-length>\0` + bytes) and matched the pinned tree. No terminal newlines or other content normalization needed to be restored in these 49 responses. The cache writer would restore a terminal LF only if that exact change recovered the pinned blob hash. `byte_audit.json` records this per file. `source_closure.json` records the root traversal and complete list. A Git hash match authenticates bytes, not mathematical correctness.

This lead actually read 15 files, 6,472 lines: README; root; modules; preamble; introduction; proof-routes figure; neutral setup; all neutral sections 02–07; neutral-bridge; all bibliography metadata. Every other file has an empty `read_ranges` in this base receipt. Their download is not counted as reading. Separate reader continuations are required before full-body promotion.

## Claims read with their premises

The neutral exterior theorem (`build/sections/neutral/setup.tex:92`) assumes a smooth connected orientable three-dimensional exterior, complete with the compact nonempty boundary included, one asymptotically flat end, DEC, `g-delta=O_2(r^-q)`, `K=O_1(r^(-1-q))`, `q>1/2`, integrable constraint densities, finite ADM charges, and weakly future outer trapped boundary toward the end. Initially assuming `E>|P|`, it asserts

\[
\sqrt{E^2-|P|^2}\ge\sqrt{a_g(\partial\Omega)/(16\pi)}.
\]

The area is a **full obstacle perimeter infimum**, allowing disconnected cuts and coincidence with the original boundary, and filling bounded complementary pockets. No extension of constraints across the boundary is assumed. These cut conventions are part of the theorem and are used repeatedly.

The bridge (`neutral-bridge.tex:99`) removes the timelike assumption: the positive conformal family raises energy linearly, leaves momentum fixed, and never lowers the positive enclosing infimum. If the original energy were at or below `|P|`, applying the already established timelike inequality just past the null crossing and taking the family parameter down to that crossing would contradict positive area. The complete-data theorem (`neutral-bridge.tex:42`) allows finitely many ends and disconnected admissible boundary/selected exterior. It restricts to the selected-end component, removes other tails with spheres whose normal into the retained exterior has strictly negative expansion, and uses the **inequality** between the retained one-end area class and the complete-data class; it never identifies distinct infima without proof.

The introduction also states the scope of the delegated modules. Neutral equality concerns the original connected one-ended exterior under extra MOTS, outer area-minimizing and precise outermostness hypotheses. Charged one-end data use electric plus magnetic charge and electromagnetic momentum; the advertised result is an upper-area inequality `r<=m+sqrt(m^2-Q^2)` with `m>=Q`, and its usual mass-form consequence requires `r>=Q`. There is no charged equality classification. Four-dimensional equality is a connected one-ended horizon theorem. Local AdS is for a fixed TT seed and specified branch on a small interval depending on both choices. These are overview coverage here, not substitutes for full module reading.

## Neutral proof architecture

1. **Strict direction and controlled rest end.** The exact conformal formulas and regularized drift-Poisson problem produce a positive function with a negative finite outward flux (`02-exterior-reduction.tex:81–312`). The family preserves DEC and boundary signs for every finite nonnegative parameter, increases energy by an explicit positive linear term, leaves momentum fixed, and bounds all enclosing areas multiplicatively. The end reduction constructs a matching-charge boosted Schwarzschild plane, blends far away, projects away finite moments of the linear constraint error, applies compact inverses, bends the exact vacuum end to rest, and repairs its vanishing DEC deficit by an explicit radial conformal ODE (`02:409–982`). Final strictification is chosen separately for each approximant; no unsupported uniform higher-derivative estimate of the original end is used.
2. **Compactness of enclosing minimizers.** A fixed competitor and uniform end geometry give an `r^2` lower bound on the portion of any minimizer reaching radius `r`, forcing all minimizers into a common compact set (`02:314–401`). Local metric convergence then controls the full infimum from both directions. Obstacle regularity and outward smooth approximation preserve coincidence conventions.
3. **Prepared geometry.** Shifted black and white total trapped-region theorems provide disjoint stable frontiers at negative, right-continuous thresholds. A zero stability eigenvalue is handled with an augmented constant-expansion implicit-function theorem and maximality (`03-trapped-domains.tex:208–287`). Arbitrary weak exterior-support sets admit smooth larger-threshold enclosures, without assuming reach or support curvature bounds (`03:295–555`). The geometry, threshold gaps, collars, offset and positive margin weight are fixed before deformation parameters.
4. **Coupled deformation and algebra.** Unknowns `f,Z` determine `t` globally through `Z=t+log(1+l(t)^2|Df|^2)/8`. The comparison metric is `hat(g)=exp(4t)(g+l^2 df^2)`. A stationary Lorentz calculation gives its scalar identity with a divergence and a coercive quadratic form (`04-deformation.tex:97–311`). Completing actual squares gives polynomial losses in `n`, retaining the degenerate axial normal coefficient. A prescribed inner flux makes the comparison boundary strictly mean concave. A four-stage homotopy removes drift, adds directional collar corrections, interpolates boundary heights to zero, then scales the final source/flux to zero. Its last compact map has output identically zero.
5. **A priori order.** The floor `t>=-epsilon` is excluded at contact by direct curvature/flux algebra plus a polynomial penalty, without an upper `Z` or Hessian bound (`04:541–812`). Relaxed limits, monotone sigmoid compositions and weak-support enclosures exclude wrong boundary layers (`05-apriori.tex:30–212`). Signed collar slopes and a polynomial weighted trace estimate precede graph Sobolev comparisons. The factor `tau/ell=ell^(1/2)` absorbs polynomial boundary losses before arbitrary fixed-`n` graph/Moser constants are introduced (`05:346–536`). Fixed-`n` Moser iteration, maximum principles and the implicit relation then bound `Z,t,Df`.
6. **Regularity and actual existence.** An axial gradient identity avoids differentiating measurable axial coefficients. Cordes estimates, a gradient drop/affine-flatness alternative, a fixed-axis correction, boundary reflection and Campanato iteration give gradient Hölder and Hessian Morrey estimates (`06-elliptic-existence.tex:19–387`). An exponential-deficit oscillation argument handles the second equation's quadratic growth with the Hessian-square Morrey measure. Schauder/interpolation then gives classical estimates. A logarithmic two-point modulus proves the scalar trace solution operator exists uniquely for every bounded `Z` input **without a floor assumption** (`06:654–818`); only afterward is the compact coupled map defined and degree applied on a moving admissible open set (`06:821–999`).
7. **End, flux and comparison.** Exhaustion is taken at fixed `n`; exponential graph decay and a scalar reparametrization show `Z,t=O_2(r^-1)` and a finite flux (`07-mass-limit.tex:28–137`). The curvature identity, fixed margins and exact boundary condition give `R_hat>=0`, completeness and strict mean concavity. The exact mass formula is `E_hat=E-F_n/(8*pi)`. Polynomial collar traces absorb the boundary loss, while the floor penalty costs `poly(n)*(ell+ell^2/tau)` (`07:225–318`). With `ell=exp(-epsilon*n)` and `tau=ell^(3/2)`, this tends to zero. A minimizing minimal enclosure satisfies the external Riemannian Penrose hypotheses, and graph domination plus the floor gives its area at least `exp(-4epsilon)*A_*`. Limits are taken in the stated order: exhaustion at fixed parameters, then `n->infinity`, then `epsilon->0`, then the separate rest-end approximation.

## Conditional construction ideas and missing arithmetic bridges

The project's target remains the BA/non-BA behavior of `q^2*|pi-p/q|`: a uniform positive lower bound for all integer `p,q>0` establishes BA; a sequence tending to zero establishes non-BA. Neither conclusion is obtained here. The following ideas are conditional research interfaces only.

### Finite obstruction projection on a common basis

**Source anchors:** `02-exterior-reduction.tex:518–597,719–773`.

The compact linear inverse works only after affine and Killing-field moments vanish. The paper lists the adjoint obstruction kernels explicitly, constructs positive-definite Gram matrices from a fixed bump and basis, subtracts the corresponding moments exactly, and then applies a support-preserving inverse. The matched translation charges and integrable-source first-moment estimates make the correction small enough for the nonlinear DEC repair.

**Possible arithmetic use:** before designing an approximation family, isolate every finite obstruction functional and construct one common correction basis. A candidate exact row-rank argument over `Q(Z)` should identify an actual nonzero minor, not merely a real span at individual parameter values.

**Required premises:** common rational/integer coordinates; exact compatibility identities; a controlled denominator-clearing lattice; a nonvanishing minor valid across the necessary specialization range; conditioning or valuation bounds sufficient for every quantitative step.

**Missing bridge:** the source uses smooth real bumps, and inverse Gram entries are real. Its positive definiteness does not yield integral corrections, specialization-safe row rank over `Q(Z)`, a common integer basis, or any nonzero determinant involving pi. Pointwise rank can coexist with degenerating minors. Quadratic occupancy does not repair those omissions.

### Exact identity before analytic smallness

**Source anchors:** `03-trapped-domains.tex:325–442`.

The parallel-neighborhood estimate cannot tolerate an error `z*|Hess(phi)|`, because the support Hessian is uncontrolled. The author solves a Jacobi endpoint problem and chooses a local unit vector field whose endpoint map has differential **exactly** equal to parallel transport. The support Hessian cancels without estimation; only a bounded second-jet `O(z)` error survives.

**Possible arithmetic use:** make the determinant in the exact arithmetic nonvanishing lemma literally the same denominator-cleared integer as the determinant in the analytic smallness estimate. Exact basis and normalization bookkeeping must precede estimates.

**Required premises:** an equality of algebraic expressions on one common basis, a proven nonzero integer witness, and an analytic bound below one for that same witness.

**Missing bridge:** no pi-dependent determinant is constructed here. Near equality of two expressions can leave an error multiplied by an uncontrolled minor or denominator. The geometry illustrates why a formal resemblance between witnesses is insufficient; it supplies no arithmetic cancellation theorem.

### Quantitative firewall between polynomial and arbitrary constants

**Source anchors:** `04-deformation.tex:26–32,271–311,748–812`; `05-apriori.tex:13–28,427–536`; `07-mass-limit.tex:225–318`.

The source makes reciprocal-polynomial coercivity explicit by completing squares. It uses `ell=exp(-epsilon*n)`, `tau=ell^(3/2)` so exponential smallness defeats proved polynomial losses. Arbitrary `C(n)` from fixed-parameter regularity is used only after the relevant polynomial absorption and is excluded from the final flux-loss estimate.

**Possible arithmetic use:** audit denominator growth, sizes of basis-change matrices, minors, occupancy multiplicities and analytic remainder bounds in one parameter. An exponential approximation estimate is useful only if all arithmetic costs entering the *same integer* have proved compatible growth.

**Required premises:** explicit uniform growth estimates in a common parameter, including denominator clearing and specialization; an independent nonzero integral witness; a sufficient relation between parameter and rational denominator `q`.

**Missing bridge:** no polynomial loss bound for a pi construction follows from the paper. An unspecified fixed-parameter constant may overwhelm the useful error. A quadratic occupancy bound does not by itself give rank, determinant integrality, or the necessary determinant-size inequality.

### Fixed competitor class before limiting optimization

**Source anchors:** `neutral/setup.tex:59–87`; `02-exterior-reduction.tex:314–401`; `neutral-bridge.tex:133–198`.

The paper counts all obstacle-coincident boundary pieces and forces **all** minimizing enclosures into one fixed compact set before moving the metric to a limit. In the many-end passage it proves the correct inequality between distinct cut classes instead of declaring their infima equal.

**Possible arithmetic use:** retain a single candidate class and determinant normalization through auxiliary transformations. If admissible rational/integer bases change, prove the needed inclusion/comparison and show how nonvanishing survives.

**Required premises:** an exact common arithmetic domain or proven class inclusions, tight height/denominator controls, and preservation of the same witness.

**Missing bridge:** perimeter compactness is not compactness of integer bases; real convergence does not preserve nonzero integer minors under unbounded denominator clearing. No common basis or `Q(Z)` row-rank statement is supplied.

### Define the construction globally before restricting to a good set

**Source anchors:** `06-elliptic-existence.tex:654–818,897–991`.

The scalar solve is defined on every bounded input set without assuming the floor. The source then defines a compact coupled map and excludes fixed points on the floor and norm boundaries before invoking degree on the moving admissible set.

**Possible arithmetic use:** prevent circular constructions that define an auxiliary map only for already-good rational approximations. A real common domain and an independently proved boundary-exclusion lemma may be useful if a genuine arithmetic encoding can be found.

**Required premises:** well-defined continuous map on a common domain, compactness and degree hypotheses, and a separate theorem converting its output into rational/integer approximation data.

**Missing bridge:** topological degree supplies real PDE solutions, not rational/integer outputs or exact determinant nonvanishing. No such conversion has been identified for pi.

### Do not infer original-data rigidity from a numerical deformation

**Source anchors:** `00-introduction.tex:139–163`; `figures/proof-routes.tex:7–20`; `neutral-bridge.tex:201–205`.

The numerical path controls inequalities and limits. The equality route begins afresh with original data under additional horizon assumptions. Its existence is described in the overview; the detailed rigidity modules belong to delegated readers.

**Possible arithmetic use:** separate analytic smallness from arithmetic nonvanishing. A deformation of a witness may establish a bound while destroying the original nonzero minor; one cannot transfer the latter backward without an exact identity.

**Missing bridge:** this is a scope safeguard, not a new pi theorem. Both exact nonzero integer determinant and same-integer analytic smallness remain explicit open obligations.

## Limits of this receipt

No cited dependency's proof was read or verified solely by reading its bibliography entry. The newly stated axial regularity, natural-growth oscillation, weak-support smoothing, boundary reflection, compact inverse and degree arguments were read as source arguments; no correctness certification is implied. No search for later validation or external theorem review was part of this reading assignment. Numerical, equality, charged, many-end and local-branch statements retain their distinct premises. There is no publishing step.
