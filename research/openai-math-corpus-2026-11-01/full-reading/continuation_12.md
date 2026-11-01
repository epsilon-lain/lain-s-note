# 第46项续读：局部AdS模块

七节逐行读完，共2974行。Complete source reading of the seven AdS module sections, including full analytic expansion/flux justification, CMC sphere-jet recursion, TT square and complete kernel, quartic null-transport/tangency obstruction, inline schematic caption and exact radial conclusion. Exact SHA and line counts match the lead closure. Other position46 modules are covered by separate receipts; external cited papers and independent correctness are not certified.

## Classify the full null space of a nonnegative quadratic deficit before testing its first surviving higher-order term.

前提：A continuous nonnegative rotation-invariant quadratic form on the full real TT seed space; finite-angular approximations preserving TT and the closed kernel; a sharp TT boundary square; actual controlled expansions for a given branch.

机制：The quadratic kernel is exactly four-dimensional (radial plus degree-one). The degree-one directions are genuine null directions even with mixed seeds meeting the boundary sign. They cannot be discarded by the sign assumption. A fourth-order square and a local tangency obstruction force c4>0 for every nonradial null direction; radial seeds instead have an exact conserved-mass identity.

候选接口：Speculative interface: if a determinant-cost deficit has a nonnegative leading form, identify every degenerate family and expand the SAME certificate on that family to its first surviving order. Separate exact symmetry/equality families from obstructed directions. This may be more effective than seeking a uniform positive quadratic gap.

缺口：No such differentiable deficit on integer common-basis configurations is established. A zero quadratic coefficient does not by itself imply a finite-dimensional kernel or positivity at fourth order. Expansions must preserve exact rank, integer normalization and full occupancy, with usable remainder estimates.

- `build/sections/ads/sections/04-quadratic.tex:16` — `ads:lem:tt-square`
- `build/sections/ads/sections/04-quadratic.tex:212` — `ads:cor:kernel`
- `build/sections/ads/sections/05-quartic.tex:510` — `ads:quart:square`
- `build/sections/ads/sections/06-conclusion.tex:5` — `ads:prop:radial`

## Convert a boundary/asymptotic observable to a bounded exact integral functional before passing to approximation or Taylor limits.

前提：Each actual member has the stated C^{2,α}_τ decay with τ>3/2 and finite mass flux. The chosen weight β satisfies 3/2<β<min(τ,√3). The Green identity is first taken at fixed member; source products have weight2β and the functional is bounded in the three stated norms.

机制：An exact Green representation computes all four mass components from u,F,B. Uniform weighted function/source/boundary remainders then yield actual mass expansions without asserting uniform convergence of discarded derivative flux errors. Finite angular coefficients are also identified with their own flux before taking angular approximation limits.

候选接口：Speculative interface: for determinant evaluation, seek an exact coefficient/integral representation continuous in norms whose approximation errors are controlled. Keep evaluation at the actual π parameter and integer normalization in the identity before taking asymptotic limits.

缺口：The arithmetic determinant observable and its normalization are not represented by the source Green operator. Uniform functional bounds must be proved for the actual row space and increasing dimension; a fixed-dimension continuity assertion gives no needed height bound.

- `build/sections/ads/sections/02-analysis.tex:266` — `ads:lem:green`
- `build/sections/ads/sections/02-analysis.tex:433` — `ads:prop:expansions`

## Eliminate auxiliary local generators to obtain a global obstruction even when the generators do not agree on overlaps.

前提：Equality in a positive quartic integral forces mean-free lapse and tracefree second form to vanish through second order. Round-jet decomposition yields local Lie-derivative-plus-trace representatives on every compact radial interval. A separated decaying degree-one mode solves a fixed homogeneous ODE.

机制：Project the necessary tangency onto degree-two scalar,gradient,Hessian modes; solve for three local generator coefficients and eliminate them. The remaining identity (V′−NV/r)²=6m0 V²/r³ holds globally independently of choices. For z=V/r, its integrable relative derivative bound and decay at infinity force z≡0 by backward Gronwall.

候选接口：Speculative interface: use exact elimination of local row-reduction choices to find an invariant polynomial identity for a putative extremizer, avoiding unjustified global gluing of locally selected bases. Combine the invariant with a genuine boundary/asymptotic restriction.

缺口：Local compatible bases need not exist, and elimination can lose the nonzero minor or introduce uncontrolled denominators. There is no π-specific integrable differential constraint here; the source obstruction cannot be transplanted merely by naming a deficit.

- `build/sections/ads/sections/05-quartic.tex:758` — `ads:lem:dipole-obstruction`
- `build/sections/ads/sections/05-quartic.tex:812` — `ads:quart:dipole-identity`

## 范围边界

- Reading is not independent theorem verification.
- The theorem is local for each fixed TT seed and prescribed decaying branch; no branch existence for every seed or uniform neighborhood is asserted.
- Finite Taylor geometry is a jet identity, not an actual global CMC foliation or Einstein extension; boundary displacements of either sign are evaluated at the level of justified jets.
- No arithmetic transfer or π BA/nonBA result is proved.
