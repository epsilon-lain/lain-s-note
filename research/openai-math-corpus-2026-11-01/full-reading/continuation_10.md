# 第46项续读：带电模块与纯电多端延拓

阅读范围：2112行；12个实际读过的源码文件。Complete actual reading of the charged module and its electric finite-end extension, not the rest of the large manuscript. All 12 receipt files match the pinned source SHA and line totals. Charged final numerical theorem is imported from companion ChargedOriginalData2026, not independently proved here. Charge preservation concerns total charges, not individual component charges; mass-lower-bound inversion requires r>=|Q|. No global Lorentz slicing or limiting metric as n→∞ is claimed.

## Preserve the observable exactly while eliminating only a direction proved invisible to it.

前提：Two closed flux forms with a fixed total charge vector, a constant SO(2) rotation, and a pointwise graph-norm square completion. The intended observable is the total charge across the full enclosing cut, not separate charges on every component.

机制：Rotate total charges to (Q,0), keep the first rotated form and discard the second zero-total-charge form, then inverse-rotate. Both original total charges and closedness survive exactly, while the discarded square is nonnegative. This is sufficient for the Riemannian comparison and does not preserve all local fields or component charges.

候选接口：Speculative interface: compress row/column data only after proving the discarded directions lie in the kernel of the actual determinant/coefficient observable. Distinguish preservation of one full common constraint from preservation of each marginal component. Track the positive form under the same transformation.

缺口：A real orthogonal rotation generally destroys integer lattices and coefficient heights. A common integer determinant need not factor through a total-charge-like linear observable. A usable analogue requires rational or unimodular maps, an exact factorization identity, and explicit denominator and occupancy costs.

- `build/sections/charged/sections/06-comparison.tex:134` — `ch:cmp:charged-comparison`
- `build/sections/charged/sections/06-comparison.tex:186` — `ch:cmp:field-square`
- `build/sections/charged/sections/06-comparison.tex:204` — `ch:cmp:new-flux-forms`

## Keep error estimates uniform only in parameters that are actually sent to infinity at that stage.

前提：Strict prepared data, collars and weights fixed first; fixed ε, then n; exhaustion R→∞ at fixed n; polynomial penalty/coercivity and an explicit integrated flux bound independent of possibly large fixed-n classical regularity constants.

机制：The energy error is bounded by C(1+n)^N(ℓ+ℓ²/τ0), with ℓ=e^(-εn) and τ0=ℓ^(3/2), so polynomial losses are dominated by exponential smallness. Only numerical bounds are passed through n→∞, ε↓0 and the final data approximation; no n-limit metric is asserted.

候选接口：Speculative interface: record separately polynomial losses in rank/differentiation and arbitrary constants attached to a fixed configuration. A determinant contradiction needs a final smallness estimate in which only controlled losses occur; avoid accidentally inserting a fixed-n existence bound into a growing-n asymptotic.

缺口：The project has not proved a comparable arithmetic flux/trace estimate independent of all uncontrolled constants. An exponential symbol is useful only with justified parameter order and a preserved same integer witness. No source PDE existence theorem supplies that arithmetic estimate.

- `build/sections/charged/sections/05-elliptic.tex:1` — `ch:ell:section`
- `build/sections/charged/sections/07-mass.tex:77` — `ch:mass:vanishing-error`

## Compare whole admissible competitor classes directly, so neither infimum attainment nor a geometric limiting object is required.

前提：A chosen end, full smooth enclosing cuts with coincident obstacle frontier counted, finitely many ends, conserved flux and the monotone branch r>=|Q|.

机制：Removing distant tails creates strictly trapped inner spheres. Every retained one-ended full cut is also a designated-end cut of the original manifold, so a_R>=a0 directly. Monotonicity of (r+Q²/r)/2 on the stated branch transfers the bound without taking the tail radii to infinity.

候选接口：Speculative interface: prove inclusion of actual feasible common-basis/minor certificate sets under a reduction before comparing optimal full costs. A direct set comparison can avoid relying on a limit or relaxed optimum being attained by an integer configuration.

缺口：Deleting supports can remove every full-rank common basis. Inclusion and witness preservation must be proved for the real candidate sets, and the quadratic occupancy objective must be monotone in the appropriate direction. Separate per-block optima cannot replace one common optimum.

- `build/sections/charged-extensions.tex:124` — `chext:thm:charged-electric-finite-ends`
- `build/sections/charged-extensions.tex:207` — `chext:eq:charged-area-under-tail-removal`

