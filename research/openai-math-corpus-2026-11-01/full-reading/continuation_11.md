# 顺序46：四维解析模块续读11

上游固定版本：`openai/math@4eaee87548250b2682affda7d05d947c05554bae`。本续读逐行读完下列8文件，共5339行、240022字节；所有Git blob SHA及行数与缓存source_closure.json相符。阅读不等于独立验算，也未核验外部引文。

| 源文件 | SHA | 实际完整阅读区间 |
|---|---|---|
| `build/sections/four/setup.tex` | `aba8103cf4256bc89a7c5dd0f36c28ed1a2117dd` | 1–252 |
| `build/sections/four/sections/01-end.tex` | `acf7ad73829c0756f49a9dabb2275cff30ef2ceb` | 1–280, 281–560, 561–830, 831–1030 |
| `build/sections/four/sections/02-geometry.tex` | `034961cbc943748393c074d9da5ccce90c870ee3` | 1–290, 291–633 |
| `build/sections/four/sections/03-system.tex` | `78abc670bd9d28e4ae69ea0ad00d1b5aa3440f1c` | 1–285, 286–570, 571–785 |
| `build/sections/four/sections/04-bounds.tex` | `5c4fcb79019dec1ae5d022609548a887800c1dba` | 1–290, 291–560, 561–760 |
| `build/sections/four/sections/05-regularity.tex` | `32b7bc368e811e80f52ab6d7727e135e5a6a797f` | 1–300, 301–615, 616–975 |
| `build/sections/four/sections/06-existence.tex` | `39bb7775ee80b21d5c0b76b71d1140baf973c069` | 1–300, 301–600, 601–875 |
| `build/sections/four/figures/two-frontiers.tex` | `84ecc9b2d3430c85d47789c1141aead04e18c782` | 1–29 |

## 稿件主张与证明结构

The manuscript claims E_e≥sqrt(|P_e|²+(1/4)(A_e/ω3)^(4/3)) for smooth orientable4-dimensional initial data, arbitrary K, DEC, integrable densities, finitely many AF ends with q>1 and the common future weakly trapped boundary sign. Setup also states an original-data rigidity equivalence for the connected outermost MOTS and outer area-minimizing horizon class; that proof lies in separately assigned rigidity files.

- Reduce each end to strict one-ended prepared rest data, preserving DEC at vanishing energy cost and comparing every entire smooth cut by a proper compression; a positive energy shell treats non-timelike charges.
- Construct total regions for K and −K, choose right-continuity thresholds before collar amplitudes, and budget the trace allowance from strict DEC slack.
- Use f,Z,t graph/conformal variables with three exactly unit transverse eigenvalues. Retain −Fτ in the scalar identity, prove polynomial coercivity and coarse-height floor exclusion, then sharp height bounds and polynomial weighted traces.
- Get bounded range, fixed-n axial scalar regularity and a measure correction for the second equation, then alternate scalar Schauder estimates. Solve scalar trial problems on the whole trial space and use a compact frozen return map with moving-section degree.
- Exhaust at fixed n, prove end decay and finite ADM flux, estimate energy with polynomial losses only, invoke the separately proved full-perimeter enclosure and cited Bray–Lee inequality, then send n→∞, ε↓0 and undo earlier reductions.

## 条件性构造线索

### 1. Use an exact first-order identification so an uncontrolled higher-order quantity never enters the error.

**前提：** A compact closed set has the expansion inequality for every smooth exterior support; smooth background metric and tensor have uniform compact-neighborhood bounds; offset distance z is below the injectivity scale. No support-curvature bound is assumed.

**机制：** Solve the Jacobi two-endpoint problem J_X(0)=X, J_X(z)=P_zX. It gives a unit field and a local map F_z with dF_z(x)=P_z exactly and |∇dF_z(x)|=O(z). The Hessian chain rule then contracts the unknown support Hessian with precisely the same tangent metric at both contacts, leaving only O(z) background errors. Approximate metric alignment would multiply the uncontrolled Hessian.

**可能桥梁：** Conditional arithmetic interface: seek an exact rational or unimodular coordinate identity preserving the very same integer determinant. Align all row blocks and the one selected common basis exactly before estimating residual coefficients; this may prevent a height perturbation from multiplying an uncontrolled occupancy or derivative term.

**缺失步骤：** No such arithmetic transport is constructed. Real parallel transport need not preserve lattices, integer determinants, column heights or nonadditive quadratic occupancy. Required: exact row rank over Q(Z), one common basis, explicit lattice/denominator and full-cost control, and analytic smallness for that same nonzero integer witness. Separate small minors or separate block bases are insufficient.

**源码锚点：** `build/sections/four/sections/02-geometry.tex:279` (`four:geo:support-transport`); `build/sections/four/sections/02-geometry.tex:349` (`four:geo:exact-isometry`); `build/sections/four/sections/02-geometry.tex:373` (`four:geo:mean-curvature-transport`)

### 2. Cancel a rough coefficient structurally before estimating derivatives or asymptotic size.

**前提：** In four dimensions A_z=I−(1−χ̃)e_z⊗e_z has three transverse eigenvalues exactly one and axial eigenvalue χ̃≥χ0>0; e_z is the gradient direction, gradient and scalar source G are bounded, background is smooth. χ̃ and G may be merely measurable.

**机制：** A_z∇(|P|²/2)−GP=(H_z−Δz g)P cancels the measurable axial factor. Divergence cancels third derivatives and yields |H_z|²−(Δz)²+Ric(P,P), coercive up to C|G|². A rank-one Frobenius defect below one gives L² absorption and some p0>2. Gradient-drop/affine alternatives and a fixed-axis correction produce Hölder gradient and Hessian Morrey bounds without derivatives of χ̃ or G.

**可能桥梁：** Conditional construction search: find an exact syzygy or antisymmetrized combination of candidate rows eliminating an uncontrolled scale-dependent coefficient before estimating it. If all blocks share one rational basis and the residual unknown direction has a quantified coercive gap, an analytic residual bound might avoid a modulus or derivative bound for that coefficient.

**缺失步骤：** The PDE axial decomposition does not establish rational row independence. Needed: an explicit exact row identity in Q(Z), a common basis, a quantitative arithmetic coercive gap, and complete transformed column-height plus quadratic occupancy cost. No regularity of an unspecified π-dependent coefficient may be assumed. Cancellation alone does not force the same integer determinant below one.

**源码锚点：** `build/sections/four/sections/05-regularity.tex:16` (`four:thm:axial-regularity`); `build/sections/four/sections/05-regularity.tex:82` (`four:reg:defect`); `build/sections/four/sections/05-regularity.tex:145` (`four:reg:flux-identity`); `build/sections/four/sections/05-regularity.tex:148` (`four:reg:flux-divergence`); `build/sections/four/sections/05-regularity.tex:279` (`four:reg:fixed-axis-operator`)

### 3. Separate existence constants from the final asymptotic budget and avoid using the sharp estimate circularly to obtain admissibility.

**前提：** Prepared geometry, strict DEC slack, thresholds, collars and ε are fixed before n; ℓ=e^(−εn), η_n=ℓ^(3/2); the floor penalty and weighted boundary estimates have polynomial losses Π_n. Fixed-n regularity constants may be arbitrarily large. Exhaustion R→∞ is at fixed n.

**机制：** First-floor exclusion pays for unrestricted −Fτ using only a coarse height bound and polynomial penalty. Later height separation permits strict-slack absorption in the final scalar identity. Final energy error uses only Π_n(ℓ+ℓ²/η_n)→0; arbitrary fixed-n existence/Schauder constants never enter it. No limiting metric as n→∞ is required.

**可能桥梁：** Conditional arithmetic workflow: establish an admissible family and a common nonzero integer witness using coarse exact information first. Then derive sharp analytic determinant smallness with a dependency graph tracking rank, derivative order, transformed column heights and the full nonadditive occupancy objective. Only controlled losses may enter the growing-parameter limit.

**缺失步骤：** No arithmetic analogue of the polynomial weighted-flux bound is proved. A formal exponential factor cannot dominate an arbitrary hidden constant depending on the same growing parameter. The π target still needs actual common witnesses, certified nonvanishing, controlled full cost and a smallness bound in the stated order of limits.

**源码锚点：** `build/sections/four/setup.tex:243` (`four:rem:constants`); `build/sections/four/sections/03-system.tex:695` (`four:prop:floor-exclusion`); `build/sections/four/sections/06-existence.tex:522` (`four:prop:mass-comparison`); `build/sections/four/sections/06-existence.tex:610` (`four:ex:penalty-absorption`)

### 4. Compare every full feasible object under one proper map without requiring attainment of an infimum.

**前提：** The rest-end replacement agrees with original data in a growing core, has a fixed positive metric lower bound remotely, and admits a proper radial compression fixing the compact boundary. Every cut counts the entire intrinsic frontier, including disconnected and coincident pieces.

**机制：** A constructed diffeomorphism Φ_R satisfies Φ_R* g≤(1+ε_R)g_R globally, ε_R→0. It sends every new admissible outer domain to an original admissible outer domain. Taking infima directly gives A*(g_R)≥(1−o(1))A*(g), without compactness or attainment of a minimizing sequence.

**可能桥梁：** Conditional common-basis comparison: map every feasible full-rank configuration to an admissible comparison configuration and prove a single inequality for its complete column-height plus quadratic occupancy cost. This could avoid replacing a common integer configuration with a relaxed optimum or independently minimized block costs.

**缺失步骤：** A reduction can delete the only common full-rank basis, and quadratic occupancy is not additive across blocks. A valid arithmetic map must preserve the actual common witness, row rank over Q(Z), integer nonvanishing and the full cost with explicit constants. A proper geometric map provides none of these algebraic facts automatically.

**源码锚点：** `build/sections/four/sections/01-end.tex:436` (`four:prop:rest-end`); `build/sections/four/sections/01-end.tex:773` (`four:end:compression-metric`)

## 范围限制

- Not independent proof verification. External citations, including total-region, regularity and Riemannian comparison inputs, were not checked externally.
- The minimal enclosure imports the separately assigned four/rigidity/01-enclosures.tex proof, which this reader has not read.
- Four-dimensional equality/rigidity was read only as stated in setup; its proof body is outside this receipt.
- No arithmetic construction, common-basis certificate, row rank over Q(Z), or full height-plus-quadratic-occupancy bound is proved here.
- Whether π is badly approximable and the target behavior of q²|π−p/q| remain unresolved. These are conditional interfaces, not transfer theorems or relevance filters.

Actual untruncated reading of exactly8 four-dimensional analytic-module files:5339 source lines,240022 bytes. Git blob SHAs and line totals were checked against source_closure.json. This receipt does not claim other position46 bodies were read. Setup includes an equality theorem statement, but its rigidity proof and the full-perimeter enclosure proof are assigned elsewhere. Actual reading is not independent validation of the manuscript or its cited dependencies.
