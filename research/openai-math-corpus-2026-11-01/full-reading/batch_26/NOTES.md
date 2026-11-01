# 第54项：严格逆一次幂界的正文首读

上游固定 `4eaee87548250b2682affda7d05d947c05554bae`。完整阅读11个文本源文件、1281行；这是阅读记录，并非独立正确性认证。

主张与架构：The manuscript claims universal C and ε>0 with M_{-1}[f′](r)<=C(1-r)^(-1/4+ε) for all normalized univalent disk maps and r>=1/2; hence B_b(-1)<1/4 and the stated Kraetzer formula fails at -1. It gives no numerical ε and no conclusion at +1.

- Compactness/distortion yield a single positive dyadic operator L whose spectral growth β controls all disk means.
- Assume β>=1/4; resolvent and affine averages give an exact inverse-first-derivative weighted scaling law. Differentiating gives two logarithmic derivative moments.
- For k=(β+2)/3 the target potential y^k/|F-ξ| has compact positive superlevels almost everywhere. Critical points are fibers of G=F-iyF′/k, with saddle/maximum sign determined by its normalized Jacobian.
- An index count dominates saddle counts by finite maxima above every positive threshold. Generic absence of ties upgrades equal-rank pairing to strict increase, measurably and equivariantly.
- Change of variables matches the cubed critical value exactly to transported weights. Expanding-box limits show E[J]<=0, with equality only when J(i)=0 almost surely. Weighted covariance handles the fixed target.
- The moment formula excludes β>1/4. At equality, the law and a countable dense source set force J≡0; the nonvanishing pre-Schwarzian then obeys an impossible harmonic lower barrier. β<1/4 yields one uniform exponent bound.

## 按原条件保留的构造机制

### Use strict improvement only on positive witness mass, then classify equality separately; a common pointwise gap is unnecessary.

前提：A measurable injective same-fiber pairing, exact target weights, finite signed mass, equivariance, and a proved passage from finite boxes to the invariant law. The ratio c is strictly greater than one on the positive-Jacobian domain; generic bad fibers carry zero positive weighted mass.

机制：Area change of variables converts the cube of a critical value into exactly the partner weight. The comparison gives b<=A<=d. If b=d, the nonnegative integral of J_+(c-1) vanishes and hence b=d=0, even when inf(c-1)=0. A separate harmonic barrier rules out the resulting global zero Jacobian.

候选接口：Speculative interface: in an extremal law on admissible row/common-basis configurations, seek a pairing which preserves the same nonzero minor and makes the full height-plus-occupancy score strictly better wherever a certified defect occurs. Separate equality rigidity from quantitative gain. This could avoid demanding a uniform improvement for every configuration.

缺口：No such invariant arithmetic law or injective score-improving common-basis move is known here. Improvement of an averaged surrogate cannot replace a nonzero and small estimate for one normalized integer determinant. The equality classification and preservation of actual row/column constraints must be proved afresh.

- `preprints/A-strict-inverse-first-power-bound-for-univalent-functions-September-24-2026/build/sections/04-strict-pairing.tex:211` — `prop:pair`
- `preprints/A-strict-inverse-first-power-bound-for-univalent-functions-September-24-2026/build/sections/05-strict-transport.tex:55` — `prop:strict-transport`
- `preprints/A-strict-inverse-first-power-bound-for-univalent-functions-September-24-2026/build/sections/06-equality-barrier.tex:8` — `lem:no-zero-J`

### Promote target-genericity to a fixed distinguished target through weighted covariance, rather than assuming a generic target includes the chosen point.

前提：Bad targets have planar measure zero for almost every map, the local change-of-variables formula includes J_+, and rebasing covariance and a strictly positive scaling weight apply to the same Borel event. A finite spatial box has positive dx dy/y mass.

机制：Although G_F(i) is always the single target -i/k, spatial integration kills bad targets with positive Jacobian; the affine law transports that zero integral to the base point. The conclusion is weighted goodness E[J_+(i)1_bad]=0, not unweighted almost-sure goodness.

候选接口：Speculative interface: when generic nonzero minors over Q(Z) are used at the fixed transcendental point 2πi, explicitly identify the algebraic or weighted transport which carries the generic certificate to the chosen point. Keep the transported weight and null event attached to the same witness.

缺口：The project already has transcendence-based rank transfer for polynomial entries, but a generic measurable or analytic exceptional set need not avoid 2πi. The source covariance is geometric, not an available arithmetic action. It supplies neither a height estimate nor an automatic fixed-point principle.

- `preprints/A-strict-inverse-first-power-bound-for-univalent-functions-September-24-2026/build/sections/05-strict-transport.tex:12` — `eq:area-formula`
- `preprints/A-strict-inverse-first-power-bound-for-univalent-functions-September-24-2026/build/sections/05-strict-transport.tex:27` — `lem:base-good`

### Preserve the exact growth parameter during compactification before attempting a strict endpoint contradiction.

前提：Compact normalized map family; bounded positive operators with a multiplicative affine action; spectral growth rho>1; endpoint cancellation whose error decays as rho^-n; finite positive scale-average mass.

机制：A resolvent eigenmeasure is horizontally averaged and then logarithmically scale-averaged. The resulting law retains y^-β for β=log2 rho itself. Differentiation gives exact moments, and E[J]=(4β-1)(β+2)/(9k²) exposes the equality threshold.

候选接口：Speculative interface: if row configurations admit a positive transfer operator respecting true rank witnesses and full costs, carry the actual growth exponent into a compact extremal law. Replacing it by an arbitrarily larger exponent too early could destroy the equality obstruction needed for a strict result.

缺口：A normalized compact state space and a positive action preserving a nonzero integer determinant have not been constructed. Integer heights and quadratic occupancy are not known to be multiplicative cocycles. Even an exact spectral law would still need a π-specific obstruction and a return to a finite witness.

- `preprints/A-strict-inverse-first-power-bound-for-univalent-functions-September-24-2026/build/sections/01-compact-model.tex:106` — `eq:rho`
- `preprints/A-strict-inverse-first-power-bound-for-univalent-functions-September-24-2026/build/sections/02-weighted-law.tex:11` — `prop:law`
- `preprints/A-strict-inverse-first-power-bound-for-univalent-functions-September-24-2026/build/sections/02-weighted-law.tex:108` — `lem:moments`
- `preprints/A-strict-inverse-first-power-bound-for-univalent-functions-September-24-2026/build/sections/04-strict-pairing.tex:18` — `eq:Jacobian`

## 量词和边界

- Complete source first reading is not independent theorem verification; cited classical estimates and external papers were not independently re-proved or read here.
- The ε and C are existential, not numerically evaluated. No uniform pointwise pairing gap is asserted.
- The fixed-base conclusion is weighted by J_+; unweighted target genericity cannot be substituted. The expanding-box limit order S, thenY, thenN is part of the proof.
- No automatic transfer of geometric invariant laws, Jacobian signs or harmonic barriers to arithmetic determinants is established; no π BA/nonBA conclusion follows.

Read every line of the pinned README, build/main.tex, all seven included sections, recursively included rank-pairing TikZ source/caption and full bibliography (11 exact SHA-verified text files, 1281 lines). The only remaining inputs are standard TeX glyph-to-Unicode infrastructure; there is no mathematical appendix or further include in the complete directory tree. Figure is schematic with illustrative gaps; no PDF page/layout verification or external cited-proof validation is claimed.
