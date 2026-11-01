# Continuation 07 — Uniform bounds for planar polynomial limit cycles

Pinned source: `openai/math@4eaee87548250b2682affda7d05d947c05554bae`.

**Assigned parts fully read: 5,595 lines.** All lines of sections 03–06 and the three assigned figure sources were read in contiguous output chunks. Seven cached Git blob hashes match the pinned source records. This is a reading record, not independent proof verification. Whole-paper closure must be merged with the parent reader.

| Source | Full read range | Git blob SHA |
|---|---:|---|
| `build/sections/03-separation.tex` | 1–1428 | `32c6ea9d1301e6e348a10b9769d3407ae2e54d14` |
| `build/sections/04-calculus.tex` | 1–1144 | `5046d23f5d0af4795d858c3ceacdc68a7b11b85d` |
| `build/sections/05-elimination.tex` | 1–1147 | `f5e785e410ec7fbd15aa6a7a32f36daa459593b7` |
| `build/sections/06-additive.tex` | 1–1748 | `949baad6174a517369f2906a1e793f881a9a0dfa` |
| `build/figures/coupled-bands.tex` | 1–50 | `0c2c1fcdb2b7b1047d5a653ddfd5f3c0384d3b6b` |
| `build/figures/matching-arguments.tex` | 1–46 | `6c2e615ac8e0f1ff57bd70d7785406571aebdb10` |
| `build/figures/retained-charges.tex` | 1–32 | `9bb12dde52d6c214b21b39a28cd8f3878bc4fa2f` |

No assigned ranges remain unread. The JSON records the exact tool-output chunk ranges. `batch_02/records.json` was not changed.

## Proof mechanisms

### 1. Zero formal subtrees are upgraded before earlier exponential multipliers are applied.

Couple only finitely many selected prefixes through fixed good joins. Correct later rates to absolute o(1) holomorphic proxies by area Cauchy transforms and propagate the next-band growth through the replaced inner wedge. For a zero subtree, first obtain |G(u)|<=C_N exp(-N Z_n(u)) separately for every fixed N. Multiply by exp(H Lhat_n-C Zhat_n); reserve half of every negative source cost, exp(Hx-cU)<=exp(K_H)exp(-cU/2), so the Cauchy-correction norm has an H-independent prefactor. A two-stage maximum principle first proves boundedness using the preliminary real decay, then treats the real ray as interior and removes the arbitrary C_N,u_N. Source-ratio estimates choose H~exp(g)/gprime at bounded slopes, giving Hx_a~V_(n-1); large slopes use a central double-exponential bound followed by a wide-strip harmonic weight. This proves exp(-c V_(n-1)) suppression, which absorbs the finitely many growing exp((a-b)Z_(n-1)) factors when passing earlier zero children. At n=1 there are no sources and H->infinity yields exact zero.

Sources: `build/sections/03-separation.tex:267-500`; `build/sections/03-separation.tex:502-1058`; `build/sections/03-separation.tex:1061-1327`.

Assumptions and limits:
- Requires the complete packet domain, fringe, source, safe-column, independent-test and controlled-family hypotheses; arbitrary real flatness does not imply exact zero.
- The full coupled functions and earlier joins stay fixed while H and preliminary real orders vary. The common full-forward-domain source envelope is essential.
- No convergence of the simultaneous full exponential series and no uniformity over all packet labels or all budgets is claimed.

### 2. Differentiated closure uses prefix-finite coefficient rules, cleanup, and a fixed coherence order distinct from working accuracy.

Finite unions, convolution and strict-positive semigroups of iterated left-finite supports give finite coefficient calculations at successive prefixes. Dependent powers are converted exactly into later exponential shifts. Extracted current W_j=e^(-Z_(j,j)) is never frozen at the next background; only later rates are translated. Differentiation preserves strict positivity, while derivative induction plus zero-subtree suppression justifies deleting negative local terms. Positive expressions have bounded ordinary derivatives and exponentially small free derivatives; real local boundedness proves positivity. Regular implicit closure first solves the real equation, then uses bounded real-symmetric tanh interpolation preserving a finite lower jet while attaining each complex trial target. A broad Taylor tube fixes sheet coherence at a stack-dependent order r_d; a separate arbitrarily higher working order N controls residuals. Piecewise Taylor integration across seams gives O(|v|^N), and a real Jacobian contraction constructs unique smooth band roots. Exact implicit equations, rather than differentiation of a pointwise asymptotic statement, supply derivative/seam/source bounds and the actual positive transition series.

Sources: `build/sections/04-calculus.tex:18-241`; `build/sections/04-calculus.tex:243-652`; `build/sections/04-calculus.tex:654-1144`.

Assumptions and limits:
- Primitive independent derivatives, anti-CR defects, retained-input independence, and controlled-family uniformity are hypotheses verified separately for each primitive library.
- Finite Taylor matching controls geometry but gives no exponential anti-CR bound by itself; chain rules require holomorphic inputs or already controlled sheet defects.
- A fixed complex ordinary disk is not assumed because dependent exponentials can rotate rapidly there. Sheet coherence is local on fixed lifted determinations, not global monodromy independence.

### 3. Ordinary elimination combines exact zeros, visible normalization and negligible-center approximation with fixed conditioning powers.

Promote a vanishing ordinary coordinate as x=epsilon rho with rho=h^nu e^(-B), preserving old free inputs and retaining the replacement degree of freedom. Moderate monomial logarithmic sensitivities are only polylogarithmic in 1/rho, so current-background feedback rho*d_x log rho=o(1) gives a regular equation. Reverse Taylor replay approximates center recipe values by finite expressions at current backgrounds. The completed-map Schur complement includes record feedback K=R_xi-R_q Y_q^(-1)Y_xi: the promotion factor is L_xi+L_R K, and the physical complement is w/(1-wL_w) times that factor. An induction in ordinary dimension and directional vanishing order has three simultaneous outputs: exact-zero chart, visible F=tau U, or full-dimensional |F|<=r^M with arbitrarily accurate records. At a derivative center w=d_z^(m-1)G=0, the nonnegative H=sum c_l^(2L/(m-l)) controls all lower coefficients. A visible center and w/s bounded yield a degree-m limiting polynomial with absent degree-(m-1) term and nonzero lower part, hence smaller zero multiplicity. A negligible center is lifted at sufficiently high accuracy only after map/inverse/neighborhood exponents are fixed. Common-domain lemmas keep these exponents independent of later record accuracy. Retestability then turns a recoverable surviving free-input perturbation into a second zero, contradicting absolute isolation.

Sources: `build/sections/05-elimination.tex:41-590`; `build/sections/05-elimination.tex:591-1029`; `build/sections/05-elimination.tex:1031-1147`.

Assumptions and limits:
- Restricting to a derivative center alone does not preserve a zero; the visible/negligible alternatives and exact hit correction repair this issue.
- Moderation applies only to specified two-sided moderate monomials. Much smaller positive expressions are kept whole; their large factors cannot be separated casually.
- Absolute isolated-zero finiteness requires every escaping sequence to admit the stated analytic or retestable-flag description, and all coefficient parameters, auxiliary coordinates, strict domains and graph ties remain part of the system. It is not automatically a uniform fiber bound.

### 4. Actual additive passages retain endpoint data that cancel lateral ambiguities.

Use the exact action A(w)=w+delta*q*log(w/q), A(p)=Q and P=exp(-(Q-q)). The actual passage is first continued along a polygon in the A-plane with increasing real part and bounded integral slow/fast profiles, before auxiliary sectorial normalizations are chosen. In model A, a graph plus convergent fiber linearization gives the exact transfer T(w_out, P exp(I(out)-I(in)) T^(-1)(w_in,v_in)); keep its actual inner endpoint kernels. The example V(p)+P(v_in-V(q)) shows why adding K exp(-w) to a lateral particular solution cancels only when both endpoints use the same determination. In model B remove only nonzero imbalance Delta=b_f-a, while balanced charge a=b_f=j>=1 produces the surviving P^j drift. Endpoint fast profiles and finite Taylor-Picard transfer entries construct finite charge coefficients with error C_K |Q|^(C_K)|P|^(K+1)+C_K exp(-c_K|Q|).

Sources: `build/sections/06-additive.tex:26-189`; `build/sections/06-additive.tex:191-283`; `build/sections/06-additive.tex:1065-1207`.

Assumptions and limits:
- Inputs are independent passage arguments before physical graph ties; Q/q->infinity, q>=q_*>1, small amplitudes on larger analytic boxes, rational profile powers, and a fixed discrete positive charge lattice are required.
- Both endpoints must retain consistent lateral information before any retirement. Independently replacing exact inner kernels by unrelated formal series can erase the cancellation.
- P is extracted as the original exact affine exponential exp(-n(Q-q)); it is not Taylor-altered through endpoint ratio corrections.

### 5. A finite charge request is a compatible quotient of convergent fixed-side tails, with balanced survivors and explicit Gevrey losses.

Retain min(a,b_f)<=K, a finite union of fixed-side tails rather than a finite total-degree triangle. At each fixed left multi-degree, solve a finite opposite-degree list then the high tail; there Delta is comparable to opposite degree and the Volterra inverse gains C/N. The pure-side multiplier has no constant fast term, so convolution raises opposite degree and gives a common fast radius independent of the fixed side/K, without uniform coefficient bounds in that side. Summable reservoirs reserve state, profile and spatial margins before slow order grows. Truncated norms at theta=epsilon/(J+1) control hD_X and h^2 d_h and give C^(j+1)j! slow coefficients. A coefficientwise weak forcing norm divided by 1+|n| is mapped to the strong tail norm because (1+|n|)/|Delta|<=C; this avoids a hidden fast-radius loss in optimal truncation. Normalized triangular primitives invert rational c+d_L, integrate resonances with zero constant, and keep fixed log degree. Their weighted residual identities and inward integral estimates control inner kernels uniformly in Q/q. In the omitted quadrant both sides have charge>K; selecting bounded subcharges gives P^(K+1) suppression with explicit polynomial Q losses.

Sources: `build/sections/06-additive.tex:285-727`; `build/sections/06-additive.tex:729-1063`; `build/sections/06-additive.tex:1095-1171`; `build/figures/retained-charges.tex:1-32`.

Assumptions and limits:
- The inverse bounds use triangular charge equations and positive generator weights. A common radius is not convergence of the unrestricted two-sided mixed series.
- Every coefficient/word/jet request is finite, and constants can depend on it. No uniformity as charge K or the fixed left degree tends to infinity is asserted.
- The ordering sum outside coefficient suprema in the weighted norm is part of the proof; replacing it by a supremum of the weak sum is not justified.

### 6. Numerical anchors and optimal cutoffs are retired at owning transitions, producing canonical germs without declaring flat errors zero.

Write Q=Q_0(1+Theta_Q), q=q_0(1+Theta_q); let s=sigma(i) be the first index of 1/Q_0, r the first index of log(Q_0/q_0), and p_0=min(l,r)<=s. If l<=r, replace inner kernels by optimal q expansions at l, remove outer 1/Q at s and inner 1/q at sigma(l). If r<l (including bounded q), then r=s: Taylor-remove positive profiles and replace finite zero-profile jets by canonical infinite-anchor jets precisely there; arbitrary Q^(-N) accuracy is enough at this index but would not suffice prematurely. Rectangular Gevrey bounds C^(a+b+1)a!b! preserve independent cutoff precisions. Bin smoothing has only simple source costs exp(Re f_i), exp(Re f_l) while their log formulas remain retained; four angular-size cases show U/|Z_j|->infinity, paying raw, derivative and partition losses. All columns are safe. Formal nonzero-diagonal recursions, normalized resonance constants, compatible finite-charge quotients and canonical infinity jets identify terminal germs under permitted flag refinements and tiny old-free perturbations.

Sources: `build/sections/06-additive.tex:1209-1434`; `build/sections/06-additive.tex:1436-1748`.

Assumptions and limits:
- Finite and infinite primitive normalizations agree to arbitrary fixed Q-power accuracy, not asserted exponential Q-accuracy; the retirement order is chosen to match that exact strength.
- Full-forward controlled-family uniformity is derived from common finite box/jet/rate/gap/slack data; eventual estimates on unrelated paths would not suffice.
- The conclusion is local for selected lifted lateral determinations, with old active free nodes preserved. It is not global history or monodromy invariance.

## Construction ideas for the π determinant project

These are proposed interfaces, not deductions about π. Each still requires a symbolic rank proof, a nonzero integer specialization and a cost bound for the same determinant.

### 1. Use an exact multiscale leading-coefficient certificate for surviving determinant terms, and keep formal rank separate from analytic separation.

For a proposed integer determinant Delta_N built from costed columns and candidate row transversals over Q(Z), expand the exact determinant into signature coefficients before taking absolute values. Put its genuinely occurring scales into a finite lexicographic flag, reduce exact dependencies, and exhibit a first nonzero coefficient by a symbolic minor/wedge or an explicit uncancelled monomial. Record its row-rank certificate and the same selected columns. Only then use scale domination to compare all later signatures and numerical errors. A tree is useful when an infinite lexicographic initial segment would make a flat truncation misleading: work at successive prefixes and pass only finitely many siblings for each requested budget.

Source mechanism: Zero-subtree suppression and finite leading-prefix isolation in section 03; exact current-symbol bookkeeping and cleaned positive trees in section 04.

Required assumptions:
- Actual determinant terms admit an exact finite or iterated left-finite multiscale expansion with controlled coefficients.
- The leading coefficient is proved nonzero over Q(Z) using permitted rows, signature caps and integer normalizations; it must not be inferred merely from formal distinctness.
- If infinite asymptotic tails are used, a relevant separation theorem with complex-domain/source hypotheses is independently proved; numerical smallness to all real orders does not suffice.

**Missing bridge:** Prove an explicit project-specific nonzero minor and surviving coefficient, a controlled integer specialization, and a bound for that exact determinant. The manuscript supplies no arithmetic rank theorem and no pi approximation estimate.

Potential use: A way to isolate where cancellation destroys all purported surviving transversals, while preventing a rank argument for one matrix from being combined with a smallness bound for another.

### 2. Fix the structural determinant cost and conditioning exponents before asking for higher cancellation accuracy.

Design a finite row/column construction with a frozen chronological list of permitted operations, a symbolic rank witness, signature-cap ledger and explicit denominator/height generators. Compute the actual forward/inverse losses and column cost, including quadratic occupancy, on this fixed structure. Then test whether cancellation/approximation order can increase by refining coefficient formulas only, leaving those structural exponents fixed. If yes, select accuracy last so the analytic determinant upper bound beats all fixed losses and falls below 1, while an independently certified nonzero integer determinant has absolute value at least 1. Any parameter correction must use the full derivative including induced record/coefficient feedback, and preserve a power-sized available neighborhood.

Source mechanism: Common domains and deferred-floor elimination in section 05, especially the fixed chronology/prefix invariant and complete Schur feedback formulas.

Required assumptions:
- Increasing accuracy does not add columns/rows, enlarge signature capacities, increase quadratic occupancy exponents or change the selected rank certificate.
- All row operations are unimodular or their exact integer scaling/denominators are included in the same cost ledger.
- The refined actual maps/coefficients have uniform power-domain and derivative bounds; a finite formula having larger denominator multiplicities is harmless only if an actual-map comparison proves the uniform bounds.

**Missing bridge:** Exhibit a concrete arbitrary-accuracy cancellation mechanism with cost exponent independent of its accuracy, and show it acts on the same integer minor. None is supplied here for BA/nonBA pi; the manuscript pattern is an audit criterion, not a completed construction.

Potential use: Detect circular choices where a later approximation order secretly raises the determinant dimension, denominator height or occupancy cost enough to destroy the final inequality.

### 3. Try a two-sided signature-charge filtration that retains all cheap-side survivors rather than truncating by total complexity.

If determinant terms have two additive nonnegative counters a,b measuring opposing families of expensive column choices, retain the union a<=K or b<=K and treat each fixed-side fiber exactly. Look for exact permitted row/column identities that eliminate an imbalance while isolating balanced a=b terms as candidate surviving minors. Bound only the omitted quadrant a>K and b>K through a multiplicative small parameter, with every polynomial/denominator factor explicit. A compatible sequence of such finite quotients could preserve a row-rank witness which total-degree truncation accidentally discards. Test this first on a finite symbolic determinant and its actual signature occupancy table.

Source mechanism: The min(a,b_f)<=K fixed-side quotient, balanced drift, compatible budgets and omitted-quadrant suppression in section 06.

Required assumptions:
- The two counters are truly additive under determinant coefficient products and the omitted quadrant forces a quantifiable factor that pays its costs.
- Fixed-side fibers have finite exact sums or proved convergent common-radius bounds; the ODE Volterra tail inverse cannot be assumed for arbitrary determinants.
- Balanced signatures really contain a nonzero allowed row transversal after cancellations, and their capacity/quadratic occupancy cost is measured for the selected minor.

**Missing bridge:** Construct the actual signature grading and exact elimination identity, prove compatible coefficient budgets and a nonzero balanced minor, then perform integer-height and pi-error accounting. A grading alone neither proves row rank nor guarantees that any term survives.

Potential use: A sharper truncation strategy for zero-term analysis when one side can be arbitrarily large but the other side is bounded, while explicitly exposing the missing balanced-survivor proof.

## Scope limits

- This records complete reading of assigned parts, not independent validation of the manuscript's claimed resolution of uniform planar polynomial limit-cycle bounds.
- Sections 01-02 and 07-11 are outside this continuation and are being handled by the parent or other readers; no full-paper completion claim is made here.
- The separation and implicit-calculus conclusions have strong complex-domain, source, finite-request and controlled-family hypotheses. They cannot be imported as general quasianalyticity for arbitrary determinant asymptotics.
- Elimination keeps all coefficient parameters and auxiliary variables variable; absolute isolated-zero finiteness is not itself a uniform bound for a fixed-parameter fiber. Later counting results were not read here.
- No mechanism proves symbolic rank over Q(Z), low-cost surviving row transversals, integer specialization, denominator/height bounds, or BA/nonBA conclusions for pi without new project-specific arguments.
