# 41 — A counterexample to Hadwiger's conjecture：09–13续读收据

固定提交 `4eaee87548250b2682affda7d05d947c05554bae`；原记录 `batch_22/records.json`、`batch_22/NOTES.md`。本收据仅新增五个原先未读的body文件，不重复领取其他文件的阅读。每一行均已按下表实际分块读完，所有源正文输出未截断。

**新增实际阅读：141663 bytes / 2939 lines。** 五个缓存文件逐一重新计算Git blob SHA与总行数，全部匹配原记录；没有修改缓存。原收据的5126行与本收据并集为20个text blobs、370905 bytes、8065 lines。仅在内存中运行 `reading.py.merge_continuation` 与 `validate_record`，得到 `full_body_read`、8065/8065行、无unread files；没有运行build、改台账或发布。

完整源码阅读不等于独立证明正确。主定理及后续推导在这里均记作稿件主张，未独立验证；PDF binary没有读取、编译或逐页视觉核验。

## 文件与实际阅读块

| 文件 | 实际块（inclusive） | Git blob SHA |
|---|---|---|
| build/sections/09-histograms.tex | 1–190；191–400；401–620；621–800；801–928 | `f38c13fc50a33a751722218c6271db9bd2710d90` |
| build/sections/10-product-tests.tex | 1–200；201–382 | `47a433dc0913018c06d08e09249144777152035e` |
| build/sections/11-common-limit.tex | 1–200；201–400；401–582 | `ebe1c6b2dc5eea6639502f3920d5304bd3604522` |
| build/sections/12-status-span.tex | 1–200；201–400；401–567 | `d5478451087ee1caac50b3ce3f39f08df15260cb` |
| build/sections/13-rare-overlap.tex | 1–180；181–360；361–480 | `26c5da67411fcf66eb37f93635ad1023cea22ce7` |

## 正文结构与完成路径

09-histograms先明确三种对象：common key tests、common input-and-key tests、orientation/input flags。所有filter形成subprobability，不除以filter概率。原始leaf仅需exponential joint cap；many-query Fourier按跨draw coefficient矩形rank计数，以method-of-types压过leaf cap，导出bounded coarse entropy；Pinsker refinement给all-filter L1 nets，small-set control给uniform integrability。nominal independence的失败在逐orientation下先删除，不能作为raw additive error再乘exponential cap。

同节的mixed-law结论另用subexponential marginal caps。一endpoint variance先用于任意固定external datum；paired typicality再以linear-in-n oversampling把第二endpoint异常压成2^(-Ω(n²))，因此joint cap不能抵消它。binary mixing对允许的unary weights统一成立。base-status positivity只对common global key cells成立；通过GL_n变换、quadratic polar rank及martingale difference平均，保持J和beta独立于cell mass，较小cell只改变onset。不是所有tests共享一个exceptional set。

10-product-tests的weak product theorem同时提供reference weak testing与actual-query comparison。finite twisted regularity用partition energy increment；terminal remainder经两次Cauchy–Schwarz变成four-copy common test，消掉任意orientation-dependent separate weights。Gram characters分为group内部、同侧copy交互、四quadrants：同侧非零项用Walsh rank，余项用twisted cut norm。递归最多r−1层，tolerance按已知系数与项数逐层分配，不反向改construction parameters。误差是weak testing而非L1；仅允许unary input weights、指定Gram/L,R constraints和orientation-only filters。

11-common-limit的option collection是固定有限数。成功query格式是orientation-only flag × unary statuses × specified binary/Gram indicators。每leaf自己的all-filter nets先固定；covariance operator trace bounds再提供统一deterministic projection，使pair-dependent选family member仍可控。tuple signature lists闭包加入projection cells与weak-product unary factors；显式保存conditional laws，避免把disintegration与weak convergence擅自互换。unflagged arrays在对leaf-pair record平均后仍iid；conditional on该record的两个orientation draws独立。

该节先用truncation与uniform integrability把vanishing overlap传到full limiting densities，再用actual和reference两部分weak-product comparison识别density：G_Z(u)=κ∫∏h^θ_t,status(ξ_l(u)) dα_Z(θ)。reference下unary projections独立，所以若所有finite options都vanish，almost every passed usable orientation pair必有某matched position的unary status overlap为零。

12-status-span把scalar recipe写成四edge spans的coupled system。每base vector都有fresh generic single-generator lift，projection surjective且kernel dimension≤2。dual obstruction必须使用非零coupled weights。no-prediction branch若obstructable pairs质量≥.1，conditioned-on-signature independent agreeing bits必须deterministic；generic positivity与共同surviving label使opposite affine offset唯一，只剩≤32对offsets；再用cylinder approximations及strict-error open events返回finite actual global-key predictors，违反no-prediction。不能逐pair重选互不相同的predictor后叫共同预测。

predicted branch将exact prediction代入dual relation，由base positivity迫使系数，并以prepared phase compatibility令所有target evaluations为零。no-prediction下passed-table质量≥.99与nonobstructable质量>.9相交> .89，不假设这两事件独立。最后删除≤100 labels能缩span就删除；四个kernel每个维数≤2，至多8次降维。稳定后≤3-generator fresh kernel words张成kernel，故一base lift加两个kernel corrections，总长≤7；每endpoint两incident lists总共≤14 labels，保留每个generator的positive-overlap certificate。由有限options取得subsequential positive constant overlap。

13-rare-overlap独立处理compatible leaf-pair质量Ω(N^-4)的exact zero-class，不能用positive-limit路径代替。exact prediction使每cross edge只需一个generic atom且无额外p-status filter。对任意预先指定的大reference set，假设polynomial bad mass并restrict-normalize后mixed caps仍成立；common-signature density identification配合h≥beta迫使fixed positive hit mass，矛盾，得到每固定p的o(n^-p) exceptional estimate。该restriction不保证每个normalized old leaf仍有原cap，证明也不用这种断言。

随后在原prepared B-leaves用all-filter L1 nets，先取uniform integrability threshold T与t,δ0，再small-set tolerance eta、net accuracy epsilon，最后才是bounded family size L。每B-family在independent A-draw之前已固定；具体table/leaf-pair filter只能从其中选择member。union/Markov删除o(N^-4) leaf-pair mass，再以A small-set control去掉W中B-density低于t/4的一小部分。剩余每compatible pair overlap≥t A0/8，总期望≥c0 t A0 N^-4/16。

最终两支带全部retained-law成本进入先前collision proposition，给four-hole probability ≥2^(-(k_max+.03)N-o(N))。k_max=56(g−1)，与100g的margin为44g+55.97>0，稿件据此完成uniform raw supersaturation contradiction，再引用先前raw-to-graph得到arbitrarily large m=2^(C0 g N)、alpha≤2、cm<m/100图。这里只报告实际读到的完成路径，并未独立证明这个counterexample。

## 有条件构造接口

### 1. Separate leaf compactness from mixed-law typicality, then build common projections before adaptive witness selection.

**精确前提：** Normalized original leaves have the exponential joint cap; the mixed law additionally retains subexponential marginal caps. Each own-leaf density family has bounded L1 nets, filters are unnormalized, and the two leaves are independent before pair-dependent usability and admissibility choices. Truncation height is fixed before projection accuracy.

**机制：** Many-query Fourier bounds give bounded coarse entropy and all-filter nets on each leaf. A covariance operator built from each own-leaf net has bounded trace. Quantizing its large eigenfunctions produces a deterministic common partition. Independence controls the sum of squared projection errors for every pair of net representatives, so either selected member may depend on both leaves without rebuilding its net using the opposite leaf. Uniform integrability then passes truncated zero products to full densities.

**可能桥接：** Conditional interface for Q(Z): if candidate coefficient or minor data admit bounded-complexity families and an independently sampled pair representation, construct one common finite observable before selecting compatible candidates. This could retain pair-dependent constraints while making limiting overlap claims reviewable.

**缺失步骤：** No arithmetic candidate family, covariance trace bound, independent representation, or common-basis preservation is established here. Quadratic occupancy is a coupled cost and cannot be replaced by sums of marginal costs. The construction does not provide the same integer determinant with both nonzero and analytic smallness.

**源码锚点：** `build/sections/09-histograms.tex:215` (lem:histogram-many-queries)；`build/sections/09-histograms.tex:344` (proof of lem:histogram-control)；`build/sections/11-common-limit.tex:132` (lem:limit-global-projection)；`build/sections/11-common-limit.tex:375` (lem:limit-zero-product)。

### 2. Use a four-copy test to transfer weak product approximation uniformly over allowed adaptive unary weights.

**精确前提：** Bounded common global tuple-key test, distinct labels at each endpoint, finite available Gram twists, raw frame comparisons, mixed joint and marginal caps, and successful queries containing only unary status weights, specified Gram/L,R bits, and orientation-only filters. The reference and actual-query comparisons are both needed.

**机制：** Finite twisted regularity gives a bounded structured sum plus small twisted-cut remainders. Two Cauchy-Schwarz steps remove arbitrary orientation-dependent unary weights and produce a common four-corner test. Typicality compares this test with its raw mean; a Gram-character classification bounds either nonzero same-side-copy interactions by Walsh rank or remaining quadrant interactions by the twisted cut norm. Recursion ends at unary factors; on the Gram-zero reference all twists equal one.

**可能桥接：** Conditional interface: for an explicitly defined arithmetic witness functional, seek a small common test algebra whose error is controlled after adaptive but separately acting row/column choices. A multi-copy moment can expose dependence lost by inspecting unary marginals.

**缺失步骤：** This is weak testing error, not L1 approximation, and it allows no arbitrary new joint input filter. Determinant integrality, height growth, full row rank of Q(Z), common bases, analytic smallness, and nonlinear occupancy need their own simultaneous estimates. F2 Walsh rank is not an estimate for the real absolute value of an integer determinant.

**源码锚点：** `build/sections/10-product-tests.tex:19` (thm:weak-product)；`build/sections/10-product-tests.tex:140` (lem:product-terminal-remainder)；`build/sections/10-product-tests.tex:170` (eq:product-four-corner)；`build/sections/10-product-tests.tex:366` (rem:weak-product-filter-class)；`build/sections/11-common-limit.tex:452` (lem:limit-density)。

### 3. Convert frequently recurring dual obstructions into one shared predictor, rather than permitting a new predictor for each local pair.

**精确前提：** Binary edge generator spans surject onto the full base status space; generic role/basis assignments have positive density on every common-key signature; bounded exclusion sets leave a common generic label; only finitely many coefficient and status formats are possible. Prediction must be a global function of actual keys common across retained units.

**机制：** A scalar obstruction has nonzero coupled weights. Conditional independence of the two status draws at the same signature forces agreeing bits to be deterministic. Generic positivity makes a deterministic opposite-side affine expression unique even as exclusions vary. At most 32 pairs of offsets remain after fixing one opposite array. Positive pair mass therefore yields a fixed shared offset on positive unit mass. Cylinder approximations and strict-error events transfer the limiting function back to genuine finite-key predictors.

**可能桥接：** Conditional rank interface: classify recurring annihilator obstructions to Q(Z) by a bounded family of common functionals before concluding that rank deficiency is structural. This may distinguish an actual shared obstruction from unrelated pair-specific failures.

**缺失步骤：** The source works over F2 with fixed finite metadata. An arithmetic analogue needs uniform coefficient-height and denominator bounds, a valid finite classification of annihilators, and proof that the same candidate common basis survives. Frequency alone gives neither an integer nonzero minor nor its analytic smallness, and no claim about π BA/nonBA follows.

**源码锚点：** `build/sections/12-status-span.tex:100` (lem:status-duality)；`build/sections/12-status-span.tex:145` (lem:status-global-prediction)；`build/sections/12-status-span.tex:193` (bounded and unique opposite offsets)；`build/sections/12-status-span.tex:230` (finite-key transfer of common prediction)。

### 4. Stabilize spans under bounded support deletion, then lift a target with one generator and at most two short kernel corrections.

**精确前提：** For each of four binary edge spans V, the base projection is onto and its kernel has dimension at most two. Every base vector admits a single-generator generic lift with fresh labels. Mandatory exclusions and prediction exceptions fit a fixed label budget, and scalar solvability survives the stated further exclusions.

**机制：** Whenever deleting at most 100 labels shrinks an edge span, make the deletion. Since only the kernel dimensions can fall, at most eight dimension losses occur over four edges. In the stabilized spans, label-disjoint words of length at most three span the kernel: compare two lifts through a fresh third lift, and enforce additivity with lifts of b,c,b+c. One base lift plus at most two three-generator kernel words realizes any target in a nonempty list of at most seven, retaining positive-overlap certificates and distinct labels.

**可能桥接：** Conditional construction interface for a common-basis search: after proving a target projection is onto, isolate a small defect kernel and realize corrections through fresh support rather than allowing support conflicts to accumulate. A deletion-stable certificate set can make a bounded support budget explicit.

**缺失步骤：** No surjective integer-lattice projection or bounded defect kernel for Q(Z) is supplied. The F2 identity b+c+(b+c)=0 and bounded word lengths do not bound integer coefficients, denominators, lattice indices, column heights, or nonadditive quadratic occupancy. Any analogue must retain the actual determinant witness and common basis.

**源码锚点：** `build/sections/12-status-span.tex:51` (eq:status-full-projection)；`build/sections/12-status-span.tex:415` (stable spans and fresh short lists)；`build/sections/12-status-span.tex:425` (lem:fresh-lists)；`build/sections/12-status-span.tex:460` (eq:status-short-kernel)。

### 5. Protect an inverse-polynomial compatibility event by proving smaller exceptional errors before restricting to it.

**精确前提：** The rare branch has exact prediction p_i=T(lambda_i,.)+a, pinned compatibility of original leaf pairs of mass at least c0 N^-4, separate admissibility filters of mass at least a0, and generic base queries with no additional p-status filter. Mixed caps survive polynomial-mass law restriction; original prepared leaves separately retain their histogram caps.

**机制：** If a deterministic reference set of fixed positive measure had polynomially many bad orientations, restrict to them; the mixed caps survive and limit density identification with unary densities at least beta contradicts the proposed small hit mass. This gives a superpolynomial exceptional probability uniformly per prescribed set. On each B-leaf, all-filter nets yield finitely many large superlevel sets before the independent A-leaf is drawn. Apply the hit bound to that predetermined family, discard o(N^-4) pair mass, then use A small-set control to remove the part where the selected B-density is low. The remaining expected overlap is at least c0*t*A0*N^-4/16.

**可能桥接：** Conditional arithmetic interface: if full-rank/common-basis candidates occupy only an inverse-polynomial fraction, first obtain an exceptional-set estimate below that fraction and a bounded opposite-side candidate family. Then a coupled height or analytic estimate may be intersected without incorrectly treating conditioning as independence.

**缺失步骤：** No such arithmetic large-set hit estimate, original-family compactness, or exact prediction identity has been proved. A rare feasible determinant must also satisfy the analytic estimate, column-height bound and true quadratic occupancy bound simultaneously. The source cannot supply these by analogy; π BA/nonBA and q²|π-p/q| remain unresolved.

**源码锚点：** `build/sections/13-rare-overlap.tex:73` (lem:rare-large-set-hit)；`build/sections/13-rare-overlap.tex:216` (rem:rare-restriction-caps)；`build/sections/13-rare-overlap.tex:235` (lem:rare-superlevel-nets)；`build/sections/13-rare-overlap.tex:385` (eq:rare-independent-net-error)；`build/sections/13-rare-overlap.tex:427` (eq:rare-overlap-bound)。

## 限制

- Full source coverage is a report of reading, not independent verification of the Hadwiger counterexample claim, raw supersaturation, imported lemmas, or cited dependencies. No PDF compilation, visual page check, formal proof verification, or external literature validation was performed.
- Leaf joint caps do not imply mixed marginal caps or the converse. Polynomial restriction preserves mixed caps but need not preserve caps for each normalized old leaf. Rare-overlap nets are used on original prepared leaves only.
- Common-test typicality has a test-dependent exceptional set, not simultaneous control over every test. Generic base-status positivity applies to common global key cells, not arbitrary input-and-key cells. Countable limiting cylinder inequalities do not assert finite-n simultaneous positivity for all measurable cells.
- Weak product comparison is weak testing error rather than L1 error. Its filter class is unary orientation-dependent input weights, specified binary Gram/L,R conditions, and orientation-only filters. An unprescribed joint input predicate is not permitted.
- Conditional orientation laws are stored through the limit; unflagged base arrays are iid only after averaging the leaf-pair record. Leaves are independent before pair-dependent restrictions. No conditional independence after a joint compatibility event is assumed.
- All algebra, status spans, and rank arguments here are over F2. Positivity, probability, entropy, and Hilbert-space density controls are real. These distinctions cannot be erased in an integer determinant application.
- The rare branch needs its own superpolynomial exceptional control because Ω(N^-4) compatibility can disappear in a positive-mass common limit. The positive-overlap branch alone cannot cover it.
- No conclusion on π BA/nonBA or q²|π-p/q| is obtained. An arithmetic transfer must prove nonzero and analytic smallness for the same integer determinant, together with Q(Z) row rank, a candidate common basis, column-height control and nonadditive quadratic occupancy. Every construction idea states unproved bridge requirements.

π的BA/nonBA及q²|π−p/q|仍未解决。这里没有从F2 rank、probabilistic overlap或finite compactness自动转移到算术结论。必须对**同一个integer determinant**同时给出nonzero与real analytic smallness，还要在同一构造里满足Q(Z) row rank、candidate common bases、columnheight及nonadditive quadratic occupancy。
