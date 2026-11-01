# 第 25 组：队列 47、49 正文首读

固定上游：`openai/math@4eaee87548250b2682affda7d05d947c05554bae`。按全库队列阅读，无相关性筛选。本组只提交 47、49；48 的全部本地源码已缓存并逐一校验 SHA，但本阅读分支没有读其正文，也没有为其创建完成或部分阅读记录。48 交由新的阅读分支处理。

本组完整阅读 21 个本地文本文件，共 4,002 行、161,682 字节。正文首读完成只说明相应主 TeX、递归数学输入、图文字和参考文献已读；不等于独立证明验证，未逐项核验引用文献，也未渲染 PDF。

| 队列 | 稿件 | 本地文本 | 已读行数 | 字节 | 根文件 |
|---|---|---:|---:|---:|---|
| 47 | An irrational-trace counterexample to reduced Baum–Connes | 11 | 1996 | 84649 | `paper.tex` |
| 49 | Lipschitz Equivalent Separable Banach Spaces Need Not Be Linearly Isomorphic | 10 | 2006 | 77033 | `main.tex` |

## 47. An irrational-trace counterexample to reduced Baum–Connes

The paper claims a finitely generated discrete group G_sur=V semidirect Gamma whose ordinary reduced group C*-algebra contains a projection b_0 with irrational canonical trace sum_i p^(-k_i), k_i=floor(n_i/2) increasing with gaps tending to infinity. Its integral K_0 class, and its rationalization, lie outside the coefficient-free proper-action reduced Baum–Connes assembly image, so the paper claims a counterexample to that conjecture for countable discrete groups.

### 实际阅读与构造架构

- README 1–18, paper.tex 1–29 and preamble 1–43: PDF identity, abstract, actual main input graph, macros, theorem environments and references.
- Introduction 1–195: coefficient-free reduced assembly, full proper classifying space, Lueck rational trace image deduction, prior irrational von Neumann kernels versus reduced C*-projections, and complete construction overview.
- Geometry 1–374: deterministic labelled Cayley graphs; distinct-copy argument; singular plane diagrams with adjustable sides and contour incidences; all face/bridge/leaf surgeries; arc count; avoidance and simplicity; universal coarse filling/thinness and root-based simultaneous overlap ball.
- Data 1–166: choices in dependency order, expansion from Laplacian gap, girth/tree ball estimates, Osajda labelling hypotheses and embedded relators, parameter table, and spanning/independence of the finite-support copy cycle spaces.
- Module 1–313: port decoder and prototype; prescribed-word codes with simultaneous primal and dual distance by finite-field counting; quotient relation module; exact intersection using farthest-copy support; finite generation; compact dual; exact uniform single-copy Haar marginals without independence.
- Spectrum 1–499: faithful evaluation regular representations and trace integral; enabled graph operator; valid components; invalid-copy deficit; algebraic cycle-coordinate elimination; localization of projected cycles to a tree union; exterior cycle rank, separated attachments, disjoint charge sets; coarea/Cauchy–Schwarz global gap; continuous functional-calculus projection.
- Figure 1–11 and its full surrounding caption: overlap tree, exterior components and attachment charging. Trace 1–67: exact number of copies through identity, mutually exclusive validity events, trace series and base-p nonperiodicity, assembly obstruction and rationalization. Bibliography 1–281 read in full.

按原稿的证明链：

- Use Lueck’s cited trace theorem after scalar extension to the finite-subgroup-denominator ring Lambda_G to force every integral assembly trace into Q, with the full proper-action direct-limit domain explicitly compared.
- Choose universal coned-geometry constants before lambda; select fixed-degree high-girth expanders, then one finite reduced label alphabet with globally unique long paths and isometrically embedded graphical relators. Minimal adjustable-side plane diagrams give short arcs, universal linear fillings, avoidance, and thin triangles.
- Relative to a root and fixed root geodesic, every overlap of a copy with all copies no farther from the root lies in one intrinsic induced tree ball of radius <=girth/16 and <eta times its vertex count. The same ordering proves an internal algebraic direct sum of finite copy cycle spaces.
- Over one sufficiently large prime field choose dimension floor(n_i/2) codes containing the exact port pattern with both primal and dual nonzero supports >eta n_i. A farthest-copy uncancelled coordinate proves W_0 intersect F_p^P=C_P-perp for the translation-invariant sum of all dual-code relations.
- Duality gives exactly uniform C_P marginals on each copy; no independence is asserted. The semidirect product is finitely generated. Clopen port projections define d minus the enabled adjacency operator in the ordinary reduced algebra via Fourier transformation.
- Valid patterns are isolated disjoint expander components. Invalid copies have a linear degree deficit. To globalize it, eliminate algebraic copy cycle coordinates in farthest-first order, retain newly introduced edges only inside the common overlap tree, bound attachment multiplicity by girth, and charge each rank to disjoint exterior vertex sets.
- The resulting finite-set deficit plus coarea and Cauchy–Schwarz bounds the entire remaining graph below by c²/(8d). Together with the expander gap rho, all evaluation representations have a common gap epsilon=min(rho,c²/(8d)); faithfulness and continuous functional calculus produce b_0 inside C_r*(G_sur).
- Exactly n_i indexed copies of type i pass through identity; mutually exclusive validity events each have measure p^(-k_i), and projection diagonal 1/n_i cancels n_i. Sparse base-p digits yield irrationality; the rational assembly trace restriction excludes the class, also after rationalization.

### 保留的条件构造思路

1. **Farthest-first elimination with one simultaneous overlap container can certify exact local relations after imposing infinitely many global relations.**

   源锚点：build/geometry.tex, Lemma lem:overlap and equation eq:overlap-radius；build/data.tex, Lemma lem:cycle-sum；build/module.tex, Lemma lem:module-intersection and equation eq:module-intersection

   前提：Every finite participating family admits an order in which all intersections of the current block with remaining blocks lie in one designated small intrinsic ball. For cycles that ball is an induced tree; for module relations every nonzero dual codeword has support larger than the ball. All sums are finite-support algebraic sums, and the chosen root/ball works simultaneously for every earlier-or-equal-distance block.

   候选接口：For the current Q(Z) row-rank/common-basis problem, seek structured row blocks and a finite-support elimination order with an exact shared-coordinate container that cannot support a nonzero relation. This can turn local independence into an exact global rank certificate under explicit premises.

   缺失桥梁：Pairwise small overlaps do not establish the simultaneous container. The source works over F_2 and F_p, whereas the target needs characteristic-zero rank and one integer determinant. A lift preserving exact lattice relations and height, plus a quantitative bound on that same nonzero determinant, is absent.

2. **Convert overlapping local deficits to a global estimate by charging rank only to disjoint exterior resources.**

   源锚点：build/spectrum.tex, Lemma lem:global-deficit, eq:cycle-elimination, eq:projected-cycle-support and eq:exterior-rank-bound；build/figures/cycle-charging.tex, full figure and caption

   前提：The finite cycle space decomposes as an internal direct sum of copy cycle spaces. After farthest-first coordinate elimination, every introduced edge outside the enabled graph lies in the induced overlap tree. An attachment pair creates a long simple cycle; girth packs its endpoints. Local invalid-copy degree deficits hold uniformly, and charged exterior vertex sets are pairwise disjoint.

   候选接口：The current column-height plus nonadditive quadratic-occupancy task could search for a rank/dimension surrogate whose successive contributions can be assigned to disjoint exterior coordinate sets, instead of adding overlapping local bounds. The exact tree-attachment identity b_1(D)=sum_C(b_1(C)+max(0,l_C-1)) illustrates how to retain an interaction correction rather than discard it.

   缺失桥梁：No target occupancy functional, elimination order or exact interaction formula has been constructed. Cycle rank over F_2 is not a quadratic arithmetic occupancy estimate; one must prove a target-specific packing mechanism and translate its constants into the height of the same integer determinant.

3. **Enforce one prescribed finite-pattern law with both primal and dual code distances, then separately globalize spectral isolation.**

   源锚点：build/module.tex, Lemmas lem:codes and lem:marginals；build/spectrum.tex, Lemmas lem:local-deficit and lem:spectral-gap; Proposition prop:projection；build/trace.tex, Proposition prop:irrational-trace

   前提：A fixed full-support word belongs to a half-dimensional code; primal distance creates discrepancies and dual distance prevents additional local relations. Exact quotient intersection gives Haar marginals p^(-k_i). A uniform gap holds in every faithful evaluation representation on both valid components and their entire complement.

   候选接口：Distant construction template: combine exact local statistical constraints with a separately proved global analytic separation, keeping local law and global gap as distinct obligations. Sparse prescribed probabilities then encode a nonperiodic scalar trace.

   缺失桥梁：The resulting scalar is an engineered sparse base-p series, not pi. It provides neither BA nor non-BA behavior of pi and no estimate on q²|pi-p/q|. A uniform spectral gap is not an integer lower bound; no map from the trace construction to the target approximation determinants is supplied.

### 阅读边界

- Full source reading is not independent theorem validation. The asserted coefficient-free Baum–Connes counterexample and every novel diagram, simultaneous-overlap, module and global-rank argument remain claims of the manuscript.
- External inputs Lueck trace theorem and formulation comparison, LPS expanders, Osajda labelling/embedding, and coarse-filling hyperbolicity theorem were identified but not independently checked against their cited sources. The local bibliography was read, not the cited works.
- The simultaneous overlap lemma, exact module intersection and disjoint exterior cycle-rank charging are essential; pairwise small cancellation or separate expander gaps alone do not establish the conclusion.
- All copy marginals are single-copy statements; independence across copies is neither proved nor needed. Algebraic finite-support cycle projections are not edge restrictions and can introduce edges, which the proof explicitly localizes to the overlap tree.
- The uniform bound is on all evaluation representations and the entire remaining enabled graph, not only the selected expanders. The paper uses continuous C*-functional calculus only after that step.
- No consequence for BA/non-BA of pi or q²|pi-p/q| has been established. No characteristic-zero row-rank theorem, common basis, target quadratic occupancy bound or integer determinant simultaneously nonzero and analytically below one is supplied.
- All local mathematical input, figure source/caption and bibliography are covered. PDF rendering/byte comparison and recursively reading external references are outside this source-reading receipt.

README and all ten local build text blobs (main, preamble, six body sections, TikZ figure and bibliography) were actually read through their full line extents. A truncated initial concatenation lost the tail of data.tex; it was explicitly recovered by an untruncated read of data.tex 115–166, overlapping the observed table and covering the entire direct-sum proof. All other coverage used complete visible chunks; spectrum was read 1–260 and 261–499. Exact cached bytes match the pinned recursive-tree blob SHAs. Source mapping and local closure are complete; reading coverage does not certify theorem correctness.

## 49. Lipschitz Equivalent Separable Banach Spaces Need Not Be Linearly Isomorphic

The paper claims separable real Banach spaces X,Y and a global bijection Psi with (4/21)||s-t|| <= ||Psi(s)-Psi(t)|| <= (76/25)||s-t||, while X contains a linearly isometric c_0(ell_2) and Y has no closed linearly isomorphic copy of c_0(ell_2). Thus they are Lipschitz equivalent but not linearly isomorphic. It asserts neither reflexivity nor a complex-scalar analogue; exclusion of ordinary scalar c_0 from Y is not asserted.

### 实际阅读与构造架构

- README 1–18 and main.tex 1–60: title/date, abstract, complete package/macro declarations, actual six-section input graph and bibliography.
- Introduction 1–183: exact constants and c_0(ell_2) obstruction, historical scope, full construction overview and distinction between metric equivalence and linear structure.
- Free graph 1–357: evaluation isometry, separability, same-norm McShane extension and free-space linearization; precise external compact-reduction and weak sequential-completeness inputs; closed weighted graph space; no-c_0 in ell_1(E); Hilbert-block gliding-hump obstruction using complete continuity.
- Triangular 1–178: weighted graph criterion; B(e,v)=e+delta(qe,v); both explicit inverse identities, sum-norm constants, coordinatewise c_0 membership and global two-point estimates; isometric c_0(V) copy and no-section warning.
- Localization 1–395: full shrinking-ball cutoff/annulus/Schur lemma; diagonal subsequence of local Lipschitz constants; finite exceptional-point cardinal bound L²/eta²; disjoint cutoff control; compact reduction and finite-cover remainder bound; exact ordered limits without uniform bounded-set tail assumptions.
- Hilbert 1–578: frozen even/odd slot isometry and inverse; weights; stage dimensions and separating bump coordinates; shared radial schedule with infinite first integral and finite weighted square integral; telescoping derivative budget; local slot support separation; finite-dimensional fiber proper covering argument; two-point inverse bounds; exact inverse-tail cancellation, compact-head diagonal limit and surjectivity.
- Assembly 1–52: applies all hypotheses, q norm <=51/25, Psi constants 76/25 and inverse 21/4, and linear c_0(ell_2) obstruction. Figure 1–17 and caption: perpendicular/even/odd slot routing. Bibliography 1–168 read in full.

按原稿的证明链：

- For a complete superreflexive Hilbert metric domain M, take the Lipschitz-free space E=F(M). Use the stated compact-reduction and weak sequential-completeness theorems with their actual hypotheses; establish evaluation/linearization and weighted graph-space facts directly.
- Given E weakly sequentially complete and q completely continuous, show Z_q has no c_0(ell_2): each Hilbert block sends an orthonormal sequence to coordinatewise norm-null c_0(U) outputs; a disjoint-support subsequence has bounded sums against sqrt(N) input norm. Pick one vector per block with tiny output, forcing a scalar c_0 embedding in ell_1(E), which a second gliding-hump argument excludes.
- Assuming onto bi-Lipschitz h:U direct-sum ell_2→U, linearize D=h-pi_U as q. The triangular B(e,v)=e+delta(qe,v) satisfies qB=hQ_1; recover x=h^{-1}(qz) and B^{-1}(z)=(z-delta(x),pi_Vx). Apply coordinatewise to graph spaces; both c_0 membership and weighted sum norms are controlled. This is not a Lipschitz section of q.
- Prove local finite-family orthogonality of high output coordinate subspaces makes q completely continuous. Scalar tests from a hypothetical nonzero weakly null image have only finitely many persistent high local-Lipschitz points by an orthogonal-square budget. Shrinking-ball tests pair negligibly via disjoint annuli and Schur; compact reduction controls the remainder with the specified order of sequence/radius/error limits.
- Construct h as an input-dependent frozen onto isometric slot rearrangement: perpendicular components stay, old slots route to even blocks, additional ell_2 slots route to odd blocks. Smooth stage-k bump coordinates in tanh coordinates eventually separate neighborhoods of every finite family.
- Set a_n²=r_n²-r_(n+1)² and rotate stages according to gamma(r)=c/(r log(1/r)), c=1/100: infinite integral near zero permits arbitrarily many stages while integral 2r gamma²=c² bounds summed weighted derivatives. This gives finite approximant derivatives within 4c of frozen isometries and hence local inverses.
- For finite approximants, the deviation from a fixed isometry lives in a finite output head. Norm preservation makes each finite-dimensional affine fiber proper; covering of Euclidean space proves global bijectivity before integrating inverse derivative bounds. Passing to pointwise limits preserves two-point bounds.
- Surjectivity of the limiting h requires a separate exact inverse-tail identity: preimage tail squared equals the output tail squared minus a finite nonnegative slot-coefficient sum. Uniform tails and bounded finite heads give a norm-convergent preimage subsequence. This yields onto h with constants 24/25 and 26/25 and local orthogonality.
- Assemble E,q,Q_1,X,Y; bound ||q||<=51/25. The weighted criterion gives Psi Lipschitz <=76/25 and inverse <=21/4. X visibly contains c_0(ell_2); complete continuity plus E weak sequential completeness excludes it from Y, giving the claimed linear nonisomorphism.

### 保留的条件构造思路

1. **Allocate infinitely many separation stages using a divergent first-integral but finite quadratic weighted derivative budget.**

   源锚点：build/sections/05-hilbert.tex, eq:gamma, eq:gamma-integrals, eq:angle-bounds and eq:slot-budget；build/sections/05-hilbert.tex, eq:radii-identities and Lemma lem:slot-locality

   前提：Nonnegative tail radii r_n decrease to zero, and a_n²=r_n²-r_(n+1)². Slot rotations have derivative <=gamma(r_n), where gamma decreases, integral gamma diverges near zero and integral 2r gamma² is finite. Orthogonal finite stages have finite derivative constants and may start at arbitrarily small radii; separating coordinates distinguish every finite input family.

   候选接口：Distant multiscale construction interface: a finite square budget can pay for infinitely many successively finer separation requirements. If a target occupancy/height construction has a telescoping resource and quadratic motion cost, search for a scheduling function separating total progress from total squared cost.

   缺失桥梁：The target determinant problem has no established tail-radius identity, smooth rotation freedom or orthogonal stage encoding. These operations need not preserve an integer lattice, rational dependence or height. A target-specific discretization and bound for the exact same nonzero determinant are missing; no conclusion on q²|pi-p/q| follows.

2. **Use orthogonal-square occupancy to isolate finitely many bad locations before compact reduction.**

   源锚点：build/sections/04-localization.tex, Theorem thm:localized-linearization and eq:exceptional-cardinality；build/sections/04-localization.tex, Lemma lem:shrinking-tests; eq:compact-remainder and eq:final-local-bound

   前提：A complete separable metric domain, finite-dimensional Hilbert output heads, and for every finite distinct family a common tail cutoff with mutually orthogonal coordinate subspaces on neighborhoods. Weakly null unit tests then spend at least (eta/L)² mass per persistent bad point, so there are at most L²/eta². Uniformly Lipschitz shrinking-ball tests and compact reduction control the exceptional and remaining pieces in a fixed limit order.

   候选接口：For nonadditive quadratic occupancy, this suggests defining a normalized square-mass resource and showing each distinct persistent obstruction consumes a definite amount. It can turn infinitely many local tests into finitely many exceptional pieces plus a compact remainder, conditional on actual orthogonality and tail-head separation.

   缺失桥梁：No target analogue of local orthogonality, weakly null tests or compact reduction has been proved. This analytic finiteness bound does not imply a common rational row basis or integer determinant nonvanishing. Quantitative arithmetic height and analytic smallness of the same determinant remain separate obligations.

3. **Prove surjectivity of a limit by an exact inverse-tail cancellation rather than pointwise convergence of bijections alone.**

   源锚点：build/sections/05-hilbert.tex, Lemma lem:frozen-inverse and eq:frozen-inverse；build/sections/05-hilbert.tex, eq:inverse-tail and completion of Theorem thm:hilbert-map

   前提：Each finite approximant is globally invertible with a uniform upper Lipschitz bound. Frozen inverse coordinates satisfy an exact square-norm formula, and even/odd routed slots partition the output tail. The remaining term is a finite nonnegative subtraction. Finite output heads are finite dimensional; fixed-output tails vanish uniformly over approximant inverse images.

   候选接口：For an infinite family of candidate witnesses or bases, seek an exact conservation identity that bounds every approximant witness tail by one fixed target tail. Together with finite-head compactness, this can produce a convergent actual witness and prevent mass from escaping.

   缺失桥梁：Pointwise full rank or the existence of separate bases does not itself give such a conservation identity or finite-head compactness. Rational/integer minors are discontinuous at rank loss and their denominators may grow; a target-specific lattice/height bound and stability of the chosen common minor are needed.

4. **Keep the exact intertwining identity and the obstruction property as distinct construction interfaces.**

   源锚点：build/sections/03-triangular.tex, Proposition prop:criterion, eq:diagram and eq:inverse-diagram；build/sections/03-triangular.tex, Remark rem:no-section；build/sections/02-free-graph.tex, Lemma lem:graph-obstruction

   前提：A genuinely onto bi-Lipschitz h with h(0)=0, its exact free-space linearization q of h-pi_U, and complete continuity of q on a weakly sequentially complete E. The triangular B has an explicit inverse derived from h^{-1}; graph-space membership includes norm-null Q images, not merely bounded images.

   候选接口：This is a useful discipline for proposed rank/analytic constructions: certify the exact identity linking the chosen object to the target, then prove the required obstruction or compactness property separately. Metric equivalence can coexist with failure of linear isomorphism, so approximate metric control cannot replace algebraic rank evidence.

   缺失桥梁：No target h,q or graph-space model exists here. The source identities are nonlinear and do not give a section of q, a common row basis or determinant integrality; any such inference would require a new exact algebraic bridge.

### 阅读边界

- The claimed solution of the separable real Banach Lipschitz-isomorphism problem is reported as stated, not independently validated. Complete source reading does not certify the localized-linearization theorem, moving-slot construction or final counterexample.
- External compact-reduction and weak sequential-completeness theorems and cited version-specific locator statements were identified but not independently read/verified. The comparison with the corpus Tingley manuscript is contextual and not used to establish the construction here.
- The obstruction is c_0(ell_2), not scalar c_0. The graph criterion needs complete continuity, not compactness, and the proof does not claim all operators into c_0(U) are compact.
- Local orthogonality is for neighborhoods of arbitrary finite distinct families after a shared tail cutoff. Neither a global partition nor uniform tail decay on bounded sets is claimed. The compact-reduction argument preserves the stated order of limits.
- Pointwise norm preservation and derivative bounds alone do not prove global invertibility. The paper separately proves finite-dimensional fiber properness/covering, then an exact uniform inverse-tail identity for the limiting surjectivity.
- Nonlinear identities qB=hQ_1 and Q_1B^{-1}=h^{-1}q do not provide a Lipschitz right inverse of q. The cited separable quotient-splitting corollary is not automatically applicable.
- No direct BA/non-BA conclusion for pi, approximation estimate q²|pi-p/q|, Q(Z) rank theorem, candidate common basis or lattice determinant simultaneously nonzero and below one is derived. All proposed arithmetic interfaces remain conditional and missing bridges are explicit.
- All mathematical local sections, macros, complete figure source/caption and bibliography are covered. PDF rendering/byte comparison and recursively reading references are outside this receipt.

All ten mapped local text blobs were actually read through their full line extents in untruncated visible chunks, including all six section proofs, main preamble, figure source and full bibliography. Hilbert section was read 1–225 and 226–578; localization 1–395. Exact cached bytes match the pinned untruncated recursive-tree Git blob SHAs. Mapping and local closure are complete. This is source-body reading, not independent theorem validation.

## 与当前 π 路线的边界

两篇都保留跨领域的远端构造线索；没有把相似的“秩”“平方预算”“谱隙”词语当作现成迁移。当前问题仍须在同一个整数行列式上同时证明非零与解析小量；Q(Z) 行秩、候选共同基、列高及非加性二次占据各有独立义务。有限域循环秩或 Hilbert 正交平方预算不会自动补齐这些义务。π 的 BA / non-BA 和 q²|π−p/q| 在本组阅读后仍未解决。
