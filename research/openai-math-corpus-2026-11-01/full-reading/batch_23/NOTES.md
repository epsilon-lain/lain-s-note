# Round 03 · Batch 23 · 全文首读（队列43、45）

Pinned upstream commit: `4eaee87548250b2682affda7d05d947c05554bae`。无relevance筛选。全部文字源码已实际逐块阅读；download、SHA核验不代替阅读。完整首读不等于独立证明验证。

## 43. A Linear Clock for Random Walk on Tree-Weighted Planar Maps（family 211）

状态：`full_body_read`；11 files，4400/4400 lines，208959 exact Git-blob bytes。Mapping/closure complete；unread files为空。

对n-edge rooted planar maps按spanning-tree数量加权，stationary degree measure起始的rate-one incident-half-edge walk X^n，在空间deterministic a_n缩放和时间nt下，作为conditional path probability law与metric-measure space联合收敛到ordinary unit-area sqrt2-LQG sphere上的Liouville Brownian motion；continuum Dirichlet form=(1/2)∫|∇u|²，故exact clock multiplier c=1，through all n，任意fixed finite horizon。结论使用companion同ensemble的metric/contour/topological limits；另得finite-law small-mass centered Green bound，任意rho∈(0,1/2)，||G_nf||∞≤A_n mu_n(|f|)^rho、sup E A_n<∞，同时所有bounded centered f。

### 证明架构

- 在pinned companion contour reconstruction inputs之上，先把全部exact primal/dual incidences保留下来；finite instruction graphs对non-tree edges按matched occurrence而非端点去重。triple contacts以middle run实际discrete minimum精确cap。固定countable compactified energy/EL outputs在抽subsequence之前加入；input-only change-of-density与bridge law给conditional kernels，forest offset基变换unimodular，避免把旧kernel取forced-offset diagonal。
- complete buffered prescriptions和strict selector sandwich把抽象outputs等同actual local tests。field-local kernels在disjoint buffers conditional独立；GFF point germs因point zero capacity trivial，于Lebesgue-a.e.点任何array germ确定。quantum unit dyadic transforms和random cone zoom在joint extraction中保留，Gaussian reference比较只用fixed density truncation。
- 若annular EL quantile b_m→0，record lows保留所有backward offsets的barriers；all-point annuli构造charges macroscopic displacement且b_m||g||²≤area。small traffic flows由crossing-count紧性给rectifiable limit，average arclength≪area。tangent crossings成为local germ events，reflection/rotation给a.e. perpendicular primal-dual pair；每对paths必须遇paired edge，1≤〈j,j*〉≤||j||||j*||=O(b_m)矛盾。于是constant energy ring cutoffs。
- dual cutoff coarea给separating primal circuits与oscillation bound，最大值原理给harmonic compactness。uniform relaxed variational envelope有lower/recovery和bounded-energy gluing；quadraticity/measure property给A(z)，neighbor grid circuits必须相交从而Sobolev coercivity，不要求face regularity。point germs/coordinate symmetries得同deterministic scalar A=A*=lambda I。
- annular condenser capacities趋lambda C0，rotated discrete current的dual potential period恰为C_m，period energy lower bound给lambda≤1。另从canonical angular cochain a_m=(1/2π)∫darg构造alpha_m=a_m+dq_m，q_m是local recovery errors的partition gluing；period逐finite resolution精确保留为1，sharp weighted recovery energy≤lambda/C0。与同condenser current的signed pairing literal |Σdu_m alpha_m|=1，Cauchy给lambda≥1。故exact energy normalization。
- fixed-area sphere interior excursion density只依赖increment block且uniform L2；加入全部array outputs仍保留同tilt。interval interiors覆盖除tour endpoint，point zero-capacity removability与target constant-nearpoint approximation补端点，并与同contour metric/dense projection/quantum area coupling保持一致。
- finite graph centered inverse经rank-one det(Hhat+zqq^T)展开：coefficient source determinant为q(S)，Cauchy–Binet给two-forest sum，再按F=T\e的b(F) extensions重组。Mullin tree cuts小侧tour intervals mass=ell/N、b_e=1+dual distance；condition schedule后primal/dual Dyck independent。negative distance L2∼M^-1/2乘nested primal cut count L2∼sqrt M给one-scale bound，grid与source min(M/N,m_f) interpolation求dyadic sum，产生map-measurable A_n。
- Green bound使小quantum-mass source potential统一小，subtract后resolvent harmonic所以uniform compactness；variational minimization识别continuous resolvents，uniform approximate identity构造bounded-generator spatial separators，optional stopping给uniform exit与Aldous tightness。stationary nested resolvent multivariate Laplace transforms识别conditional path law，exact degree identity nL_n=−(1/2)mu_n^-1 H_n与lambda=1给time n。

### 有条件构造线索（含远域线索）

- **先构造具有精确整数period的cochain，再仅加single-valued gradient修正以降低能量；同一finite object保留literal非零pairing，之后才施解析Cauchy upper bounds。**
  - Source: build/sections/04-energy.tex:645-722 (en:angular-recovery; alpha_m=a_m+dq_m)；build/sections/04-energy.tex:724-755 (en:reciprocity; exact signed pairing=1)；build/sections/02-inputs.tex:150-216 (inp:topology; integral winding/degree sign)
  - 前提：canonical a_m来自projected integer winding，overlay mesh与cell closedness；q_m全局single-valued、local errors uniform→0，partition weights bounded energy且sum1；cochain closed区域必须覆盖同一condenser全部nonzero increments；primal/dual incidences和orientation exact。
  - 可能接口：π方向最有结构针对性的条件设计：找一个integer-valued topological/algebraic functional I(D)=1，允许的analytic basis corrections落在I的kernel，从而同一个candidate minor/linear form保持非零而优化height或smallness。把nonzero证据先固定，再估计同一对象的column norms，类似这里pairing=1与energy bound共同约束lambda。
  - 缺桥：正文的pairing是cohomological current/cochain pairing，非Q(Z) row determinant；没有π、integer coefficients、q-denominators、Q(Z) rank或quadratic occupancy界。需新证候选rows的canonical integral class、correction kernel、共同basis、height代价和π-analytic pairing完全一致；这一步未构造。
- **rank-one Laplacian determinant的同一系数，以Cauchy–Binet转成非负two-forest squares，再polarize并精确计cut multiplicities，得到inverse bilinear identity。**
  - Source: build/sections/06-green.tex:53-130 (gr:forest; gr:rank-one through gr:two-forest-bilinear)；build/sections/06-green.tex:147-193 (gr:tour-cut; b_e=1+D_e*)
  - 前提：finite connected multigraph、integer incidence matrix去一row、positive unit conductances；q/r总和零；同一个reduced H=BB^T和tau=det H；Cauchy–Binet source determinant经leaf elimination为±q(S)，每forest的tree extensions数b(F)必须保留。
  - 可能接口：检视candidate common bases时可尝试把Gram/minor determinant表示成同一integer incidence matrix的forest/minor sum，利用positive squares给rank与nonzero的独立certificate，再以相同expansion组织column height/occupancy估计。也可把normalizing denominator tau在algebraic与analytic两侧显式清除，避免换掉determinant。
  - 缺桥：π approximation rows目前没有graph incidence/total-unimodularity或forest parameterization，source vectors未必centered；positive Gram determinant只给Euclidean rank，不能自行升级为Q(Z)rank、uniform integer lower bound及π smallness。需要同一matrix的representation和定量规模控制，本文不提供。
- **连续目标的variational convergence不足以控制稀疏质量source；单独构造uniform small-mass inverse estimate，再把inhomogeneous对象拆为small potential+harmonic part。**
  - Source: build/sections/06-green.tex:324-483 (one-scale bound, source interpolation, rho<1/2 dyadic sum)；build/sections/07-dynamics.tex:134-213 (uniform resolvent convergence; localized q_n)
  - 前提：all bounded centered sources同时受同一个map-measurable A_n控制；fixed-size law的conditional Dyck independence、one-scale second moment与small-side mass interpolation；A_n tight必须保留在同一joint extraction，最终pathwise sup A_n有限。
  - 可能接口：若column height/degree与quadratic occupancy只控制典型columns，另设对exceptional/sparse coefficient support的inverse or operator-norm estimate，再作small-error+structured-rank分解；强调same object、uniform-over-source quantifier，避免只对每fixed row的估计替代所有rows的rank保障。
  - 缺桥：mass、Green operator与Q(Z) rows无现成对应；sparse support不等于height小，期望/概率控制不能直接推出deterministic π分母序列。需显式candidate operator与uniform arithmetic sparse-error bound；仅作为远域construction clue保留。

### 边界

- 完整source首读并不独立验证此random walk theorem或其companion metric/conformal papers；numerical-array conditional transfer、nondegeneration local germ步骤和all-point annular lemma等作为稿件论证阅读，未重证外部dependencies。
- 该稿不讨论π或Diophantine approximation；以上候选接口都需新增arithmetic/common-basis/height/analytic桥，不推出π的BA/nonBA，q²|π−p/q|目标仍未证明。
- exact period=1是本文finite cohomology pairing identity，不是已获得同一个nonzero integer determinant；all-minors source denominator和finite-map probability estimates不能无桥转移至Q(Z) row rank、column heights或quadratic occupancy。
- source closure：README目标linear-clock-random-walk-tree-weighted-planar-maps.pdf对应build/paper.tex；main依次input introduction、inputs、nondegeneration、energy、sphere、green、dynamics、bibliography；introduction再input history。全部11文字files与inline figures已读，无独立appendix或外置figure文件。

## 45. An infinite finitely presented periodic group（family 247）

状态：`full_body_read`；9 files，3880/3880 lines，175187 exact Git-blob bytes。Mapping/closure complete；unread files为空。

构造同一个无限有限呈示graded unital F2-algebra R=F2⊕I，其中每个有限矩阵M_d(I)的每个元素都nilpotent（无共同index）；G=St_12(R)是普通有限呈示、无限、periodic的group，且property(T)、nonamenable。I作为nonunital algebra普通有限呈示、无限维、nil且Jacobson radical但非nilpotent；R有限呈示、无限维且algebraic。所有elementwise性质由有限relations导出，不额外添加variety identities；不宣称bounded group exponent。

### 证明架构

- 仅level0有限homogeneous binomial/zero word relations。上层letter编码为q_h core+q_h palette+stop，bounded local reversible verifier以addressable accepting tableau检查下一层rule；Kleene fixed point使table-membership program自模拟，D先定、h0后定。clean blocks的ownership/committed slots确保任意替换只模拟真正upper relations，非零可反复编码成任意高degree。
- arbitrary words不假定从clean可达：完整forged full-head states必须可rollback至acquisition或显式zero；非zero迫使release。语法与known-endpoint unary zeros限制headless fragments；letter-recognized cuts给唯一bounded-end token decomposition，isolated token normalization至≤k个upper letters；所有长非零word含long distinct-color eligible palette run。
- 在characteristic2，m>d²个numerical d×d matrix的全排列和是alternating multilinear并为零。按color-erased word选首long run，固定跨端chunks、自由排列完整chunks，所有words有同一algebra image，segmentation与identifier multiplicity保留，故每个orbit矩阵权和零。cut reduction及prefix automaton将fixed test推至更大palette，D>10(c+1)(k+1)保证palette增长快过dimension²，有限层后反向传vanishing。
- 对任意T∈M_l(I)以finite nonempty word expansion，用K=F2(z)给每chunk乘central z；eventual word-length cancellation先得P_n(z)=0，scalar extension injective后提[z^a]，T^a=Σ_{n=a}^{ab}[z^a]P_n=0。nonzero degrees无限给R无限，positive-degree presentation给I有限呈示和上述nil/radical/algebraic后果。
- 有限群呈示用length≤2 roots和ring relators，spare coordinates构造所有word roots，验证split/index independence、all commuting roots、involutions及empty-word products，再由additivity得到全部Steinberg relations和two-sided ideal closure。R无限使E_n无限；matrix-nil和finite GL_n(F2) quotient使GL_n(R) periodic，尚须J_n torsion。
- frames的finite labeled independence complexes给homology/simple connectivity，small-coordinate full matrix presentation与faithful rectangular strips把matrix action提升为Aut(St)，所以J_n central。Q-central extensions在locally finite unipotent subgroups唯一分裂，root relations提升；自然识别Hom(J_n,Q)与H2(E_n;Q) dual。unordered frame orientation的transposition在Q上为−1，杀l=2..4 stabilizers；double-complex消元保留标准inclusion，n≥12使任意rational hom extend。
- grading R→R[t]将kernel word连接u(1)=u和ground-field u(0)。N大于显式word support bound时，cyclic与truncated shift substitutions之差全部收集于disjoint-endpoint faithful strip，matrix image=1迫使同一Steinberg word equality。N取2-power，cyclic conjugate upper-unipotent给stabilized u^N，truncated给field-kernel torsion。Q-dual extension降回rank12，abelian rational detection给J_12 torsion，最后lift group periodic。

### 有条件构造线索（含远域线索）

- **有限局部规则必须同时控制clean起点和任意伪造中间态；以可逆rollback或显式zero消除后者。**
  - Source: build/sections/finite-rules.tex:472-616 (context,verification,rollback)；build/sections/finite-rules.tex:840-1047 (faithful simulation/arbitrary-word token normalization)
  - 前提：需显式finite alphabet/local syntax、transaction ownership、固定payload的可逆tour、在所有可表示register tuples上保留edge，semantic zero不能删去rollback edge；known-endpoint unary tests须独立于缺失neighbor。
  - 可能接口：若未来Q(Z)候选row由递归/automaton生成，可借此分离：合法seed生成的row非零/independent性，与所有可能path/partial state都保持integrality、height、同一determinant identity。把证明不变量编译成局部certificate并处理forged state，防止只验证规范构造路径却漏掉正规化后的任意states。
  - 缺桥：本文characteristic2 word algebra和unbounded simulation scale没有给integer rational-approximation rows、Q(Z) row rank、column height或q²|π−p/q|估计；需要新构造使局部certificate控制实际同一非零整数determinant及解析smallness。当前只是结构方法，未转移。
- **用额外central indeterminate分离不同计数，先在field extension证明整体vanishing，再凭injectivity逐系数返回原对象。**
  - Source: build/sections/cancellation.tex:446-535, especially ca:length-polynomial/ca:power-coefficient
  - 前提：每个fixed test有限支撑、chunk labels非空且有统一finite length bound b；R[z]→F2(z)⊗R injective；eventual vanishing适用于同一个含z test，而不是每个z单独取变化threshold。
  - 可能接口：在候选common-basis determinant construction中可用形式variables分别追踪column degree、coefficient height contribution和quadratic occupancy，要求同一generating identity逐coefficient取出的项仍对应同一个integer minor。此法提醒将factor count和word/degree count分开，避免从总和小或总和零误推每层。
  - 缺桥：这里所得是char2 exact zero和nilpotence，通常会消灭nonzero minor，不能直接供π的integer lower bound。须另证所取coefficient非零、统一Q(Z)rank及专门analytic bound；未给出任何这类桥。
- **不同realizations的误差先收集到一个map-faithful子对象，再以共同image identity证明source层exact equality。**
  - Source: build/sections/kernel-torsion.tex:179-267 (triangular faithfulness/evaluation)；build/sections/kernel-torsion.tex:274-395 (wrap-removal; N>2(m+1)δ)；build/sections/stability.tex:696-711 and kernel-torsion.tex:468-481
  - 前提：correction source/target bands在全过程不交；每次conjugation只扩张δ、无opposite roots、有限m；faithful strip而非仅matrix quotient；rational dual extension精确保留standard inclusions，且kernel已central/abelian。
  - 可能接口：对于同一integer determinant与analytic linear form的两种basis计算，先找一个受控error module，证明evaluation map在该module injective，从而获得exact same-object identity；类似rational-dual延拓可作为common-basis/rank transport需要的自然性检查。
  - 缺桥：Steinberg root strip injectivity与Q(Z) polynomial-row evaluation没有既成对应；rational-dual extension不提供定量denominators/height或determinant lower bound。需为具体π rows构造faithful error module并证定量bounds；不假定成立。

### 边界

- 这是pinned source的完整首读记录，非对声称解决finitely presented Burnside/nil问题的独立证明验证；引用定理、local rule table exhaustiveness、frame/homology chain arguments未独立重建。
- 全文无π、Diophantine approximation或badly approximable主张；上述跨域ideas只保留条件接口与缺桥，不能推出π是否BA或nonBA，q²|π−p/q|目标仍未证明。
- characteristic2 permutation cancellation依赖exact zero，与需要同一nonzero integer determinant的目标方向可能冲突；matrix tests任意dimension≠Q(Z) row rank/height/occupancy theorem。
- Source closure按README目标paper.pdf对应build/paper.tex及其全部inputs映射；included figure同时完整阅读，正文inline correction-strip figure也在kernel-torsion区间内阅读。无未读text dependencies；PDF binary不作为独立body来源重复计数。

π的BA/nonBA与`q²|π−p/q|`目标仍未证明。当前共同审查对象仍是同一个非零整数determinant与analytic smallness是否确实重合、Q(Z) row rank、candidate common bases、column height加quadratic occupancy；上述稿件没有自动完成任何迁移。
