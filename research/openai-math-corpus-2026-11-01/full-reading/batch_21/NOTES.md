# Batch 21：队列40、42、44的完整源码首读

固定上游commit：`4eaee87548250b2682affda7d05d947c05554bae`。这三篇按全库队列阅读，没有relevance筛选。完整覆盖只指实际读完列出的正文、图文字和bibliography，不等于独立验证论文主张。

| 顺序 | 稿件 | 文件 | 已读源码行 | bytes |
|---|---|---:|---:|---:|
| 40 | A smooth three-torus diffeomorphism with simple Lebesgue spectrum | 9 | 2371 | 102856 |
| 42 | A Complete Local Domain without a Small Cohen–Macaulay Module | 4 | 1266 | 53763 |
| 44 | Rigidity of the Turing degrees | 10 | 899 | 35670 |

全部23个文本blobs均以git blob SHA校验；原始字节存在source_context/fullread_round03_21/{position}/。apply_patch补入的唯一末尾LF仅在删去恰好一个LF恢复目标SHA时才删除。详细byte_audit.json记录修正；所有source_files read_ranges覆盖1..N，source_mapping与closure均complete，unread_files为空。

## 40. A smooth three-torus diffeomorphism with simple Lebesgue spectrum

README链接根paper.pdf；目录唯一manuscript entrypoint为build/main.tex，其title/date/abstract相符。main48–53依次输入全部六个section，57使用build/references.bib；六section无其他外部稿件include，无appendix或figure input。全部仓库文本文件均列入本receipt。

存在保持标准归一化体积的C∞微分同胚T:T³→T³及实f∈L²₀，使{f∘Tⁿ:n∈ℤ}成为整个L²₀的正交归一基，故Koopman限制具有multiplicity one的Lebesgue spectrum且T ergodic。进一步主张mixing、zero entropy和全Lyapunov exponents为零；all-order mixing另依赖同collection的multiple-mixing论文。

证明结构：

- 以S(θ,u,a)=(θ+φ(u),u,a)和good arc上的有限disjoint angular packets为每次有限问题的精确谱模型。单packet密度为(mφ mod1)_*(|A|²du da)，不同color的cyclic spaces正交；谱乘子作用和density-norm estimate供后续误差比较。
- 令H′_N=N(c+φ+Δ_N)、σ_N=N⁻¹+N⁻²J，用C_N(θ,u,a)=(u,a,θ−H_N(u))定义新twist。固定c、J、cutoff和物理chart以后N→∞同时控制map及inverse的任意C^k proximity。stationary phase定位j/N≈m(c+φ)，给squared coefficient的summed o(1)，先在mass level删除rounding，再以K(s,Nt−j)计算fine temporal density。
- 大整数c使不同old indices的stationary images精确分离；spare coordinate提供可指定jitter law。fast-variable averaging只对有限固定慢权重给积分比较，不从TV推出pointwise fast-diagonal比较。neutral refresh的宽线性jitter通过Riemann-sum error将每color在modP各residue的谱密度置为p_r/P，同时保持物理向量近似。
- 构造严格正q0，Fourier support在小band且实轴tail为1/(1+R²)。在R×ℤ/P的joint law中加入小Schwartz perturbation γ C_P(h)g(R)。P由有限小素数构成，筛掉所有higher surviving harmonics；γ依赖normalized-index band与d但不依赖未来scaleB。c可依实际有限indices增长，N最后选，neutral-refresh后各color获得1+2λRe(d_r(z)e(−Nz))。
- 每block用实际输入p_r、tags b_r和β=∑b_rp_r/p选择固定smooth q，控制∫√p|β−q|²<δ。先定d、P、λ和n使δ/2≤nλ²≤δ，再anchor并做n个refresh-transmit pairs。有限Laurent phase tests选快速增长L_i，把product-Haar下的有限多项式moment比较转到实际一维diagonal；不假设整个phase diagonal随机独立。
- 产品model中total density ratio R的一阶cross contribution为2nλ²|β−q|²，故Hellinger loss≤300δ²+ν；D=q+(nλ)⁻¹∑e(L_i z)在每color tilt下均值恰b_r、variance≤2/δ，对previous bounded estimator的cross term逐z消失。K块后prediction error²≤5/(Kδ)，Hellinger loss≤400Kδ²，可使Kδ大而Kδ²小。
- 将任意observable的θ、u、a三个方向全部转为非零packets；small perturbation F+tH处理ratio zeros，再用rectangular constant tags近似H/F。两次neutral passages把原支撑分离转换成纯angular-index分离且moduli独立于spare coordinate。prediction后以p^−1/2 whitening，并在{p=0}补flat vector，恢复exact density1且L2改变量由∫(√p−1)²给出；prediction改写为新generator的L2谱乘子再近似成finite Laurent polynomial。
- 重复dense smooth tests，各次保留旧finite-polynomial predictions，令good-arc complement、map/inverse C^k errors、generator L2 errors可求和。limit是volume-preserving C∞ diffeomorphism，flat density在极限给bilateral orthonormal orbit，重复tests使其cyclic span为整个centered space。后果分别用Riemann–Lebesgue、外部multiple-mixing、Rokhlin entropy与Pesin公式。

构造候选（均有条件，不作π结论）：

1. **以当前weighted mean消掉一阶总损耗，保留各color信号的二阶可积预算**

   锚点：build/sections/04-blocks.tex:121–148 (eq:block-centering)；build/sections/04-blocks.tex:290–369 (eq:block-second-moment-bound, eq:block-hellinger-step)；build/sections/04-blocks.tex:372–505 (eq:block-estimator-moments, eq:block-total-prediction)。

   机制：把每color偏置写成a_r=b_r−q，以当前真实densities的β作center；总损耗中一阶项变成|β−q|²，单块代价O(δ²)，可累积的信息budget为Kδ，prediction error²为O((Kδ)⁻¹)。

   前提：必须有非负可归一化weights、同一个作用于全部colors的center、可计算或可控制的二阶moment及cross-term cancellation，且每block可从实际输入重新center；单纯形式上平均误差不够。

   可能接口：对候选common basis下的列height/二次occupancy代价，尝试先把columns分成保留tag的组，建立total analytic mass或Gram quantity的weighted mean，再考察centered modifications是否能将总解析损耗降为二次而仍积累可区分的独立列信息。Kδ与Kδ²分开列预算，避免把occupancy当可加。

   缺失桥梁：稿件仅在positive spectral densities和可任意小smooth perturbation中实现这些moment identities。尚无π线性形式的正权表示、保持整数/同一nonzero determinant的合法centered操作，也未把prediction信息转换成Q(Z) row rank或列height优势。没有推出π的BA或nonBA。

2. **有限测试决定参数下界，实际复杂对象生成后才选最终scale，消除循环的uniformity要求**

   锚点：build/sections/02-passages.tex:54–128 (lem:passage)；build/sections/04-blocks.tex:64–108 (lem:phase-tests)；build/sections/04-blocks.tex:150–288 (anchoring and block-error-propagation)；build/sections/05-completion.tex:346–392 (strong Koopman continuity and retention)。

   机制：每stage只保留有限测试；先固定block length、coefficients、cutoffs和误差放大常数，再用最终自由整数同时满足所有lower bounds。对旧结论保留有限open conditions，末尾diagonalize。

   前提：最终参数必须真正自由且所有要求同向成为有限lower bounds；有限测试系数在最终scale选择前固定；后续扰动对已有证书有已证明continuity或quantitative stability。

   可能接口：对Q(Z)上的多个候选row spaces或共同basis，先保存有限子式/有限valuation或degree测试，再安排一个新的整数specialization同时避开有限零点并满足analytic-smallness阈值；每轮明确区分已经固定的符号数据和依赖新参数的actual coefficients。

   缺失桥梁：π目标同时包含height upper bounds、integrality/divisibility和analytic decay，这些未必都是lower-bound条件；specialization可能增加列height且旧nonzero minors未必在同一basis下保留。需要精确共同格、避零与height的同时估计。

限制：

- 完整源码首读不等于独立验证Banach问题的解答。所有stationary-phase、signal、block估计按文中陈述理解；未独立复核所有analytic uniformity或外部引文。
- all-order mixing显式依赖另篇MultipleMixing2026；Rokhlin finite-multiplicity entropy及Pesin公式属于外部定理，读到参考文献不计作已读那些文献正文。
- construction_ideas是有前提的迁移候选；物理L2 closeness和Hellinger loss不能自动替换为整数行列式的非零性、Q(Z)满秩、height控制或π近似结论。

## 42. A Complete Local Domain without a Small Cohen–Macaulay Module

README链接根paper.pdf；唯一article entrypoint是build/paper.tex，title/date与README一致。paper629输入surface.tex，696–697使用references.bib；surface中的TikZ branch-incidence figure内嵌，无额外figure/appendix文件。所有仓库文本文件均列入receipt。

存在含ℂ且residue field为ℂ的三维complete Noetherian normal local domain，没有任何非零finitely generated maximal Cohen–Macaulay module。更一般，smooth integral projective complex surface X和ample globally generated H的completed section cone若有small CM module，则K_X²−2c₂(X)≤6H²。文中构造n=7 cover使H²=1715、K²−2c₂=13377>10290，且未完成vertex local ring也无small CM module。

证明结构：

- Global generation选三个无共同零点sections，定义有限f:X→P²且f*O(1)=H、degree h=H²。full section ringT有限于S=ℂ[x0,x1,x2]，completionR=T⊗Sℂ[[x0,x1,x2]]=∏T_j，最低非零homogeneous term证明domain。total line bundleW0=Tot O_X(−H)经regular local base change构成proper regular model W，off exceptional X与punctured cone同构，不需T在degree1生成。
- Γ(W,O)=R且W integral regular，故R normal。若M small CM，parameters x0,x1,x2是regular sequence，Auslander–Buchsbaum使M在A上free rankN；punctured cone regular，M局部free。E=(p*M)**的restrictionP到X torsion-free positive rank；g*E在base blowup V的D外为trivial rankN。
- HN factor slopes按c1·H/(rank h)规范化。取最低slope quotientG，沿X elementary modification使restriction有0→G(H)→P′→B→0。若width>1，maximum slope不增而average slope每次至少增加1/r，故有限步使所有factor slopes在长度≤1区间；punctured bundle不变。
- Divisor lattice lemma以I^bL⊂L′⊂I^−bL清分母并twist到L′⊂L。Q=L/L′有限I-filtration的Q_j是O_D(j)^N的quotient，故0≤n_j≤N及c1(Q_j)≥j n_j。derived restriction给[F]=[O_D^N]+∑([Q_j(1)]−[Q_j])，于是c1(F)=∑n_j，2ch2(F)≥∑(2j+1)n_j；quadratic occupancy identity推v(F)=2ch2(F)/N−(c1(F)/N)²≥0。用于F=f*P且N=rh。
- 对每torsion-free slope-semistable G，比较X与P²的twisted Euler characteristics得normalized c1/ch2关系；Bogomolov上界2u ch2(G)≤c1(G)²，再令C=c1(G)/u−K/2用Hodge index，得v(f*G)≤1/4−(K²−2c2)/(12h)。HN mixture identity追加factor means的weighted variance；width≤1令该variance≤1/4，故0≤v(f*P)≤1/2−(K²−2c2)/(12h)。
- Blowup P²的四点complete quadrangle，十条branch curves的SNC divisorΔ=−2K_B，K_B²=5。五个line-ratio n-th roots经valuation independence给degree d=n^5；每node处两branch valuation vectors independent，因此étale局部为(t1,t2)↦(t1^n,t2^n)，使normalization smooth且ramification exact。
- Canonical divisorK_X=(1−2/n)π*(−K_B)与branch strata Euler characteristics(2,−10,15)给K²−2c2=d(1−10/n²)。ramification divisors构造J及E，six line sections使J globally generated，three opposite-pair products使2J−E globally generated，H=3J−E且nH=π*(−K_B) ample，H²=5n³。n=7违反obstruction，排除完成环所有finite CM modules；faithfully flat completion保regular sequence，进一步排除affine-local modules。

构造候选（均有条件，不作π结论）：

1. **有容量约束的height layers给真正二次occupancy下界，保留层间cross terms**

   锚点：build/paper.tex:467–523 (lem:lattice; eq:layers, eq:Tor)；build/paper.tex:510–522 (quadratic occupancy identity)。

   机制：对有限层rank n_j∈[0,N]，N∑(2j+1)n_j−(∑n_j)²=∑n_j(N−n_j)+2∑_{j>ℓ}n_j(N−n_ℓ)≥0。权重2j+1是j²到(j+1)²的离散差；逐层容量而非独立可加成本产生总占用平方。

   前提：需要同一ambient rank N、共同filtration、每层确为被N个generator控制的quotient，及可把valuation/height或degree成本写成该层weighted sum；不能混用来自不同basis的n_j。

   可能接口：为π项目候选common bases上的column heights建立shared flag，令n_j计超过第j层的独立columns或leading terms。若同一非零整数determinant的degree/height可从这些layers得到∑(2j+1)n_j，则把总occupancy S=∑n_j与最小quadratic cost S²/N同时纳入analytic-smallness比较；该式显式处理非加性的跨层项。

   缺失桥梁：尚未构造Q(Z) row spaces的共同integral filtration，未证明层ranks给目标determinant的实际height成本或收益，也未把capacity bound和π线性形式解析衰减连接。这里只提供条件性离散不等式，不能由sheaf Chern数字直接推算π或合并不同basis的最优height。

2. **保持外部对象不变，用离散递增势能把内部斜率宽度压到一个单位**

   锚点：build/paper.tex:402–443 (lem:balance, eq:transform)；build/paper.tex:605–626 (eq:mixture and variance≤1/4)。

   机制：每次对最低slope quotient作elementary modification，把它提高1；max slope不增，而average每次至少增1/r，因此有限终止于width≤1。最终用mixture variance而非假设整体semistable，控制不同pieces合并时额外quadratic cost。

   前提：必须有合法操作精确保留所需外部bundle/格数据；固定rank、离散step下界和被初始max控制的potential；mixture identity须在同一normalization下成立。

   可能接口：若Q(Z) row modules可在保持目标子式/整数determinant与解析evaluation的条件下作elementary lattice modifications，可尝试平衡各column的valuation或degree widths，并对不同row blocks的mean height保留variance correction，寻找common basis的finite normalization procedure。

   缺失桥梁：文中操作保持punctured sheaf，不等于整数unimodular row/column operation；提高低slope可能增加目标height，解析smallness未必不变。需要证明对应合法操作、同一determinant保持和variance-to-height公式，才能应用。

3. **局部valuation向量独立性作为构造前的精确rank证书**

   锚点：build/surface.tex:135–217 (lem:cover-smooth)。

   机制：先以五条divisor valuations证radicands独立；再在每node检验两branch vectors independent，才将generator equations可逆地化为coordinate power maps，区分whole generic algebra与subextension。

   前提：具有明确valuation homomorphisms、同一finite-field coefficient space及可逆generator change；局部rank证书覆盖所有可能singular strata。

   可能接口：当候选common basis需同时满足多个local places的independence时，可将local valuation矩阵列为逐点证书；显式区分全row module的rank与某个subextension/submatrix rank，防止用局部可逆化冒充全局共同basis。

   缺失桥梁：ℱ_n上的cover rank并不推出Q(Z)上的row rank、整数determinant非零或archimedean估计；需独立构造valuation-compatible module和global lifting。

限制：

- 源码首读覆盖完整；没有独立验证small Cohen–Macaulay conjecture反例，尤其未独立证明extension restriction torsion-freeness、divisor lattice lemma或所有外部surface-theoretic dependencies。
- Bogomolov、HN、Hodge index、RR、Stacks tags、normalization与ampleness引用均按文中使用理解。bibliography全读不计作这些外部文献正文已读。
- π迁移仅为条件性接口；sheaf slopes、Chern variance和finite-field valuation ranks不能直接当作π线性形式的height/rank/analytic证据。

## 44. Rigidity of the Turing degrees

README链接根paper.pdf；build/main.tex与title/date相符，6输入glyph-mappings，45–49依次输入five sections，50输入references.tex；recovery139输入figures/bounded-recurrence。glyph-mappings1的glyphtounicode为TeX distribution字体Unicode映射，无仓库manuscript文字，不列为未读稿件文件。无appendix或其它repository text sources。

在ZFC中，full Turing-degree ordering的任意order automorphism都是identity，无definability或regularity假设。中间结论：对Borel且countable fibers、满足F(x+y)≤_T F(x)⊕F(y)的F:ℝ→2^ℕ，每irrational t∈(0,1)在comeager多(u,v)上可从四个指定rational-affine arguments的F-values恢复C(t)。另以Slaman–Woodin外部equivalence推出full second-order arithmetic的无参数biinterpretability。

证明结构：

- 用strict rational cuts C(x)表示ordinary reals，点态计算rational affine combinations；rational output时其cut computable，irrational时比较终止。每noncomputable degree由某irrational t∈(0,1)代表。Order automorphisms保finite joins和least degree。
- 调用Slaman–Woodin 2005 manuscript的all-input arithmetic total representative Fhat，Borel graph criterion使它Borel；F=Fhat∘C的fibers注入单个countable Turing degree。degree-join preservation使F(x+y)≤_T F(x)⊕F(y)，不要求uniform reduction。
- 对least successful addition-program index e(x,y)，correctness是Borel，故e是Borel。利用Borel functions在comeager restriction上的continuity，取F|D、e|E连续；对countable H=ℚ+ℚt要求所有仿射shift位于D/E，因相关二维affine maps为homeomorphisms而得comeager G_t。
- 固定(u,v)∈G_t后，e在两个中心的relative continuity给固定程序p,q及小ρ，使两次addition在所有H内小offset上有效。countable fibers使u两侧可选F-values不同的endpoints，从F|D连续性选一bit与rational δ，使低端strip强制up-step δt，高端strip强制down-step δ(t−1)。
- recurrence s_{n+1}=s_n+δ(t−r_n)、s0=0由endpoint bit确保s_n∈[b−,b+]⊂(−ρ,ρ)。four oracle tuple由F(u)、F(u−v)、两个increment values组成；递归构造program index P_n，先调用p到v neighborhood，再用q返回u neighborhood，从tuple计算branch bits r_n；不把t,u,v或数值s_n当额外oracle。
- telescoping给s_N=δ(Nt−∑r_n)，故A_N=N⁻¹∑r_n满足|t−A_N|<1/(Nδ)。δ rational使误差界effective，irrationality保证每rational query最终分离，tuple计算C(t)。δ、程序indices和finite constants允许依(u,v,t)非uniform地存在。
- 各tuple argument的degree≤d(t)∨d(u)∨d(v)，因此recovery后apply π^−1得π^−1(d(t))≤d(t)∨d(u)∨d(v)在comeager pairs成立。relative category avoidance对预先固定的Y代表π^−1(d(t))排除Y不≤C(t)，去掉auxiliary degrees。对inverse automorphism重复并用antisymmetry得identity。biinterpretability再依Slaman–Woodin外部equivalence。

构造候选（均有条件，不作π结论）：

1. **bounded internal state加可计算branch sequence提供显式discrepancy账本**

   锚点：build/sections/recovery.tex:90–138 (eq:recurrence, eq:invariant)；build/sections/recovery.tex:142–205 (program composition, eq:error)；build/figures/bounded-recurrence.tex:21–25。

   机制：只需在区间两端保证branch把state推回内部，中间不要求单调；telescoping产生N t−∑r_n=s_N/δ。branch bits可由固定四oracle逐步计算而无需计算数值state，boundedness给O(1/N) approximation error。

   前提：需要真实bounded-state invariant、固定可执行branch rule、rational δ及不偷偷使用目标实数的新oracle；branch程序沿全部实际states合法。

   可能接口：对π−3或其它参数构造integer-valued counts时，可将每步rounding error放入精确state ledger，尝试让多个rows共用branch更新，并以state confinement显式控制误差，区分symbolic degree/recoverability与数值逼近速度。

   缺失桥梁：原构造只给|π−p_N/N|≤C/N，因此N²误差至多O(N)，不触及BA下界或nonBA的liminf=0。要得到1/N²级或更强，需要额外small return s_N=O(1/N)或special subsequence，且同时保证同一个integer determinant非零与height控制；文中没有提供。π可计算，Turing-degree结论本身对π为trivial。

2. **用countable-shift稳定域把点态存在的操作局部固定，再以avoidance去掉auxiliary witnesses**

   锚点：build/sections/recovery.tex:40–88 (least-index Borel map and eq:all-shifts)；build/sections/preliminaries.tex:94–144 (lem:avoidance)；build/sections/rigidity.tex:28–47 (fixed Y before pairs vary)。

   机制：Borel least-witness selector只在comeager restriction连续；将countable future shifts全部拉回相交，得到沿整条实际recurrence固定的两个程序。最后固定目标Y，再对全部可能program indices取countable union的meager obstruction，从comeager成功去掉auxiliary degrees。

   前提：需要精确Borel/Polish结构、countable family的future operations、pullback preserves meagerness的open或homeomorphism条件，以及独立的relative avoidance theorem；target必须在auxiliary parameters变化前固定。

   可能接口：若common-basis候选可编码为countably many symbolic operations，而每一步所需rank/minor certificate的有效参数域具有可证明regularity，可借此思路先选择一个对全部预定有限/可数更新稳定的域，再讨论auxiliary choices的消去。对determinant项目特别提醒必须先固定所比较的同一整数对象，再变化auxiliary parameters。

   缺失桥梁：category genericity不等于integer parameters中有点，更不提供effective height、uniform scale或analytic smallness；Q(Z)满秩和common-basis selection没有现成的category-avoidance对应物。文中非uniform existence不能直接作为可实施的π辅助行构造。

限制：

- 这是完整源码首读，不是独立证明Turing-degree rigidity。外部all-input arithmetic representation是关键前提，稿件特别区分generic-input版本；未独立核验2005未出版manuscript的准确scope。
- 未独立核验biinterpretability equivalence及其它外部references；读完references.tex不等于读完引用文献。
- four-value recovery是pointwise/nonuniform Turing reduction，无运行时间、height或Diophantine quality估计；bounded recurrence只给O(1/N)数值误差。不能由computability rigidity断言π的BA/nonBA或行列式秩性质。

π项目的未完成边界：以上没有建立q²|π−p/q|的BA或nonBA行为，没有把不同整数对象的非零性与解析小量混作同一determinant，也没有从局部rank或category genericity推得Q(Z)上的common basis。42的occupancy式可直接作为条件性离散工具保存；40的Kδ/Kδ²机制和44的bounded-state/error账本保留为较远构造来源。
