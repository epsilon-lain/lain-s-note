# Batch 20: queue positions 37–39

Fixed upstream commit: `4eaee87548250b2682affda7d05d947c05554bae`. All three manuscripts were read in ascending queue order, without title filtering. Reading is not independent proof verification.

Closure: 35 text blobs, all exact Git blob SHA and byte-length verified; every source line was actually read in untruncated tool outputs. Position 38 also has one external figure PDF: exact SHA verified, all text extracted, then its single rendered page visually read. Main paper PDFs were not page-by-page compared to TeX.

## 37. Brennan's conjecture and sharp inverse-square integral means

Status: `full_body_read`; complete source/input closure. 12 text files, 2620 lines.

对整个 normalized schlicht class S 及每个 ε>0，有与 f 无关的 Cε 使 M_{−2}[f′](r)≤Cε(1−r)^{−1−ε}（1/2≤r<1）；B_S(−2)=1，且 B_S(t)=|t|−1 对所有 t≤−2。任意单连通平面域 W（球面边界至少两点）到单位圆盘的 conformal bijection φ 满足 ∫W|φ′|^s<∞，4/3<s<4，等价于任意 univalent disk map 的 ∫D|f′|^t<∞，−2<t<2/3。Koebe 给两端发散，另推出文中指定 homogeneous Sobolev pullback 范围。

Proof architecture:

1. 用 F(i)=0,F′(i)=1 的 compact half-plane class K、Q=1/F′ 与 affine rerooting T_z 构造正算子 A_zΦ=|Q_F(z)|²Φ(T_zF)，L=(A_{−1/2+i/2}+A_{1/2+i/2})/2；cocycle 使 L^n 成为 dyadic grid sums，ρ=lim||L^n||^{1/n}。紧集 distortion 把整圆 inverse-square means 控制于 ||L^n||，故只需 ρ≤2。
2. 反设 ρ>2，β=log₂ρ。正 resolvent 在 ρ 处 divergence 给 probability eigenmeasure μL=ρμ；grid endpoint cancellation 与 ρ>2 给 horizontal period2。先水平平均、再对 v^β μ̄A_iv 作一个 log-scale period 平均，得对全部 affine z 和 nonnegative Borel Φ 的精确 weighted law E[|Q_F(z)|²Φ(T_zF)]=y^{−β}EΦ。微分 constant-test identity 得 E a=−β/2、E b=0、E|q|²=β(β+1)/4。
3. k=(β+3)/4>1，Vξ=y^k/|F−ξ|。affine inverse z# 的 Jacobian y^{−3} 与 conformal area change 得 E∫S0 y^{β+1}|F|^{−4}≤256π；submean 给无穷远逃逸，integration by parts 给 compact target 的 horizontal preimage length 一致界。每个 dyadic witness cell 的 target radius O(h^k)、个数 O(h^{−1})，面积 O(h^{2k−1}) 可求和，故 a.e. target 的正 superlevel compact。
4. G_F=F−iyF′/k 的 normalized J 满足 k²J=k(k−1)+(2k−1)a+|q|²，故同一 P 下 k²EJ(i)=k(k−1)>0。critical points 正是 G_F=ξ，Hess logV 的 trace<0、det=−k⁴J/y⁴，所以 J>0 是 saddle、J<0 是 maximum。regular target 加 compact superlevels 使各正阈值上 critical points 有限；Jordan domains/index sum Cλ=Mλ−Sλ+B 且 Cλ≥B，给每阈值 M≥S，pole 不计作 finite maximum。
5. 分别按 V 降序、(x,y) lex tie-break 给 saddle/maxima 排名，同 rank 配对得 injective、V 不减；local inverse branches 加 countable activation 证明 Borel，正 affine scale 保序给 equivariance。area formula 把 source mass ωJ± 变为 fiber sums y^{β−1}|Q|⁴=(V/k)⁴，pairing 给正/负质量比较。
6. 先用 area formula 与 weighted law 把 target-a.e. goodness 转成 root fixed target −i/k 的零 exceptional positive-J weight。截断 partner displacement 为 D_N，在 B_R,Y 与扩张 B′ 比较，期望 source density恰为 dxdy/y；固定 N、R=Y² 后 Y→∞，再 N→∞，得同一 J 同一 P 的 EJ(i)≤0，与上式相矛盾。由 ρ≤2 得 uniform radial estimate，Hölder 给 negative interior area，weighted derivative area + function-power bound 给 positive t<2/3；Koebe 验证 sharpness。

Construction ideas and speculative π interfaces:

- **用全部阈值上的累计优势，把两个存在性列表配成同一个不降权证书。** Premises: 两组有限或可数对象 A、B 带正 score；每个 score≥h 的子集有限，并对每个 h>0 有 #B_{≥h}≥#A_{≥h}；对象带在允许变换下保序的确定性 tie-break。 Mechanism: 两组分别降序排列，同 rank 配对；任何 rank j 的 A 对象上方有 j+1 个 A，累计优势保证至少 j+1 个 B 具有不小的 score，因此存在 injective score-nondecreasing assignment。论文还用 nonsingular local branches证明其 Borel，正 affine scale保持 lex顺序。 π interface: 若 π 项目能在同一个 Q(Z) 上，把 formal nonzero minor/common-basis 证书列成 A，把满足列高/occupancy预算的整数评价证书列成 B，并证明各阈值的累计优势，就可避免‘存在非零对象’与‘另有小对象’分属不同基；rank配对是明确的同对象选择机制。 Missing: 目前没有这些两列表、累计优势或允许变换保序的实现。分别有 row-rank 与低高元素不推出累计优势，更不自动给 common bases。配对本身也不保存整性、determinant 非零、columnheight 或 quadratic occupancy；这些须在 score定义和对象关系中证明。
- **把假定的超额增长转成精确协变概率律，再让一个可积 observable 产生相反符号。** Premises: compact metric 状态空间上的 bounded positive transfer operator；weights 满足精确 multiplicative cocycle；离散采样存在明确 endpoint cancellation，且假设 spectral growth严格大于 boundary-error growth；可做水平/尺度 averaging。几何比较还须有同一 law 下的可积 observable 与 negligible boundary cost。 Mechanism: normalized resolvent eigenmeasure→有限period→连续horizontal average→logscale average，得到所有 basepoints 的 weighted covariance；moment computation 给 EJ>0，equivariant mass comparison给相同 EJ≤0。显式 root-null transfer阻止 a.e. target 的论断被误用在固定 target。 π interface: 可探索将 Q(Z) 的 bounded-normalized support/column configurations作状态，bad approximation假设生成增长算子；若能让 determinant 或 rank-defect 的同一加权统计同时受代数非零与 occupancy transport约束，就有比逐个 minor 追踪更全局的矛盾框架。 Missing: 整数列操作未给出 compact scale-covariant action，determinant zeros也没有本文的 critical-index geometry；尚无候选算子、cocycle、可积同一observable 或边界代价估计。这里只是构造方向，不能把 expectation sign矛盾当成已得到的非零小整数矛盾。

Limits: 全文源码首读，不是独立证明认证；未核验外部引文的定理适用条件，未编译或逐页比对主 PDF。 结论保留 arbitrary ε loss；B_S(−2)=1 不声称 inverse-square area 在端点可积。a.e. targets、root weighted null set 和全 K uniform bounds 是不同陈述，记录时未互换。 几何 signed-Jacobian contradiction 用同一 P、同一 J、先固定 N 后做 expanding-box 极限；任何拟议 π 移植需重新证明这些兼容条件。 未证明 π 的 BA/non-BA、未建立实际 Q(Z) 的新 row-rank/common-basis certificate，也未产生受控整系数非零小 determinant。

已逐行读完 README、main、全部 8 个 section TeX、rank-pairing figure TeX 及完整 bibliography；所有文本缓存 bytes 与固定提交树的 Git blob SHA 精确一致。全部仓库正文/图文/参考文献输入闭包完成；没有附录或额外正文文件。主 paper.pdf 未逐页核对。

## 38. Symplectic Balls in Symmetric Polar Products

Status: `full_body_read`; complete source/input closure. 10 text files, 1589 lines.

对任意 n≥2 与任意 origin-symmetric convex body K⊂R^n，U_K=intK×intK° 的 Gromov width 为4；每个 0<c<4 都有 smooth symplectic ball embedding B^{2n}(c)→U_K，无须 K 光滑、strictly convex、unconditional 或 product。由此宣称所有维度 n≥1 的 symmetric Mahler inequality vol(K)vol(K°)≥4^n/n!，常数sharp；n≥3 另给满足严格 volume及pair-capacity条件的 finite ball packings。

Proof architecture:

1. 建立 holomorphic ball principle：f在Ω中唯一共同零点0、|f|=O(|z|^k)、每个 τ=|f|² 的 s<1 closed sublevel compact，则 (V={τ<1},dd^c(|z|²+τ)) 含任意 c<πk 的 standard ball。对 χ(logτ) 取斜率 b，a<μ<bk，regularized max与 μlog(|z|²+δ) gluing；先固定r1再选δ，近零径向form含πa球。
2. τ与glued potential差支撑在一个固定compact sublevel C；Moser vector field也支撑于C，保证非紧V中的完整flow。径向 R(z)=sqrt(1+μ/(|z|²+δ))z 满足 R*ω0=modified ω，tΨ′(t)严格递增给smooth inverse，再由 inverse Moser transport回原form。
3. 显式 F(w)=8/π² Σ(−1)^l w^{2l+1}/(2l+1)²，wF′=(8/π²)arctan w；right semicircle angular imaginary derivative A>0。固定height的radius参数给 ∂v/∂r=|wF′|²/(rA)>0，证明global biholomorphism及中心水平切片，tips为±i；T′=(8/π²)atanh(r)/r→∞。
4. J_k=k²∫0^v |g|^{2k−2}|g′|² 的radius表达 J_k/k=k∫r0^R r^{2k−1}/A dr。主项为πR^{2k}/4；tip bound 1/A≤π/2+C sqrt((1−r)/(T(r)−|t|))，compact part≤C*k*r*^{2k−2}，outer moving-endpoint integral经 H=k(1−ra) 控制≤1+π/e。先r*→1选cutoff、再k大，得全D一致εk→0。
5. 固定有限strip body K={|b_j·x|≤1}，b_j span R^n，取 f_j=g(b_j·z)^k。唯一common zero与compact sublevels满足ball principle。q=x,p=−u−ΣJ_k(b_j·u,b_j·x)b_j 的pullback正是原ω；各J随v非减给 (p2−p1)·(u2−u1)≤−|u2−u1|²，即使V fiber disconnected仍global injective。
6. 同一phase map的 h_K(p)/k≤C_K/k+(π/4)τ+mεk；m,K先固定再k→∞，得到 image⊂intK×S intK°，S=(1+η)πk/4。对c<4先选η使Sc<πk，源ball dilation√S与target momentum /S的pullback factors互消，给真正symplectic embedding入U_K。
7. 从K°有限net得到K⊂K′⊂(1+α)K，(K′)°⊂K°；先fix finite K′再k与η，position /λ与source dilation√λ再互消，覆盖所有K。supporting conjugate pair把U_K放入area4 square cylinder，square-area-preserving disk map+nonsqueezing给上界4。球volume preservation→Mahler，cube/crosspolytope给sharpness；packing corollary调用文中引用的外部ball-packing theorem。

Construction ideas and speculative π interfaces:

- **把移动端点的误差压成无量纲参数，并严格固定‘有限约束先、指数后’的量词顺序。** Premises: 有非负kernel k r^{2k−1}sqrt((1−r)/(r−ra))，0≤ra<1；主系数的excess可由cutoff r*处的 C*/sqrt(L*) 控制，L*→∞。constraint count m 与几何常数C须在k变化前固定。 Mechanism: H=k(1−ra)、y=k(1−r)把所有随k移动的lower endpoint同时变为 ∫0^H e^{−y}sqrt(y/(H−y))dy≤1+π/e；先选r*使outer error小，再选k使compact exponential error小，得到对height/r0/R全部一致的εk→0，再使 C/k+mεk 小于固定slack。 π interface: 对π项目中的columnheight或quadratic occupancy估计，应尝试将随截断变化的最坏support边界写成无量纲拥挤参数，给统一端点界；若先选有限候选columns/constraints再放大degree，可防止 mε_k 因候选数同时增长而失控。 Missing: 这里没有Q(Z)的kernel表示或可用εk；若实际minor大小、column数m也随k增长，固定m结论不足，必须证明速率及m(k)εk→0。该方法本身不说明formal row-rank或common-basis存在。
- **用同一映射的强单调性获得全局唯一性，即使可行域的单纤维不连通。** Premises: 固定x，所有允许u满足每个b_j·u处于同一个可积分的interval；各 J_j(v,t_j) 在完整interval上非减。设 p(u)=−u−ΣJ_j(b_j·u,t_j)b_j。 Mechanism: 对任意u1,u2直接相减得 (p(u2)−p(u1))·Δu=−|Δu|²−ΣΔJ_j Δv_j≤−|Δu|²；不必在可行fiber内连接两点。加上local非退化pullback给global smooth embedding。 π interface: 可探索给受限column系数选择器加入可分的monotone penalties，使同一support配置下的校准参数有唯一解；这样在非零证书与高度优化间移动时，参数选择不因可行域分支而换对象。 Missing: 尚无π项目的合法monotone parametrization、存在性或整点恢复。连续injective map不提供整系数列、可逆基变换或低height整数评价；penalties若改变Q(Z)还需重新证明rank及same-object条件。
- **所有增益与成本在同一个显式构造上做互逆缩放核算。** Premises: 一个构造同时提供capacity<πk的ball与其image momentum scale S=(1+η)πk/4；源dilation pullback是Sω0，目标momentum division pullback是S^{−1}ω0；Sc<πk。 Mechanism: L_S∘Φ_k∘E∘D_sqrtS 这一个复合的pullback精确为ω0，且同一image满足polar bound。近似任意K时也用source dilation与position division相互补偿，并保留严格slack。 π interface: 可把π determinant构造的degree/columnheight/denominator成本记录为同一复合操作的scale ledger：任何用来增强非零minor的列变换，也必须在它实际输出的同一矩阵上算高度和小量，而不是比较两个独立最优构造。 Missing: symplectic reciprocal real scaling不对应整数determinant的不变性；column缩放会引入整性和denominator代价，必须显式清分母、计算determinant scale及quadratic occupancy。没有由此得到新数论不等式。

Limits: 全文源码首读与figure文字/视觉核对，不是独立验证；未核验外部nonsqueezing、Moser背景或packing theorem，未编译/逐页比对主PDF。 每个c<4由一次有限构造取得，不声称capacity4端点embedding；无equality-case classification。主结果限定origin-symmetric K的特定polar product，不推广到任意phase-space convex domain。 全D一致slice bound与固定finite constraint list是实质条件，不能把逐点o(k)或增长m(k)下的误差当成同一结论。 所有π接口是推测性构造建议；未提供Q(Z) actual row-rank、common bases、整数非零小determinant，或可直接使用的columnheight/quadratic occupancy估计。

已逐行读完README、main、6个sections、figure wrapper及references；唯一外部figure PDF已按精确Git blob下载，提取全部文字并render视觉核对所有标签。无附录或额外body。所有35文本blob中的本篇10个均SHA精确核验；主paper.pdf未逐页核对。

## 39. The Factor-Two Hardness Threshold for Vertex Cover

Status: `full_body_read`; complete source/input closure. 13 text files, 2121 lines.

对每个固定整数m≥4，存在从3SAT到explicit finite simple unweighted graph G的deterministic polynomial-time reduction：SAT⇒τ(G)<(1/2+1/m)|V|；UNSAT⇒τ(G)>(1−1/m)|V|。故每个固定1≤α<2的Vertex Cover approximation均NP-hard，除非P=NP不存在固定2−ε因子。论文唯一外部hardness input是perfect-completeness Label Cover value1 vs任意固定σ的gap，不假设Unique Games。

Proof architecture:

1. 每个unordered position pair取独立Label Cover occurrence，j<k分别只得到left/right endpoint question；query不暴露occurrence、map或opposite endpoint。valid local label满足其scope所有internal constraints；compatible selection至多每query一label、overlap一致且满足combined union scope的全部internal constraints。其incidence sums的绝对maximum定义hereditary compatibility norm。
2. graph vertex枚举所有pair seeds与d个grid-weight blocks；S_v为各query增量之和，distinct v,w的edge iff ||S_v+S_w||*≤2t。只需最多2d active blocks，枚举≤(1+r)^{2d} local choices即可exact norm；edge test化为integer 2mζ≤p√d。固定参数使n=M^{binom(d,2)}p^{rd}显式polynomial。satisfying global labeling的compatible functional>{t}选出independent set，midpoint-grid tail给density>1/2−1/m。
3. 若independent set≥n/m，discard最多1个norm≤t vertex；余者sum-vectors A*与其negative相隔>2t。一次固定 F=(dist*(z,−A*)−dist*(z,A*))/2，odd且1-Lipschitz；fresh uniform continuous weights rounding error≤d/p<t/2，得 E_s f=0 每seedtuple、Ef²≥νd。单weight/seed/保留opposite questions的whole-star oscillations分别≤2,4,2。
4. 对weights作orthogonal support decomposition并给每非空support指定j；每seed再按constant/own-question/residual分解。ordinary budgets为Σ|S|energy≤4d、Σ|T|energy≤8d²；star D_k保留opposite questions，仅k的query能变，故||f−D_kf||²≤4，charging invisible incident seed residual给Σ|N|energy≤4d。不要求同seed两endpoint projections commute。
5. 随机h-batch J并冻结outside weights/noninternal seeds。某component在指定j的own-question mean存活需 S∩J={j}、N∩J=∅、T避开J\{j}的internalpairs；conditional probabilities p1=(h−1)/(d−1), p2=(h−1)(h−2)/((d−1)(d−2))给平均variance≥(h/d)(νd−16h−32h²)≥νh/2。取一restriction使ΣEVar_s H_j≥νh/4；剩余seeds/weights仍原product law。
6. 固定batch querytuple时，F∘Σincrements是finite piecewise affine；choose true gradient off proper hyperplanes，exceptions setg=0，始终g∈conv{±compatible incidence vectors}，每block l1≤1。own-question conditional mean导数为bar g_j a.e.，dimension-free cube variance≤2E||bar g_j||²，使总energy≥βh。thresholdβ/2给own-input list≤ell=2/β；outside-list squaredmass≤β/2，conditional Jensen把retained mean mass送回同一个fullg，然后full compatibility hull给其pointwise mass≤Z，故EZ≥βh/2=2。
7. 固定一pair的otherseeds及全部weights，唯一remaining occurrence仍uniform。fixed listmaps在j,k只见自身endpoint，因此每pair ranks定义整个source的一次global labeling，与sampled occurrence无关；source value≤σ给pairhit≤ell²σ，pointwise Z≤1+Σpairhit推出EZ≤1+binom(h,2)ell²σ<3/2。该上界重新平均到同一restricted law，与EZ≥2冲突。先按m定ν,β,ell,h,d=2^38m^10,σ=2^(−46)m^(−12)，再请求alphabets，避免alphabet dependence；cover gap ratio2−6/(m+2)给任意fixedα<2的rational threshold。

Construction ideas and speculative π interfaces:

- **把两端共享依赖的二次总能量，经保持对端信息的whole-star重采样降到线性预算。** Premises: 独立pair seeds c_e与独立position weight blocks；f每seedtuple的weight mean为0、Ef²≥νd；单weight/seed oscillation≤2/4；同时重采样k所有incident seeds但保留每个opposite position收到的question时，f oscillation≤2。需h满足16h+32h²≤νd/2。 Mechanism: weight-support component指定j后，invisible residual on{j,k}归账到star centered k。D_k只改一个position的可见query，所以residual energy总和≤Σk||f−D_kf||²≤4d，替代逐pairO(d²)。随机h restriction的linear inclusion与pair inclusion概率分别h/d与h²/d²，保留own-question variance≥νh/4。 π interface: 若Q(Z)的列/条件能用独立pair parameters表示，且一组共享参数的重采样在保留其他列的可见数据时只改变一列，可尝试把quadratic occupancy的‘所有pair’代价分为星残差与真正pair项，获得只对不可见依赖计费的selection lemma。 Missing: 尚无实际Q(Z)的independent-seed representation或whole-star oscillation bound；rank/determinant非线性跳变不受本文Lipschitz variance控制。pairwise低占用不自动满足star条件，且variance存活不等于一个integer common basis存活。
- **局部平均不必共享一个全局证书；先用条件Jensen把质量送回同一full witness，再应用共同可行集合的凸包。** Premises: 有一个full vector g∈conv{±v_I:I为共同compatible selections}，每block l1≤1；local information G_j给bar g_j=E[g_j|G_j]，总EΣ||bar g_j||²≥βh；list T_j只由G_j决定，thresholdβ/2。 Mechanism: T_j大小≤2/β；off-list squared energy≤β/2。对list mass先用 conditional Jensen EΣ_T|g|≥Σ_T|bar g|≥Σ_T|bar g|²，再对同一fullg用signed-incidence hull约束ΣjΣ_Tj|g_j|≤最大一次compatible selection可命中的list数Z。由此EZ≥βh/2，而不声称拼起的conditional means本身属于原hull。 π interface: 对π的common bases问题，可尝试令v_I编码同一个Q(Z)上同时满足非零与height/support约束的可行基，构造full certificate vector后才做各列局部平均/threshold；把局部质量回送fullvector再提取一份共同basis，能直接针对‘rank证书与小量证书不是同一个对象’的缺口。 Missing: 需要真正的共同basis族、其凸包表示、至少βh的energy来源和一个已合法的fullcertificate；本文并不为任意两个matroids提供common basis，也不把formal rowrank转成该hull。若仅各自有bases再对local均值threshold，不能套用。整性与determinant height还须独立控制。
- **选择函数在进一步条件化时保持固定，禁止让certificate偷看被测试的共享seed。** Premises: fixed total local list maps只依赖其own endpoint questions与weights；所有其他seeds/weights固定后，remainingpair occurrence依旧uniform；每个rank pair定义对所有source variables的endpoint-only total labeling。 Mechanism: 后续conditioning仅改变map输入的分布，不重新计算conditional kernels或conditional gradients。于是rank-pairdecoder是一个sampled occurrence之前已固定的global assignment，source promise才可用于其成功概率；最后把uniform upper bound平均回原law，与原law的energy lower bound比较。 π interface: π项目若用按列/按support的局部选择器，应在最后评价/取minor前冻结选择规则和公共Q(Z)，并明确哪些参数被隐藏；证明非零与小量时须在同一剩余law或同一个整数评价上闭合，不能在不同条件化后重新优化各自的basis。 Missing: 这是证书依赖审计原则，没有给出π的decoder或probability gap；若选择器确需见完整shared parameters，本文bound失效。也不能从a.e. derivative identity推出每个特定integer input的rank或smallness。

Limits: 全文源码首读，不是独立认证NP-hardness theorem；未核验PCP/parallel repetition等外部source theorem适配，未运行reduction或形式化验证其analytic inequalities。 polytime只对每个固定m成立，d、alphabets、r与exponent可以极大；不声称统一实用算法或m作为input时的同一复杂度。 own conditional gradients一般不同时在full compatibility hull内；文章先Jensen回到fullg才作hull比较。后续fixweights/seed不得重算kernels，否则pair decoding不再对应fixedglobalassignment。 π桥接均明确是待证构造假设；没有新Q(Z) rowrank/common-basis、integer determinant非零小量矛盾，或可用columnheight/quadratic occupancy结论。

已逐行读完README、main、全部9个sections（含history）、star figure所有标签/caption与完整references；本篇所有13文本blob均与固定提交tree的Git blob SHA一致。输入闭包完成，无附录或额外body。主paper.pdf未逐页核对；reading不等于独立verification。

## Cross-paper observations

The strongest common theme is closure on one witness. Position 37 compares the identical normalized Jacobian under the identical weighted law, after explicitly transferring exceptional-target null sets to the root. Position 38 checks both capacity gain and momentum cost on one explicit composite map, with reciprocal form factors and finite constraints fixed before the power. Position 39 sends local conditional-gradient mass back to one full compatible gradient before invoking its hull; its local means are not falsely assumed to share that hull.

For the π BA/non-BA project, none of these papers supplies an actual common basis of Q(Z), a new row-rank proof, an integer nonzero-and-small determinant, or a columnheight/quadratic occupancy bound. The ideas above are conditional construction programs with their missing premises stated. Distant geometric and complexity mechanisms were retained rather than filtered away.
