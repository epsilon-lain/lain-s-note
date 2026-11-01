# Batch 07：完整源码首读记录

固定上游：openai/math@4eaee87548250b2682affda7d05d947c05554bae。本文目录 tree bcd86ed1f932263ae0387e15469e1966c7c2e94c，recursive preflight 未截断。只完成队列位置19；位置30、32仅做文件目录检查，没有正文已读量，交回root另行分配。本批records目前只登记19，避免把未读paper作为已读计数。

## 19. Algebraic Kuga–Satake correspondences and Hodge conjectures on a K3 quadratic locus

完整读过 README（18行）及自包含 build/paper.tex（5391行）。所有源码bytes均以Git blob SHA校验：README 080907511932a1e3a66cd137613669ce64081126，main 20ec0749de32f6e97cdac6e34fcd0cf2d49bfe6e。唯一TeX中没有input/include外部正文，没有外部bibliography命令；参考文献全部内嵌，无附录。README链接paper.pdf，main title与其对应；未编译、未逐页核对PDF。这是完整body首读，不是独立证明认证。

### 主张与结构

论文主张：若projective complex K3的实际transcendental cup-product quadratic space能rationally isometrically embed into U_Q^⊕2⊥〈−1〉⁴，则所有指定full even-Clifford Kuga–Satake tensors都代数，并有所有self-powers、所有次数的rational HC/GHC。此locus包含所有ample P=U⊕D8(−1)⊕D4(−1) polarized K3，含Picard jumps。还给一般criterion：任意projective K3只要一个exact standard even-Clifford tensor代数，就有全幂HC/GHC。后续把结论转到Hilbert schemes、指定twisted sheaf/Bridgeland moduli和相关birational hyperkähler varieties。

一般criterion先行而且不依赖后半的具体K3 family。关键是保存exact j_w(v)(a)=vaw：alternating compositions在1处评价为c_k(α)w^k，PBW与w非isotropic给所有exterior powers的injectivity。取algebraic adjoint、用Hodge polarization的nondegenerate restriction，再以Cayley–Hamilton polynomial inverse建return R_k。这既使exterior Hodge classes由divisor返回成cycles，又以support dimension给whole exterior coniveau≥k−1。TR case的relative volume tensors和CM case的standard–dual pairings经Zarhin group+FFT组装成普通HC。GHC另需最高权carrier argument：按partition columns放入full exterior carriers，TR/CM多embeddings的宽度取最大值，Galois rationality控制每一factor；ordinary HC algebraizes真实rational quotient，再用return的codimension支持整个substructure。

几何半部从abelian-surface double cover X、resolved inversion quotient Y建rank16 anti-invariant lattice V；剥去八个algebraic square−4 directions得到rank8 U。其Y-form经local quadratic classification为V_P，而resolved pullback U→U_X使form乘2。logarithmic tangent deformation和theta restriction把period differential化成4×4 tangency map。一个rank3 pencil direction A′的kernel是skew gradient；residue pairing=4证明A_d(kernel)不落在im A′，故det(A_d+λA′)的λ³系数非零，period rank6。

曲线源通过Gauss fibers的4组opposite points建coherent half sums。assignment torsor在每个puncture有fixed local selection，translation monodromy因此trivial；两个parity covers给curves。signed relations的排除很具体：3-term relations组成zero-diagonal、off-diagonal±1的4×4 matrix，模2是J−I而平方为I，故odd determinant；4-term relations则由local singleton monodromy与connectedness排除。addition map总degree8，四条constructed branches各multiplicity≥2耗尽valuation degree，从而pullback branch divisor为2E。剩余étale double-cover character在product fundamental group上分成两个factor，分别cover后得到product of curves dominates X。

为了把U的given Hodge comparison真做成algebraic map，保留大integral V以管控−2 roots，嵌入K3 lattice并通过period-stack canonical contraction下沉到原t-line。文中对finite Weyl inertia、polarizing line scalar和Hodge-line character均有单独处理。t=s+u²后Hodge line是O(1)，canonical degree为−1，得到K3 total threefold uniruled；X side也是conic threefold而uniruled。

两次square substitutions是最需要精确tracking的桥：第一步H³可通过curves algebraically split，故given H³ Hodge comparison由cycles实现；固定models、good G、relative cycle C后，再令s=b+v²。conic在t=b的两components差形成cylinder d，normal/self/cross intersections给t(d)d=−4id。通过Thom naturality与affine Leray的H¹(R³)层，保证C的degree3 agreement正好足以比较cylinder。扩大open时的boundary Hodge structures在b变化前已固定；simple nonconstant U_b对固定H_F的Hom除可数点为0，消去error，最后−1/4 transpose cylinder返回exact surface map。

curve injection仍不是所需tensor。full SO(2,6) seed处，weight-one effectivity把common Hodge Lie algebra的非trivial constituents限制为两half-spins；full C+中各出现8份。opposite-spin products各含standard U一次，same-spin不含U；两orientation加rank-one multiplicity maps span全部128维Hom。用rational Hodge-Hom base change下沉，不要求complex half-spin projectors单独rational；flip在H¹⊗H¹上须用负的geometric interchange graph。

every-period步骤在target local ample-P ball中固定whole rank8 U_P和full tensor action，不能只看holomorphic line或随jump变动的T_S。proper Hilbert cycle parameters上ambient Betti classes连singular/nonreduced parameters都locally constant；ordered generic seed逃开任何proposed countable proper images迫使一个proper image覆盖whole ball。Picard jump以后才做exact双槽compression：p⊗((E_T^flat)^−1 i^dual E_U^flat)，得到tensor of pji。最终用rational orthogonal density、Witt extension与Buskin达到全部rational embedding locus。

### 可以尝试的构造接口

records.json中完整列出六项假设、机制、潜在桥接与缺口，最值得跟进的是：

1. Corank-one pencil certificate：若Q(Z)的某个真实candidate minor可写成A+λB、rank B=r−1，则用left-kernel/right-kernel pairing ℓ(Ak)而非全排列展开来证明formal determinant非零。必须另证实际低height实现与最终整数评价成本。
2. Mod-prime certificate：尝试选actual integer/polynomial minor使其mod2（或另prime）图样可逆；论文4×4 signed pattern是一个明确例子。仍须同一列集合同时满足另一rank条件，不能把single-matroid证书误作common basis。
3. Exact operator exterior evaluation：寻找Q(Z)的faithful algebraic realization，令wedge determinant有可追踪evaluation identity。abstract nonzero embedding不够，polynomial inverse/retraction也可能引入巨大分母。
4. Shared carriers：按同一column index装载不同signature blocks，某些carrier count可取maximum widths。但quadratic occupancy下合并也会增加平方成本，需要在真实supports与rank constraints中重算。
5. Freeze before parameter choice：先冻结minor、transformations、denominators和forbidden polynomials，再同时分析nonzero与smallness。不能将Baire或uncountability用于获得小整数点；可数proper sets完全可能覆盖所有integers。
6. Ambient exactness then compression：先在固定ambient object追踪full map，再用explicit inclusion/retraction和dual normalization压缩。integer determinant目标还必须算p、dual map的分母、height和保共同基。

上述均是跨领域启发，未给pi项目新增已证明lemma、common basis、actual Q(Z)满行秩证书，或small nonzero integer determinant。主结果范围、Y/X的form normalization、exact Clifford map与普通embedding的区别，在后续使用时不可删去。

### 阅读覆盖

无截断阅读块为1–650、651–1300、1301–2000、2001–2700、2701–3400、3401–4100、4101–4750、4751–5100、5101–5391。最后641行首次包装输出有截断，已按后两块重读，故无阅读缺口。没有运行Lean或逐篇核验外部参考文献。

