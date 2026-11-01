# Batch 18 — queue 35完整正文首读

固定来源：`openai/math@4eaee87548250b2682affda7d05d947c05554bae`，目录tree `963c61e950f08f4f55576651f7995711596635b8`。论文为 **The Selmer converse for elliptic curves at every prime**，2026-10-07版本。actual main是`build/main.tex`。这是全文源码阅读记录，不是独立证明核验。

已完整读取README、main、preamble、全部11个section/appendix输入、嵌入proof-roadmap图和完整bibliography。共16个text blobs、334,261 bytes、7,160行；每个文件字节的Git blob SHA与size都对pinned recursive tree重算通过。所有正文分块顺序读取、输出未截断。`records.json`逐文件列出完整read ranges和实际chunk ranges；`byte_audit.json`保留核对结果。唯一未纳入源码阅读的blob为`main.pdf`；没有PDF阅读或渲染核验主张。

| 文件 | 全读行数 |
|---|---:|
| README / main / preamble | 21 / 33 / 82 |
| introduction / roadmap图 | 269 / 36 |
| arithmetic reduction | 257 |
| tame deformation | 868 |
| theta construction | 687 |
| constant terms | 722 |
| positive Fourier–Jacobi coefficient | 553 |
| cusp lifting | 634 |
| characteristic-zero Galois parameters | 717 |
| extension extraction | 609 |
| Sylvester cube sums | 118 |
| nonstandard-transfer appendix | 484 |
| references.bib | 1,070 |

## 稿件的精确主张与依赖边界

对任意椭圆曲线A₀/ℚ、任意素数p，full p-power Selmer group若corank为r∈{0,1}，则analytic rank与Mordell–Weil rank同为r，并且整个Sha有限。没有reduction、CM、残余表示、torsion或isogeny限制；不计算refined BSD leading coefficient。应用是在additive prime 3处处理所有ℓ≡4,7,8 mod 9的prime cube sums。

p=2的converse与用于选辅助twists的analytic density theorem均被明确引用自同collection的另一稿件。本次未回读该外部稿件，亦未验证Gross–Zagier–Kolyvagin、theta/period formulas、PEL geometry、Labesse、RACSD Galois与BHTY nonvanishing等引用证明。读全bibliography仅表示读了文献记录、版本和locator。以下架构均描述稿件如何论证。

## 证明结构

1. 奇p下用analytic twist density在指定positive-density congruence classes选D,D′，构造K/F/L。两条K-Selmer line中，辅助twist line由已知non-torsion point生成并在p非零局部化；原line若也nonzero localization，one-sided Greenberg base space为零。唯一strict情形有两维Greenberg base〈a,h〉而Selmer line a处处locally zero。只需证conductor-one Heegner trace y非torsion。
2. 选moving split qᵢ及ring-class cyclic p-quotients，固定处character平凡、p-local degree的p-part有界；Frobenius趋γ，ε(γ)=1、V(γ)无eigenvalue1、C(γ)非root of unity。定义固定Artin ring R_b和bounded admissible finite-stage cochains；identity必须在**所有元素序列**上成立，才可在uniform precision下specialize和除t。
3. strict时first-order lift有cup-product obstruction。修正cup primitive的Heisenberg commutator只改变central primitive value，保留V、所有characters、radicals及Euler bounds；可选γ使obstruction为非退化skew 2×2矩阵。exact coefficient sequence把任意高阶class挤到last t-layer，length≤2；non-strict为0。故b*=3或1总比upper length大。
4. y torsion杀Heegner weighted Kummer trace的constant coefficient；strict时conjugation反转t，而Selmer line invariant，杀第一次t-division的specialization；第二次division落在locally zero line。uniform local logarithm比较得到logs在t^{b*}前趋零。
5. 构造两个real factors均U(3,1)的holomorphic theta family。Fock权(m,m,1;−1)及literal monomial-dual/divided-power lattice始终一致，split pivot tests给两个normalized parabolic units。unary moments经common Igusa square、compatible cusp norms和least-denominator saturation给uniform bounds。
6. binary constant-term scalar S满足清分母bilinear period identity A S²=B P₀²P₁²。P₀在m→0时经inverse Katz primitive/path integration变成Heegner log，包含任意bad reduction。p处gamma factors按相同orientation抵消；q处保留完整raw factor而不除growing group order。先在group algebra清分母，再固定阶Artin invert；从S²的2b*阶小与coefficients bounded逐偶次归纳得S的b*阶小。
7. 独立构造positive Fourier–Jacobi coefficient，其augmentation valuation有统一上界。两个unary periods分别验证root sign，外部nonvanishing选固定finite characters，在compatible degree discs同时避period和moving Euler factor零点；固定处Cramer-rule spanning以bounded denominators施加额外local controls。
8. target PEL ordinary几何不要求A₀ ordinary。actual J-level的canonical branch只标记multiplicative flag/frame，boundary为finite étale covers。positive Jacobi base change、Hasse延拓和Serre vanishing lift到genuine cusp forms。Serre–Tate trace在同一integral lattice给normalized units，bounded augmentation content允许取消所有Hecke relations，得到**full Artin ring** congruence。
9. 一个Labesse parameter同时保持full infinity character和所有split local data；逐real-factor证明individual holomorphic pseudocoefficient sign。Arthur blocks用实际shifted-copy direct sum构造4维Galois representation。两个Hodge–Newton contacts给local flags；resultants、direct-sum monodromy分析和Steinberg parahoric rank bound给actual inertia identities。bounded allocation filtration使operator nilpotents只损固定precision。
10. uniform PI-algebra module generators保证limiting corner functional仍是admissible cochain。2+1+1 corner algebra先消去cyclotomic extension path；common-kernel条件改成corner minor identities，在quotient及nilpotents下存活。D_fy若在任何ambient field消失，则stable y-line的actual conjugate weight sum3和residual ε⁻²的sum4在固定local units上矛盾。故Fitt₀(D_fy)=0，base change后dual length≥b*，并注入同一个Greenberg module，得到0<1或2<3的矛盾。
11. Sylvester应用保留paired-isogeny descents中的Sha kernels a,b，得到full corank≤(1−a−b)+(a+b)=1；随后negative root sign排除rank0。

appendix另外处理isogenous tori的组合normalized centralizer factor：保留function-field degree lattice index I(T)，只有split degree-one place才令I=1；不能把rational torus correspondence误认成integral isomorphism。其change-of-characteristic段在同一definable参数空间、compatible Lie differential frames和共同discriminant normalization下把paired orbital-integral identity化为constructible function零性。正文Labesse simple-trace路线明确不使用weighted fundamental lemma；appendix不验证其外部positive-characteristic输入。

## 对当前π / Q(Z)项目的推测性构造接口

完整的七条idea、precise hypotheses、anchors和missing steps见`records.json`，没有按题目距离排除想法。

- **Artin length而非仅specialization rank。** 若能给Q(Z)候选row transversals建立一个exact fixed-order deformation，first-order obstruction可能记录rank degeneration的持久程度。缺口是实际module、同一整数determinant及相应height；generic rank或t-leading coefficient还不够。
- **先清分母再截断，augmentation统一取消。** fixed b下bounded coefficients和nonzero bounded-valuation常数项给bounded inverse，取消成本至多bC。可用来设计清分母与column-height账本，但p-adic precision不会自动成为π所需Archimedean smallness。
- **literal dual monomial lattice。** polynomial/derivative构行时可检查factorials是否仅由basis convention引入；须同时核算basis-change determinant和integer lattice index，不能把成本搬到另一处后当作节省。
- **将局部subspace证书转成minor identities。** 在projection、quotient或column compression前保留具体pivot/Plücker证书，可能更适合row-transversal审计。缺口仍是允许行的实际surviving minor、primitive index和cost。
- **Fitting/support到总length。** 多个短torsion summands可共同贡献所需length，不必先找free rank-one summand。对occupancy/matching obstruction是远接口；Fitting零本身不给指定integer minor或size bound。
- **约束kernel中的commutator。** 寻找保持column occupancy与低阶moments的noncommuting operations，仅改变higher obstruction。必须实际算integral operations及height；类比不保证能修复determinant。
- **rational equivalence下的组合不变量。** row-rank保持不等于integer determinant或covolume保持；记录index、saturation、normalization的最终组合，避免算术下界和解析上界落在不同对象上。

本批不推出新的π determinant构造、BA/nonBA判断或算术矛盾。全文首读已完成；独立核验和所有上述迁移的关键构造尚未完成。
