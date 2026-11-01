# Batch 12 — positions 24–25, complete source-first body reading

固定上游：`openai/math@4eaee87548250b2682affda7d05d947c05554bae`。两个条目均未作主题相关性排除，均读到正文结尾及全部递归输入。此处 `full_body_read` 表示完整源码首读，**不表示独立正确性认证**。依 root 后续调度，position 27 留给 fresh-context 阅读者；本批不为未开始的 27 写 partial 占位。

| Position | 题目 | 状态 | 文字文件 | 源码行 | 源码字节 |
|---|---|---|---:|---:|---:|
| 24 | A nine-dimensional counterexample to Borsuk's covering assertion | full_body_read | 14 | 1,533 | 64,332 |
| 25 | Lech's multiplicity conjecture | full_body_read | 9 | 3,238 | 139,130 |

所有缓存位于 `/workspace/scratch/66870cb2f937/source_context/round02_group6/{24,25}/`。每个文件均用 Git `blob <byte_count>\0<bytes>` SHA-1 与预检树条目核对；apply_patch 产生的额外末尾 LF 仅在移除后精确重现原 SHA 时删除。字节数与 tree size 一致。`records.json` 包含逐文件路径、SHA、行数和实际阅读范围；`byte_audit.json` 汇总。两个 root 的递归 input/include/bibliography 图均闭合，全部 tree 中非 PDF 文件都已读，无附录或未读 body。

## 24 — 九维 Borsuk covering 反例

稿件主张完整紧 projector image `X={uuᵀ: ||u||=1}⊂Sym₄,tr=1≅R⁹` 的直径为 √2，不能由十个严格小直径集合覆盖；不确定最小反例维数或 `b(X)` 精确值。几何恒等式 `||Pₓ−Pᵧ||²=2−2⟨u,v⟩²` 把直径对变为正交线。严格间隙可把任意覆盖扩大成仍避正交对的开覆盖，smooth partition of unity 因而给非负、和为一、正交 supports 不交的 admissible map。

核心从 map 而非 triangulation 推导 support complex 结构。PSD 积分延拓 `E(Q)` 精确保 trace 和 range 上所有标签的并；`F(A)=E(A₊)−E(A₋)` 因正交谱空间的标签不交而没有抵消，`||F(A)||₁=tr|A|`。odd sphere map 给标签数下界 `m≥k(k+1)/2`。等号时，稿件用规则值 level set、affine-chart 平凡化的 orthogonal-complement sphere bundle 和 projective double-cover class 比较奇偶 preimage counts，逼出所有 maximal supports 恰有 k 标签，facets 之和成为 F₂ cycle。

六标签系统为互补 triples 恰选一、每 pair 恰含于两 triangles；点 link 为五圈，五标签 restriction 唯一决定全系统，transposition 不能保系统。十标签 complex 的 facet complement 强迫该系统，exchange 闭路排除三重 triangle 延拓并把 edge links 限制为三圈/四圈。partner 交换保留 blocks，剩余标签计数只容许 block sizes `3,2,2`；一个外部 pair 的 triangle 遂仅能有一个 tetrahedral 延拓，与要求恰两个矛盾。末节以 `I_s+(1/4)11ᵀ` 的 Gram 矩阵添加 diameter-distance 点，推广到每个 d≥9。

可迁移构造记录三条（详见 JSON 的 precise hypotheses、missing steps 和精确 anchors）：

- `extension.tex:25,112`：通过正锥输出 support 分离得到精确无抵消 norm identity；探索 Gram/PSD lifting 是否可把 pi 的 quadratic column occupancy 与 trace quantity 联系起来。缺的是实际 admissible partition 以及与同一整数 determinant 的定量联系。
- `topology.tex:88,262,290`：阈值等号处由 odd degree 和 cover class 推出所有局部 maximal supports 的纯度及 cycle coefficients。Q(Z) 候选若无真实 manifold/odd-map 结构，不能只凭 supports 仿用结论。
- `six-labels.tex:44; ten-labels.tex:41,73,110`：局部唯一恢复 + exchange holonomy 排除 multiplicity，并产生不变 blocks。可尝试给 row-transversals 定义 transport/forbidden permutation 后导出 occupancy 限制；普通 basis exchange 不具备这些刚性，仍需保真实形式秩。

关键界限：完整源码首读未验证 topology equality argument，也未把外部 triangulation vertex bounds 误作此 support complex 的直接定理。严格小直径是平滑化约化的真实假设。

## 25 — Lech multiplicity

稿件主张任意 nonzero Noetherian local rings 的 flat local map 满足 `e(R)≤e(S)`，不限维数、residue field、characteristic。另称 characteristic-p complete local domain 上任意 short finite-free complex 的 Dutta multiplicity 至少 `e(R)`。

统一核心：可能非 Noetherian 的 `C`，`I=(z₁,…,z_h)`，additive length `λ`；不仅要 additive，还要求固定 I-power annihilation 下 zero-length sums/colimits 稳定、每 parameter-power quotient length 精确等于 `a^h μ`，positive Koszul homology length zero。对长度≤h、off V(I) locally contractible、alternating rank zero、entries∈I^s 的 free complex，有 `λH₀≥μs^h(b₀−A_hΣbᵢ/s)`。常数只依 h 和 ranks，不依 ring/length/finite extension。

blowup 让 I invertible；graded comparison modulo length zero 后，不能直接用有限 thickening 的 Euler（alternating-rank-zero 时它们全为零），而用大截断到固定小截断的稳定 homology images 来比较原 complex 与 blowup pullback。将 differential 除适当 I-powers 后，Eisenbud–Schreyer product-of-projective-lines 构造的 root bundle 控制 cohomological degrees。选择 `∏(x+i)` 导数根附近的整数 roots，使每个 layer sum 主项相同；alternating rank identity 恰隔离 b₀。

char p 的 perfection normalized length 通过 trace-dual conductor 压低 quotient kernel 为 `O(q⁻¹)`，识别为 Hilbert–Kunz，并用 Roberts estimate 控制 Koszul homology。Cohen factorization 给 `T=R⊗_A U`、`e(T)=e(R)`、perfect kernel J 与 `F=P⊗K(x)`。Frobenius 不改变 ranks 而加深 entries；componentwise lower bound、associativity 及 Ma perfect-ideal identity 合成 char p 结论。

mixed char 借 parameter roots 和 Heitmann integral extension 加深 entries。扩张 rank r 可随 q 增长，但 `q^h μ=r e(B)` 精确抵消。固定 q 的数据先取 parameter-power N-limit，prime filtration 给不含扩张的 upper bound，之后才 q→∞。generic regularity 使每个 top filtration prime 仅属于一个 component，故没有 coefficient 损失。residue char0 则 finite-free/common-coefficient-field reduction 后，把 algebra/locality/freeness 和单向 multiplicity bounds 编为有限整系数方程；source presentation minors 消失 + reduction relation 给下界，target monic initial relations 给上界，Artin/étale specialization 转到 char p。

可迁移构造记录四条：

- `complex.tex:468,548,603`：导数根把一系列 layer integrals 的主项精确平衡；对 column-cost/二次 occupancy 可寻找同类 antiderivative 权重，但必须先有真实层模型与正性。
- `complex.tex:314,361,408`：稳定 images 恢复原对象，避免有限截断数值掩盖真实存活信息。Q(Z) 若用 jets/specializations，需要 compatible maps、有限 filtration 和最终稳定证书。
- `mixed.tex:255,279,404,414`：精确归一化扩张成本、误差常数对扩张独立、先内层后外层取极限以及 coefficient-one incidence。pi 的 denominator/height/列数必须独立重新计算，不能借用 multiplicity rank formula。
- `charzero.tex:82,159,276,324`：有限 polynomial certificate 仅保最终矛盾需要的单向 bounds，monic leading terms 在特化时保持。对目标 minor 非零需额外 inverse witness，且存在性/Artin approximation 不保证低高度整数点或 analytic smallness。

界限：CLMST、Heitmann、Cohen–Gabber、Roberts/Hochster–Huneke、Ma、Artin 等按稿件给定版本和假设记录，未逐篇独立回读验证。文中 uniformity 不包括所有 homotopy exponents/truncation levels；mixed case 也没有交换 N 与 q 极限。没有得到 pi 的最终同一整数 determinant 非零且绝对值小于一的证书。
