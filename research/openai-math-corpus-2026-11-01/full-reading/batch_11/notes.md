# Batch 11 — positions 23 and 34

Pinned upstream: `openai/math@4eaee87548250b2682affda7d05d947c05554bae`.

两篇均按 actual main 的全部正文输入闭包逐段完整阅读，无相关性筛选。`full_body_read` 只表示阅读完成，不表示独立核验定理正确。下载与阅读严格分开。

| position | manuscript | files read | lines read | status |
|---|---|---:|---:|---|
| 23 | Two limit cycles for quintic Liénard systems | 11 | 3265 | full_body_read |
| 34 | Nonuniqueness with local conservation for the hard-sphere Boltzmann equation | 18 | 6197 | full_body_read |

已读包括 README、actual main、全部 body、preamble（34）、嵌入图和bibliography；无未读附录或正文。缓存目录为 `/workspace/scratch/66870cb2f937/source_context/round02_group5/{23,34}`。所有文本逐个核对原Git blob SHA、size和line count。`apply_patch`多加的terminal LF仅在删除一个LF后精确恢复pinned SHA时删除；无正文修改。完整审计见同目录 `byte_audit.json`。

35的全部文本也已精确缓存并audit，但未阅读其正文、未收入本batch的records；交由fresh-context reader。35的main/README曾用于映射探查，不构成body completion。

## 23: proof map and reusable mechanisms

主张：任意实系数primitive F次数≤5，经典 `ẋ=y−F(x), ẏ=−x` 至多两条几何不同的孤立周期轨道，且可达；涵盖退化孤立环与零区间的正确排除。

1. 两半轨道用 `u=x²/2` 化为 `u_y=φ±(u)−y`。半宽r是峰值的全局参数；轴端点的全部宽导数通过逆支方程延拓到非光滑profile边界。轨道与同一开宽区间上 `Δ=M+−M−` 的孤立根对应。
2. `H(λ,κ,r)` 的二次profile拟合是全局坐标。先固定κ拟合M，再由 `∂κ(H_r|H=M)=G/P>0` 和独立properness证明完整范围。这里局部Jacobian非零与全局可达明确分开。
3. 最重步骤是固定模型参数的 `(R/G)_r≤0`。W的特征输运、端点代数、Riccati的正规2×2系统、共同projective标签和Schwarzian给首零排除。标量消元方程在y=0奇异，但同一fundamental pair跨越它，因此两端对象一致。
4. 四个e,c符号情形只有一个允许两根；该情形 `Δ′=AΔ+Bq, B>0, q nondecreasing`，正积分因子后先降、平、升。其他情形给每个零处正导数，至多一根。
5. 下界使用显式small-ε quintic的精确能量商，在s=1,2有两个单根，从而两条双曲环。

对当前π项目最有价值的迁移是：Schur补式条件拟合与全局范围分开；同一对象跨奇异坐标的精确归一化；把复杂匹配差的零点计数变为单调障碍符号。它们都还没有给actual Q(Z)存活minor或整数高度预算。

## 34: proof map and reusable mechanisms

主张：三维periodic hard-sphere Boltzmann存在同一速度有界、质量/能量/absolute entropy有限初值的两条global entropy solutions；均strong L¹ continuous，满足全部五个局部collision-invariant balances和任意polynomial-weight的Q±可积性。

1. √2方向的nonzero integer norm给periodic tube separation。geometric annuli+fine opposite stripes的cold jets质量为O(s²)，密度为O(1/s)，velocity width为s¹⁰。bath接近Maxwellian但短轴周围挖hole；finite colored sum始终解同一cutoff equation。
2. 单cold color由free core+small Gaussian barrier控制，active tubes在 `L_s=s^(−1/2)` 前不交；到期转入regular part。cap-selection以连续seed参数选达到预设horizon-dependent bound的解，零seed给dormant family。
3. 新global continuation使用small entropy的spatial L¹ defect、free-flight velocity averaging、两个分离速度球、small-mass gain，把正frequency与最终nonlinear gain常数固定在bootstrap B之前。指数衰减使退休增量可求和。参数顺序是实质步骤。
4. 同位置二次collision products由平方包络控制：disjoint support消除cross terms，每尺度 `∫s²/(t+s)² dt=O(s)`。velocity averaging+Vitali识别products，强path traces与weighted flux控制传全部local laws。
5. hit dilation给atomic two-stream medium中的linear limit。cap达到sup norm不自动保留nonzero；two gains排除axis concentration，five gains的analytic rank7与robust pullback排除nonaxial concentration。
6. constant-medium branching给positive potential、Doeblin minorization、mixing和drift block；finite compact history tests把它稳定到measurable medium。每次mass增长≥3，cap仅O(t⁴)，选择 `3>(1+a)^4` 排除nonzero global dilated limit。
7. finite U让hit scale与horizon scale的ratio有正下界。dormant horizon dilation为零而cap hit dilation非零，两尺度Jacobian精确比较迫使两条global solutions不同。

当前项目的具体跨领域候选：

- **Ambient rank与admissible selection分开。** §10在所有collision times相等的ambient点展示rank7，再用analytic minor零集null；witness无需满足spatial masks。bounded masks留在密度中，C¹ pullback仍成立。用于Q(Z)时必须另证实际row-transversal支持不完全落入坏簇；a.e.rank没有整数高度或最小minor下界。
- **二次预算先于弱极限。** 低mass/低线性cost不足以控制同位置products；真实disjointness或明确overlap预算才令平方成本可加。这与quadratic occupancy设计直接相通。
- **有限寿命替代永久分离。** 只在所需height/time window内要求每block横截；过后以衰减转入remainder。√2的norm是已知算术输入，不能替代π所缺的integer determinant非零性。
- **同一量的可重复增长对全局cap。** 有严格余量的finite local blocks可在粗糙环境中反复放大，最后比较同一个M(t)的两个指数。迁移到π仍须把同一个非零integer determinant及清分母、column-cost、occupancy始终绑定。

## Scope limits

未编译PDF、未作逐页视觉核对、未逐条独立复算长符号证明或回读引用。34中嵌套figure input按recursive tree的唯一suffix匹配到 `build/paper/figures/*.tex`；已完整读取，未测试TeX search path。两篇的跨领域想法都标明假设、机制、missing step和源码锚点，均不构成π BA/nonBA结论。
