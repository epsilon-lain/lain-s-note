# 第 01 批：正文首读记录（队列 1–3）

固定提交：`4eaee87548250b2682affda7d05d947c05554bae`。本批包含 001（Milne rationality）、032（CM Hodge）、071（removable boundaries）。每篇均先读 README 与 main，再读取 main 的全部正文输入及图、宏、文献表；recursive trees 均 `truncated=false`。没有按与数论题目的表面相关性排除稿件。

三篇的正文及实际输入闭包均已读完，分别记录 10／14／12 个文本文件，3086／4672／2066 行，共 9824 行。行数包含 README、preamble、图源码及文献表；不是 PDF 页数。没有附录或未解析的 TeX include。记录在 records.json 中逐文件列原路径、blob SHA、总行数及完整区间。PDF 二进制未逐页比对；引用文献和 Lean 未逐一验证。因此“首读完成”不表示论文定理已独立证明。

## 001：精确有理归一化与最大障碍

作者主张：任意好约化 abelian variety 的有理 Hodge 类与特殊纤维上的互补除子乘积，在各个 prime-to-p 以及 crystalline 实现中的迹是同一有理数，含 p=2。结合 032 才进一步得到一个共同有理代数 cycle；这个加强并非主有理性证明的输入。

正文的关键链条是：共同有理 Lefschetz 空间的完美配对 → 在幂的特殊纤维上经一个 quasi-isogeny 输运 → 用有理交数比校准所有局部 similitude → CM eigenspace 平衡约化 → 平衡 Weil 类的“正交为零／存在可逆自伴反线性 endomorphism”二分 → 在后一分支构造保持额外除子的 admissible filtration → 另一个辅助 lift → 除子多项式的系数提取 → 同一有理恒等式在所有实现共同特殊化。

可迁移的具体工具：

- **有界次数同时避零。** 对固定域 K1 上总次数≤b 的所有非零多项式，取 N>b，B=bΣN^i，代入 Xi=θ^(N^i)，并令 θ 的域次数大于 B。原文在离散赋值域用 θ^q=p、q>max(B,e) 实现。不同 monomials 不碰撞，故一个点同时保所有 bounded-degree minors 非零。可用于先做形式 rank/nonvanishing，再选参数点；但代数次数和幂指数很大，尚无低高度整数 specialization。
- **总缺陷的超模最大值。** f(P)=dim(P∩L)−tN(P) 有有限值，且 f(P)+f(Q)≤f(P∩Q)+f(P+Q)。最大值对象在交与和下封闭，故唯一最大障碍被所有保持 f 的 automorphisms 固定。数论模型若能把容量缺陷写为 modular cost−matroid rank，并证明支持家族是相应的格，可只检验稳定障碍。加入 quadratic occupancy 或 transversal 限制后需要重新证明闭包与超模性。
- **同一对象的归一化。** λv=deg(a*θB·θ*^(D−1))/deg(θ*^D) 是同一个有理 λ。这是当前“同一个非零 determinant 的算术下界与解析上界”路线的精确结构提醒：任何基变换或投影都要有一个共同有理清分母因子，不能在不同局部缩放对象上分别估计。
- **系数提取而非逐 embedding 下降。** P(x)=(Σxτ²Dτ)^s 的纯系数给 determinant lines；整个有理空间通过 scalar extension 下降。若用有限差分/Vandermonde 迁移到整数 determinants，评价点高度与逆矩阵分母必须显式计价。

独立有效性边界：积分模型、special-point、admissible period、uniformization、p=2 条款及同态提升都按源码引用读取，未回读原文逐项认证。此类输入不应被本批阅读状态隐藏。

## 032：局部交换、保非零拼接与固定总输入

作者主张：所有复 CM abelian varieties、任意维数和余维的 rational Hodge conjecture。正文的主证明不使用同 001 的交叉特殊化后果。

主要构造链条：所有 Galois 共轭平衡的 embedding-line tensors → 同 column sums 的 sign-row 表之间的 two-row switches → 四因子 switch → surface 上非零四重 period → 四个 divisor kernels 把 functional 变为精确 pure tensor。analytic 部分在 compact two-ball quotient 上构造 theta one-forms，以 rank-two double alternation 得 wedge；保持总 rational space、finite input、quotient 固定，只在 compact-place decompositions 中移动，用非零 self-pairing 的连续性得到 mixed period。arithmetic 部分以 full tuple moduli 的 ordinary-neighbor 计数得到 integral Frobenius cubic，再由 inverse Satake parameter、central finite order、degree-one Weil weight 与 filtration-preserving generalized Hecke projectors 抽出所需 rational CM sources。最后回同一 projective surface，再逐步 contraction 回原 variety，下降有理系数。

与 signature capacities 的直接结构联系：

- **固定 column marginals 的 fiber 连通性。** 两张 ±1 表若每列 plus 个数相同，可逐列配对两种 mismatch，仅交换两行的相反 entries；每一步满足四行关系 a+b=a′+b′。这不要求整个高阶关系一次生成。对 Q(Z) 可尝试定义同 signature 的 assignment fiber；但禁止边、真实有限行池与 rank 条件可能断开 fiber，必须逐步检查 surviving transversal。
- **辅助 slots 只能通过已证非零的 dual pairing 消去。** 原文端点张量拼接会乘明确的 sign 与 Πqi，qi 是由 ample polarization 保证非零的 pairing。一般 determinant 的 Cauchy–Binet 展开并无同样纯度，不能省略 cancellation 问题。若迁移，要构造 rational dual bases，并计算 pairing 分母。
- **先固定总输入，再改变分解。** 非零 mixed period 来自同一总 theta tensor 的连续变动；只有选择 rational endpoint 后才改坐标并有限分解。对参数 determinant 的启发是先固定 D(t) 和解析界，找非零 seed，量化保非零邻域，再选低高度 rational t。源码这里只给存在，不控制分母，故不能直接进入当前指数竞争。
- **projector 必须同时保两种结构。** 纯 slope 0／1 只有在 projector 保 filtration 且与 Frobenius 交换时，才能用互补分块的 weak admissibility 得纯 Hodge type。文中给出不保 filtration 的 slope-projector 反例。矩阵问题里对应要问：解析有利的投影是否也保整数格／可控分母？

已完整读到 theta covariance 的正负与对偶、metaplectic/circle splitting 区别、finite stabilizers、inverse Hecke convention、generic-neighbor 的 q/q² 重数、q-power Frobenius 方向、generalized eigenspaces 的 nilpotent 部分和所有 arithmetic consequences。这里仍是作者论证的完整首读，不是逐节点 audit 通过。

## 071：保留有限证书、固定能量下界与平方和相消

本文精确主张是 removable boundary ⇒ circle-domain rigidity，任意补集分量基数；不是本文自身的 general existence theorem，也不是双向 He–Schramm 等价。引言明确区分这些断言。

主链条：任意连通度的 component correspondence／collapse → 引用伴随 Koebe 稿件的 retained-test construction → 补洞 jump error 和 contours 确认弱极限的 exact values → 近极大 Laurent coefficient 与 conductivity deformation → pair-null level-set compatibility → 给 finite lists 加 separating tests 并保存 small-energy barriers → 极限若有 nonpoint component，其轴投影含固定区间 J，强制 energy≥|J|/(2L)，与保存的 k^−2 上界矛盾 → dust scalar Sobolev extension → masks 处理 disk/point mismatch → common-point traces 在典型轴线投影长度零 → slice flux 下界与两侧 area 上界相消，得到非负平方和≤0，逼 rigidity。

可迁移的四点：

- **旧证书的 strict margins。** 第 k 阶段保存 tolerance tk；以后误差按 2^−n min(t1,…,t_(n−1)) 递减。保非零/保边界值只需有限阶段精确，未偷用 infinite-list model。构造越来越大的 determinant family 时也可保存已有 nonzero minors 的 quantitative slack；但 rational rounding 与分母成本另算。
- **稳定障碍给固定下界。** 一个 nonpoint limiting component 的投影 J 不随后续阶段消失；每条相应 slice 在长度≤2L 内必须从1变0，由 Cauchy–Schwarz付固定能量。数论迁移需要找到不随参数规模塌缩的 rank/occupancy obstruction，并把它转为同一对象上的定量预算下界。
- **三条相容不等式产生平方和。** 同一 disk、同一 f 的 source area、target area、flux 合起来精确给 ∫(|f′|−1)²+πΣ(rj−sj)²≤0。可以尝试把 column-cost 与 quadratic occupancy 改写成两个 quadratic budgets 与一个 bilinear lower bound。但只有从实际 determinant 构造证明三式才有内容，形式配方不补 nonvanishing。
- **重叠错误可由 l2 core packing 综合。** balls 内若有相对面积有正下界的互不相交 cores，则 maximal theorem 给 ||Σ(aj/dj)1Bj||2≲||(aj)||l2，无需 balls disjoint。离散版本需要真实的 averaging/maximal-operator bound 与 disjoint cores，不能直接把 l1损失改写成l2。

关键边界：071 的 retained-test transfer/selection 是伴随《Koebe’s Circle-Domain Conjecture》Theorem 7.5／Corollary 7.6 的外部输入。本批没有声称那篇也已读完。它不是未解析的 TeX include，所以不妨碍本稿正文闭包完成，但妨碍把本稿阅读等同于独立证明验证。

## 对当前 pi 路线的可检验后续

这三篇没有给 pi 的 BA 或 nonBA 结论。先保留以下具体研究接口，而不按学科否定其价值：

1. 为同一 signature 下的 surviving Q(Z) assignments 定义允许的 two-row moves，证明 fiber 是否连通；不连通时记录真正 invariant。
2. 把“容量−rank＋occupancy”写成明确集合函数，检查交并闭包及超模方向；只在成立时使用 canonical 最大坏集合。
3. 把目标 integer determinant 编码为 bounded-degree polynomial，先独立证明形式非零；再单独优化低高度 specialization，而不把 Zariski genericity 当作算术成本已解决。
4. 每次加入辅助 slots、投影或基变换，维护共同清分母、非零 scalar、解析范数和 column/occupancy 成本的同一张账。目标仍是同一个非零整数 determinant：下界至少1，上界严格小于1。


## 源码缓存精确字节审计

对全部 36 个缓存文本重新计算 Git blob SHA。初始每个文件均仅多出一个末尾 LF；且只有删去这一个字节后才精确匹配记录中的固定 blob SHA，故逐个按此条件修复。修复后 36/36 SHA 全匹配。原记录的 total_lines 与完整 read_ranges 已按原始内容计算，分别 3086／4672／2066 行，无需改数。正文内容未改变，不重复计算首读覆盖；审计明细见 byte_audit.json。
