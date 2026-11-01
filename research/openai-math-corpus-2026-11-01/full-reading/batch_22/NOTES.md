# 41 — A counterexample to Hadwiger's conjecture：精确续读点

固定上游提交 `4eaee87548250b2682affda7d05d947c05554bae`，目录树 `c9ca7036ef4365d1a74b260de5581508a75ca78f`。这是无相关性筛选队列位置41，family157。状态是 **partial_body_read**，不是full。

全部20个text blobs已fetch，并按原始Git blob SHA及tree size核验一致；源码缓存根为 `/workspace/scratch/66870cb2f937/source_context/fullread_round03_22/41`。总370905bytes、8065lines。其中15个文件、229242bytes、5126lines由本读者完整分块阅读，所有输出未截断。只下载、hash核验的内容没有算作阅读。

README的PDF是paper.pdf，actual source root是build/main.tex；main输入preamble及15个section/appendix文件，intro另嵌入01-history，bibliography为build/sources/references.bib。intro103–125的TikZ four-hole图所有节点、边、文字、caption已读；没有其他figure input或未知include。PDF binary本身没有读、没有编译或逐页视觉核验。

## 已读

| 文件 | 实际阅读块（inclusive lines） |
|---|---|
| README.md | 1–18 |
| build/main.tex | 1–41 |
| build/preamble.tex | 1–70 |
| build/sections/01-history.tex | 1–103 |
| build/sections/01-introduction.tex | 1–232 |
| build/sections/02-geometry.tex | 1–321 |
| build/sections/03-distributions.tex | 1–250 |
| build/sections/03-frame-laws.tex | 1–250；251–482 |
| build/sections/04-moments.tex | 1–230；231–418 |
| build/sections/05-realization.tex | 1–240；241–450；451–623 |
| build/sections/06-collision.tex | 1–291 |
| build/sections/07-phases.tex | 1–260；261–530；531–791 |
| build/sections/08-preparation.tex | 1–250；251–500；501–716 |
| build/sections/14-parameters.tex | 1–260；261–505 |
| build/sources/references.bib | 1–265 |

## 必须继续阅读

以下五个body文件共141663bytes、2939lines。全部已在上述缓存中，SHA及size均正确，但本读者没有看其正文，不能由intro/appendix的roadmap补算：

| 文件 | 全部未读范围 | 固定blob SHA |
|---|---|---|
| build/sections/09-histograms.tex | 1–928 | f38c13fc50a33a751722218c6271db9bd2710d90 |
| build/sections/10-product-tests.tex | 1–382 | 47a433dc0913018c06d08e09249144777152035e |
| build/sections/11-common-limit.tex | 1–582 | ebe1c6b2dc5eea6639502f3920d5304bd3604522 |
| build/sections/12-status-span.tex | 1–567 | d5478451087ee1caac50b3ce3f39f08df15260cb |
| build/sections/13-rare-overlap.tex | 1–480 | 26c5da67411fcf66eb37f93635ad1023cea22ce7 |

新读者应按实际分块全读，写续读receipt；仅在行覆盖并集完整且无未映射input、无unread时，才能把full与closure设真。读取完成≠独立证明验证。请避免重新fetch或改动缓存换版本。

## 已读部分的数学结构

主定理声称存在arbitrarily large m点图，alpha≤2且connected matching cm<m/100；connected是每对matching edges均touch。一般小branch counting给h≤(m+4cm+2)/3，故h<26m/75+2/3<m/2≤chi_f≤chi。普通Hadwiger与fractional-coloring weakening均在稿件声称范围。它不声称cm=o(m)，不解决较小正线性系数的弱猜想。所有代数在F2，概率、密度、熵为real。

holes由shared tensor、两个全域affine-gradient contraction identities和opposite a-bits定义。triangle parity把六个cross-gradient等式求和，bilinear项由symmetry取消而留下0=1。sample graph在positions上取hole complement，即便raw types重复也无independent triple。

raw supersaturation假设每个unit law有两个marginal caps Mmu及joint cap2^(DN)mu²，内部两个endpoints可任意依赖，两个units独立。已读distribution→graph证明用relative-entropy minimizer作为进度函数；删除conflict-neighborhood时熵升高。无可行law时terminal cut给small endpoint set及small pair set，fingerprints数exp(o(m))与sample tail控制large connected matching。cannot replace joint cap by marginals。

exact pin peeling把law保留为实际权重的leaves，取得modpins fresh-image minentropy；leaf仅有joint cap，mixed law仍保有marginal cap，不能互换。cut-profile true intersections低rank，canonical sparse representation由strict-majority kernel归一。selector moment低rank decomposition用flat quotient上的commuting self-adjoint idempotents；projected pin支持排除保证selected label rays的nominal independence。注意projected pins不同于individual pin intersections。

gradient realization先冻结pin/key Grams，放宽响应为两个derivative项加quotient pure bilinear forms。annihilator分解为effective-space profiles及selected atoms，有限scalar recipe令目标在此annihilator上为零。explicit tester routing先提取O(r0) rank，余项用保留pins/keys/target row-column forms的blockwise compression压缩，再以private channel dot pairing factorization实现。所需实际cross Gram scalar tests仅O(n)，不是O(n²)。

collision proposition把accepted key subprobability的overlap变成真正four-hole事件。小record cells在固定parameter pair时prune再average，不能把不同opposite parameters的survivors全部求交。key-vector equality主成本2^(−kN)，scalar bits<.01N、record logcount<.005N；quotient character rank t利用cell minentropy给Walsh2^(−.45tN)，无需独立Gram bits。必须跟同一orientation/key/witness对象保留所有条件。

phase alternative用no-cover rank-growth moment或bounded tiny cover。tiny-cover kernel F(A,B)=u_B,sum(Delta_A) rank≤d0N；若F(A,A)=0，则(1+F)⊙(1+F)^T在incompatible subset上为I，其rank平方界给compatible pairs Ω(N^−4)。singleton phases exponential小，Fourier inversion保留four-bit条件。额外support/cover evaluation pins使phase event由leaves固定，relative injection适用于rare accepting families，不能joint condition后再声称units独立。

common prediction给p_i(x)=T(lambda_i,x)+alpha_i a(x)+f_i(key)且f_i对全部retained units相同。zero classes先用degree2 Boolean support gap exactify，再用selector degree排除少量bad labels，最后奇数tag shared-only cut sum消gamma*q。Delta=0时alpha1=alpha2证明构造共同synthetic keys；只在small affine slices实现local tests，得到每tag cross-coefficient matrix rank≤2s1，却总和为I_(2gs1+1)，矛盾。common image-array test不同于unit-specific mark-dependent predicate。

preparation三支是no prediction、positive compatible table mass、zero-class rare Ω(N^−4) leaf pairs；后者选每对leaves的specific compatible injecting table，admissibility分解为两个unary orientation filters，个别filter下界a0依赖固定table数。channel tolerance在table数之前固定，所以a0不反向影响h。

appendix完整已读。构造参数链为early law/pin budgets→selector/mixers→r0,h→M0→n，mixers每n先于unit law选择。N=M0n且M0固定，O(n)不是o(N)，但系数能由M0最后选择变成指定小fraction。后续test accuracies/partition sizes只增加n onset，不能反向改J/h/pin budgets。

## 条件构造接口

records.json保留6条详尽idea，每条有精确premises、mechanism、possible bridge、missing step与source anchors：entropy law/cut机制；annihilator及保约束压缩；selector labels避projected pins；lowdegree exactification与共同test；低rank kernel保留稀有compatibility；acyclic vector/scalar/概率预算。

它们都没有建立π的BA或nonBA。π目标的q²|π−p/q|仍未证明；现有Q(Z)需求是同一个integer determinant的nonzero与real analytic smallness，以及row rank、candidate common bases、columnheight和nonadditive quadratic occupancy同时成立。F2非零minor只能在具体integer reduction正确时证integer nonzero，不能给real smallness；rational compression亦须另外追踪integer lattice index及height。远接口保留，不能用类比替代missing bridge。

## 字节与收据

`records.json`含所有20源路径、原SHA、总行数、实际read ranges、实际chunk ranges、五个unread ranges及source mapping。`byte_audit.json`逐文件记录size/SHA和总数。apply_patch写缓存时每文件多加了一个terminal LF；仅在去掉该字节后Git blob SHA恰好恢复原SHA的条件下去掉，20个均已恢复原始bytes；没有正文改写。此存储修复没有计作阅读。
