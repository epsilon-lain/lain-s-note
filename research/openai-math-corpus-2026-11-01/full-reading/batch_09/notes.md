# Batch 09 — position 21

固定提交：`4eaee87548250b2682affda7d05d947c05554bae`。作者主张与证明路线按源码首读记录；不作独立认证。

## 阅读范围

完整首读固定提交中的 README、主文件、摘要、八个主文节（00–06、09）、两个完整附录（07、08）、正文嵌入 proof-map 图、306 行参考文献表，以及 verification/README.md，共16个文本文件5558行233552字节。所有缓存按 Git blob SHA 与 tree size 核实；主文递归输入闭包全部映射并读完。main.tex 的 glyphtounicode-cmex 为 TeX 发行版的字形映射辅助输入，不是仓库数学正文。full_body_read 仅表示源码首读，不表示独立证明核验；未编译/逐页比对PDF、未运行有限计算程序、未验证外部引用。29和33仅为后续读者缓存，不在本批阅读记录内。

| Source | Lines read | Git blob SHA |
| --- | ---: | --- |
| README.md | 1–32 | c325604eb268d9fc5bb4939a891c703ba2a424fe |
| build/abstract.tex | 1–5 | f4e2ade72e1338638625be194b03301477dce6f4 |
| build/figures/proof-map.tex | 1–42 | 48c7bb7d8758aa6af9a3576aad56dbac8a45d377 |
| build/main.tex | 1–57 | 2e74f2aff1a505f414a784d348499edb221f5ab2 |
| build/references.bib | 1–306 | 1616df468f4118e100dc4ce53ea30a3db1851d7d |
| build/sections/00-introduction.tex | 1–201 | 1d4b1e83a4c3fba5e5fdaf53cfbc1d90343bae52 |
| build/sections/01-cones.tex | 1–212 | 2d8ec5007ff5c42500d916efef804478eb901d97 |
| build/sections/02-projections.tex | 1–333 | 244e740205d8cac2fd35de12118ec1167ceb6350 |
| build/sections/03-entropy.tex | 1–303 | c44fe12e8bf903bad78a613165a1e45d2b452e5d |
| build/sections/04-quadratic.tex | 1–629 | dbd08cfa51c2099162530060aaafd89b0a9e7f21 |
| build/sections/05-layers.tex | 1–621 | 8014ea471d8f684fdffa6374ea0f1281613285dc |
| build/sections/06-linear-equality.tex | 1–473 | 3db21d724e60b14c9ef99b51e5a2bf9adb1e7757 |
| build/sections/07-scalar.tex | 1–1012 | 056bae738fa56b94f314099cbb3b53d720517b39 |
| build/sections/08-segments.tex | 1–1196 | 6d41f4ffe977032d41bf54d7588b1d21084a0010 |
| build/sections/09-consequences.tex | 1–88 | 13093b433ccf0faa9108d6e48aee0423b73d1250 |
| verification/README.md | 1–48 | 045bf48992f2e7da2393cab90104d25e1bddd08d |

`main.tex`的所有仓库内数学输入、两附录、`00-introduction.tex`嵌入的图及bibliography均已闭合。读取`02-projections.tex`时一次工具输出中断在307附近，已以300–320行完整补读；最终coverage为1–333，无缺口。

## 主张原文范围

稿件主张：每个 n≥1 的一般凸体 K⊂R^n 都满足 |K| |(K-s(K))°|≥(n+1)^(n+1)/(n!)²，等号当且仅当 K 为单纯形；任意内部中心的对应乘积也有同一下界。无中心对称/边界光滑假设。另由引用的全维几何到函数蕴涵得任意proper l.s.c.凸势φ、0<∫e^-φ<∞时(∫e^-φ)(∫e^-φ*)≥e^n，并给出sharp例，但未分类全部函数等号。对两full-dimensional、essentially continuous的log-concave概率密度，所引等价给H(η1)+H(η2)≤-3n+T(ν1,ν2)。这里记录作者原文主张，不据首读认证成立。

## 证明结构

1. 在m=n+1维将K提升为solid pointed cone C，取正对偶D；chi_C(V)chi_D(U)=(n!)²/m^m·|K||(K-z0)°|，故证明U·V=m时Laplace乘积≥1即可。

2. 选Σ=I+T与Z=Σ^(1/2)G，对每个z最小化coercive Gaussian-smoothed projection potential得唯一bias ξ_z，使EΠ_C(Z+ξ_z)=a(z)U。Moreau正对偶分解给X_z−Y_z=Z+ξ_z及orthogonality，构造P_z、layer积分A、其degree-one部分L=Σ_i G_i M_i和余项R。

3. 联立坐标变换e^Q与covariance feedback。对Cbar=E mathsfC(L)、Kbar=E mathsfK(L)，以Tnew=(-Cbar−sqrt(Cbar²−4λKbar))/(2λ)更新；scalar bounds给spectral box与discriminant gap，Q球的uniform inward估计配合Brouwer得到EA=0及λT²+Cbar∘T+Kbar=0，全程不假定commutation。

4. 将Π_C、Π_D按正weights在层上积分，先截断与Gaussian平滑证明injectivity及log determinant integrability，再换元/Jensen。log-Jacobian surplus以非嵌套正kernel N/8保留；profile differential identities把成本整理成log Laplace product/m≥−mathcalE。

5. Gaussian inverse-generator covariance与两层scalar kernel比较给δ[l]、δ[l²]的精确恒等式；对d=F−|u|²且F=c+2kq组合，分出linear-field、Hermite余项和spectral threshold defects。构造β≥0并保留strict layer remainder，奇偶degree gap与weighted PSD blocks控制余项，不要求P_z nested。

6. 对L在自身spectral basis用divided-difference公式与Gaussian integration by parts，引入ν(dx,dy)权重Σ_l |(M_l)ij|²。归一化quadratic matrix equation消掉S−I项；Young/Gram估计与commutator控制之后，所有ν积分被严格scalar segment variance inequality吸收。得到mathcalE≤−.14τT²−∫(W−W0)dβ≤0。

7. 等号迫使T=0、ν支撑在diagonal、β=0。Gaussian orthogonality把ν[(x−y)²]=0变成全部确定M_i两两commute；固定basis中P_z等于L的谱阈值a.e.。选一个层后projection cross weak derivatives为0，C分裂成非退化half-lines的积，即orthant；原cone的height-one截面为simplex。反方向直接计算simplex polar volume。

8. Appendix07以stable Gaussian-tail formulas避免指数抵消，在152节点用精确有理外包，先证radicand再用q导数；analytic strip与16-node interpolation把节点推到连续区间，x≥64的导数tail bounds闭合实线。Appendix08按h与tail位置穷尽所有segments：短段按h²归一化，中段保留q variance并用cubic/degree32、30 Bernstein证书，宽段与远尾以strict budgets闭合，所有共享边界均覆盖。

## 可迁移构造机制

### 1. 联立选取规范与误差度量，使后续取消发生在同一对象上。

**精确前提：** 固定有限维solid pointed正对偶cone及内点U,V且U·V=m；T在[t_-,t_+]的spectral box、Σ=I+T>0；Cbar=E mathsfC(L)、Kbar=E mathsfK(L)连续；scalar输入给Δ=Cbar²−4λKbar≥m_*²I及matrix square sandwich；另有所有T统一的tr(Q EA)<0球面向内估计。

**机制：** 在Q球×T box上共同做Brouwer fixed point。取Tnew=(-Cbar−sqrtΔ)/(2λ)，Loewner–Heinz保box；fixed point给EA=0及λT²+Cbar∘T+Kbar=0。展开时mixed products在Jordan product中取消，无需Cbar与sqrtΔ commute；后续(S−I)项被精确消去。

**跨域/π方案桥接：** pi方案若可将候选transversal的scaling/gauge与column-cost的正定权矩阵联立选择，可先设定一个产生精确抵消的matrix equation，再证明有界盒内存在解，避免分别优化rank和height后无法对齐。

**尚缺步骤：** 当前Q(Z) rank与row transversals为代数/离散问题；Brouwer输出实参数不保证整系数、低高度或保留同一个非零整数determinant。需证明离散可实现性、分母/缩放成本及所需向内估计，不能把连续归一化直接当作determinant witness。

**源码锚点：** build/sections/02-projections.tex:175 (prop:normalization; equations(8))；build/sections/02-projections.tex:213 (covariance update)；build/sections/06-linear-equality.tex:169 (equation(25) cancellation)。

### 2. 对不嵌套的收缩矩阵族建立正阈值缺陷测度，再把误差与commutator分开计费。

**精确前提：** P_z为可测实对称PSD contractions，0≤P_z≤I，层尾部有充分矩以使积分收敛；B_*为具有全部有限矩的对称矩阵场；Θ_z=1_{(-∞,z]}(B_*)，E_z=P_z−Θ_z。测试函数及导数polynomial growth；trace权Λ为确定对称矩阵。

**机制：** B_*的谱basis中(E_z)ii(x_i−z)≥0，故β^{B_*}为正测度。无需P_z nested或与Θ_z commute即有||M_f||²_Λ=2e_Λ(f,f)−Nt_Λ(f,f)。PSD同号blocks保留trace正性；weighted layer的剩余项通过Cauchy–Schwarz计入||[B_*,T]||，把不能直接取正的权重误差显式收费。

**跨域/π方案桥接：** 可以尝试把row/column选择的Gram projector或其凸松弛与容量阈值作比较，将occupancy缺陷写成正测度，再把交互项单独计为commutator/混合代价。这一框架尤其提醒：非嵌套候选族也可保留全族信息，不能把各列独立估计。

**尚缺步骤：** Q(Z)本身不是PSD contraction；必须构造保留rank与transversal约束的Gram/projector表示。β控制的是谱阈值偏差，不自动等价于列二次occupancy或整数determinant非零。

**源码锚点：** build/sections/04-quadratic.tex:278 (spectral thresholds and positive beta)；build/sections/04-quadratic.tex:324 (lem:threshold-defect)；build/sections/05-layers.tex:135 (lem:layers-diagonal)；build/sections/05-layers.tex:395 (lem:layers-variable)。

### 3. 严格pairwise不等式把全局等号压缩成commuting结构；不要求统一strict margin。

**精确前提：** L=Σ_l G_l M_l，其中M_l为确定实对称矩阵，G_l独立standard Gaussian；正对称测度ν由m^-1 EΣ_lij |(M_l)ij|² δ_(xi,xj)定义；scalar deficit在x≠y严格为负、diagonal为0且可积；同时另有strict positive layer integrand与β≥0。

**机制：** 积分达到零迫使ν离开diagonal的质量为0；ν[(x−y)²]=m^-1Σ_lk||[M_k,M_l]||²_HS。于是M_l两两commute，可在一个固定basis同时diagonalize。β=0再将P_z压成同一阈值，最终由projection弱导数的分离推出product cone。

**跨域/π方案桥接：** 若determinant construction中某个cost下界可精确表示为非负pair interactions之和，达到极小值便可分类候选支持或block decomposition；这种先给正缺陷、后读出等号结构的做法可用于识别导致rank塌陷的唯一退化型。

**尚缺步骤：** 现有结论仅为exact equality rigidity；strict但无统一margin不能直接给near-extremizer稳定性，更不产生目标row transversal。要用在低高度非零minor上，需单独证明可积pair measure的表示与定量离退化距离。

**源码锚点：** build/sections/06-linear-equality.tex:17 (spectral pair measure nu)；build/sections/06-linear-equality.tex:341 (strictness and support)；build/sections/06-linear-equality.tex:374 (commutators and orthant rigidity)。

### 4. 有限精确证书、解析插值、无限尾部三层闭合连续参数范围。

**精确前提：** profiles在|Im t|≤.45邻域analytic且已独立给modulus bounds；δ=1/40、16-node stencil，节点由有理外包给值/second differences；Cauchy在内strip半宽.36给order16 remainder<10^-5；外部x≥64另有显式导数tail estimates。所有除法/square root使用前先证明分母/radicand下界。

**机制：** 中央chord误差≤(4387/8192)M+10^-5<.536M+10^-5；对radicand先做d插值，后计算q，避免循环。节点阈值积分靠convexity与exact piecewise-affine integration变成真实积分界。中宽segment的degree32、30多项式转Bernstein basis，全部正系数证整个闭区间；tail解析预算覆盖其余实线。

**跨域/π方案桥接：** 若pi的determinant小量或column-cost预算依赖少量连续归一化参数，可把探索得到的数值参数变成可复核的精确证书：有理节点外包、明确interpolation remainder、边界/尾部独立控制，尤其适合二次occupancy预算的参数区间验证。

**尚缺步骤：** 有限证书仅认证给定显式不等式，不能替代形式非零minor的证明。此首读未实际重算表格、运行附属Python或核验analytic strip的每个外部基础结论；移植必须为新函数重新给branch/nonvanishing与导数界。

**源码锚点：** build/sections/07-scalar.tex:138 (lem:scalar-certificates)；build/sections/07-scalar.tex:272 (lem:scalar-strip)；build/sections/07-scalar.tex:378 (lem:scalar-interpolation)；build/sections/07-scalar.tex:596 (lem:scalar-tail)；build/sections/08-segments.tex:980 (exact polynomial certificates)。

### 5. 保留辅助profile方差，以三次endpoint displacement代价吸收二次端点误差。

**精确前提：** t区间[l,y]长度2h，物理x均匀平均在dt下密度≥v_h/(2h)，v_h=h coth h−h>0；q连续且|dq/dt|≤L=.193；端点差|q_+−q_-|=2hz，z≥0。无需q monotone。

**机制：** 对非负ψ，用FTC得∫ψ(q(t))dt≥L^-1|∫_(q_-)^(q_+)ψ(v)dv|。对ψ=(q−bar q)²，端点值区间上的最小积分为(2hz)³/12，故V_q/h²≥v_h z³/(3L)。将b_a拆成V_q加均值到chord误差，保留负V_q后，剩余成本成为可在z≥0显式求最大值的cubic。

**跨域/π方案桥接：** 在有平滑插值的离散选择问题中，若column-cost的变化有Lipschitz控制且endpoint代表不同occupancy配置，可尝试以方差的三次下界惩罚大跨度选择，补偿独立quadratic bounds丢失的相关性；也适用于transport/energy的一维过渡预算。

**尚缺步骤：** 离散row transversals未给有正密度且保rank的插值路径。variance只给解析代价关系，不能在pi场景中自动保证整数性、低高度或非零determinant。

**源码锚点：** build/sections/08-segments.tex:847 (nonnegative endpoints and q variance)；build/sections/08-segments.tex:905 (cubic remaining variance cost)。

## 限制与后续

- 这是作者稿件源码首读，未独立认证一般Mahler主张，未把仓库发布状态视为同行审查或定理成立的证明。
- 全部正文与两个附录已读；finite rational tables、Bernstein系数和tail constants被作为原文证明步骤记录，未重新执行或独立重算。
- verification/README.md已读，但十个programs、runner和requirements仅有tree映射，未读其程序体、未运行；它们是正文之外的计算复现附件。
- 未编译TeX、未逐页核对paper.pdf，未逐项检查参考文献的准确性或其假设。主文数学input闭包完整；TeX发行版字形辅助输入未另读。
- 一般凸体常数与中心对称类4^n/n!的更强常数是不同结论，本文只引用对称伴随稿而不在此证明；functional e^n后果依赖引用的全维几何到函数蕴涵，不能从固定维结果或simplex equality直接转移函数等号分类。
- 可迁移构造思想保留精确前提。当前未给出针对pi方案的低成本Q(Z) rank/minor、row transversal、整数determinant的可用新定理。

本批只提交position21的阅读记录。29/33的源码已预取并按Git blob SHA核实，位于`scratch/source_context/round02_group3/{29,33}`；预取本身不算正文阅读。Root随后独立完成33，29由后续读者继续，本批不重复计数。

