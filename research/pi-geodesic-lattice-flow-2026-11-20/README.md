# π 的连续格流与低冗余有理逼近

本目录把用户最初的“一维连续参数—二维整数格点—有理数”想法改写为连续的 unimodular lattice flow。

- [连续格流的定义与坏逼近等价式](CONTINUOUS_LATTICE_FLOW_FOR_PI.md)：显式曲线 \(A_\pi(t)\in SL_2(\mathbb R)\)，整数向量映射 \((q,p)\mapsto(qe^{-t},(p-\pi q)e^t)\)，以及最短向量与 \(q^2|p/q-\pi|\) 的精确关系。
- [最短向量下包络实验](LOWER_ENVELOPE_EXPERIMENT.md)：对 \(q\le200000\) 的有限候选集做 80 位精度下包络搜索，列出活跃向量与切换时间，并对照连分数收敛分数。
- [实验程序](pi_lattice_flow_experiment.py)：使用 mpmath 和单调凸包算法复现。
- [实验结果 JSON](pi_lattice_flow_experiment.json)：保存活跃向量、误差、切换时间及连分数对照。

**范围说明：** 连续格流与坏逼近之间的等价关系是经典的丢番图逼近/齐次动力系统表述。有限下包络结果仅针对 \(q\le200000\)，不证明无限活跃序列，也不证明 \(\pi\) 是或不是坏逼近数。

- [活跃向量与连分数定理](ACTIVE_VECTORS_AND_CONTINUED_FRACTIONS.md)：证明格流下包络中的每个活跃分母必为 \(\|q\pi\|\) 的严格纪录，因此必是正则连分数收敛分母；再给出哪些收敛分数留在下凸包上的精确判据。

- [活跃向量与连分数定理](ACTIVE_VECTORS_AND_CONTINUED_FRACTIONS.md)：证明每个活跃分母必是 \(\|q\pi\|\) 的严格纪录，因此必为正则连分数收敛分母；并给出收敛分数的下凸包生存判据。
- [仅用收敛分数的下包络实验](CONTINUED_FRACTION_LOWER_HULL.md)：推导完整商给出的精确误差表达式，并将全分母扫描降为收敛分数扫描。
- [收敛分数下包络复核程序](continued_fraction_hull_check.py) 与 [有限对照结果](continued_fraction_hull_check.json)：对 \(q\le200000\) 比较全扫描和收敛分数扫描。

- [全局活跃判据与部分商曲率](GLOBAL_ACTIVE_CRITERION_PARTIAL_QUOTIENTS.md)：给出每个收敛分母的精确活跃时间区间 \((L_n,U_n)\)，并推导相邻三收敛点的局部曲率不等式。明确局部条件必要但不充分，避免把相邻斜率测试误当成全局凸包判据。
- [部分商曲率复核程序](partial_quotient_curvature_check.py)：以高精度数值检查局部公式，并与有限下包络结果对照。

- [部分商局部凸性与全局加权平均判据](PARTIAL_QUOTIENT_LOCAL_GLOBAL_TEST.md)：证明若下一部分商 \(a_{n+1}\ge2\)，该收敛点自动通过局部凸性；若 \(a_{n+1}=1\)，局部条件精确化为 \((r_n+z_n)^2+2r_nz_n(1+r_n+z_n)>1\)。同时把全局活跃条件改写成左右所有加权割线斜率的比较，强调局部通过不保证全局活跃。

- [流式下凸包与加权 PAVA](STREAMING_LOWER_HULL_AND_PAVA.md)：将收敛分母点按顺序输入单调栈；当相邻下降斜率不再严格递减时弹出中间点。证明每点至多入栈/出栈一次，并把整形过程识别为带权 pool-adjacent-violators 算法。
- [π 的百万分母流式实验程序](pi_lattice_stack_experiment.py) 与 [JSON 输出](pi_lattice_stack_experiment.json)：处理所有 \(q\le10^6\) 的收敛分母，列出最终有限前缀下包络、被弹出的收敛点和连续切换时间。

- [部分商 1 的精确局部阈值](SHARP_LOCAL_FILTER_PARTIAL_QUOTIENT_ONE.md)：将局部凸性化简为 \(\alpha_{n+2}< (\beta_n+2)/(\beta_n-1)\)，其中 \(\alpha_{n+2}\) 是前向完全商、\(\beta_n=q_n/q_{n-1}\) 是反向有限商；若 \(a_{n+1}\ge2\)，局部凸性自动成立。
- [阈值复核程序](partial_quotient_one_threshold_check.py) 与 [100 位数值结果](partial_quotient_one_threshold_check.json)：检查 \(\pi\) 的部分商为 1 的收敛点，包括 106、33102、33215、66317、265381。

- [有限未来否决证书定理](FINITE_FUTURE_VETO_CERTIFICATE.md)：证明对每个固定收敛点，未来否决阈值 \(F_n\) 的下确界必在有限未来指标处达到；给出由 \((q_{m+1}^2-q_n^2)/E_n^2\ge U_m\) 触发的可验证停止条件。区分逐点有限前瞻与尚未证明的统一有界前瞻。
- [π 的未来否决停止实验](future_veto_stopping_experiment.py) 与 [数值记录](future_veto_stopping_experiment.json)：展示 \(q=99532\) 的候选经两步未来检查后通过，以及 \(q=265381\) 被 \(q=364913\) 否决。

- [由两个部分商控制的显式前瞻上界](EXPLICIT_LOOKAHEAD_FROM_PARTIAL_QUOTIENTS.md)：证明若 \(a_{n+2}\ge2\)，最多检查两个未来收敛分数；若 \(a_{n+2}=1\)，令 \(A=a_{n+3}\)，则最多检查 \(k(A)-1\) 个，其中 \(k(A)=\min\{k:F_k\ge(A+2)/\sqrt{2A+3}\}\)，故前瞻长度为 \(O(\log A)\)。
- [π 的前瞻上界实验](explicit_lookahead_pi_check.py) 与 [20 个候选的 JSON 记录](explicit_lookahead_pi_check.json)：比较理论上界和实际停止长度；数值仅作有限诊断。

- [双边选择器的有限传递矩阵公式](FINITE_TRANSFER_MATRIX_TWO_SIDED_SELECTOR.md)：将过去/未来阈值同时除以 \((q_n/E_n)^2\)，用 \(M(a)=\begin{psmallmatrix}a&1\\1&0\end{psmallmatrix}\) 与 \(N(a)=\begin{psmallmatrix}-a&1\\1&0\end{psmallmatrix}\) 的有限乘积计算全部归一化比值；未来长度由显式前瞻定理控制，剩余唯一的尾部参数是完全商 \(z_n=1/\alpha_{n+2}\)。
- [双边传递矩阵复核程序](two_sided_transfer_matrix_pi_check.py) 与 [π 的有限测试数据](two_sided_transfer_matrix_pi_check.json)：将矩阵公式与直接收敛分数比值对照，检查前 25 个候选。

该双边选择器还可用有限连分数柱区间包住唯一尾部参数 \(z_n=1/\alpha_{n+2}\)，并通过外向舍入区间运算认证选择器符号；见上述传递矩阵笔记第 7 节。

- [椭圆形式群与 adelic 格流桥接](ELLIPTIC_FORMAL_GROUP_ADELIC_LATTICE_BRIDGE.md)：从 Pasten–Vidaux Lemma 4.1 推出同时实数逼近与有限素数处 \(p\)-adic 同余的加强引理；将同余条件识别为有限指数子格；并证明固定同余子格的坏逼近常数只是 \(M\kappa((\pi-k)/M)\)，所以朴素局部筛选不会解决原问题。明确指出真正缺少的是椭圆见证的定量高度/分母控制。

- [任意阶椭圆形式群消项](ELLIPTIC_HIGH_ORDER_FORMAL_CANCELLATION.md)：构造整数权重 \(\lambda_{m,j}=(-1)^j\binom{m+1}{j-1}\)，证明 \(\sum_j\lambda_{m,j}[j]_F(T)=T+O(T^{m+1})\)。应用于椭圆曲线，得到任意阶实数逼近，同时保留所有指定 \(p\mid M\) 处的同余。附有粗略高度预算，并明确尚无足够的分母控制来改善 \(\pi\) 的逼近。
- [权重的精确整数验证程序](verify_elliptic_high_order_weights.py) 与 [\(m=1,\ldots,30\) 的检查结果](elliptic_high_order_weight_checks.json)：验证插值恒等式、乘法形式群的首个非零误差系数及权重范数。
