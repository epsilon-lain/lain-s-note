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
