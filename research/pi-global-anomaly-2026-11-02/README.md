# 全局异常引理：精确预算与内部区间约化

2026-11-02。集中攻击上半区间第一层异常总量，不继续提高上一轮的弱指数保证。本轮核实了全局预算约化，并证明一个统一的边缘质量界；**仍未证明 S(n)≤n/6 对所有充分大的 n 成立，也未排除当前递推族的全部成功子序列。**

## 固定对象与待证目标

整数列 A、B 均满足

\[
X_{n+1}=(2n+1)X_n+n^2X_{n-1},\quad
A_0=A_1=1,\quad B_0=0,\quad B_1=1.
\]

定义 g_n=gcd(A_n,4B_n)，p_n=4B_n/g_n，q_n=A_n/g_n，

\[
Q_n=q_n^2|\pi-p_n/q_n|,
\quad e_p(n)=v_p(g_n)-v_p(n!),
\]
\[
E_n=\prod_{p\text{ odd}}p^{\max(e_p(n),0)},\qquad
D_n=\prod_{p\text{ odd}}p^{\max(-e_p(n),0)}.
\]

研究的第一层异常量是

\[
S(n)=\sum_{\substack{n/2<p<n\\p\text{ odd prime},\ p\mid A_{n-p}}}\ln p.
\]

这里没有把第一层异常视为已经获得正超额。目标仍是一个逐指标的全局上界；平均界或密度为一的结论不能排除极稀疏指标。

## 1. 预算约化成立，而且有精确整数证书

令 H_n 为上半区间全部奇素数的乘积，P_n=exp(S(n)) 为异常素数乘积，L_n 为 lcm(1,…,n−1) 的奇数部分。对每个 n≥2，证明

\[
\boxed{E_n\operatorname{odd}(n)H_n^2
\mid D_nP_n^2L_n.}
\]

其中 odd(n) 表示正整数 n 的奇数部分。这是各素数赋值的不等式，不使用素数定理；它保留了所有常规亏损以及区间端点。取对数后再使用经典素数定理，才得到

\[
\boxed{\ln(E_n/D_n)\le2S(n)+o(n).}
\]

如果另用 T(n) 统计上半区间真正达到 p²∣g_n 的素数权重，还可加强为 ln(E_n/D_n)≤S(n)+T(n)+o(n)，其中 0≤T(n)≤S(n)。这项加强没有改变本轮优先研究第一层 S(n) 的目标。[证明、严格有限 Q 下界及端点审核](global-s-budget/GLOBAL_FIRST_LAYER_BUDGET.md)。

因此，S(n)≤n/6 的假设确实足以给出

\[
Q_n\ge\exp\bigl[(\ln2-2/3+o(1))n\bigr]\longrightarrow\infty.
\]

更一般地，只要证明 limsup S(n)/n<ln2/4≈0.1732868 即足够；不必拘泥于恰好 1/6。这个充分条件本身尚未得到。

## 2. 可以统一消掉两端的小指标质量

经典反射关系允许将每个异常折叠到

\[
k(n,p)=\min\{n-p,\ 2p-n-1\}\le(n-1)/3,
\quad p\mid A_{n-p}\iff p\mid A_{k(n,p)}.
\]

令 C_k 是整数 2^k A_k/k! 的奇数部分。现有精确二进赋值给出

\[
\boxed{C_k\le(2+\sqrt2)^k.}
\]

对任意整数 K≥0，只统计折叠指标 k(n,p)≤K 的异常权重，记为 S_edge(n,K)。这些不同素数的乘积整除 ∏_(k=0)^K C_k，于是

\[
\boxed{S_{\rm edge}(n,K)
\le\frac{\ln(2+\sqrt2)}2K(K+1).}
\]

这是真正覆盖所有 n 的结论。取 K=⌊n^(1/3)⌋，便有

\[
S(n)=S_{\rm interior}(n)+O(n^{2/3}),
\]

其中 S_interior 只统计 n^(1/3)<k(n,p)≤(n−1)/3 的异常。两端这部分不会贡献线性规模；剩下的统一问题已经明确集中到内部区域。若能证明 S_interior(n)≤n/6+o(n)，同样足以排除本族，但这项内部界仍未证明。

还证明异常必须满足 k(n,p)>ln(n/2)/ln(2+√2)，从而两个端点附近确有无异常的对数宽度条带。[完整反射、乘积证明与适用范围](global-s-reflection/REFLECTION_AND_INTERIOR_REDUCTION.md)。

## 3. 实验专门审查聚集结构

完整精确扫描扩展至 1≤n≤100000，另以不使用反射的不同实现独立重算至 10000。整数证书确认有限范围内均有 S(n)<n/6，唯一最大比值仍为 S(21)/21=ln17/21≈0.134915。

100000 范围内共有 7026 个异常事件、6800 个非空指标；至多三个异常素数同时出现。六个三素数聚集的指标为 2772、9559、29158、29259、40219、73219。对这六处的全部 18 个素数再检查模 p²，均得到 v_p(g_n)=1，即相对 n! 的超额为零；最大比值处 n=21、p=17 也不是正超额。

这些结果描述了最危险的有限聚集，**不推断三素数是全局上限，也不推断无穷多个指标仍满足 n/6**。关于反射图和中点素数的额外诊断同样保留有限范围。[实验说明](global-s-experiments/EXACT_ANOMALY_CLUSTERS.md)与[精确数据](global-s-experiments/global_s_cluster_checks.json)。

## 尚缺的机制

仅有反射、单个模数的零点稀疏和相邻不同时为零，仍不足以约束不同素数在同一 n 上命中。独立审核给出了一个抽象零集模型：每个素数至多两个零点、满足反射且没有相邻零点，仍可在稀疏指标产生约 n/2 的对角权重。这个模型不满足原递推，**不是本族引理的反例**；它准确说明缺少的是跨模数的算术约束，不能把单素数零点计数当成统一引理的证明。[审核第 5 节](global-s-review/INDEPENDENT_GLOBAL_S_REVIEW.md)。

下一处具体断点是内部区域的逐指标集中上界。需要从原递推本身导出跨素数的乘积或排斥约束；目前没有得到这样的定理，也不能用独立随机模型替代。

## 复现与出处

- [来源与前置证明](SOURCES.md)。经典 Schur／Holt 对称和素数定理均有明确出处，不将其重建称为新的一般结果。
- [独立审核](global-s-review/INDEPENDENT_GLOBAL_S_REVIEW.md)是另一研究代理重新推导并核算，不是外部同行评审。
- [运行记录及文件清单](CHECKS.json)。有限核验与全指标证明分开记录。

从本目录运行，均仅用 Python 标准库：

```bash
python -B global-s-budget/verify_global_s_budget.py
python -B global-s-reflection/verify_reflection_interior.py
python -B global-s-experiments/scan_global_s_clusters.py --limit 100000 --direct-check-limit 10000
python -B global-s-review/verify_independent_global_s.py
```

本轮没有提高保守逼近指数，没有更换递推族，没有获得新的 π 有理逼近纪录。即使将来证明统一引理，也只排除这一族，不会单独决定 π 的坏逼近性。
