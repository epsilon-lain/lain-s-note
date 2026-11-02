# 可用素数池：中点的零密度与平均零点数目标

2026-11-02。接续 [实际初值筛与可用素数池](../pi-initial-value-sieve-2026-11-02/README.md)，固定同一个递推族。**尚未证明可用素数池的最终容量上界，也没有构造新的成功逼近。** 本轮严格解决的是平均零点数约化中的中点修正项，并核查了观测约 0.4 的文献解释。

## 1. 研究对象及已经成立的充分条件

\[
U_0=1,\quad U_1=2,\qquad
(n+1)U_{n+1}=2(2n+1)U_n+4nU_{n-1},
\]
\[
\mathcal Z=\{p\text{ 奇素数}:\exists\,0\le r<p,\quad p\mid U_r\},\qquad
\vartheta_{\mathcal Z}(x)=\sum_{p\le x,\,p\in\mathcal Z}\ln p.
\]

此前的约分逼近记为 \(p_n/q_n\)，质量为 \(Q_n=q_n^2|\pi-p_n/q_n|\)。已有多区间亏损定理给出

\[
\ln Q_n\ge(2+\ln2)n-4W(n)+o(n),\qquad W(n)\le\vartheta_{\mathcal Z}(n).
\]

因此用户提出的 \(\vartheta_{\mathcal Z}(x)\le((2+\ln2)/4-\eta)x\) 确为排除本族全部质量有界子序列的充分条件。这里 \(\eta>0\) 固定，必须对所有充分大的 \(x\) 成立。前次十万项有限证书不替代此统一条件，本轮没有重跑或延长那次搜索。

## 2. 新的严格步骤：中点总权重没有主项

对每个奇素数 \(p\)，定义

\[
z_p=\#\{0\le r<p:p\mid U_r\},\qquad
m_p=\mathbf1_{p\mid U_{(p-1)/2}}.
\]

Holt 反射把非中点零点配成两两一组，所以

\[
z_p=2k_p+m_p,\quad k_p\in\mathbb Z_{\ge0},\qquad
\boxed{\mathbf1_{p\in\mathcal Z}\le\frac{z_p+m_p}{2}.}
\]

中点不能逐素数直接删除：例如 \(p=31\) 时恰有一个零点。

但是它对应于固定椭圆曲线

\[
E:\ y^2=x(x^2+2x+2)
\]

的 Frobenius 迹为零。这条曲线没有复乘；一个短证书是它在 \(p=3,5\) 的普通约化分别有迹 \(2,-2\)，Frobenius 二次域分别为 \(\mathbb Q(\sqrt{-2})\)、\(\mathbb Q(i)\)。若原曲线有复乘，复乘域必须同时嵌入并等于这两个普通约化的有理自同态代数，矛盾。

因此可调用 Serre 的非复乘迹零素数定理，无条件得到

\[
\boxed{\sum_{p\le x,\,p\text{ 奇}}m_p\ln p=o(x).}
\]

这一步无需假设 GRH 或 Lang–Trotter 猜想。特别地，具有奇数个数位零点的素数集合也有零对数密度，因为 \(z_p\) 为奇数当且仅当 \(m_p=1\)。这只处理不成对的中点，不控制其余成对零点。

完整证明及定理使用范围：[PAIRED_ZERO_MEAN_AND_MIDPOINT.md](prime-pool-structure/PAIRED_ZERO_MEAN_AND_MIDPOINT.md)。

## 3. 下一步可以精确改写为一个第一矩上界

定义加权零点数

\[
T(x)=\sum_{p\le x,\,p\text{ 奇}}z_p\ln p.
\]

上述结论严格推出

\[
\boxed{\vartheta_{\mathcal Z}(x)\le\tfrac12T(x)+o(x).}
\]

所以一个明确的充分目标是：存在固定 \(\eta>0\)，使

\[
\boxed{T(x)\le\left(1+\frac{\ln2}{2}-\eta\right)x
=(1.3465735902\ldots-\eta)x}
\]

对所有充分大的 \(x\) 成立。届时

\[
\ln Q_n\ge(2+\ln2)n-2T(n)+o(n)
\ge2\eta n+o(n)\longrightarrow+\infty.
\]

**第一矩上界仍然未证明。** 它把“至少有一个零点”的并集问题转成加法统计，但还需要控制零点跨素数出现的频率。反射本身不提供这样的平均估计。仅使用反射且没有中点时，系数 \(1/2\) 也已是最好的：一个非中点零点对就使不等式取等。

## 4. 现有数据与约 0.4 的含义

重新汇总前次已经存在的精确分类表，没有新增素数分类：

| \(x\) | \(T(x)/x\) | \(\vartheta_{\mathcal Z}(x)/x\) | 有限配对上包络 \((T+\sum m_p\ln p)/(2x)\) |
|---:|---:|---:|---:|
| 1000 | 0.993576 | 0.391316 | 0.507268 |
| 5000 | 1.068361 | 0.407169 | 0.537062 |
| 10000 | 1.060914 | 0.407388 | 0.532777 |
| 20000 | 1.018619 | 0.399696 | 0.511203 |

这些是有限描述值。它们未证明 \(T(x)\) 有极限、未证明最终系数小于 1.346573，也不提供有效的渐近起点。

Kohen (2026) 的 Remark 1，第 7 页，已经对通常中央三项式序列提出相关概率启发：反射后约 \(p/2\) 次零点试验，如果把它们当作独立、每次概率 \(1/p\)，预测活跃比例为 \(1-e^{-1/2}\approx0.393469\)。本序列数据与这一启发相容，属于推测性比较；它不是针对本序列的密度定理。论文固定素数后的指标密度结论也不等于本题变化素数的容量结论。

文献范围核查：[PRIMARY_SOURCE_SCOPE.md](prime-pool-literature/PRIMARY_SOURCE_SCOPE.md)。

## 5. 可复现范围及保留的缺口

- 主推导程序逐素数检查至 2000，含反射和独立椭圆曲线点计数；实际有限计数验证不替代 Serre 定理。
- 独立审计另在有限域中展开二次多项式的幂计算序列，并检查平均约化、非复乘证书和字符矩公式。详见 [独立审查](prime-pool-audit/INDEPENDENT_PRIME_POOL_AUDIT.md)。
- 对前次 2261 行分类表检查完整反射、奇偶数与零点计数；文件 Git blob 摘要必须匹配前次提交。对数只用于表格展示。

在本目录运行：

```bash
python prime-pool-structure/verify_paired_zero_mean.py
python prime-pool-audit/verify_prime_pool_independent.py
python summarize_existing_pool.py
```

第三个入口读取相邻的前次研究目录中的 `active-primes/active_prime_roots.csv`。文件摘要及各检查命令见 [CHECKS.json](CHECKS.json)。

下一处未闭合目标就是上述 \(T(x)\) 的统一第一矩上界，或者能直接控制 \(\vartheta_{\mathcal Z}(x)\) 的额外算术结构。现有局部对称、固定素数 Lucas 数位规律以及概率独立模型均不能替代它。即使此路线完成，也只排除当前逼近族，不证明 \(\pi\) 坏逼近。
