# gcd 爆破：额外约分反例、二进赋值与有限排除

2026-11-02。接续上一轮的整数线性形式与实际约分阈值，直接攻击稀疏的大公因子能否发生。本轮得到全指标的素数赋值约束、一个明确的算术反例，以及覆盖 n≤20000 的完整有限排除。没有构造出质量趋零的无限子序列，也没有证明这一个族在所有 n 上都会失败。

## 固定对象

仍令两列整数分别满足

\[
X_{n+1}=(2n+1)X_n+n^2X_{n-1},\qquad n\ge1,
\]

初值为 \(A_0=A_1=1\)、\(B_0=0,B_1=1\)。定义

\[
g_n=\gcd(A_n,4B_n),\quad
p_n=4B_n/g_n,\quad q_n=A_n/g_n,
\quad Q_n=q_n^2|\pi-p_n/q_n|.
\]

这里的 p、q 已完全约分。上一轮已经证明 \(Q_n\sim(\sqrt2/n)(n!/g_n)^2\)，以及以下所有 n≥1 都成立的严格界：

\[
\frac{2}{2n+1}\left(\frac{n!}{g_n}\right)^2
<Q_n<\frac{4}{2n+1}\left(\frac{n!}{g_n}\right)^2.
\]

本轮的目标是进攻 g 本身，而非只延长 π 的小数位。

## 1. 一个试图封死奇素数约分的捷径，被明确反例否定

本轮尝试的简单上界是：g_n 的奇数部分是否总能整除 n! 的奇数部分？答案是否定的。第一个反例出现在

\[
\boxed{n=72,\quad p=17,\quad
v_{17}(A_{72})=v_{17}(B_{72})=5,
\quad v_{17}(72!)=4.}
\]

因此 \(v_{17}(g_{72})=5\)，比阶乘多出一个 17。这里 \(v_p\) 表示素数 p 的指数，允许作用于非零有理数。这是**对该简便上界的反例**，不是对 π 坏逼近性的反例；前一轮没有把这个新尝试的上界当作已证明结果。

精确整数及独立模运算证书见[奇素数分块分析](gcd-prime-blocks/ODD_PRIME_BLOCKS_AND_EXCESS.md)和[反例与奇数超额门槛](gcd-valuation-theory/ODD_EXCESS_BARRIER.md)。重复 base-17 零数位也未形成一个持续有效的无限反例族。

## 2. 2 的幂次可以对所有指标精确求出

对每个整数 m≥0，证明了

\[
\boxed{
\begin{aligned}
v_2(A_{4m})&=2m, &v_2(A_{4m+1})&=2m,\\
v_2(A_{4m+2})&=2m+2,
&v_2(A_{4m+3})&=2m+3+v_2(m+1).
\end{aligned}}
\]

偶数指标通过反向二项式和式控制；较难的奇数指标通过一个形式倍角恒等式与五项多项式证书完成。证明中使用的逆二项式级数在 2-adic 意义下取值为零，与它的实数求和值不同；两种极限没有混用。完整论证见[二进赋值证明](gcd-valuation-theory/two-adic-proof/TWO_ADIC_NUMERATOR_VALUATION.md)。

令 u_n 为有理数 g_n/n! 约分后分子的奇数部分。它只计算超过阶乘的正奇素数增益，并不是 g_n 的整个奇数部分。由于 g_n 整除 A_n，上述公式给出

\[
\frac{g_n}{n!}\le8n^2\,2^{-n/2}u_n.
\]

与严格质量下界结合，得到统一的全指标定理：

\[
\boxed{Q_n>\frac{2^n}{96n^5u_n^2}\qquad(n\ge1).}
\]

因此，任何质量有界的候选都要求 u_n 至少具有 \(2^{n/2}\) 的指数主尺度，只能差一个多项式因子。若能证明 u_n 在所有充分大 n 上满足某个低于这个临界指数的上界，就能排除本族内的一切质量趋零子序列。**这一步所需的奇数部分总增长界尚未得到。**

这条结果将“需要很大的 g”细化成了可区分的机制：2 的幂次有确定的指数损失，必须由大量或高次的异常奇素数约分补偿。

## 3. 奇素数异常必须通过具体的数位与邻近检验

设 \(\ell_n=A_n/n!\)、\(w_n=B_n/n!\)，对奇素数 p 定义超额

\[
e_p(n)=v_p(g_n)-v_p(n!).
\]

由精确有理 Wronskian 及 w 的分母界，证明对 n≥2：

\[
\boxed{e_p(n)\le\min\left\{
\lfloor\log_p(n-1)\rfloor-v_p(n),
\lfloor\log_p(n+1)\rfloor-v_p(n+1)
\right\}.}
\]

这里的整数对数表示满足 p^h≤对应整数的最大 h。若 e_p(n)>0，则更强的是 \(v_p(\ell_n)\) 本身也不超过这个上限。超过上限的素数可以从超额候选中整项删除。

经典 Schur 型数位分解还给出 \(\ell_{bp+r}\equiv\ell_b\ell_r\pmod p\)。结合 w 的生成函数，设 \(n=bp^h+r\)，其中 \(h=\lfloor\log_p n\rfloor\)、\(1\le b<p\)、\(0\le r<p^h\)，则

\[
p^h w_n\equiv(-1)^{h(p-1)/2}\ell_r w_b\pmod p.
\]

右端非零时，直接得出 \(e_p(n)=-h\)：此时约分相对阶乘有亏损。出现正超额必须先消去这个主极点，再通过更高次的模运算检验。数位零点本身不足以认证正超额，也不能把零数位的个数直接当作赋值。

证明、素数幂转移分块及古典结果的出处见[奇素数分块](gcd-prime-blocks/ODD_PRIME_BLOCKS_AND_EXCESS.md)。更强的邻近极点条件及无需因数分解的整素数删除算法见[局部极点筛法](gcd-integrality-obstruction/LOCAL_POLES_AND_EXCESS_SIEVE.md)。

## 4. 真正执行的有限爆破

新增筛法只需计算整数 \(U_n=2^n\ell_n\) 与一个最小公倍数，便可得出 Q_n 的严格下界。它不需要 B_n 的巨大整数，也不需要任何 π 的近似值。仍未排除的少数指标，再用实际 gcd 处理。

| 检查 | 严格结果 |
|---|---|
| 逐项计算 A、B、n!、g，覆盖 1≤n≤5000 | 本族最佳指标为 n=2，对应 3/1；所有其他指标有 Q_n>9/14 |
| U／lcm 必要条件筛选，覆盖 1≤n≤20000 | 19995 个指标被直接排除于 Q_n≤1/5，剩余 2、3、6、7、15 |
| 对五个剩余指标作实际约分 | 只有 n=2 未被排除，且 1/10<Q_2<1/5；故 3/1 是完整检查范围内的唯一最佳 |
| 另一套独立实现核查小范围 | 模素数公式、异常反例、筛法与 2-adic 多项式证书均通过 |

扫描至 5000 时，观察到 u_n 最大为 2677，出现在 n=4762。我们没有将这个有限现象外推为 u_n≤n 或多项式增长的定理。

这些范围都是**这一固定构造的指标范围**，不是对分母 q≤20000 的全部有理数穷举。筛法的保留项也不是成功候选证书，必须继续经过精确第二阶段。完整证明、算法和记录见[有限排除说明](gcd-blast-search/SIEVE_THEOREM_AND_FINITE_EXCLUSION.md)。

## 尚未闭合的部分

现在缺口可以精确表述：能否统一控制各个异常奇素数同时出现时的乘积，使 u_n 的增长低于 \(2^{n/2}\) 的临界尺度；或者反过来构造并认证真正超过该尺度、且保留全部质量收益的稀疏指标族。

本轮没有证明这两种情形中的任何一种。单独排除这一族也不能证明 π 坏逼近，因为后者要求覆盖所有既约有理数。本轮未获得新的 p/q 纪录。

## 证明、来源与复现

- [独立 gcd 理论审查](gcd-integrality-obstruction/INDEPENDENT_GCD_THEORY_REVIEW.md)。
- [第二份独立审查](gcd-audit/INDEPENDENT_GCD_REVIEW.md)，另有独立标准库验证。
- [来源核对与适用范围](SOURCES.md)。本轮重新识别出该递推属于经典反正切连分数及 Legendre 构造，不把重建既有族称为新发现。
- [运行记录与文件校验和](CHECKS.json)。所有新增复现入口仅使用 Python 标准库；全指标结论依赖证明，有限扫描不承担外推。

从本目录运行：

```bash
python -B gcd-valuation-theory/two-adic-proof/verify_two_adic_valuation.py
python -B gcd-valuation-theory/verify_odd_excess_counterexamples.py
python -B gcd-prime-blocks/verify_odd_prime_blocks.py
python -B gcd-integrality-obstruction/verify_local_poles_sieve.py
python -B gcd-blast-search/scan_gcd.py --limit 5000
python -B gcd-blast-search/sieve_gcd.py --limit 20000
python -B gcd-audit/verify_independent_gcd.py
```

“独立审查”指其他研究代理重新推导及运行核验，不是外部同行评审。上轮的[精确整数线性形式与联合质量公式](../pi-joint-arithmetic-2026-11-02/README.md)仍是本轮的明确输入。
