# 构造接续：规定相消、提取周期、同时核算分母

2026-11-02。接续用户的新构造检查点 `db1290c827479ff85beaeca8606e8258d72f0c69`，遵循 [构造优先的当前方向](../CURRENT_DIRECTION.md)。**本轮给出了明确乘子生成规则、全规模非零周期与误差收敛证明，以及实际既约分母的联合质量上界；没有证明质量趋零，也没有恢复主攻旧素数池排除问题。**

固定用户的参数族

\[
I_m(H)=\int_0^1\frac{x^{2m}(1-x^2)^{2m}H(x^2)}{(1+x^2)^{4m+1}}\,dx,
\quad m\ge1,\quad \deg H\le m+1.
\]

对有理多项式 \(H\)，精确原函数给出有理数 \(\alpha,\beta\)，使 \(I_m(H)=\beta+\alpha\pi/4\)。当 \(\alpha\ne0\) 时，取完全约分的 \(p/q=-4\beta/\alpha\)，于是

\[
Q=q^2|\pi-p/q|=\frac{4q^2}{|\alpha|}|I_m(H)|.
\]

有理 \(H\) 可统一清分母放回整系数原函数格；整体缩放不改变上述既约分数及质量。

## 1. 原核的解析上界可以立即加强

换元 \(x=\tan\theta\) 后，恰好有

\[
I_m(H)=16^{-m}\int_0^{\pi/4}
\sin^{2m}(4\theta)H(\tan^2\theta)\,d\theta,
\]
\[
I_m(1)=\frac{\pi}{4\,64^m}\binom{2m}{m}.
\]

因此，令 \(h=\|H\|_1\)，得到更强的完全有理上界

\[
\boxed{Q<\frac{22}{7}\frac{q^2h}{|\alpha|}
\frac{\binom{2m}{m}}{64^m}.}
\]

解析因子由旧的 \((2/27)^m\) 改善为 \(O(16^{-m}/\sqrt m)\)，相应充分目标为

\[
\boxed{q_m^2\|H_m\|_1/|\alpha_m|=o(16^m\sqrt m).}
\]

这个高度目标尚未达到。对用户已有 \(m=6\) 候选，新范数上界比旧范数上界改善约 15.6365 倍，但仍约为 \(1.57\times10^7\)，远大于有符号积分证书给出的实际质量约 0.638。因此范数改善真实存在，却仍不能替代有符号相消。

证明：[精确核与相消搬运](constructive-cancellation/RODRIGUES_CANCELLATION_AND_KERNEL_NORM.md)。[既有候选的界比较](existing_candidate_bound_comparison.json) 只读取前次固定数据，没有新增候选搜索。

## 2. 逐个消掉多项式矩，会精确变成参数搬运

写 \(t=x^2\)，令 \(w(t)=t^{a-1/2}(1-t)^b\)。若规定

\[
\int_0^1 w(t)t^jH(t)\,dt=0\qquad(0\le j<k),
\]

全部解恰好是

\[
H=w^{-1}\frac{d^k}{dt^k}
\bigl[t^{a+k-1/2}(1-t)^{b+k}P(t)\bigr],
\qquad\deg P\le c-a-b-k.
\]

分部积分后，原积分等于参数 \((a+k,b+k,c+k)\) 的核乘以 \((c)_k\)，两个周期同时乘同一个有理数，所以产生**完全相同的既约分数**。这对整套相消解空间成立，不只是一个例子。

在用户族中，把全部 \(m+1\) 个矩消掉，得到 \((2m+1,3m+1,5m+2)\) 的常乘子核，耗尽了剩余自由度。我们给出了其全规模非零周期证明和包含分母成本的质量上界；该上界没有衰减。这里没有据此证明实际质量必然失败，也没有把这个结果改成主攻旧族排除。

这说明“精确相消”本身必须做等价性审查。完整规则、整数坐标和成本见 [显式周期与账式](constructive-arithmetic/EXPLICIT_PERIOD_EXTRACTION_AND_COST.md)。

## 3. 改为一次消掉完整截断积分，得到明确生成规则

保持 \(a=m,b=2m,c=4m+1\)，取固定线性截断规模

\[
K_m=13m,
\]
\[
S_K(t)=\left(\frac23\right)^c
\sum_{k=0}^{K}(-1)^k\binom{c+k-1}{k}
\left(\frac{2t-1}{3}\right)^k.
\]

定义精确有理矩

\[
\mu_j=\int_0^1x^{2m+2j}(1-x^2)^{2m}S_K(x^2)\,dx,
\quad j=0,1,
\]
\[
\boxed{H_m(t)=\mu_1-\mu_0t.}
\]

它们可以直接用有理数有限求和计算，不输入 \(\pi\)、连分数或目标分数。规则自动保证

\[
\boxed{\int_0^1x^{2m}(1-x^2)^{2m}H_m(x^2)S_K(x^2)\,dx=0.}
\]

这里只有一个完整截断相消条件，不要求每个低阶矩分别为零，故不属于上一节的最大 Rodrigues 特化。这里没有声称它是文献中全新的逼近方法。

两个基础周期有闭式。设 \(J_j=I_m(t^j)=\beta_j+\alpha_j\pi/4\)，则

\[
\alpha_0=\frac{\binom{2m}{m}}{64^m},\quad\beta_0=0,
\qquad
\alpha_1=\frac{6m+1}{2m-1}\alpha_0,
\quad\beta_1=-\frac1{(2m-1)\binom{4m}{2m}}.
\]

因此 \(\Delta=\alpha_0\beta_1-\alpha_1\beta_0\ne0\)。构造周期为

\[
\alpha_m=\mu_1\alpha_0-\mu_0\alpha_1,\qquad
\beta_m=\mu_1\beta_0-\mu_0\beta_1.
\]

证明给出 \(K=13m\) 对每个 \(m\ge1\) 都保证 \(|\alpha_m|\ge|\Delta|/2>0\)，不靠有限数据外推。

更直接地，令 \(L_K\) 为 \(\mu_0,\mu_1\) 分母的最小公倍数，\(u_j=L_K\mu_j\in\mathbb Z\)。构造出的分数就是下列整数比的完全约分：

\[
\boxed{\frac{p_m}{q_m}=
\frac{2^{6m+2}u_0}
{\binom{4m}{2m}\binom{2m}{m}
[(6m+1)u_0-(2m-1)u_1]}.}
\]

这是一条全规模的明确生成规则，没有逐次的 LLL 系数选择。证明见 [完整截断生成器](constructive-review/COMPLETE_TRUNCATION_GENERATOR.md)，基础周期的独立推导见 [闭式周期与截断规模](constructive-arithmetic/CLOSED_BASE_PERIODS_AND_K13.md)。

## 4. 误差和既约分母的联合上界

令

\[
M_j=\int_0^1x^{2m+2j}(1-x^2)^{2m}\,dx
=\frac{2^{2m}(2m)!}{\prod_{r=0}^{2m}(2m+2j+2r+1)},
\]
\[
\varepsilon_K=\left(\frac23\right)^c
\frac{\binom{c+K}{K+1}}{3^{K+1}}
\left(1-\frac{c+K+1}{3(K+2)}\right)^{-1}.
\]

则有包含实际完全约分分母的严格账式

\[
\boxed{Q_m\le\frac{4q_m^2}{|\alpha_m|}\,
\varepsilon_K\bigl(|\mu_1|M_0+|\mu_0|M_1\bigr).}
\]

系数分母也没有略去：

\[
L_K\mid3^{c+K}\operatorname{odd}\operatorname{lcm}
(1,\ldots,2(3m+K+1)+1).
\]

闭式周期还给出更直接的分母上界

\[
\boxed{q_m\le\frac32\binom{2m}{m}L_K.}
\]

结合周期清分母可得到全规模的先验联合上界。它目前不能证明 \(Q_m\to0\)：用现有通用分母估计替换实际分母后，得到的上界增长。这个事实限制的是当前证书的能力，不是对实际质量的全局失败证明。

另一方面，闭式周期和 \(K=13m\) 确实足以证明通常误差收敛：

\[
\boxed{\left|\pi-\frac{p_m}{q_m}\right|
\le\frac{128}{3}(4m^2-1)
\left(\frac{2^{33}}{3^{23}}\right)^m\longrightarrow0.}
\]

由 \(\pi\) 无理，进一步有 \(q_m\to\infty\)，没有断言分母单调。但通常误差趋零仍弱于用户要求的 \(q_m^2\) 归一化质量趋零。

## 5. 当前该保留的构造任务

本轮已完成“明确生成规则＋误差和分母的联合上界”，但**尚未完成其中最关键的净收益：质量上界趋零。** 更高阶逐矩相消会被参数搬运吸收；完整截断规则则暴露了其有理矩分母成本。后续应针对显式整数对 \((u_0,u_1)\) 的约分、截断多项式的整数化，或保留乘子自由度的核改造，寻找真实的分母节省。不能仅以增大 \(K\) 宣称收益，也不能用排除旧素数池来替代这一构造目标。

原来候选的 \(Q\approx0.638066\) 是有符号积分认证的有限表现；本轮没有产生新的 \(\pi\) 逼近纪录，也没有证明 \(\pi\) 的坏逼近性或非坏逼近性。

## 复现和独立审查

全部本轮验证只使用 Python 标准库；没有调用 LLL 或新的大范围分数搜索。入口为：

```bash
python constructive-cancellation/verify_cancellation_identity.py
python constructive-arithmetic/explicit_periods.py
python constructive-review/audit_rodrigues.py
python constructive-review/audit_complete_truncation.py
python compare_existing_bounds.py
```

最后一个入口读取相邻的前次研究目录中的固定候选 JSON。检查摘要见 [CHECKS.json](CHECKS.json)，独立证明审阅见 [INDEPENDENT_CONSTRUCTIVE_AUDIT.md](constructive-review/INDEPENDENT_CONSTRUCTIVE_AUDIT.md)。[来源与依赖](SOURCES.md) 区分直接推导和已有文献中的 Rodrigues 表示。
