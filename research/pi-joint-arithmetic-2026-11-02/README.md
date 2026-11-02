# 从近根几何转到误差与约分后分母的联合控制

2026-11-02。本轮读取并核查用户提供的 `pi_next_work_20261102.zip`，随后沿固定参数的反正切 Padé 积分构造实际整数线性形式。仓库核查基准是 `befa7cf7a515a18c3bd00ef17a4c6f53dd9ae6df`。下列新推导是相对本项目检查点的推进，不声明文献首创或外部同行评审。

**主要结果：** 对一个明确的正积分／Legendre–反正切族，误差与真正约分后的分母满足一个精确渐近公式。该族能否产生质量趋零的子序列，等价于一个完全由整数递推定义的最大公因子增长问题。我们没有证明这个增长条件发生，也没有决定 π 的坏逼近性。

## 1. 接收材料中的结论，经重新核查成立

输入归档 SHA-256 为 `45f0f7c529c2002b46214b6158821bdded09995771690f1de5ba47e7efddd5f5`。原文件及其清单完整保存在 [imported-next-work](imported-next-work/PI_NEXT_WORK_20261102.md)，没有改写输入。

- 对原双频率辅助多项式，足够大的 n 下，最近于 \(\omega=2\pi i\) 的根唯一、简单、位于虚轴。全阶混合端点系数界补上了此前缺少的高阶余项控制，支持原文的 Newton 相对误差 \(O(n^{-2})\)。
- 对偶数 n 的两个完整反射因子，最大系数高度满足 \(\log H=(n^2/4)\log n+O(n^2)\)。独立整数行列式及有限域程序重新认证 n=4、6、8、10、12 的十个因子及其有理不可约性。这里没有推断一般 n 的不可约性。
- 指数 Padé 对照中的完整分母与高度成本计算成立。一个需补充的范围是其式 (13) 的直接推导要求 b≥2，即既约 p/q 的 q≥3；q=1、2 应另行处理。[审查记录](reviews/INDEPENDENT_HEIGHT_PADE_REVIEW.md)明确保留了这一限定。

审查还发现原全阶递推可以多保留一个因子：

\[
B_{n,2}=\frac{K_nK_{n-1}}{2n-2}.
\]

它把缩放圆盘内的高阶余项加强至 \(O(n^{-3})\)，从而得到单独证明的推论

\[
\boxed{\zeta_n-\omega
=-\frac{P_n(\omega)}{P_n'(\omega)}
\bigl(1+O(n^{-3})\bigr).}
\]

该加强及其假设见 [Newton 审查第 4 节](reviews/INDEPENDENT_NEWTON_REVIEW.md)。这仍是同一辅助多项式的解析定位，没有降低其算术成本。

## 2. 明确输出整数线性形式

令 n 为非负整数，定义正积分

\[
I_n=\int_0^1\frac{x^n(1-x)^n}{(1+x^2)^{n+1}}\,dx>0.
\]

定义两列整数 \(A_n,B_n\)，初值为

\[
A_0=A_1=1,\qquad B_0=0,\quad B_1=1,
\]

并让每一列分别满足

\[
X_{n+1}=(2n+1)X_n+n^2X_{n-1},\qquad n\ge1.
\]

取最大公因子及既约分数

\[
g_n=\gcd(A_n,4B_n),\qquad
p_n=4B_n/g_n,\quad q_n=A_n/g_n>0.
\]

则 \(\gcd(p_n,q_n)=1\)，而且有精确恒等式

\[
\boxed{q_n\pi-p_n
=\frac{4(-2)^n n!}{g_n}I_n.}
\]

这给出的就是以 π 为一次项的整数线性形式。正积分保证它不为零，并固定误差的交替符号。分数 \(p_n/q_n\) 收敛到 π，分母趋于无穷；这不代表线性形式本身或归一化质量已经趋于零。

一个重要的结构核对是：此正积分族恰好等于固定参数 a=1 的 Legendre–反正切 Padé 族。它们是同一组分数的两种表示，不能计作两种独立的逼近机制。

## 3. 把误差与约分收益放进同一个公式

令

\[
\mathcal Q_n=q_n^2\left|\pi-\frac{p_n}{q_n}\right|.
\]

对上面的确定族，解析积分的 Laplace 估计与精确约分恒等式共同给出

\[
\boxed{\mathcal Q_n=
\frac{\sqrt2}{n}\left(\frac{n!}{g_n}\right)^2
\bigl(1+O(n^{-1})\bigr).}
\]

这里已经使用真正的既约分母，未将一个方便的公共分母当成 q。其解析常数来自两个互补的驻点：令 \(\lambda=1+\sqrt2\)、\(\rho=(\sqrt2-1)/2\)，则 \(2\rho\lambda=1\)。积分带来的 \(\rho^n\) 与 Legendre 系数的 \(\lambda^n\) 在质量表达式中精确配合，最后留下 \(\sqrt2/n\) 和实际公因子。

因此，在本族内部，对任何 \(n_j\to\infty\)，

\[
\boxed{\mathcal Q_{n_j}\to0
\quad\Longleftrightarrow\quad
\frac{n_j!}{g_{n_j}\sqrt{n_j}}\to0.}
\]

等价地，令 \(h_n=g_n\sqrt n/n!\)（n≥1），则本族存在质量趋零的子序列，当且仅当 \(\limsup h_n=\infty\)。这个待证明条件只包含整数递推、最大公因子和阶乘。若能证明它，就能证明 π 非坏逼近；π 非坏逼近本身并不保证这一个指定族能够实现它。

完整推导见 [正积分与 gcd 临界条件](joint-arctan-integral/DIRECT_INTEGRAL_GCD_THRESHOLD.md)，以及[独立联合算术审查](reviews/INDEPENDENT_JOINT_ARITHMETIC_REVIEW.md)。

## 4. 同时证明了一条限制：好指标必须越来越稀疏

整数递推的 Wronskian 精确等于

\[
A_nB_{n-1}-A_{n-1}B_n=(-1)^n((n-1)!)^2.
\]

因此

\[
\boxed{g_ng_{n-1}\mid4((n-1)!)^2.}
\]

这个整除关系已经约束相邻的约分收益。另从 Legendre 正误差公式和相邻既约分数的非零整数交叉差，直接得到全部 n≥1 的界

\[
\mathcal Q_n\mathcal Q_{n-1}>
\frac{n^2}{4(4n^2-1)}>\frac1{16}.
\]

利用精确误差渐近，还能控制任意固定指标差 k≥1：

\[
\liminf_{n\to\infty}\mathcal Q_n\mathcal Q_{n-k}
\ge
\frac{\lambda^{-2k}}
{\bigl(1-(-1)^k\lambda^{-2k}\bigr)^2}>0.
\]

k=1、2 的右端分别为 1/8、1/32。如果某子序列 \(\mathcal Q_{n_j}\to0\)，它必然满足 \(n_{j+1}-n_j\to\infty\)。这些结论排除了逐项使用完整序列、或使用任何最终有界间隔子序列达到质量趋零；它们仍允许某个非常稀疏的成功子序列。

一般固定有理参数 a>0 的误差公式、全 n 的一阶及二阶邻近约束，见 [Legendre 联合成本](joint-log-pade/LEGENDRE_ARCTAN_JOINT_COST.md)。这些约束针对指定构造，不能用来断言所有有理数都满足相同下界。

## 5. 实际检查到什么

有限计算对每个候选执行约分，并用严格有理区间认证误差。没有使用浮点 π 来判定好坏。

| 构造与检查范围 | 结果 |
|---|---|
| 本正积分族，1≤n≤100 | 最小质量在 n=2，分数为 3/1，质量处于 (0.14159265358, 0.14159265359)；这只是这一有限族的最小值 |
| 同族 n=100 | 既约分母有 90 位，质量约为 \(1.323\times10^{103}\)；绝对误差缩小没有抵消分母平方成本 |
| 同指标 Machin 组合，1≤n≤100 | 未产生改进；n=100 的既约分母有 426 位，质量约为 \(5.988\times10^{650}\) |

Machin 组合是在求和后重新约分的，未假设两个分母一定互素。上面的 Legendre 邻近定理没有被移用于这个组合。有限结果不用于推断未来所有 n，也不构成新的有理逼近纪录。

## 6. 下一步的具体算术任务

本轮把一条路线收敛到了可以单独进攻的目标：研究上述确定递推的 \(g_n\)，特别是其素数赋值和稀疏指标上的大公因子机制，尝试证明或排除 \(g_n\sqrt n/n!\) 无界。若改用共同分母的多点／多重 Padé，也应保留同样的验收方式：输出实际 p、q，完整约分，并控制 \(q|q\pi-p|\)。

本族若最终失败，只能排除本族；证明 π 坏逼近仍要求覆盖全部既约有理数。本轮没有完成该全局问题，也没有完成本族的无限子序列问题。

## 证明依赖、来源和复现

原近根推导依赖此前的[周期非零下界](../pi-period-mechanism-2026-11-02/raw-root-mechanism/HANKEL_PERIOD_MAGNITUDE.md)和[实际端点迹渐近](../pi-endpoint-sensitivity-2026-11-02/period-derivative-scale/PERIOD_LOG_DERIVATIVE.md)。本轮联合算术推导是独立的一维整数／积分论证，不依赖原稿的高维插值定理或 μ(π)=2 的宣称。

来源核对及输入小范围修订说明见 [SOURCE_AND_SCOPE.md](SOURCE_AND_SCOPE.md)。本轮新增的复现入口仅使用 Python 标准库：

```bash
python -B reviews/verify_newton_exact.py
python -B reviews/finite-height/independent_finite_certificates.py
python -B reviews/independent_pade_checks.py
python -B joint-log-pade/verify_legendre_arctan_joint.py
python -B joint-arctan-integral/verify_direct_integral.py
python -B reviews/independent_joint_arithmetic_checks.py
```

原压缩包中的高精度根图例及符号程序仍按原样保留，部分需要 mpmath 或 SymPy；它们不承担本轮解析渐近结论的证明。独立复核指另一研究代理重新检查推导，不是外部同行评审。完整文件校验和与实际运行结果见 [CHECKS.json](CHECKS.json)。
