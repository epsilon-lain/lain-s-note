# 从素数提升到实际分母的误差指数

2026-11-02。接续[素数赋值与平方提升](../pi-gcd-prime-adic-2026-11-02/README.md)，本轮把研究目标统一写成

\[
\varepsilon_n=|\pi-p_n/q_n|\le C/q_n^c,
\qquad C>0,\quad c>0,
\]

其中 p、q 始终完全约分，C 是与 n 无关的常数；取 C=e 就是自然常数 e/q^c，取任意小的固定 C 则是另一个误差预算。这两种读法都适用。得到一个全指标二进定理、精确的指数换算，以及成功指标必须满足的大素数联合提升密度。尚未构造出质量趋零的无限子序列。

## 1. 原先剩下的二进猜测可以证明

仍令 A、B 满足 X_(n+1)=(2n+1)X_n+n²X_(n−1)，初值 A_0=A_1=1、B_0=0、B_1=1，并置

\[
g_n=\gcd(A_n,4B_n),\quad p_n=4B_n/g_n,
\quad q_n=A_n/g_n,\quad R_n=g_n/n!.
\]

本轮证明

\[
\boxed{v_2(g_n)=v_2(A_n)\ (n\ne3),\qquad
v_2(g_3)=v_2(A_3)-1.}
\]

所以除 q_3=6 外，所有约分后的 q_n 都是奇数。证明从 r_4=160/51 出发，利用相邻差 r_n−r_(n−1) 在 n≥5 时具有非负二进赋值，不需要猜测 B_n 的赋值。[完整证明](exponent-audit/TWO_ADIC_GCD_SATURATION.md)。

## 2. 可以精确换算为 C/q^c

记 λ=1+√2、a=ln λ、b=(ln2)/2。由此前已经证明的正积分误差和渐近式，目标等价于精确阈值

\[
R_n\ge T_n(C,c):=\frac{A_n}{n!}(\varepsilon_n/C)^{1/c},
\]

而该阈值满足

\[
T_n(C,c)\sim
\frac{\lambda^{1/2-1/c}(4\pi/C)^{1/c}}{\sqrt{2\pi\sqrt2}}
n^{-1/2}\lambda^{(1-2/c)n}.
\]

令 E_n 为超过 n! 的正奇素数乘积，D_n 为约分后仍然缺少的奇素数乘积，即

\[
R_n=2^{d_n}E_n/D_n,
\quad d_n=v_2(g_n)-v_2(n!),
\quad J_n=\ln(E_n/D_n).
\]

现在 d_n=−n/2+O(ln n) 已完全确定。若沿所选子序列 J_n/n→j<a+b，则这族分数的实际误差指数满足

\[
\boxed{\frac{-\ln\varepsilon_n}{\ln q_n}
\longrightarrow\frac{2a}{a+b-j}.}
\]

无须假设增长率存在的 limsup／liminf 版本、所有 n 的必要和充分不等式，见[指数阈值](exponent-error-map/REDUCED_DENOMINATOR_EXPONENT_MAP.md)与[净素数预算](exponent-prime-budget/NET_PRIME_BUDGET_AND_ERROR_EXPONENT.md)。这里只看 E_n 不够；D_n 的亏损必须同时扣除。

此外，D_n 整除 odd(lcm(1,…,n))。结合精确二进赋值及经典素数定理，得到一个较弱但实际成立的保证：对任何固定 C>0 及

\[
0<c<\frac{2a}{a+b+1}=0.7911979206\ldots,
\]

均有 ε_n q_n^c→0，因此 ε_n≤C/q_n^c 最终成立；例如 ε_n≤e/q_n^(3/4) 最终成立。这是本族的保守保证，不是最优指数，也不是新的有理逼近纪录。等号处没有结论。

## 3. 成功需要大量大素数在同一个指标提升

定义

\[
\Theta_c=a+b-2a/c.
\]

如果 Θ_c>0，在任何满足 ε_n≤C/q_n^c 的无穷成功子序列上，令 N_large(n) 为满足 p>√(n−1) 且 v_p(g_n)>v_p(n!) 的不同奇素数个数，则必须有

\[
\boxed{N_{\rm large}(n)\ge
(\Theta_c-o(1))\,n/\ln n.}
\]

原因是每个素数的正超额最多贡献 ln(n−1)，小素数的总容量只有 O(√n ln n)=o(n)；大素数每个又最多增加一个因子。因而目标 c=2 要求至少约 0.3466 n/ln n 个大素数在同一 n 同时提升。单个提升例子或有限个固定素数的组合不能达到这个规模。

一个无需计算 A_50000、B_50000 的精确整数证书给出：若 n=50000 满足 ε_n≤e/q_n²，则至少需要 **1601 个正超额奇素数，其中至少 1561 个大于 √49999**。这是必要数量，不是断言该指标实际存在这些素数；计算使用 e<3，所以还留有余量。[证明与证书](exponent-prime-budget/NET_PRIME_BUDGET_AND_ERROR_EXPONENT.md)。

## 临界处仍然缺什么

在 c=2，写 h_n=√n g_n/n!，则

\[
q_n^2\varepsilon_n\sim\sqrt2/h_n^2.
\]

固定的 e/q² 只要求一个固定质量门槛；要证明 π 非坏逼近，仍须找到 q_n→∞ 且 h_n→∞ 的子序列。仅仅让指数趋近 2，或刚好补齐线性的奇素数亏损，都没有完成这一步。

下一步可直接研究“大素数同时提升的总质量”，证明它在所有大指标不足以支付 Θ_c n，或构造达到该密度且保留净收益 E_n/D_n 的无限指标族。目前两者均未得到；即使排除本族，也不能推出 π 坏逼近。

## 审核、出处与复现

- [独立指数审核](exponent-audit/INDEPENDENT_EXPONENT_REVIEW.md)；[输入平方提升审核](exponent-audit/INDEPENDENT_PRIME_LIFT_REVIEW.md)。这里的独立审核是其他研究代理重推与核算，不是外部同行评审。
- [来源和已有框架对应](SOURCES.md)。增长率换算为逼近指数的通用框架已有文献；本轮不把这一框架或经典 Padé 族称作新发现。
- [运行记录及文件校验和](CHECKS.json)。精确有限检查不承担任何无穷外推。

从本目录运行，均仅用 Python 标准库：

```bash
python -B exponent-error-map/verify_exponent_error_map.py
python -B exponent-prime-budget/verify_net_prime_budget.py
python -B exponent-audit/verify_prime_lift_and_saturation.py
python -B exponent-audit/verify_independent_exponent.py
```

原压缩包的 20 项清单核验可另向第三个脚本提供 `--source ORIGINAL_EXTRACTED_DIRECTORY`。已有的 n≤50000 部分素数扫描本轮没有重跑，覆盖范围仍以原报告为准。
