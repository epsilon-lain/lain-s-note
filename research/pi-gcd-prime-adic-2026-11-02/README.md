# π 递推逼近的大公因子提升与累计成本

2026 年 11 月 2 日。研究基准为 [epsilon-lain/lain-s-note 的 dff048f7 提交](https://github.com/epsilon-lain/lain-s-note/tree/dff048f73da83e63978e1c01576d38ea71af1546/research/pi-joint-arithmetic-2026-11-02)。本轮沿既定整数递推推进，证明了全指标二进赋值、奇素数的数位判据及大素数平方提升条件。它们识别出真实的大公因子增益，同时量化了增益必须达到的总规模。

本轮没有构造出质量趋零的无限子序列，也没有证明这一族的所有稀疏子序列都失败。下列新推导是相对该仓库检查点的推进，不声明文献首创。独立复核是其他研究代理逐步检查证明，不是外部同行评审。

## 固定对象和成功条件

两列整数分别由

\[
A_0=A_1=1,\qquad B_0=0,\quad B_1=1,
\]

\[
X_{n+1}=(2n+1)X_n+n^2X_{n-1}\quad(n\ge1)
\]

定义。沿用

\[
g_n=\gcd(A_n,4B_n),\quad
p_n=4B_n/g_n,\quad q_n=A_n/g_n,
\quad h_n=\frac{g_n\sqrt n}{n!}\quad(n\ge1).
\]

原检查点已证明

\[
\mathcal Q_n=q_n^2\left|\pi-\frac{p_n}{q_n}\right|
=\frac{\sqrt2}{h_n^2}(1+O(n^{-1})).
\]

因此本族能否有质量趋零的子序列，取决于 \(h_n\) 是否无界。原来的固定间隔限制仍成立；本轮研究的是稀疏指标上 \(g_n\) 可能增长的机制。

## 一个真实的额外公因子机制

用 \(v_p(m)\) 表示素数 \(p\) 在非零整数 \(m\) 的分解中出现的次数。对奇素数 \(p\)，当 \(p<n<2p\) 时写 \(n=p+r\)，\(0<r<p\)。定义整数辅助序列

\[
D_0=0,\quad D_1=2,
\]

\[
D_{r+1}=(2r+1)D_r+r^2D_{r-1}+2A_r+2rA_{r-1}.
\]

本轮证明了完整的三种情形：

\[
\boxed{v_p(g_{p+r})=
\begin{cases}
0,&p\nmid A_r,\\
1,&p\mid A_r,\ A_r/p+D_r\not\equiv0\pmod p,\\
2,&p\mid A_r,\ A_r/p+D_r\equiv0\pmod p.
\end{cases}}
\]

这不是渐近猜测。证明使用移位连续分式多项式 \(C_r(u)\)：\(C_r(0)=A_r\)、\(C_r'(0)=D_r\)，且原递推在第 \(p\) 步的耦合项为 \(p^2\)。因此对两列解都有

\[
X_{p+r}\equiv C_r(p)X_p\pmod{p^2},
\quad C_r(p)\equiv A_r+pD_r\pmod{p^2}.
\]

\(B_p\) 在模 \(p\) 意义下非零；Wronskian 又保证该范围的公因子赋值至多为 2，于是得到上述精确分类。

例如 \(p=73,r=34\) 时，

\[
A_{34}/73\equiv37,\qquad D_{34}\equiv36\pmod{73}.
\]

两项相加恰为零，所以

\[
v_{73}(g_{107})=2>1=v_{73}(107!).
\]

同样，\(p=131,r=35\) 的两个剩余分别为 104、27，也产生额外提升。更早还有 \(v_{17}(g_{72})=5>4=v_{17}(72!)\)。因此，假设 \(g_n\) 的奇数部分总整除 \(n!\) 会错误地丢掉真实机会。

证明与精确验证见 `gcd_blocks/PRIME_BLOCKS_AND_LIFTING.md`，独立复核见 `gcd_audit/INDEPENDENT_GCD_BLOCK_AUDIT.md`。

## 每个指标都存在指数规模的二进亏损

对全部整数 \(k\ge0\)，本轮证明

\[
\boxed{\begin{aligned}
v_2(A_{4k})&=2k,&v_2(A_{4k+1})&=2k,\\
v_2(A_{4k+2})&=2k+2,&
v_2(A_{4k+3})&=2k+3+v_2(k+1).
\end{aligned}}
\]

最后一类的额外相消需要单独处理。证明先将 Legendre 系数倒序写成二进收敛的二项式级数，再建立在 \(\mathbb Q_2\) 中的恒等式

\[
\sum_{\ell\ge0}\frac{\ell!}{(2\ell+1)!!}=0.
\]

该恒等式通过可逐项积分的几何级数及 \(\mathbb Q_2(i)\) 中收敛的对数证明；它不把正实数级数的值认作零。剩下的赋值判断归结为一个有限多项式模 4 恒等式。全证明及独立复核分别在 `gcd_arithmetic/TWO_ADIC_VALUATION.md` 和 `gcd_audit2/INDEPENDENT_TWO_ADIC_AUDIT.md`。

定义奇素数的正超额乘积

\[
E_n=\prod_{p\text{ 为奇素数}}
p^{\max(0,v_p(g_n)-v_p(n!))}.
\]

\(E_n\) 只计算比阶乘多出的部分，不包括阶乘本身已有的约分。由 \(g_n\mid A_n\) 和上述精确公式，严格得到

\[
\boxed{h_n\le2(n+1)^2\sqrt n\,2^{-n/2}E_n.}
\]

所以，任何成功子序列都必须满足

\[
\boxed{\ln E_n\ge\frac{\ln2}{2}n-O(\ln n).}
\]

单个额外素因子的增益可以是真的，却仍远小于所需的累计补偿。若能证明 \(E_n=\exp(o(n))\)，就能排除本族全部成功子序列；该估计目前尚未建立。反过来，达到这个必要规模也不自动充分，因为奇素数中的负亏损尚需一并抵消。

## 一般奇素数的赋值能用数位筛选

令 \(a_n=A_n/n!\)、\(b_n=B_n/n!\)。精确生成函数给出

\[
\sum a_nz^n=(1-2z-z^2)^{-1/2},\quad
b_n=\sum_{j=1}^n\frac{a_{n-j}a_{j-1}}j.
\]

\(a_n\in\mathbb Z[1/2]\)，且对奇素数 \(p\) 有 Lucas 型关系

\[
a_{mp+r}\equiv a_ma_r\pmod p\quad(0\le r<p).
\]

若 \(n=\sum_{j=0}^e d_jp^j\)，\(d_e\ne0\)，则

\[
\boxed{p^eb_n\equiv
(-1)^{e(p-1)/2}b_{d_e}\prod_{j<e}a_{d_j}\pmod p.}
\]

右侧非零时，准确得到

\[
v_p(g_n)=v_p(n!)-e.
\]

这既确定常规亏损，也给出大公因子候选必须满足的异常数位配置。右侧为零只表示需要继续提升，不能直接视为正超额。

另外，对任意奇素数和 \(n\ge2\)，归一化 Wronskian 还给出统一上限

\[
v_p(g_n)-v_p(n!)
\le\lfloor\log_p(n-1)\rfloor-v_p(n).
\]

令 \(L_n=\operatorname{lcm}(1,\ldots,n-1)\)，则

\[
\boxed{E_n\mid\operatorname{oddpart}\!\left(
\frac{L_n}{\gcd(L_n,n)}\right).}
\]

这里 \(\operatorname{oddpart}\) 表示除去全部因子 2。特别地，额外素数都小于 \(n\)，超过 \(\sqrt{n-1}\) 的素数在 \(E_n\) 中至多出现一次。该最小公倍数界仍允许指数规模的 \(E_n\)，尚不能与前述必要条件形成矛盾。

完整证明见 `gcd_arithmetic/ODD_PRIME_DIGIT_CRITERION.md`。

## 精确计算的范围

一般定理由上述证明建立；有限计算用于核对恒等式、分支与整数常数。

| 检查 | 范围或结果 |
|---|---|
| 精确二进赋值公式 | \(1\le n\le10000\)，全部通过 |
| 一般素数数位同余 | 35614 项精确核对，通过 |
| 素数平方提升分类 | 4180 项精确核对，通过 |
| 完整 gcd 与累计成本 | \(1\le n\le5000\)，精确检查通过 |
| 上述完整范围内 \(E_n>1\) 的指标 | 41 个 |
| 上述完整范围内最大的 \(E_n\) | \(E_{4762}=2677\) |

对 \(3\le n\le5000\)，还用整数算术认证

\[
4(n!)^2>(2n+1)g_n^2.
\]

结合原检查点的严格 Legendre 误差下界，这意味着这些指标全部满足 \(\mathcal Q_n>1/2\)。因此它们没有改善此前 \(n=2\) 给出的分数 \(3/1\)。这个结论是有限区间认证，不推断所有更大指标。

主计算入口为 `gcd_research/verify_global_budget.py`，其 JSON 包含全部 41 个超额指标。浮点对数只用于诊断展示，任何整除、赋值和不等式认证均使用精确整数或有理数。

还进行了针对组合提升的精确搜索：覆盖 \(n\le50000\) 与全部奇素数 \(p\le500\)。先用数位零点筛掉不可能产生正超额的素数，再模素数高幂检查其余 43 个素数，共发现 115 次单素数超额，没有同指标双素数超额。对这 115 个候选重新计算完整整数与 gcd，也未发现其他超额素数。因此在此有限范围内，没有找到至少包含一个 \(p\le500\) 的组合提升；两个素数都大于 500 的情形不在这项排除范围内。脚本、完整候选及精确结果见 `gcd_combined/`。

## 下一步的具体构造目标

本轮找到的额外提升机制是真实的，下一步需要在同一指标上累计它们。可优先开展两项工作：

1. 结合数位判据与移位连续分式的模素数幂条件，构造同时满足多个提升条件的指标。必须核算全部奇素数亏损及二进亏损；单个条件的提升次数不构成最终质量证书。
2. 研究 \(\ln E_n\) 的统一上界或成功子序列下界。上界若低于 \((\ln2)n/2\) 一个固定正比例，就能排除本族；构造方向则需要达到该指数阈值并进一步控制未计入 \(E_n\) 的亏损。

这两项分别对应严格排除与直接构造。证明本族失败仍不能推出 π 坏逼近，因为本族没有覆盖所有既约有理数。

## 文献对应与适用边界

[Ben David 等，2024，On the Connection Between Irrationality Measures and Polynomial Continued Fractions](https://link.springer.com/article/10.1007/s40598-024-00250-z) 研究的 factorial reduction 与本问题直接相关；其 Appendix D Theorem 5 可用于给出阶乘和最小公倍数形式的约分下界。但 factorial reduction 允许 \(g_n/n!\) 指数衰减，不能推出这里需要的 \(g_n\sqrt n/n!\) 无界。

[2026 年的 Legendre 赋值论文](https://link.springer.com/article/10.1007/s40993-026-00743-0) 提供数位与赋值工具，其关于特定有理输入的结论不能直接移植到 \(P_n(i)\)。相关原始文献、索引差异与应用条件见 `gcd_literature/PRIME_RESET_AND_FACTORIAL_REDUCTION.md`。

## 复现

本轮主要程序仅需 Python 标准库：

```bash
python gcd_arithmetic/verify_gcd_arithmetic.py
python gcd_blocks/verify_prime_blocks.py
python gcd_literature/verify_prime_reset.py
python gcd_research/verify_global_budget.py 5000
python gcd_combined/scan_combined_surplus.py
python gcd_combined/certify_known_candidates.py
```

独立复核中的有限多项式展开另外使用 SymPy。完整的收敛、赋值与无限指标结论均保留在证明文件中，没有用扫描结果替代。
