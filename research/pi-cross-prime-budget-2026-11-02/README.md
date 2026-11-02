# 多素数区间的亏损证书与共同整数障碍

2026-11-02。继续研究固定 Legendre–反正切 Padé 递推族。**本轮没有证明最终的 \(S(n)\le n/6\)，也没有排除本族全部稀疏成功子序列。** 得到了跨多个素数区间的精确整除证书、一个更强的有限质量下界，以及一类相邻项消元方案的严格限制。

前一检查点：[全局第一层异常](../pi-global-anomaly-2026-11-02/README.md)。本目录中的命题由下列证明给出，有限程序用于审计，不能替代统一证明。

## 1. 固定研究对象

\[
A_0=A_1=1,\quad B_0=0,\quad B_1=1,\qquad
X_{n+1}=(2n+1)X_n+n^2X_{n-1}.
\]

令

\[
g_n=\gcd(A_n,4B_n),\quad
p_n=4B_n/g_n,\quad q_n=A_n/g_n,\quad
Q_n=q_n^2|\pi-p_n/q_n|.
\]

此前已证明 \(q_n\to\infty\)、\(Q_n\sim(\sqrt2/n)(n!/g_n)^2\)。把奇素数的正超额与亏损分别记为 \(E_n,D_n\)，并记 \(J_n=\ln(E_n/D_n)\)，则

\[
\ln Q_n=n\ln2-2J_n+O(\ln n).
\]

目标仍是控制实际约分后的质量，而非只看未约分误差。

## 2. 本轮主要结果：保留几乎全部大素数的常规亏损

设

\[
U_n=2^n A_n/n!\in\mathbb Z,\qquad C_n=\operatorname{odd}(U_n).
\]

给定区间层数 \(B\ge1\)，定义

\[
t=\max\{\sqrt n,n/(B+1)\},\qquad
H=\prod_{t<p\le n,\ p\text{ 奇素数}}p,
\]
\[
P=\gcd(C_n,H),\qquad
F=\gcd\left(H,\prod_{b=1}^{B}B_b\right),\qquad
M=H/\gcd(H,PF).
\]

严格成立

\[
\boxed{M\mid D_n,\qquad H\mid D_nPF.}
\]

原因是：对 \(p\mid H\)，写 \(n=bp+r\)。若 \(p\nmid C_n\) 且 \(p\nmid B_b\)，Lucas 关系与主极点公式强制 \(e_p(n)=-1\)。因此除去共同整数 \(C_n\) 的素因子和一个短前缀的素因子，剩下的素数都给出真实亏损。这里没有假设不同素数随机或独立。

取 \(B=\lfloor n^{1/3}\rfloor\)，则 \(\ln F=o(n)\)、\(\ln H=n+o(n)\)，而低于 \(t\) 的全部可能正收益也是 \(o(n)\)。设 \(W(n)=\ln P\)，得到

\[
\boxed{J_n\le W(n)+\ln E_n-n+o(n)\le2W(n)-n+o(n),}
\]
\[
\boxed{\ln Q_n\ge(2+\ln2)n-4W(n)+o(n).}
\]

因此任何使 \(Q_n\) 有界的无穷指标子序列，都必须满足

\[
\boxed{\ln\gcd(C_n,H)\ge
\left(\frac12+\frac{\ln2}{4}-o(1)\right)n
=(0.6732867951\ldots-o(1))n.}
\]

这要求至少 \((0.6732867951\ldots-o(1))n/\ln n\) 个不同的选定素数同时整除同一个 \(C_n\)。\(H\) 的区间下端约为 \(n^{2/3}\)，上端为 \(n\)。首层异常不等于实际正超额；保留 \(E_n\) 的第一条不等式更精细。

**未解决之处很明确：尚未证明这个特定素因子乘积低于上述阈值。** 已知 \(\ln C_n=\ln(2+\sqrt2)n+O(\ln n)\approx1.227947n\)，单凭高度不足以矛盾。素数定理只能计算候选素数的总权重，不能限制其中哪些整除 \(C_n\)。

完整证明：[MULTI_BAND_DEFICIT_CERTIFICATE.md](multi-band-budget/MULTI_BAND_DEFICIT_CERTIFICATE.md)。

## 3. 已经可用的有限改进

若 \(G\) 是已经证明的 \(E_n\) 的倍数，例如 \(G=\gcd(C_n,\operatorname{odd}(\operatorname{lcm}(1,\ldots,n-1)))\)，记 \(s_n=v_2(U_n)\)，则

\[
\boxed{Q_n>\frac{2\,2^{2(n-s_n)}M^2}{(2n+1)G^2}.}
\]

同样可使用旧筛法中更小的 \(\widehat G_n\)。相对于原来只计算正收益的下界，新证书乘上精确的 \(M^2\ge1\)，所以这一有限改进无需等待新的全局猜想。它只需计算 \(U_n\) 和短前缀 \(B_1,\ldots,B_B\)，无需计算完整 \(B_n\)。

例如 \(n=21,B=2\) 时，\(H=46189,P=17,F=1,M=2717,\widehat G=17,s_n=13\)。旧下界为 \(131072/12427\)，新下界为 \(967585169408/12427\)，提高 \(7382089\) 倍。旧下界在这个例子中已经大于一；它用于展示改进因子，不冒充新排除指标。

主验证覆盖 \(2\le n\le1000\) 的 5715 个区间选择和 355049 次常规赋值检查。取 \(B=\lfloor n^{1/3}\rfloor\)，新证书证明这一有限范围内每个 \(n\ge4\) 都有 \(Q_n>1\)。这没有延伸到所有 \(n\)。独立程序另外从实际整数 gcd 复核至 \(n=500\)。

## 4. 哪种相邻项消元不能降低成本

把 \(U_{n-1},U_n,\ldots,U_{n+h}\) 视为由两个自由初值 \(y,x\) 生成的形式递推窗口。若整数多项式系数的线性组合在 \(x=0\) 时对任意 \(y\) 恒等为零，则这个组合必为

\[
G(n)U_n,\qquad G\in\mathbb Z[n].
\]

这是任意有限窗口长度的精确结论。这样的非零线性组合保留全部中央奇数因子 \(C_n\)，无法凭绝对值获得更小的共同整数。

对窗口变量总次数为 \(d\) 的多项式组合，清分母后同样保留中央因子；当 \(F_n\ne0\) 时，定量下界为

\[
\ln\operatorname{odd}(|F_n|)\ge
\ln(2+\sqrt2)n-dh\ln(n+h)-O(\ln n).
\]

所以 \(dh\ln(n+h)=o(n)\) 的这类构造也没有提供目标所需的小高度。**限定条件是仅从中央项为零、对任意邻项成立的形式恒等式；这不是所有消元、所有大窗口或利用实际初值的算术关系的失败定理。**

证明：[COMMON_WINDOW_WITNESS_BARRIER.md](cross-prime-product/COMMON_WINDOW_WITNESS_BARRIER.md)。

## 5. 辅助审计与接下来真正需要的输入

此前的边缘高度界有效，但仅剔除收缩边缘不需要那么多算术：简单计数给出 \(S_{\rm edge}(n,K)\le2(K+1)\ln n\)，普通素数定理进一步给出任意 \(K=o(n)\) 时 \(S_{\rm edge}=o(n)\)。这项简化不解决内部的跨素数集中。

精确中点可解释为固定椭圆曲线 \(y^2=x(x^2+2x+2)\) 的 Frobenius 迹为零，得到必要条件 \(p\equiv7\pmod8\)。对固定 \(n\)，中点最多对应一个素数，所以它也没有解决全局问题。观察到的更强余数类 \(31\pmod{32}\) 未证明。

详见 [EDGE_CAPACITY_AND_MIDPOINT.md](midband-arithmetic/EDGE_CAPACITY_AND_MIDPOINT.md)。

后续工作的验收条件应是下列之一，而非仅扩大有限扫描：

1. 证明 \(\limsup W(n)/n<1/2+\ln2/4\)，或者保留真实提升层，证明 \(W(n)+\ln E_n\le(1+\ln2/2-\eta)n\) 对所有充分大 \(n\) 成立，其中 \(\eta>0\)。
2. 找出并严格验证与这些界相冲突的实际递推指标结构，不能只给抽象零点模型。
3. 使用实际初值提供额外算术约束，使共同整数能够去掉已证明无关的素因子，而不只重新生成 \(U_n\) 的倍数。

即使排除整个固定族，也不推出 \(\pi\) 坏逼近。本轮没有给出 \(\pi\) 的坏逼近性、非坏逼近性或新的有理逼近纪录。

## 复现与审查

四个入口均只用 Python 标准库：

```bash
python multi-band-budget/verify_multi_band_budget.py
python cross-prime-product/verify_common_window_witness.py
python midband-arithmetic/verify_edge_midpoint.py
python crossprime-review/verify_crossprime_independent.py
```

[独立审查](crossprime-review/INDEPENDENT_CROSSPRIME_REVIEW.md)、[来源及旧命题依赖](SOURCES.md)、[检查命令与文件摘要](CHECKS.json)。
