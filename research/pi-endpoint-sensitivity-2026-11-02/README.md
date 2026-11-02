# 周期附近的真实敏感性与最近复零点

2026-11-02。这一轮继续研究同一族一维双频率原始列行列式，证明它在真实周期附近的陡峭变化，并追踪最近复零点的距离。所有渐近结论都有解析证明；有限精确计算用于核对恒等式、符号及常数，不用于外推定理。这里没有得到 π 的坏逼近性、非坏逼近性或新的有理逼近纪录。

## 这次推进到了哪里

沿用

\[
P_n(Z)=\det\bigl([z^s]e^{hz}(jZ+z)^a\bigr),
\quad j,h\in\{0,1\},\quad 0\le s,a<n,
\qquad \omega=2\pi i.
\]

它使用真实原始列；行按中心 j、再按 s 排序，列按 h、再按 a 排序。n 是构造的规模，**并不是候选有理数的分母 q**。这个稀疏一维族仍不覆盖原稿关键的高维参数范围。

前一轮已经证明归一化端点扰动的矩阵范数很大，但这没有排除行列式内部继续相消。本轮直接计算标量对数导数，得到

\[
\boxed{
\log\left|\frac{P_n'(\omega)}{P_n(\omega)}\right|
=2n\log n+\bigl[\log(16/\pi^2)-2\bigr]n+O(\log n).
}
\]

因此该阶乘尺度确实保留在行列式本身。它不只是某一个矩阵表示给出的宽松上界。

令 \(A_n=\prod_{a=0}^{n-1}a!\)、\(U_n=(2\pi)^{n^2}A_n^2/A_{2n}\)。此前已证明 \(U_n/6<|P_n(\omega)|\le U_n\)。本轮进一步得到以下三个距离尺度。

| 对象 | 精确定义 | 本轮证明的渐近尺度 |
|---|---|---|
| 整个圆盘上的小值范围 | \(\mathcal R_n=\sup\{R: |P_n(Z)|\le2U_n\text{ 对所有 }|Z-\omega|\le R\}\) | \(\log\mathcal R_n=-2n\log n+O(n)\) |
| 实际多项式最近的复零点 | \(\delta_{Z,n}=\min_{P_n(\zeta)=0}|\zeta-\omega|\) | \(\log\delta_{Z,n}=-2n\log n+O(n)\) |
| 固定 Z 的独立端点零点 | \(\delta_{Y,n}=\min_{G_n(\omega,Y)=0}|Y-1|\) | \(\log\delta_{Y,n}=-2n\log n+O(n)\) |

这里的 \(G_n(Z,Y)=\det(YET(Z)-T(Z)E)\) 保留独立端点符号 Y，满足 \(P_n(Z)=G_n(Z,1)\)。这三个对象不同；相同的对数主阶不能把它们直接认作同一条轨迹。

## 为什么会出现这么尖的变化

在矩积分表示里，对端点变量求一次导数，相当于把一个区间积分换成沿复射线 \(x=1+iv/(2\pi)\) 取值的线性泛函。其余 n−1 个变量贡献一个 2n−2 次的乘积。把这部分提出后，积分中的权重变成

\[
e^{-v}v^{2n-2},
\]

即归一化后的 Gamma 分布，主要质量位于 v 约等于 2n 的区域。它带来精确的阶乘因子

\[
K_n=
\frac{[(2n-2)!]^2(2n-1)!}
{[(n-1)!]^4(2\pi)^{2n-1}}.
\]

关键还在于相位：我们证明剩下的复数乘积在 L¹ 中趋于 −1，连靠近 v=0 的无界尾部也有单独控制。因此振荡相消无法抹去这一因子。真正的端点响应满足

\[
\operatorname{tr}(H_n(\omega)^{-1}V_n(\omega))
=i(-1)^nK_n
\frac{\chi_{n-1}(2\pi)}{\chi_n(2\pi)}
\bigl(1+O(n^{-1/2})\bigr),
\]

其中此前已经证明 \(1/6<\chi_k(2\pi)\le1\)（k≥2）。整个论证不需要中心极限定理。

这为 lain 所说的“内凹”提供了一个可检验的局部机制：在这个辅助多项式模型里，精确端点消去与快速增长的扰动响应同时存在。它在精确复数运算中也成立，不能解释成舍入误差造成的假象；这里没有假设实数之间存在物理作用力。

## 从陡峭到零点，需要补哪一步

仅有大导数还不够。本族多项式的次数满足 \(D_n\le n(n-1)\)，所以根分解给出

\[
\left|\frac{P_n'(\omega)}{P_n(\omega)}\right|
\le\frac{D_n}{\delta_{Z,n}}.
\]

另一方面，已有显式圆盘半径

\[
r_n=\frac1{16\cdot3^8n^4\,2^{8n}(2n-2)!}
\]

及中心非零下界。将多项式缩放到单位圆盘，使用一个零点的 Blaschke 因子，可得 \(\delta_{Z,n}>r_n/12\)。这两端共同证明最近复零点的对数主阶。证明没有断言最近根一定在虚轴上，也没有断言唯一、简单或 Newton 迭代收敛。

## 端点结构还有哪些约束

端点多项式满足倒数对称

\[
G_n(-Z,Y)=(-1)^nY^nG_n(Z,Y^{-1}).
\]

因此固定虚数 Z 时，Y 根按 \(Y\leftrightarrow1/\overline Y\) 配对。相对端点算子 \(M_n=H_n(\omega)^{-1}V_n(\omega)\) 满足

\[
\det M_n=P_n(\omega)^{-1},
\qquad \Re\operatorname{tr}M_n=n/2.
\]

实部被固定，并不妨碍虚部具有上面的巨大尺度。Hilbert Gram 归一化后，端点矩阵 \(\mathsf B_n\) 的全部奇异值还满足

\[
\log\sigma_j(\mathsf B_n)=2(n-j)\log n+O(n),
\quad 1\le j\le n,
\]

其中误差对 j 一致。对应的外幂界与精确行列式共同推出：固定 0<α<1 时，至少

\[
(1-\alpha)n-O_\alpha(n/\log n)
\]

个独立端点 Y 根落在 \(|Y-1|\le e^{-\alpha n\log n}\) 内，按代数重数计。奇异值结论并不自动给出各特征值的同一渐近公式，也不排除相对最大范数意义下的低秩近似。

## 与有理逼近之间仍缺什么

最近的复零点是随 n 改变的代数对象，次数和系数成本都在增长。它们不是分母受控的 p/q。此前的偶数阶最小清分母多项式，在既约有理中心 ia/b 上仍有精确的 \(b^{n^2/2}\) 求值分母成本。

所以这一轮严格解释了一个辅助构造为什么产生极窄的零点邻域，并证明单纯改好矩阵范数上界不能改变“整个复圆盘统一小值”的主尺度。它没有证明所有形状的区域、所有选列或所有方法都遇到同一障碍。接下来的算术目标仍是：构造同时保留相消、非零性，并显著降低次数或有理求值分母成本的对象。

## 证明、复核与复现

- [真实对数导数与最大小值圆盘](period-derivative-scale/PERIOD_LOG_DERIVATIVE.md)，附[独立代理复核](period-derivative-scale/INDEPENDENT_DERIVATIVE_REVIEW.md)。
- [最近 Z 根与最近 Y 根的主阶](nearest-root-scale/NEAREST_ROOT_SCALE.md)，附[独立代理复核](nearest-root-scale/INDEPENDENT_NEAREST_ROOT_REVIEW.md)。
- [端点倒数对称、精确行列式与谱约束](endpoint-pencil-structure/ENDPOINT_PENCIL_STRUCTURE.md)，附[独立代理复核](endpoint-pencil-structure/INDEPENDENT_PENCIL_REVIEW.md)。
- [完整奇异值谱、外幂成本与端点根计数](exterior-endpoint-cost/EXTERIOR_SINGULAR_PROFILE.md)，附[独立代理复核](exterior-endpoint-cost/INDEPENDENT_EXTERIOR_REVIEW.md)。
- [运行记录、校验和及范围说明](CHECKS.json)。各子目录的 `verify_*.py` 仅使用 Python 标准库，直接运行即可重建同目录的验证 JSON。

“独立代理复核”指另一研究代理重新检查论证，不是外部同行评审。核心输入来自此前的[周期双边估计](../pi-period-mechanism-2026-11-02/raw-root-mechanism/HANKEL_PERIOD_MAGNITUDE.md)与[局部 Gram 圆盘](../pi-local-gram-2026-11-02/next-local-period-sensitivity/LOCAL_PERIOD_DISC.md)。本轮不把用户提供的原稿整体结论作为已经独立认证的定理。
