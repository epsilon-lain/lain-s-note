# 两个可以闭合的问题：复采样解码阈值与迁移分形的厚度振荡

2026-10-31。核对基线：`3c82b0db8d27933928494fcc128789d86616151a`。

**范围。** 本文不再增加 P2848 的等价表述，而是完成两个独立问题：

1. 有限个带绝对误差的圆周复底数观测，能否保证恢复任意有界整数数位串的前 N 位？给出精确最坏情况阈值、达到阈值的反例及最优半径。
2. 已构造的迁移 Cantor 集 K_t 是否有一个稳定的 Minkowski 内容常数？给出否定答案、显式对数周期主项及主复维数。

这些是将已有单位根滤波、整数数位约束、裂缝计数和返回计数渐近应用到具体问题。Fourier/Cauchy 采样与分形管公式均有成熟文献，不声称基本方法或以下全部命题是文献首创；尚未做穷尽性的创新性检索。它们不是 π 的新丢番图逼近界，也不是某个已命名开放猜想的解决。

依赖的仓库记录：
- [复底数解析性与四方向滤波](../complex-base-analyticity-2026-10-31/README.md)，§§2、8。
- [迁移空间的裂缝计数](../base-transport-2026-10-31/README.md)，§§4—5。
- [返回计数的定量渐近](../pi-intrinsic-arithmetic-2026-10-31/README.md)，§4。

## A. 一个精确解码问题

### A1. 数据模型

固定正整数 A、N。未知串 a=(a_n) 满足 a_n∈{0,1,…,A}，定义

\[
F_a(z)=\sum_{n\ge1}a_nz^{-n},\qquad |z|>1.
\]

选择半径 R>1，令 ω=exp(2πi/N)，在 N 个精确位置 z_j=Rω^j 观测

\[
y_j=F_a(z_j)+e_j,\qquad |e_j|\le\eta,
\qquad 0\le j<N.
\]

误差可为任意复数，不假设随机、独立或零均值。目标只恢复 a_1,…,a_N，不恢复整个无限尾部。取样位置本身先假设精确。

记 X=R^N。对 r=1,…,N，计算

\[
B_r=\frac{R^r}{N}\sum_{j=0}^{N-1}y_j\omega^{jr}.
\tag{A1}
\]

r=N 对应离散 Fourier 的零剩余类；没有常数项 a_0。

### A2. 精确最坏情况阈值

**定理 A。** 若 X>A+1，则令

\[
T=\frac{A}{X-1},\qquad
\boxed{\eta_*(A,X)=\frac{1-T}{2X}
=\frac{X-A-1}{2X(X-1)}.}
\tag{A2}
\]

对于任意 η<η_*(A,X)，公式

\[
\boxed{\widehat a_r=
\left\lfloor\Re B_r-\frac T2+\frac12\right\rfloor,
\quad 1\le r\le N,}
\tag{A3}
\]

对所有允许的无限串和所有允许的观测误差，都精确给出 a_r。

而在 η≥η_* 时，不存在任何解码方法能对这一取样布局和整个数位类统一保证恢复前 N 位。故这是最坏情况容噪半径的精确上确界；等号不在保证范围内。

**恢复证明。** 单位根正交关系与绝对收敛给出

\[
\frac{R^r}{N}\sum_jF_a(R\omega^j)\omega^{jr}
=\sum_{k\ge0}a_{r+kN}X^{-k}
=a_r+\tau_r,
\]

其中

\[
0\le\tau_r\le\frac{A}{X-1}=T.
\]

误差的 Fourier 分量满足

\[
\left|\frac{R^r}{N}\sum_je_j\omega^{jr}\right|
\le R^r\eta\le X\eta.
\]

因此

\[
\left|\Re B_r-\frac T2-a_r\right|
\le\frac T2+X\eta<\frac12,
\]

整数舍入唯一，得到 (A3)。

**不可改进的反例。** 第一条串在 2N、3N、4N、… 位全部取 A，其余位为 0；第二条串只在第 N 位取 1，其余位为 0。两者的前 N 位不同。在每个采样点处，它们分别取同一实数

\[
F_a(R\omega^j)=\frac{A}{X(X-1)},\qquad
F_b(R\omega^j)=\frac1X.
\]

两组数据在每个坐标上的距离恰为

\[
\frac1X-\frac{A}{X(X-1)}=2\eta_*(A,X).
\]

把两组数据的中点作为观测，就与两条不同数位串都相容，且各自所需误差恰为 η_*。任何方法都无法统一区分。证毕。

这个反例允许最终循环或有限串，因定理的对象类就是所有 a_n∈{0,…,A}；它不声称是在仅限 π 合法串的更小类上仍然最优。

### A3. 半径过小时，连无噪声统一识别都不可能

若 1<X≤A+1，存在 c_k∈{0,…,A} 使 Σ_{k≥1}c_kX^{-k}=1。

当 X<A+1 时，单位的贪心 X 展开可用这些数字（整数 X≤A 时也允许首位 X）；当 X=A+1 为整数时，可用全 A 尾串。于是“仅在 N 位取 1”的串与“在 (k+1)N 位取 c_k”的串，在全部 N 个采样点的值相同，前 N 位却不同。因此 X>A+1 是这个整体数位类的可恢复性阈值。

### A4. 最优半径及 π 字母表上的简单选择

对 X>A+1 最大化 (A2)，得到唯一极大点

\[
\boxed{X_*=(A+1)+\sqrt{A(A+1)},\qquad
R_*=X_*^{1/N},}
\]

\[
\boxed{\max_X\eta_*(A,X)
=\frac1{2(\sqrt{A+1}+\sqrt A)^2}.}
\tag{A4}
\]

证明只需对 (X-A-1)/(2X(X-1)) 求导；导数分子的唯一相关零点为 X²−2(A+1)X+(A+1)=0。

对 A=3，X_*=4+2√3，最优阈值为 7/2−2√3≈0.03589838486。使用略简单的 X=8，即 R=8^{1/N}，就有

\[
\boxed{\eta<1/28\quad\Longrightarrow\quad
N\text{ 个观测保证恢复前 }N\text{ 位}.}
\tag{A5}
\]

该绝对误差阈值不随 N 缩小，是因为半径随 N 靠近单位圆，不是因为深位置不再衰减。固定 R>1 时，η_*渐近为 (1/2)R^{-N}。

**没有免费产生信息。** 这里假设能够获得 F 的带误差观测；若观测本身由已知 π 数位生成，恢复只是验证解码器，不是新的 π 数位生成算法。R=8^{1/N} 时，最坏的函数幅度 A/(R−1)=O(N)，而导数上界 A/(R−1)²=O(N²)。因此数据动态范围及采样位置精度仍有成本。小的节点误差应经导数界计入 η。该结果不支持从少量固定精度实数读取无限信息。

## B. 迁移 Cantor 集有稳定的“厚度常数”吗？

### B1. 固定几何对象

固定 π 的原生单位串 d。令

\[
\Phi(t)=\sum_{n\ge1}d_nt^{-n},\quad t>\pi,
\qquad \delta=1-\Phi(t)>0.
\]

K_t 是全部闭 π-合法数位串按底数 t 求值的像，位于 [0,Φ(t)]。已有记录证明其 Lebesgue 长度为 0、Hausdorff 维数

\[
D=\frac{\ln\pi}{\ln t}\in(0,1),
\]

且第 n 层恰有 c_n 条裂缝，每条长度为 δt^{-n}，其中

\[
c_0=1,\qquad c_n=\sum_{k=1}^n d_kc_{n-k}.
\]

记

\[
m=\sum_{n\ge1}n d_n\pi^{-n}.
\]

已证明的定量渐近为

\[
\boxed{\left|c_n-\frac{\pi^n}{m}\right|
\le\frac{2^n}{\pi-2}.}
\tag{B1}
\]

它来自余数生成函数的分解，不依赖任何关于 π 部分商的猜想。

### B2. 管长的精确公式

令 W_t(ε) 为 [0,Φ(t)] 中，到 K_t 距离小于 ε 的点集的长度。每条长度为 l 的裂缝贡献 min(2ε,l)，故

\[
\boxed{W_t(\varepsilon)
=\sum_{n\ge1}c_n\min\{2\varepsilon,\delta t^{-n}\}.}
\tag{B2}
\]

通常在整个实轴计算的 ε 邻域长度为 V_t(ε)=W_t(ε)+2ε；后面主项不受这两个外侧端点影响。

置 h=2ε，并令

\[
N=\left\lfloor\log_t\frac\delta h\right\rfloor,
\qquad \theta=\log_t\frac\delta h-N\in[0,1).
\]

当 ε 足够小时 N≥1，(B2) 可精确分为

\[
W_t(\varepsilon)
=h\sum_{n=1}^Nc_n+\delta\sum_{n>N}c_nt^{-n}.
\tag{B3}
\]

### B3. 完整结论：没有单一内容常数，有明确的周期主项

**定理 B。** 对每个固定 t>π，

\[
\boxed{\frac{W_t(\varepsilon)}{(2\varepsilon)^{1-D}}
=P_t(\theta)+O_t\big((2/\pi)^N\big),}
\tag{B4}
\]

其中

\[
\boxed{P_t(\theta)=\frac{\delta^D}{m}
\left[\frac\pi{\pi-1}\pi^{-\theta}
+\frac\pi{t-\pi}\left(\frac t\pi\right)^\theta\right],
\quad 0\le\theta\le1.}
\tag{B5}
\]

P_t(0)=P_t(1)，故可以连续地作 1 周期延拓；在 (0,1) 上 P_t''>0，所以它不是常数。于是 K_t 没有 Minkowski 内容极限，但其上下 Minkowski 内容均为有限正数。

**这里的“不具有 Minkowski 内容”不是 Lebesgue 不可测。K_t 是紧集，Lebesgue 长度已知为 0。** 问题只是其邻域长度除以维数预言的 ε^{1-D} 后，是否收敛到一个常数；答案是否定的。

**证明。** 将 c_n=π^n/m+e_n 代入 (B3)。主项中的几何级数给出

\[
\frac h m\frac{\pi^{N+1}-\pi}{\pi-1}
+\frac\delta m\frac{(\pi/t)^{N+1}}{1-\pi/t}.
\]

利用 h=δt^{-N-θ} 并除以 h^{1-D}，保留的两项恰为 (B5)。遗漏的有限下端项为 −hπ/(m(π−1))。由 (B1)，总误差的绝对值不超过

\[
\frac{h\pi}{m(\pi-1)}
+\frac h{\pi-2}(2^{N+1}-2)
+\frac\delta{\pi-2}\frac{(2/t)^{N+1}}{1-2/t}.
\tag{B6}
\]

再次除以 h^{1-D}，得到 (B4)。P_t 由两个系数为正的指数项构成，二阶导数严格正；端点同值可直接代入检验。因此它非恒定。

沿 ε_N=(δ/2)t^{-N-θ_0} 的任意固定相位子列，归一化长度趋于 P_t(θ_0)。选择 θ_0=0 和 θ_0=1/2，两个极限不同。因此单一内容常数不存在。通常的 V_t(ε)/ε^{1-D} 只是将这些极限乘以 2^{1-D}，结论相同。证毕。

**十进制例子。** t=10 时，D≈0.49714987269413385435。两条相位子列的 W/(2ε)^{1-D} 极限分别为

\[
P_{10}(0)\approx1.45297922125531,
\qquad
P_{10}(1/2)\approx1.24153992280782.
\]

数值是高精度显示；非相等由严格凸性和相同端点证明，不依赖小数判断。

周期平均仍然存在，而且

\[
\boxed{\int_0^1P_t(\theta)\,d\theta
=\frac{\delta^D}{m}
\left(\frac1{\ln\pi}+\frac1{\ln(t/\pi)}\right).}
\tag{B7}
\]

这也是相应对数尺度 Cesàro 平均的归一化内容；若改用 ε 而非 2ε 归一化，应再乘 2^{1-D}。

### B4. 复维数直接来自旧计数生成函数

把裂缝长度的 s 次幂求和，定义几何 zeta 函数

\[
\zeta_t(s)=\sum_{n\ge1}c_n(\delta t^{-n})^s,
\qquad\Re s>D.
\]

所有正数幂使用实对数定义。令 C(u)=Σ_{n≥0}c_nu^n，则

\[
\zeta_t(s)=\delta^s\big(C(t^{-s})-1\big).
\]

(B1) 保证

\[
C(u)=\frac1{m(1-\pi u)}+B(u),\qquad B\text{ 在 }|u|<1/2\text{ 全纯}.
\]

因此 ζ_t 亚纯延拓到 Re s>ln2/lnt，并在这一半平面只有如下简单极点：

\[
\boxed{s_k=\frac{\ln\pi}{\ln t}+\frac{2\pi i k}{\ln t},
\quad k\in\mathbb Z,\qquad
\operatorname{Res}_{s=s_k}\zeta_t(s)=\frac{\delta^{s_k}}{m\ln t}.}
\tag{B8}
\]

虚部的等间距对应对数尺度的周期。这是标准分形串/复维数语言在本模型中的具体计算；(B4) 的证明直接求和完成，不需要调用更强管公式定理的增长条件。

## C. 文献关系与可复现范围

### 外部背景

1. L. N. Trefethen and J. A. C. Weideman, *The Exponentially Convergent Trapezoidal Rule*, SIAM Review 56 (2014), 385–458, DOI 10.1137/130932132。§2 给出单位根采样、Taylor 系数混叠与几何尾界；§12 讨论 Cauchy 积分计算。https://people.maths.ox.ac.uk/trefethen/publication/PDF/2014_149.pdf
2. M. L. Lapidus, G. Radunović and D. Žubrinić, *Complex dimensions of fractals and meromorphic extensions of fractal zeta functions*, arXiv:1508.04784v4 (2017)。§1.2 定义 Minkowski 内容；§4 讨论对数周期管主项和复维数。https://arxiv.org/abs/1508.04784

第 A 节的明确对象类、误差模型、中心校正、达到阈值的双串反例和半径优化，都有本文的自足证明。第 B 节是仓库裂缝结构与已知渐近的直接推论。我们没有将标准 Fourier 方法、整数舍入或格点尺度振荡本身当作文献首创。

### 实际运行

`python verify_closed_problems.py` 重新用 Machin 恒等式和向外取整区间认证 512 位；没有输入预存 π 小数串或连分数。

- 4 点复采样例子全程使用高斯有理数，精确恢复 3011；同时追踪无限尾误差。
- 4,165 个分量区间检查使用 Fraction，验证恢复区间严格落在正确整数的半单位邻域内。
- 12 对精确双串数据达到 η_* 不可区分阈值，另有 4 个 X=A+1 的零噪声碰撞。
- 对返回计数验证到 512 阶的生成函数系数恒等式。
- 可选 mpmath 以 90 位精度演示 N=16、32、64 的复采样解码，以及 t=4、10、16 下两种相位、三个深度的管长趋近。这一部分明确是高精度数值演示，不是向外取整区间证书。

一般结论由上面的证明给出；有限检查只核对实现。未运行 Lean。运行摘要见 `verification_summary.json`；程序运行会生成完整 `verification_results.json`，下载验证包也保留完整结果。

**完成与未完成。** 已完成的是“所定采样模型下可恢复的精确噪声阈值”与“所定迁移分形的 Minkowski 内容及周期主项”。没有因此得到 P2848、Φ(10) 超越性或一个可商用测量装置的证明。把这两个小问题做完整，比继续给同一个大缺口改名更有实际产出。
