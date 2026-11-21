# 从零点窗口和无外力爆破的证书结构继续研究单位圆拼积

## 原始材料

- 黎曼猜想项目的 `zeroFree_of_rpow_window` 把一个窗口上界归约为零点排除。其窗口量词是每个起点后都存在一整段满足同一个上界的窗口。Alpöge 公告于 2026 年 11 月 20 日。Sutherland 的 MIT 审计已报告全库编译、逐模块回放以及两个独立内核接受完整导出。审计状态以原页为准。[Lean 判据](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab/src/Development/AutomorphicArithmetic/DedekindLFunctions/LandauTheorems/ZeroFreeOfRpowWindow.lean#L11-L14) · [Alpöge 公告](https://alpo.ge/remark.pdf) · [MIT 审计](https://math.mit.edu/~drew/grh/)
- Alpöge–Buckmaster 的 2026 年 11 月 20 日无外力 Navier–Stokes 预印本把方案分为近似自相似剖面、残差与稳定性证书、真实剖面和爆破推论。Section 7 与 Appendix B 适合借鉴证书设计。公开仓库说明数值证书和区间核查程序仍未随版本发布，补丁后的 Lean 代码在组装时未完整重编译。因此这里借用条件逻辑的结构。[论文仓库和发布状态](https://github.com/tristanbuckmaster/navier-stokes/blob/ad53419de079eb8f95f60313473374624a348546/README.md) · [Lean 与论文对照](https://github.com/tristanbuckmaster/navier-stokes/blob/ad53419de079eb8f95f60313473374624a348546/lean/README.md)
- 原拼圆笔记中的几何下界可补成双侧。Bárány–Larman 的 Theorem 4 对固定圆心、确定性整数凸包给出 $\mathrm{area}(NB\setminus\operatorname{conv}(NB\cap\mathbb Z^2))\le C N^{2/3}$。缩回单位圆盘后，$\delta_N=\pi-S_N$ 满足 $1/(50N^{4/3})<\delta_N\le C N^{-4/3}$。取固定 $C\ge1$，必要时把有限小尺度纳入常数。[Bárány–Larman 原文](https://www.math.ucdavis.edu/~deloera/MISC/LA-BIBLIO/trunk/Barany/Baranyrpdcg.pdf) · [本项目拼圆定义](../rational-disk-packing/README.zh.md) · [稀疏峰与窗口量词](../FIXED_ANGLE_WINDOWS_AND_SPARSE_PEAKS.zh.md)

## 方向 A：从一次极端逼近得到一整段网格信号

设 $p/q<\pi$ 已既约，$\Delta=\pi-p/q>0$，$E=q^2\Delta$，$M_N=2N^2 S_N$。选择二进尺度 $X$ 满足
$$
(2C/\Delta)^{3/4}\le X<2(2C/\Delta)^{3/4}.
$$
对每个 $X\le N<2X$ 且 $q\mid 2N^2$ 的网格，令 $j=M_N-2N^2p/q$。由于 $\delta_N\le\Delta/2$，有 $j\ge0$。又 $j=2N^2(\Delta-\delta_N)\le2\Delta N^2$，故
$$
0\le j\le L N^{2/3},\qquad L=4\cdot4^{4/3}C.
$$
这个整数等于从最大拼图中去掉的最小格点三角形数。特别地，所有 $q\mid N$ 的网格都合格，至少有 $X/q-1$ 个。该数量随异常程度变化，尺度满足 $X\asymp q^{3/2}E^{-3/4}$。因此即使极端分数本身很稀疏，单个极端分数也会在自己的网格区间留下多次信号。

可精确校正重复计数。取 $r(q)$ 为满足 $q\mid 2r(q)^2$ 的最小正整数。逐素数验算给出 $q\mid2N^2$ 当且仅当 $r(q)\mid N$，且 $r(q)\le q$。对二进 $X\ge1$，设 $Q=X^{2/3}$。定义 $F(X)$ 为对 $X\le N<2X$ 及整数 $0\le j\le\min(M_N,LN^{2/3})$ 求和，保留约分后分母 $q\le Q$ 的项 $(M_N-j)/(2N^2)$，每项权重为
$$
\frac{r(q)}{X}\left(\frac Qq\right)^2.
$$

**命题。** $F(X)$ 在所有二进 $X$ 上有统一上界，当且仅当 $\inf_{p/q<\pi}q^2(\pi-p/q)>0$，其中下确界取既约分数且 $q>0$。

证明上界时，窗口内的分数都落在长度 $(C+L/2)/Q^2$ 的区间中。两条不同既约分数至少相距 $1/Q^2$，因此只有固定数量的不同值。每个固定值出现至多 $X/r(q)+1$ 次。一次乘以 $r(q)/X$ 后消去了重复次数。若从下方的归一化误差下界为 $c>0$，则 $(Q/q)^2\le(C+L/2)/c$。

反向取一列 $E\to0$ 并按上式为每个分数选 $X$，则 $(Q/q)^2\ge 2C/E$，所以最终满足 $q\le Q$。整段里至少有 $X/r(q)-1$ 个允许的 $N$，除以 $X/r(q)$ 后仍有正的统一下界。于是 $F(X)$ 在这一子列发散。

**任务。** 独立重推这个命题，再从凸包整数数据和局部同余出发，尝试给 $F(X)$ 建立统一上界或给出迫使其变大的机制。可用 Pick 公式 $M_N=2I_N-B_N-2$（$I_N$ 计所有格点，$B_N$ 计边界格点）及素因子分解。若只能证明平均界或任意晚存在安静窗口，用稀疏大部分商对照验证它确实漏掉哪些波峰。一个具体的中途交付是某个说明了覆盖范围的网格家族上的统一结论及遗漏分母反例。内拼图只处理从下方逼近，外侧须另证。

**双侧版本。** 保留同一权重，把候选改成 $(M_N+j)/(2N^2)$，范围改为 $|j|\le (L+2C)N^{2/3}$，取 $\Delta=|\pi-p/q|$。对任一侧的分数，用 $j=2N^2(p/q-S_N)$，三角不等式给出 $|j|\le 2N^2(\Delta+\delta_N)$，上述窗口证明仍成立。分数集中区间的长度增加一个固定倍数，Farey 间距证明仍适用。这个带符号求和的统一上界等价于通常的双侧坏逼近。$j>0$ 是对面积做的有理外推，实际的圆内拼图对应 $j\le0$。

## 方向 B：重整化、剖面残差和跨尺度接口

令 $N=2^n$，$M_n=2\cdot4^n S_{2^n}$，$b_n=M_{n+1}-4M_n$。嵌套网格给出 $b_n\ge0$，且几何上界给出 $b_n=O(4^{n/3})$。确切的空隙可逐尺度重写为
$$
\pi-S_{2^n}=\sum_{k\ge n}\frac{b_k}{2\cdot4^{k+1}}.
$$
试用最简单的离散自相似型 $b_{n+3}=4b_n$，与几何尺度一致，但若在足够大的 $n$ 永远成立，三个余数类上的级数都成为可精确求和的有理几何级数，推出 $\pi$ 为有理数。若整数残差 $r_n=b_{n+3}-4b_n\to0$，也会落入同一矛盾。更一般地，$b_n$ 若最终满足有理常系数齐次线性递推，其生成函数为有理函数。由增长界，$1/4$ 严格落在收敛半径内，在该点求值仍将推出 $\pi$ 为有理数。这些障碍针对精确封闭的数列拟合。

**任务。** 借鉴 Section 7 与 Appendix B，设计一个有非周期残差的可验证尺度推进。候选局部剖面可取抛物线。圆顶部坐标 $x=N^{1/3}u$、$y=N-N^{-1/3}v$ 中，圆内判据精确写作 $u^2-2v+N^{-4/3}v^2\le0$。在有界盒子里把远离阈值的点确定下来，单独处理靠近抛物线的格点。说明相邻局部块如何接合、如何合并面积与模信息，以及尾误差如何在目标既约分母 $q$ 的 $q^{-2}$ 尺度上控制。若选择二进和三进网格共同推进，核对 $b_3(2N)+9b_2(N)=b_2(3N)+4b_3(N)$，其中 $b_k(N)=M(kN)-k^2M(N)$、$M(N)=2N^2S_N$，并明确在哪一步加入了格点几何信息。

请交付一项非周期残差定理、一个可验证的局部到整体证书接口，或一个说明该接口为什么失效的严格反例。先取小的整数基准独立校验；全尺度结果要注明所有统一常数和接口条件。
