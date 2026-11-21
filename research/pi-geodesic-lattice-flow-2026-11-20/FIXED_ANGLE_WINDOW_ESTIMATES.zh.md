# 固定角附近能把窗口缩多窄

沿用[主报告](FIXED_ANGLE_WINDOWS_AND_SPARSE_PEAKS.zh.md)的 $F,G,T_Q,S_Q,M_Q$，其中 $M_Q$ 是截断后的正主项和。下文在相应位置使用 [π 的无理性指数为 2](https://github.com/openai/math/blob/main/lean/docs/017.md)。归一化误差是否有统一正下界仍未决定。

## 在以 1/π 为中心的缩小区间里，正平均已经成立

令 $\alpha_0$ 为任意无理性指数为 2 的数，包括 $1/\pi$。固定整数 $H\ge1$ 和实数 $0<\kappa<2$。令 $w=Q^{-2+\kappa}$，以 $I_Q=[\alpha_0-w,\alpha_0+w]$ 为圆周上的区间。记 $\langle f\rangle_I=|I|^{-1}\int_I f(\alpha)\,d\alpha$。则
$$
\lim_{Q\to\infty}\langle F_\alpha(Q,H)\rangle_{I_Q}=\frac9{\pi^2H},\qquad
\liminf_{Q\to\infty}\langle G_\alpha(Q,H)\rangle_{I_Q}\ge\frac3{2\pi^2 H}.
$$
特别地，对每个充分大的 $Q$，在 $1/\pi$ 的这个邻域中确有超过半高的读数。$\kappa$ 可以是任意固定正数。结论在黄金比相位的缩小邻域也成立，即使黄金比中心本身在 $H\ge11$ 从不过半高。

证明先数 $I_Q$ 中的节点。置 $X=w^{-1/2}$、$R=Q\sqrt w=Q^{\kappa/2}$。把整数对送到余体积 1 的格 $(q/X,X(q\alpha_0-p))$。节点计数等于该格中位于 $R\Omega$ 的原始点数，其中 $\Omega=\{(x,y):1\le x<2,\ |y|\le x\}$，面积为 3。指数 2 保证最短非零向量 $\lambda\gg_\eta Q^{-\eta}$ 对所有 $\eta>0$ 成立。理由是 $|q|\ge X$ 时横坐标至少 1，否则用 $\|q\alpha_0\|\gg_\eta |q|^{-1-\eta}$，而 $q=0$ 直接检验。取最短向量为第一基向量，第二基向量向其投影后取余，得到直径 $O(1/\lambda)$ 的基本域。对 $R\Omega/d$ 按边界条带估计全格点数，再以 Möbius 反演筛原始点。单次边界误差为 $O(R/(d\lambda)+1/\lambda^2)$，因区域避开原点，只需 $d\ll R/\lambda$。面积主项延伸至无穷的代价 $O(R\lambda)$。从而
$$
\#(S_Q\cap I_Q)= 3R^2/\zeta(2)+O_\eta(RQ^\eta).
$$
此处对最短向量先用更小的指数，以吸收对数和 $\lambda^{-3}$ 后再记为 $Q^\eta$。取 $\eta<\kappa/2$ 即得相对误差趋零。区间端点的包含方式保留在边界估计内。

下面把节点计数用于真实核形。置 $V=Q^{\kappa/2}$、$s=V/Q^2$。把 $I_Q$ 两端各扩张或缩短 $s$ 所改变的节点数为 $O(V+1)=o(R^2)$，由 Farey 间距。位于扩张区间以外的所有核，在 $I_Q$ 上的总值一致为 $O((s^{-2}+Q^2/s)/L^2)=o(1)$，因为两侧距离至少 $s$，同侧按 $1/(4Q^2)$ 分隔。缩短区间内的节点满足 $\int_{I_Q}K_L(\alpha-x)\,d\alpha=1/L+O(1/(L^2s))$，用 $K_L(t)\le1/(4L^2\|t\|^2)$ 对尾部积分即可。$K_L^2\le K_L$ 给出同一尾误差，而其整周期积分为 $(2L^2+1)/(3L^3)$。上、下比较给出 $\langle F\rangle$ 的极限。对 $F^2$ 仅保留非负的单节点平方，得到 $\liminf\langle F^2\rangle\ge(2/3)9/(\pi^2H)$。相减即证。

这里已经能到每个固定的 $Q^{-2+\kappa}$ 宽度。若继续降到 $w\asymp Q^{-2}$，上述 $R$ 只有常数大小，计数误差无法吸收。要从区间平均推得中心的读数，还须另行控制中心附近的振荡。

## 实际缺少哪些傅里叶频率

设 $L=HQ^2$，$m_Q(d)=\sum_{Q\le dr<2Q}\mu(r)$（$r$ 取正整数），$C_Q(k)=\sum_{d\mid k,\ d<2Q}d m_Q(d)$。一次波非零频率为 $(2/L)(1-k/L) C_Q(k)\cos(2\pi k\alpha)$，$1\le k<L$。

**输入频率。** 对固定 $0<\delta<2$，在 $\mu(\alpha)=2$ 时，所有 $k\le K=\lfloor Q^{2-\delta}\rfloor$ 的这些项加在一起趋于零。取 $D=\sqrt K$。写 $k=d\ell$，对 $d\le D$ 先用有限几何和与分部求和处理 $\ell$，再用 $d|m_Q(d)|\le3Q$，得到上界
$$
\ll (HQ)^{-1}\sum_{d\le D}\|d\alpha\|^{-1}.
$$
对 $d>D$ 展开 $m_Q(d)$，固定 $r,\ell$，则 $d$ 在单一区间内，且 $d(1-d\ell/L)$ 的最大值与总变差均为 $O(Q/r)$。按 $d$ 作几何求和并对 $r$ 取绝对值得
$$
\ll (HQ)^{-1}\log(2Q)\sum_{\ell\le K/D}\|\ell\alpha\|^{-1}.
$$
对 $1\le n\le M$，点 $\{n\alpha\}$ 彼此并与零相隔 $\gg_\eta M^{-1-\eta}$，按两侧的距离排序，逆距离之和为 $O_\eta(M^{1+\eta}\log(2M))$。取 $\eta=\delta/4$，总误差为 $O_\delta(H^{-1}Q^{-\delta/4}\log^2(2Q))$，常数可依赖固定角。该误差在 $Q=2^j$ 上可求和。$F$ 一致有界，故先从 $F$ 中删掉这些非零低频再计算 $F^2-F/2$，对前 $J$ 个窗口平均的影响为 $O_{H,\delta}(1/J)$。本段并未将不同输入高频之差形成的低输出频率删去。

**固定的输出频率。** 对每个固定 $n\ne0,H\ge1$，完整二阶量满足
$$
\int_0^1G_\alpha(Q,H)e^{-2\pi i n\alpha}\,d\alpha\longrightarrow0.
$$
这一步不要求 π 的无理性指数。证明用到额外的几何结构。对每个余体积 1 的格 $\Lambda$ 定义
$$
\Psi_H(\Lambda)=\sum_{(x,y)\in\Lambda_{\rm prim},\ 1\le x<2}
 \operatorname{sinc}^2(\pi H y/x).
$$
对整数 $L$，逐周期叠加 $\operatorname{sinc}^2(\pi L(t+m))$，$m\in\mathbb Z$，得到 $K_L(t)$。这由 $\sum_m(t+m)^{-2}=\pi^2/\sin^2(\pi t)$ 及整数点的连续延拓得出。因此
$$
F_\alpha(Q,H)=\Psi_H\big(\{(q/Q,Q(q\alpha-p)):p,q\in\mathbb Z\}\big)
$$
是精确恒等式。分母窗在格空间里变成固定条带 $1\le x<2$。同一格里的不同原始向量在条带内的斜率相隔至少 $1/4$，来自整数行列式。因此 $\Psi_H$ 一致有界，丢掉 $|y/x|>A\ge1$ 的尾和为 $O(H^{-2}/A)$。每个截断和仅可能在有格向量碰到有限条边界线的格上不连续，这些格为 Haar 零测集。

[Marklof，Theorem 1](https://people.maths.bris.ac.uk/~majm/bib/fareystats.pdf) 给出相位与该伸缩格的联合测度弱收敛到均匀相位测度乘 Haar 概率。其格约定与这里相差固定的坐标旋转和参数方向，可等价地用行格 $\mathbb Z^2\begin{pmatrix}1&\alpha\\0&1\end{pmatrix}\operatorname{diag}(Q^{-1},Q)$。对有界且在极限测度下几乎处处连续的检验函数可用弱收敛。乘上 $e^{-2\pi in\alpha}$，则零均值使乘积极限为零。将 $\Psi_H$ 换成 $\Psi_H^2-\Psi_H/2$ 即得所述命题。二阶式的线性项与固定有限项的同号卷积各自趋零（由 $C_Q(k)=O_k(Q)$），故反队指出的该固定频率的长差卷积本身也趋零。此处仅逐个控制固定频率，未在 $1/\pi$ 上交换增长的频率求和与极限。

## 对几乎每个单独固定的相位，正平均确实成立

设 $\gamma_H=\int(\Psi_H^2-\Psi_H/2)\,d\mu$，其中 $\mu$ 为幺模格的 Haar 概率。则对 Lebesgue 几乎每个 $\alpha$，同时对全部正整数 $H$，
$$
\lim_{J\to\infty}J^{-1}\sum_{j=1}^JG_\alpha(2^j,H)=\gamma_H,
\qquad
\frac{3}{2\pi^2H}\le \gamma_H\le\frac{3}{2\pi^2H}+\frac{18B}{\pi^2 H^3},\quad B=8\pi^2/3.
$$
证明用[模曲面测地流的混合性](https://moore.pims.math.ca/sites/default/files/lecture-extra-files/mixing.pdf)，故行格空间上的右乘 $g=\operatorname{diag}(1/2,2)$ 遍历。对可数稠密的紧支撑连续检验函数使用 Birkhoff 定理，得 Haar 几乎处处的离散轨道对 Haar 等分布。这组典型点在右乘下三角幺幂矩阵 $v$ 时性质不变，因 $g^{-j}vg^j\to1$。右乘正对角矩阵时结论亦保持，因其与 $g$ 可交换并保持 Haar 测度。用上三角幺幂参数 $\alpha$、再乘下三角幺幂与正对角矩阵做局部坐标，Haar 在该坐标内等价于 Lebesgue 体积。由 Fubini 及上述两种不变性，几乎每个 $\alpha$ 的起始行格本身就具有等分布轨道。当检验函数有界且几乎处处连续时仍可取极限。对整数 $H$ 取可数交集即可。

计算界时先在相位上积分再用 Marklof 定理。设 $N_Q=|S_Q|$，则 $N_Q/Q^2\to9/\pi^2$，可由 $\varphi(q)=q\sum_{d\mid q}\mu(d)/d$ 求和得到。一次均值恰为 $N_Q/L$，平方中所有单项平方合计积分为 $(N_Q/L)(2/3+1/(3L^2))$。逐点从最近节点以外的背景至多 $B/H^2$ 得到
$$
0\le F^2-\sum_{x\in S_Q} K_L(\alpha-x)^2\le 2(B/H^2) F .
$$
积分并取极限给出两侧常数。主报告给出具有指数 2 而正平均仍失败的具体数，所以“几乎每个”没有决定 $1/\pi$ 是否属于典型点集。

## 截断费用的主项

令 $M_Q$ 为仅保留 $Q\le q<2Q$、$\gcd(q,a)=1$、$|\sin q|<Q^{-3/8}$ 的最近整数 $a=\operatorname{round}(q/\pi)$ 所对应的 $q^2/[Q^4(q-\pi a)^2]$ 之和。使用 $\mu(\pi)=2$，则对每个 $\varepsilon>0$，
$$
T_Q(1)-M_Q=\frac{28}{\pi^3}Q^{-5/8}+O_\varepsilon(Q^{-15/16+\varepsilon}).
$$
本结论既容许 $M_Q$ 有界，也容许它无界。

**证明。** 未截断的最近节点之和与 $T_Q(1)$ 相差 $O(Q^{-1})$，直接来自每个分母的其余等距节点估计和 $\csc^2 y=y^{-2}+O(1)$。写 $\sigma=\arcsin(Q^{-3/8})$，则被删主项是 $Q^{-4}\sum q^2/u^2$，求和遍及互素整数 $q,a$，其中 $Q\le q<2Q$，$\sigma\le u=|q-\pi a|\le\pi/2$。对 $u$ 的正负两支分别取 $h\le u<2h$ 的倍尺度段，最后一段截到 $\pi/2$。

在每段取 $X=\sqrt{Q/h}$、$R=\sqrt{Qh}$，使用余体积为一的格 $\Lambda=\{(q/X,X(q/\pi-a))\}$。该段及其轴对齐子矩形均落在直径 $O(R)$ 的框里，且与原点分开。在 $(q,u)$ 坐标下面积为 $A$ 的子矩形，变换后面积为 $A/\pi$。记 $\lambda$ 为格中非零向量的最短长度。无理性指数 2 传给倒数，给出 $\|q/\pi\|\gg_\delta |q|^{-1-\delta}$。对 $|q|\ge X$ 用横坐标，对 $1\le|q|<X$ 用纵坐标，$q=0$ 直接检查。因此 $\lambda\gg_\delta X^{-\delta}\gg_\delta Q^{-\delta}$。

取一个最短格向量和一个行列式为 1 的基补向量。将补向量沿最短向量约化，基本平行四边形的直径为 $O(1/\lambda)$，这里用到单位余体积使 $\lambda$ 有绝对上界。子矩形缩小 $d$ 倍后，用平行四边形覆盖其边界邻域，得到全部格点计数为
$$
A/(\pi d^2)+O\big(R/(d\lambda)+1/\lambda^2\big).
$$
对 $d$ 以 Möbius 权求和留下原始向量。由于缩小后的集合不含原点且每一点的长度 $\ll R/d$，只需 $d\ll R/\lambda$。将面积主项补至无穷损失 $O(R\lambda)$，其余累计 $O((R/\lambda)\log(2Q)+R/\lambda^3)$。故所有子矩形的原始格点计数均为 $A/(\pi\zeta(2))+O_\varepsilon(RQ^\varepsilon)$。

对 $s=q/Q$ 和 $v=u/h$ 作二维分部求和，权 $s^2/v^2$ 在 $[1,2]^2$ 上的一阶偏导和混合偏导均有界。每段加权和除以 $Q^4$ 后的误差为 $O_\varepsilon((Qh)^{-3/2}Q^\varepsilon)$。从 $h=\sigma$ 开始的几何级数总误差为 $O_\varepsilon(Q^{-15/16+\varepsilon})$。两支面积积分为
$$
\frac{2}{\pi\zeta(2)Q^4}\int_Q^{2Q}q^2\,dq\int_\sigma^{\pi/2}u^{-2}\,du
=\frac{14}{3\pi\zeta(2)Q}\left(\sigma^{-1}-2/\pi\right).
$$
用 $\sigma^{-1}=Q^{3/8}+O(Q^{-3/8})$ 及 $\zeta(2)=\pi^2/6$ 即得结论。
