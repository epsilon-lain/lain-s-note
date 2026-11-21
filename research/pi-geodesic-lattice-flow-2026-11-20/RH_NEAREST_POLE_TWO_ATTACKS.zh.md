# 从一个最近极点到两项 π 研究任务

目标是判别是否存在常数 $c>0$，使每个正整数 $q$ 和整数 $p$ 都满足 $q^2|\pi-p/q|\ge c$。本文建立两个可以向相反方向使用的窗口判据，证明可排除的一部分项，并给出一座有限峰的证书。两边针对 π 的关键估计均待证明。

## Lean 中实际用到的结构

[landau_pole_asym](https://github.com/SINF-KEN/rh-lean-check/blob/40e56bab/src/Development/AutomorphicArithmetic/DedekindLFunctions/LandauTheorems/LandauPoleAsym.lean#L20-L24) 从函数在孤立极点附近的分解出发。中心为 $z$，极点为 $\rho$，留数为 $c$。剥去 $c/(s-\rho)$ 后，解析余项贡献的第 $k$ 阶 Taylor 系数至多为 $M/(r')^k$。当 $d=|z-\rho|<r'$ 时，原极点的系数模为 $|c|/d^{k+1}$，两项可以按增长率分开。[landau_window](https://github.com/SINF-KEN/rh-lean-check/blob/40e56bab/src/Development/AutomorphicArithmetic/DedekindLFunctions/LandauTheorems/landauWindow.lean) 在其中一个分支隔离最近极点，与另一侧的窗口估计比较，调用 landau_window_contra 取得矛盾。对源码作了定理和调用处的核读，本次未重新编译 Lean 全库。[MIT 核验记录](https://math.mit.edu/~drew/grh/)

下面从格点斜率构造一条一实变量函数，以最近既约分数充当孤立信号。这里所有关于 π 的推论只用下述公式，无须以 RH 为假设。

## 共同几何

取 $\alpha=1/\pi$ 和二的幂 $Q=2^j$，$j\ge0$。一个分母窗口内的节点是
$$
S_Q=\{a/q: Q\le q<2Q,\ 0\le a<q,\ (a,q)=1\}\subset\mathbb R/\mathbb Z .
$$
记 $\|t\|$ 为到最近整数的距离，$\rho=\min_{x\in S_Q}\|\alpha-x\|$。取倒数保持坏逼近性，因为逼近 $\alpha$ 的 $a/q$ 满足 $a^2|\pi-q/a|=\pi(a/q)q^2|\alpha-a/q|$，其中 $\pi(a/q)\to1$，反向同理。足够近的圆周逼近可选取实直线上邻近 $\alpha$ 的代表。

两个不同节点的圆周距离至少 $1/(4Q^2)$，由交叉乘积是非零整数。半径 $r$ 的弧上至多有 $1+8Q^2r$ 个节点。将所有距离从小到大排成 $\rho_i$，即得 $\rho_i\ge (i-1)/(8Q^2)$，$i\ge2$。这个排序给了两个可显式控制的背景误差。

## A 队可以使用的极点与已排除的分母

定义随实变量 $x$ 移动的读数
$$
T_Q(x)=Q^{-4}\sum_{a/q\in S_Q}\csc^2(x-\pi a/q).
$$
在 $x=1$，抽出最近节点，再用 $|\sin \pi y|\ge 2\|y\|$ 对其余节点求和，得到
$$
0\le T_Q(1)-Q^{-4}\csc^2(\pi\rho)\le B,\qquad B=8\pi^2/3 .
$$
因此 $\sup_j T_{2^j}(1)<\infty$ 当且仅当 π 坏逼近。一枚归一化误差趋零的分数会把第一项推向无穷。反之，若 $\alpha$ 有坏逼近常数 $c$，则 $\rho\ge c/(4Q^2)$，故 $T_Q(1)\le 4/c^2+B$。

设 $\Phi_q$ 是第 $q$ 个分圆多项式。从每个线性因子的对数模求二阶导数可得
$$
T_Q(x)=-Q^{-4}\frac{d^2}{dx^2}\log\left|\prod_{Q\le q<2Q}\Phi_q(e^{2ix})\right|.
$$
定义整数 $m_Q(d)=\sum_{r\ge1,\ Q\le dr<2Q}\mu(r)$，其中 $\mu$ 是 Möbius 函数。用 Möbius 反演从所有 $a/q$ 中留下既约分数。再对余切乘法公式求导，得到 $\sum_{a=0}^{d-1}\csc^2(x-\pi a/d)=d^2\csc^2(dx)$。于是在固定点
$$
T_Q(1)=Q^{-4}\sum_{1\le d<2Q}m_Q(d)\frac{d^2}{\sin^2 d}.
$$

**已经可以丢弃的项。** 设 $E_Q$ 由 $Q^{3/4}<d<2Q$ 且 $|\sin d|<Q^{-3/8}$ 的整数构成。则
$$
T_Q(1)=Q^{-4}\sum_{d\in E_Q}m_Q(d)d^2/\sin^2 d+O(Q^{-1/4}).
$$
这里使用 [π 的无理性指数为 2](https://github.com/openai/math/blob/main/lean/docs/017.md)。对每个固定 $\eta>0$ 有 $\|h/\pi\|\gg_\eta h^{-1-\eta}$。对 $0\le h\le D$ 的圆周点 $h/\pi$，任意两点至少相隔常数倍的 $D^{-1-\eta}$。包括原点后，从两侧分别排序，得到 $\sum_{d\le D}|\sin d|^{-2}\ll_\eta D^{2+2\eta}$。又 $|m_Q(d)|\le3Q/d$。因此 $d\le D$ 的绝对贡献为 $O_\eta(D^{3+2\eta}/Q^3)$。取 $D=Q^{3/4}$ 和 $\eta=1/8$，这部分为 $O(Q^{-9/16})$。对 $|\sin d|\ge Q^{-3/8}$ 的项，直接用 $Q^{-4}(3Q)Q^{3/4}\sum_{d<2Q}d=O(Q^{-1/4})$。两部分可按不相交方式删去。

余下每项的 $r$ 都严格小于 $2Q^{1/4}$。由于间距估计，只剩 $O_\eta(Q^{5/8+\eta})$ 个候选 $d$。证明全局上界仍需要控制其中异常小的 $|\sin d|$。形式上平移可给出 $\cot d-\cot(d+h)=\sin h/[\sin d\sin(d+h)]$。固定非零整数 $h$ 时，随实变量 $d\to0$，右侧只有一阶极点。用有界系数对有限个 $1/[\sin d\sin(d+h)]$ 作线性组合，无法在极点附近逐点控制 $\csc^2d$ 的二阶极点。仅在整数点的额外算术估计仍可探索。

## B 队可以使用的有界波峰与反例提取

再取正整数 $H$，设 $L=HQ^2$，定义
$$
K_L(y)=\left|L^{-1}\sum_{t=0}^{L-1} e^{2\pi ity}\right|^2,\quad
F_\alpha(Q,H)=\sum_{x\in S_Q} K_L(\alpha-x).
$$
这是实变量 $\alpha$ 的非负周期函数。$K_L(0)=1$ 且 $K_L(y)\le\min(1,1/(4L^2\|y\|^2))$。上述距离排序给出
$$
0\le F_\alpha(Q,H)-K_L(\rho)\le B/H^2.
$$
这里抽取按距离最近的一项，其余项合计有可控上界，无须假设核在全圆周单调。若 $H\ge11$ 且 $F_\alpha(Q,H)\ge1/2$，则 $K_L(\rho)>1/4$，因此某个既约 $a/q$ 有 $q^2\|\alpha-a/q\|\le4/H$。所以只要在趋于无穷的 $H$ 上构造这种峰，π 就不是坏逼近数。$Q$ 必然趋于无穷，因为固定有限窗口在 $H\to\infty$ 时读数趋零。

更精确地，先让 $Q$ 过任意晚的二次幂再提高 $H$，有
$$
\lim_{H\to\infty}\limsup_{j\to\infty}F_\alpha(2^j,H)
=\begin{cases}0&\alpha\text{ 坏逼近},\\1&\alpha\text{ 非坏逼近}.\end{cases}
$$
在第一种情形 $F\le(4/c^2+B)/H^2$。在第二种情形，固定 $H$，取 $q^2|\alpha-a/q|\to0$ 并取包含 $q$ 的窗口，便有 $L|\alpha-a/q|\to0$，对应核值趋于 1。统一上界 $1+B/H^2$ 完成夹逼。若把两个极限的次序交换，每个无理数都给出零。

还可复制 Lean 中“持续到所有充分高阶”的性质。固定 $H$，写 $F_j=F_\alpha(2^j,H)$。令 $D_k=\sum_{j=0}^k2^{-j-1}F_j^k$，则 $\lim_k D_k^{1/k}=\sup_j F_j$。上界由权重和至多 1 得到，对固定 $j$ 保留单项给出下界。某个 $j$ 的峰在所有 $k\ge j+1$ 的根式读数中至少贡献其原高度的一半。取高次幂用于保留稀疏峰，独立的算术上界或无限峰构造仍是必要步骤。

## 用另一个公式寻找高峰

记 $N_Q=|S_Q|=\sum_{Q\le q<2Q}\varphi(q)$，$C_Q(k)=\sum_{d\mid k,\ d<2Q} d\,m_Q(d)$，$k\ge1$。有限傅里叶展开与 Möbius 反演给出
$$
F_{1/\pi}(Q,H)=\frac{N_Q}{L}+\frac2L\sum_{k=1}^{L-1}(1-k/L)C_Q(k)\cos(2k).
$$
也可从每个 $d$ 的完整根系直接求和。令 $s=L/d$，$m=\lfloor s\rfloor$，$t=s-m$。则在固定相位有较短的精确公式
$$
F_{1/\pi}(Q,H)=\sum_{d<2Q,\ d\ge1}m_Q(d)\,
\frac{(1-t)\sin^2(md)+t\sin^2((m+1)d)}{s^2\sin^2d}.
$$
它来自把长度 $s$ 的三角傅里叶权写成长度 $m$ 与 $m+1$ 的两个未归一化平方和的线性组合，然后除以 $s^2$。[Ramanujan 和的定义与 Möbius 公式](https://mathworld.wolfram.com/RamanujansSum.html)

**独立的可检验猜想。** 设 $G_\alpha(Q,H)=F_\alpha(Q,H)^2-F_\alpha(Q,H)/2$。对全体相位，平移不变性和各节点贡献的非负性给出
$$
\int_0^1 F_\alpha\,d\alpha= N_Q/L,\qquad
\int_0^1 G_\alpha\,d\alpha\ge\frac{N_Q}{L}\left(\frac16+\frac{1}{3L^2}\right).
$$
因为 $\int_0^1 K_L^2=(2L^2+1)/(3L^3)$，平方中不同节点的交叉乘积非负。**猜想的充分条件**是对无穷多个 $H\to\infty$，在固定相位 $1/\pi$，有
$$
\limsup_{J\to\infty}\frac1J\sum_{j=1}^J G_{1/\pi}(2^j,H)>0 .
$$
这将为每个这样的 $H$ 提供任意晚的 $F>1/2$。坏逼近相位如 $(\sqrt5-1)/2$ 具有同样的全相位积分，而在它自身的窗口上，当 $H$ 足够大时，全部 $G\le0$。故全相位均值不能代替这条猜想。某些非坏逼近数也可能只有极稀疏的峰，直接构造稀疏峰仍足够。任务 B 可从展开 $G$ 的非零频率入手，检验在 $Q=2^j$ 的求和中哪些固定相位项能抵消。

## 有限证书与独立检查

$113/355$ 在 $Q=256$ 的窗口内。取 $H=32$、$L=2097152$，并令 $u=L|1-113\pi/355|$。Machin 公式及交错余项给出 $u<179/1000$。于是
$$
F_{1/\pi}(256,32)\ \ge\ \left(\frac{\sin u}{L\sin(u/L)}\right)^2
\ \ge\ (1-u^2/6)^2\ >\ 989/1000.
$$
这是 $355/113$ 的倒数给出的单个有限峰。

运行 `python verify_rh_farey_filter.py`。脚本用两条整数路径核对四个小窗口的全部傅里叶系数，以有理数检查相位一阶、二阶矩等式和界，并用有理数认证上述 $0.989$。直接三角求和之间的对照另用浮点计算。无穷尺度和有限证书的逻辑分别在前面证明。

## 两段可复制的任务书

下面两队各自从一条独立算术表达式出发。记号与验收基准见上文。

**A 队，尝试证明 π 坏逼近**

```text
以“pi 不是坏逼近数”为反证假设，目标是推出矛盾。设 Q=2^j，mu 是 Mobius 函数，
m_Q(d)=sum(mu(r) for r>=1 if Q<=d*r<2*Q)。
设 E_Q 为满足 Q^(3/4)<d<2*Q 且 abs(sin(d))<Q^(-3/8) 的整数集，
正弦按弧度。目标是证明 sum(m_Q(d)*d^2/sin(d)^2 for d in E_Q)
对所有充分大 Q 都 <= C*Q^4，C 与 Q 无关。笔记证明其余项已可删去。

从分圆乘积 prod(Phi_q(exp(2*i*x)) for Q<=q<2*Q) 在 x=1 的对数模二阶导出发，
或按短余因子 r<2*Q^(1/4) 分组。找一个独立估计，写出首次用到固定角 1 的步骤。
先处理一族明确的剩余项，估计对误差和常数均须一致。
如果移动复中心，必须证明能回到实点 x=1。已有指数 2 可用，但要跟踪常数。
交付有证明的下一级估计和可复算证书。若失败，指出剩余哪类二阶极点，
给出所用不等式的具体障碍。有限计算和 Q^(4+eps) 都不足以完成目标。
```

**B 队，尝试证明 π 非坏逼近**

```text
构造使 q^2*abs(pi-p/q) 趋于零的无限见证。用 alpha=1/pi，
Q=2^j，H 为正整数，L=H*Q^2。采用笔记的 Fejer 窗口 F_alpha(Q,H)。
H>=11 且 F>=1/2 就能提取一枚既约分数，倒数变换前的归一化误差 <=4/H。

令 G=F^2-F/2。从 F 的有限 Ramanujan 傅里叶式展开 G，保留固定相位 cos(2*k)。
先复现全相位平均为正的下界。再尝试证明：对无穷多个 H，
固定 alpha=1/pi 时，G(2^j,H) 对 j=1..J 的平均在 J->infinity 的 limsup 严格为正。
也可直接构造 H_m->infinity 及 Q_m，使 F(Q_m,H_m)>=1/2，允许峰极稀疏。
用黄金比相位检验哪一步确实依赖 pi。先独立复验 F(256,32)>989/1000。
交付严格的无限机制，或一类非零频率的可证明界及其余项。
有限峰只作为证书。无需假设或证明 pi 的无理性指数大于 2。
```
