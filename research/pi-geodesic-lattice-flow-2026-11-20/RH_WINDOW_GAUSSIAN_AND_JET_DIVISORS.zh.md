# 从 RH 的窗口判据到 π 的高斯格点和与喷流子式

以用户给出的正高斯格点和及共同子式因子为起点。这里证明精确的固定相位恒等式、可用的滑窗上界判据、幂窗口移植的反例，以及原插值矩阵的一个显式共同因子。圆周率是否坏逼近仍待判定。

## RH 窗口的量词

在 [Lean 定理](https://github.com/SINF-KEN/rh-lean-check/blob/40e56bab/src/Development/AutomorphicArithmetic/DedekindLFunctions/LandauTheorems/ZeroFreeOfRpowWindow.lean#L11-L14) 中，固定 $0<\theta<1$、$\eta>0$ 和常数 $C$。若任意晚的地方都能找到起点 $Y$，使整个 $[Y,Y^{1+\eta}]$ 上的 $G_m(F,y)\le Cy^\theta$，便推出对应的零点排除。逆否命题是，若排除失败，则每个足够晚的此类窗口都有至少一点超界。

[误差预算](https://github.com/SINF-KEN/rh-lean-check/blob/40e56bab/src/Formalism/SizeWindowBrackets/GrowthBoundGivesRiemannHypothesis.lean#L109-L122) 在归一化后含有 $x/m$ 和 $m$，还要控制包含 $k$ 的项及规模因子，因此正式假设给出 $x^{1/2-\varepsilon}\le m\le x^{1/2+\varepsilon/2}$。[共同格](https://github.com/SINF-KEN/rh-lean-check/blob/40e56bab/src/Formalism/LocalToGlobal/IndefiniteForms/IndefSide.lean#L273-L289) 是高度一局部格的交，之后再与无穷处的度量一起使用。原路径出现重定向，上述源码从保留相同提交的 GitHub 镜像读取。MIT 的核验页在北京时间 2026-11-21 06:05 记录了 175187 个模块全量编译和模块级复演、标准公理检查，以及两个独立内核对完整导出的接受。原生 Lean 对整份导出的复演和整套 comparator 在该时点仍各有运行中的检查。[MIT 核验记录](https://math.mit.edu/~drew/grh/)

## 高斯峰的形状与可用的滑窗

令 $\alpha$ 为无理实数，$X\ge1$，定义
$$
L_\alpha(X)=\{(q/X,X(q\alpha-p)):p,q\in\mathbb Z\},
\qquad
\Theta_\alpha(X)=\sum_{v\in L_\alpha(X)}e^{-\pi\|v\|^2}.
$$
格的余体积为 1。若 $p/q$ 既约、$q>0$，写 $\delta=q^2|\alpha-p/q|$，$X_*=q/\sqrt{\delta}$，$u=\log(X/X_*)$，对应原始格向量的长度为 $r$，则
$$
r^2=2\delta\cosh(2u).
$$
将坐标旋转后取行列式为 1 的基 $(r,0),(r\tau,1/r)$。对第一基向量的整数倍用一维 Poisson 求和，得到
$$
r\Theta_\alpha(X)=\sum_{a,b\in\mathbb Z}
 e^{-\pi(a^2+b^2)/r^2}\,e^{2\pi iab\tau}
 =1+4e^{-\pi/r^2}+O(e^{-2\pi/r^2})\quad(0<r\le1).
$$
余项一致，因为除原点和四个平方和为 1 的点外，其余平方和至少为 2。单独保留原始向量的整数倍，再用一维 Poisson，给出 $\Theta_\alpha(X)\ge1/r$。当 $0<\delta<1$ 且 $q\delta^{-1/4}\le X\le q\delta^{-3/4}$，有 $r^2\le\sqrt\delta(1+\delta)$。这既给出用户的全窗下界 $\delta^{-1/4}/\sqrt2$，也把整个窗内的真实总和确定到相对误差 $O(e^{-\pi/(\sqrt\delta(1+\delta))})$。若下一部分商为 $A$，相应上一收敛分数的 $\delta<1/A$。取记录中的 $A=845188380288$，中心高度严格超过 $\sqrt{A/2}>650072$。[部分商记录](https://oeis.org/A033089/list)

对任意 $X>0$ 还有精确的不等式
$$
e^{-|h|}\Theta_\alpha(X)\le \Theta_\alpha(e^hX)\le e^{|h|}\Theta_\alpha(X).
$$
证明：二维余体积 1 格的对偶等于原格旋转 $90^\circ$ 的像。固定 $X$ 后设 $F(t)=\sum_v e^{-\pi t\|v\|^2}$，Poisson 给出 $F(t)=t^{-1}F(1/t)$，在 $t=1$ 求导得到 $2\pi\sum_v\|v\|^2e^{-\pi\|v\|^2}=\Theta_\alpha(X)$。沿 $\log X$ 求导则得到 $2\pi\sum_v(v_1^2-v_2^2)e^{-\pi\|v\|^2}$，绝对值不超过 $\Theta_\alpha(X)$。逐项求导在紧时间区间内成立，因为格的范数在该区间一致等价。积分对数导数即得不等式。特别地，
$$
\Theta_\alpha(X)\le 2\int_X^{2X}\Theta_\alpha(Y)\frac{dY}{Y}.
$$
因此，对每个 $X\ge1$ 的倍尺度滑窗积分给出同一上界，就足以证明坏逼近性。其逆向来自逐点上界。逐点上界与坏逼近的等价性如下：坏逼近下界让每个 $q\ne0$ 的格向量长度平方至少为 $2|q||q\alpha-p|$，$q=0$ 的非零格向量长度至少为 1。统一的最小间距结合圆盘装填使高斯和一致有界。不存在坏逼近下界时，取 $\delta\to0$，上述中心峰高发散。这里尺度限制为 $X\ge1$，因为 $X\downarrow0$ 时 $q=0$ 的一列就会抬高高斯和。

### 任意晚的好幂窗口仍然可以漏掉峰

取 $\alpha=[0;a_1,a_2,\ldots]$，除 $n_j=2^{2^j}$ 处令 $a_{n_j}=j+1$ 外，部分商均为 1。该数非坏逼近且无理性指数为 2。后一性质可由 $\log a_{n+1}=o(\log q_n)$ 与收敛分数误差公式推出。这里 $q_n$ 至少按斐波那契数增长。

若 $a_{n+1}=1$ 且 $q_n\le X\le q_{n+1}$，任一非零格向量至少有一个坐标绝对值 $\ge1/2$。当 $|q|\ge q_n$，用 $|q|/X\ge q_n/q_{n+1}\ge1/2$。当 $0<|q|<q_n$，第二类最佳逼近给出 $|q\alpha-p|\ge\|q_{n-1}\alpha\|>1/(q_n+q_{n-1})$，结合 $X\ge q_n$。$q=0$ 可直接检查。因此整个 $[q_{n_j},q_{n_{j+1}-1}]$ 的高斯和具有同一个上界。又 $\log q_{n_j}=O(n_j)$，$\log q_{n_{j+1}-1}\gg n_{j+1}$，后者除以前者趋于无穷。于是对每个固定的 $\eta>0$，存在任意晚的 $[Y,Y^{1+\eta}]$，其中 $\Theta_\alpha$ 一致有界。这个例子具有指数 2，仍会在别的尺度产生任意高峰。它说明针对 π 的证明需要补充固定相位的估计。

## 固定相位的两条波与一个加权除数和

对 $p$ 做 Poisson 求和，得到绝对收敛的恒等式
$$
X\Theta_\alpha(X)=\sum_{q,k\in\mathbb Z}
 e^{-\pi(q^2+k^2)/X^2}e^{2\pi i\alpha qk}.
$$
取 $q+k,q-k$，它们同奇偶。按奇偶拆开后得到
$$
X\Theta_\alpha(X)=|A_0(X)|^2+|A_1(X)|^2,\qquad
A_j(X)=\sum_{n\in\mathbb Z}
 e^{-2\pi(n+j/2)^2/X^2+2\pi i\alpha(n+j/2)^2}.
$$
取 $\alpha=1/\pi$，刚好回到已有波形的固定相位 $e^{2in^2}$，并出现其半整数平移。坏逼近性在取倒数后等价，因为对逼近 $\pi$ 的既约 $p/q$，有 $p^2|1/\pi-q/p|=(p/q)/\pi\cdot q^2|\pi-p/q|$，而且逆向同理。所以能对这两个固定相位和一起证明 $O(\sqrt X)$，就能完成 π 的上界。这里的平方和恒等式与此前基于四次积分的比较是两个不同的观测公式。[原波形笔记](FOUR_BUNDLE_REVIEW.zh.md)

也可按 $qk=n$ 合并，令 $b_X(n)=\sum_{d\mid n,\,d>0}e^{-\pi(d^2+(n/d)^2)/X^2}$，则
$$
X\Theta_\alpha(X)=2\sum_{k\in\mathbb Z}e^{-\pi k^2/X^2}-1
 +4\sum_{n\ge1} b_X(n)\cos(2\pi n\alpha).
$$
因此具体的算术目标是 $\sum_{n\ge1}b_X(n)\cos(2n)=O(X)$。只需上侧的 $O(X)$，因为 $\Theta\ge0$ 与轴上项 $O(X)$ 已给出统一的 $-O(X)$ 下界。使用 [OpenAI 的指数 2 定理](https://github.com/openai/math/blob/main/lean/docs/017.md) 可得 $\Theta_\pi(X)=O_\varepsilon(X^\varepsilon)$。例如对 $0<\varepsilon<1$ 分成 $|q|\ge X^{1-\varepsilon}$ 和 $|q|<X^{1-\varepsilon}$，先在前者读横坐标，再用指数 $2+\varepsilon$ 的逼近下界在后者读纵坐标。最短长度 $\gg_\varepsilon X^{-\varepsilon}$，由前节的峰公式或装填的常数区间情形得到上界。两条波的临界上界仍待证明。

### 对相位取平均的精确障碍

把 $\alpha$ 在 $[0,1]$ 上均匀平均。均值是轴上项除以 $X$，趋于 2。另一方面，
$$
\int_0^1\left(\Theta_\alpha(X)-\int_0^1\Theta_\beta(X)\,d\beta\right)^2d\alpha
=\frac{6}{\pi}\log X+O(1).
$$
由正交性，方差等于 $8X^{-2}\sum_{n\ge1}b_X(n)^2$。在平方和里将 $qk=q'k'$ 写成 $q=ad,\ q'=bd,\ k=bc,\ k'=ac$，其中 $a,b,c,d>0$ 且 $\gcd(a,b)=1$。设 $R=(a^2+b^2)^{1/2}$，得到
$$
\sum_{n\ge1}b_X(n)^2
=\sum_{\gcd(a,b)=1}\left(\sum_{c\ge1}e^{-\pi R^2c^2/X^2}\right)^2.
$$
当 $R\le X$，括号为 $X/(2R)+O(1)$，由单调函数与积分比较。平方后的总误差为 $O(X^2)$，因为 $\sum_{R\le X}R^{-1}=O(X)$。当 $R>X$，括号为 $O(e^{-\pi R^2/X^2})$，这段总和也是 $O(X^2)$。用 Möbius 反演和四分之一圆的格点计数 $N_0(T)=\pi T^2/4+O(T+1)$，互素对数是 $\pi T^2/(4\zeta(2))+O(T\log(2T))$。对 $R^{-2}$ 分部求和，得 $(\pi/(2\zeta(2)))\log X+O(1)$。代回即证。该方差来自对相位的整体统计，无法提供固定相位的统一二次矩上界。

## 原 Gaussian 整数喷流的共同因子

沿用 [临界行列式笔记](DETERMINANT_AT_THE_CRITICAL_EXPONENT.zh.md) 的清分母方式。记 $V_i=w_i/\theta$，中心序号为 $j=0,\ldots,K-1$，保留 $v_0s+\sum_i V_i\beta_i<H$ 的行。清分母后，第 $i$ 个中心增量为 $c_i=2\mathrm i\,p_i d_i/q_i\in\mathbb Z[\mathrm i]$。每列的条目是一个系数抽取
$$
[t^s u^\beta]\,P(t,jc+B(t)+u),
$$
其中 $P$ 与 $B(t)$ 均为 Gaussian 整系数多项式，并且 $B$ 与中心序号 $j$ 无关。全行子式指取全部这些行、再选同等数目的列。

固定 Gaussian 素元 $\varpi$，选 $i$ 使 $v_\varpi(c_i)=\min_\ell v_\varpi(c_\ell)$。在 $\varpi$ 处的局部环，$c_j/c_i$（此处 $j$ 为坐标）均可使用。令 $W=\max_\ell V_\ell$。留下
$$
v_0s+W\beta_i+\sum_{\ell\ne i}V_\ell\beta_\ell<H
$$
的行作为子包。在这组行上做剪切 $u_i=v_i,\ u_\ell=v_\ell+(c_\ell/c_i)v_i$。系数变换以 $\beta_i$ 为序三角且对角元为 1。因为给第 $i$ 坐标的权重为最大值 $W$，用于计算新系数的所有旧行都还在子包里。把它补成对全部行的可逆变换后，每个全行子式的局部赋值不变。新坐标下，中心仅沿第 $i$ 方向变化，剪切后的 $B(t)$ 仍与中心序号无关。

按固定 $(s,\beta_{\ne i})$ 分组，设
$$
k=\left\lceil\frac{H-v_0s-\sum_{\ell\ne i}V_\ell\beta_\ell}{W}\right\rceil>0.
$$
每列在这组的条目等于同一个一元多项式在 $0,c_i,\ldots,(K-1)c_i$ 的前 $k$ 个 Taylor 系数。先对此多项式除以首一的 $\prod_{j=0}^{K-1}(Z-jc_i)^k$ 取余，系数仍在局部环内。故这组行从一个大小 $Kk$ 的合流 Vandermonde 矩阵因子中通过。它的行列式为
$$
B_i^{k^2},\qquad B_i=\prod_{0\le j<\ell<K}c_i(\ell-j).
$$
可用首一 Newton 基 $(Z-jc_i)^a\prod_{r<j}(Z-rc_i)^k$ 验证。按 $(j,a)$ 排序的取系数矩阵是三角阵，其第 $(j,a)$ 个对角元为 $\prod_{r<j}(jc_i-rc_i)^k$，其中 $0\le a<k$。把所有分组组成块对角因子，得到每个非零全行子式的共同赋值下界
$$
v_\varpi(\Delta_{\mathrm{int}})
\ \ge\ v_\varpi(B_i)\sum_{s,\beta_{\ne i}:\,k>0} k^2.
$$

固定中心和权重，令 $M_H$ 为全部行数。单纯形求和给出
$$
\frac{\sum_{s,\beta_{\ne i}}k^2}{M_H H}\ \longrightarrow\
\frac{2V_i}{(m+2)K W^2}.
$$
分子主项为 $2H^{m+2}/[(m+2)!v_0 W^2\prod_{\ell\ne i}V_\ell]$，行数主项为 $K H^{m+1}/[(m+1)!v_0\prod_\ell V_\ell]$。例如 $K\ge2$ 时，$\varpi=1+\mathrm i$ 已给出正的主阶收益，因为全部 $c_i$ 被 2 整除。

算术端可以把上述共同因子从实际清分母的费用中减去。对这个具体证书，归一化收益总和至多为 $2\log|D|/[(m+2)K W]$，其中 $D=(\gcd_i c_i)^{K(K-1)/2}\prod_{j<\ell}(\ell-j)$。具体地，若所选子式的实际清分母成本为 $C_H$，而这次证出的共同因子为 $g_H$，则原有理子式的对数下界变为 $-C_H+\log|g_H|$。理由是 $V_i\le W$，且各局部选择正好取到 $c_i$ 的最小赋值。固定其余参数时，$\log|\gcd_i c_i|\le\log|c_i|=O(w_i)$ 对每个固定坐标成立，来自 $d_i/q_i\mid\operatorname{lcm}(1,\dots,T_i-1)$。因此分母权重高度分离时，该证书的收益允许趋零。真正的共同因子可能更大。尚未得到能跨所有候选分母抵消解析亏损的下界。

[复算脚本](verify_rh_pi_window_transfer.py) 对同一截断分别按乘积与按互素对计算方差，并用整数消元穷举一个 $8\times10$ 喷流矩阵的全行子式。有限矩阵的 gcd 为 64。方差的浮点检查仅用于对照，渐近结论来自上面的推导。
