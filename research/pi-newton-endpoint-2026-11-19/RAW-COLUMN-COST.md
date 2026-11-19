# 原始单项式选列能否完全避开最后一个逼近分数

本节只讨论非零全行子式的**名义最后次数**
$$
 C_m=\sum_{\rm columns}\alpha_m-\sum_{\rm rows}\beta_m.
$$
符号沿用本目录主文。行 $(j,s,\beta)$ 满足 $v_0s+w\cdot\beta/\theta<H$，原始单项式列 $Y^hX^\alpha$ 满足 $w_0h+w\cdot\alpha\le H$。$M(H)$ 是总行数。下文取 $m\ge2$、权重固定，并把各中心的不同性及相应插值门槛作为前提。这个判别不涉及约分后的真实分母。

记删除最后一个空间坐标后的两个体积比为
$$
 V=K\theta^{m-1},\qquad U=K(w_0/v_0)\theta^{m-1}.
$$
若 $V,U<1$，在为较低维插值重新满足下述权重门槛后，所有充分大的 $H$ 都存在 $C_m=0$ 的原始非零全行子式，其数值与 $r_m$ 无关。若 $V>1$ 或 $U>1$，固定权重下存在 $c>0$，使一切充分大的 $H$、一切原始非零全行子式均满足 $w_m C_m\ge c H M(H)$。这两种结论不给出等号边界的分类。

**充分性。** 选接近 1 的有理数 $\lambda\in(\theta,1)$，使较低维数据 $(\theta/\lambda,\lambda v_0)$ 的两个体积比仍严格小于 1。选各空间权重时，同时满足原问题与这套低维数据的逐次门槛。按 $b=\beta_m$ 分组，子包高度为 $T_b=H-w_m b/\theta>0$。只选 $\alpha_m=b$ 的列来构造对应的对角块。原矩阵在行 $\beta_m>\alpha_m$ 上为零；在 $\beta_m=\alpha_m=b$ 上取系数会完全耗尽 $X_m^b$，留下低维插值矩阵。

设加强后的低维插值可用于所有充分大的 $R$ 的倍数 $D$，$R$ 固定。取这样的 $D\ge\lambda T_b$，且 $D\le\lambda T_b+O(1)$，其中常数包括插值的起始高度。其行目标包括原权重下高度 $T_b$ 的所有低维系数。剩余列预算满足
$$
 H-w_m b-\lambda T_b=(1-\lambda)H+w_m b(\lambda/\theta-1)\ge(1-\lambda)H.
$$
因此 $D\le H-w_m b$ 对所有可行的 $b$ 同时成立。逐块选出非零方阵，沿对角块相乘即得所述原始子式。这个证明也解释了删除最后坐标的代价，行列式中对该近似的依赖完全消失。对原论文的参数 $K(w_0/v_0)\theta^m=1/2$，若 $\theta>1/2$ 且 $K\theta^m<\theta$，就处在这一侧。

**必要的层数开销。** 令 $S_k$ 为 $\beta_m<k$ 的行数，$A_k$ 为所有 $\alpha_m<k$ 的列在完整行包上的秩。这些列仅在 $S_k$ 行上有支撑。任一独立选列最多取其中 $A_k$ 列，故
$$
 C_m\ge\sum_{k\ge1}(S_k-A_k).
$$
这是逐层表示 $\sum\alpha_m-\sum\beta_m$ 得到的，每层差非负。

若 $U>1$，令 $n=m+1$、$\rho=K(w_0/v_0)\theta^m$。用所有 $\alpha_m<k$ 的列数上界 $A_k$，对权重单纯形取格点渐近，得到
$$
 \liminf_{H\to\infty}\frac{w_mC_m}{HM(H)}
 \ \ge\int_0^\theta\!
 \big[1-(1-z/\theta)^n-\rho^{-1}\{1-(1-z)^n\}\big]_+\,dz.
$$
这里 $z=w_m k/H$。具体地，$S_k=M(H)-M(H-kw_m/\theta)$，允许列总数为 $W(H)$ 时，低次数列数为 $W(H)-W(H-kw_m)$。固定正权重下的格点计数为体积主项加 $O(H^{n-1}+1)$，可用每个格点的单位立方体作内外夹逼，故求和后误差比 $HM(H)$ 低一阶。方括号内函数在零点的右导数是 $n(1/\theta-1/\rho)>0$，因为 $U=\rho/\theta>1$。积分严格为正。

若 $V>1$，选充分小的 $\varepsilon>0$，使 $K(\theta-\theta\varepsilon-\varepsilon)^{m-1}>1$。取 $L=\lfloor\varepsilon H/v_0\rfloor$ 和 $\varepsilon H/(2w_m)\le k\le\varepsilon H/w_m$ 的整数。先只看 $s<L$、$\beta_m<k$ 的小包，行数至少为
$$
 K Lk\,\#\{\beta'\in\mathbb N_0^{m-1}: w'\cdot\beta'<\theta H-\theta v_0 L-w_m k\}.
$$
固定空间列次数 $\alpha$ 后，改变 $h$ 只会乘上 $(1+t)^h\bmod t^L$，故在所看小包上这些列的秩至多为 $L$。所有 $\alpha_m<k$ 的列在此小包上的秩因此不超过
$$
 Lk\,\#\{\alpha'\in\mathbb N_0^{m-1}: w'\cdot\alpha'\le H\}.
$$
两数的差至少是 $c_1 Lk H^{m-1}$，其中 $c_1>0$ 与 $H,k$ 无关。这是同一个空间单纯形的两次格点估计，由 $\varepsilon$ 的选择保证主项之差为正。低次数列在大包 $\beta_m<k$ 上的余维 $S_k-A_k$ 至少是其在所看小包上的余维。在给定的整个 $k$ 区间求和，使用 $L\asymp H$、$k\asymp H/w_m$ 及固定权重下 $M(H)\asymp H^{m+1}$，即得 $w_m C_m\ge c HM(H)$。

称为非零的子式必有某个非零配对，因此各坐标 $C_i\ge0$。这里证明的是其中最后坐标的开销；解析展开中的真实误差阶、变换后子式的整除与代入后约分仍需分别估计。

[OpenAI 的插值定理](https://github.com/openai/math/blob/main/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/interpolation.tex) · [仓库先前的低成本选列与较强参数下的 Hermite 构造](../contradiction-rank-2026-11-01/)
