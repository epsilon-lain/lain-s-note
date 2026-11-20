# 临界指数处的行列式构造

目标是控制 $\gamma(p,q)=q^2|\pi-p/q|$ 的统一正下界。《The irrationality exponent of pi is 2》用可能存在的好逼近构造非零 Gaussian 有理行列式，移回精确周期后比较实处上界与清分母后的下界。它在每个固定 $\nu>2$ 上有功率余量。临界问题要求在 $\nu=2$ 时继续保留 $-\log\gamma$。狄利克雷逼近已经保证无穷多个 $\gamma\le1$，所以只用这个不等式无法推出矛盾。[原文](https://github.com/openai/math/blob/main/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/main.tex) · [原行列式章节](https://github.com/openai/math/blob/main/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/determinant.tex)

## 实际清分母乘子

沿用原矩阵。列由 $Y^h X^\alpha$ 给出，行由 $(j,s,\beta)$ 给出。矩阵项是把 $X_i$ 换成 $2\mathrm i j p_i/q_i+G_i(t)+u_i$、$Y$ 换成 $1+t$ 后的 $t^su^\beta$ 系数，其中 $G_i$ 是 $\log(1+t)$ 截到 $t^{T_i-1}$ 的多项式。令
$$
L_i=\operatorname{lcm}(1,\ldots,T_i-1),
\qquad d_i=\operatorname{lcm}(q_i,L_i).
$$
以 $\prod_i d_i^{\alpha_i}$ 乘第 $\alpha$ 列，以 $\prod_i d_i^{-\beta_i}$ 乘第 $\beta$ 行。非零矩阵项必有 $\alpha_i\ge\beta_i$，而乘去这个差幂后，剩下的每个 $d_i(2\mathrm i j p_i/q_i+G_i(t))$ 属于 $\mathbb Z[\mathrm i][t]$。因此对一切选中的非零方子式，
$$
\log|\Delta_H|
\ge -\sum_{\rm cols}\sum_i\alpha_i\log d_i
     +\sum_{\rm rows}\sum_i\beta_i\log d_i .
$$
这保留了原来逐项乘上公共分母时没有利用的两项信息。一项是 $q_i$ 与 $L_i$ 已共用的素数幂，另一项是各列共享同一个加权次数预算。若 $w_i=\lceil\log q_i\rceil$、$T_i=\lceil F_0w_i/v_0\rceil$，用原文的 Chebyshev 常数 $\Lambda=4\log2$，每列的额外对数代价最多为
$$
(w\alpha)\,\Lambda\left(\frac{F_0}{v_0}+\frac1{w_*}\right),
\qquad w_*=\min_i w_i .
$$
原式在同一项用的是 $H\Lambda F_0m/v_0$ 和 $H\Lambda\sum_i1/w_i$。按这里的清分母，$F_0/v_0$ 项不再乘变量个数 $m$。

## 通过秩选低代价列

固定原插值定理要求的权重及中心。从低到高扫列的加权总次数 $w_0h+w\alpha$，保留增加秩的列，直至成一个非零子式。写 $M_H$ 为行数。对每个较小的可用次数 $J\le H$，次数不超过 $J$ 的列在对应的较小行集合上满秩，因而在完整行集合上的秩至少为 $M_J$。原定理允许充分大的、充分可整除的 $J$，使用一个固定间距的子序列即可。在这一顺序下贪心基的前缀大小恰好是前缀秩。用
$M_J\sim c J^{m+1}$，从 0 到 $H$ 积分行数缺额，得到
$$
\limsup_H\frac{1}{M_H H}\sum_{\rm cols}(w_0h+w\alpha)
\le\frac{m+1}{m+2}.
$$
另外对行的单纯形积分，$\bar b_H=(M_H H)^{-1}\sum_\mathrm{rows}w\beta$ 实际趋于 $\theta m/(m+2)$。这比用 $\bar b_H\le\theta$ 精确。

复算脚本用两个有理中心搭出一个秩 39 的小矩阵。选定同一个子式后，旧清分母乘子为 5,301 位，新乘子为 972 位。脚本还直接验证后者乘以有理行列式确实等于整数矩阵的行列式。实验使用实有理中心及测试权重，检查的是代数规则。上面的证明适用于原文的 Gaussian 有理中心。[精确算术实验](verify_endpoint_minor.py)

## 同时计算两种解析节省

原行列式展开中的一项，把 $M_H$ 条翻译行分配给横向指标 $a$。写 $n_a$ 为每个指标收到的行数。原文分别估计小指标碰撞和大指标误差。可同时保留其节省。若碰撞常数是 $c=(\log2)/4$，对任何 $0<A<1$，
$$
\nu\sum_a n_a wa+c\sum_a n_a^2
\ge\nu A H M_H-\frac{\nu^2}{4c}\sum_a(AH-wa)_+^2.
$$
逐个指标配方即可证明。对右边用单纯形积分，归一化后损失为
$$
\frac{\nu^2 v_0 A^{m+2}}{2c(m+2)K\theta^m}.
$$
因此可以把低代价列及两种节省一起纳入比较。

这三项改动在原解析界下仍没有闭合临界点。可证明在充分大的 $m$ 时，沿用原解析界就算免掉额外算术代价也不够。令
$g=2A-[(m+1)+m\theta]/(m+2)>0$。要求上面的配方损失小于 $g$，给出
$v_0\ll(m+2)gK(\theta/A)^m/A^2$。原解析代价含 $100K/w_0$，而原插值参数满足 $K\theta^m w_0\le v_0$，于是
$K A^m\ll (m+2)g^2/A^2$。再用解析代价还含 $(2\log2)/v_0$，必要条件变成
$$
(m+2)^2g^4(\theta/A^2)^m\gg A^4 .
$$
但 $A\ge1/3$，且
$\theta-A^2=(1+2\theta)/(m+2)-g-(1-A)^2$。左边至多是
$e^{27}(m+2)^2g^4e^{-mg}=O(m^{-2})$，右边则有正下界。这指出仅改变清分母、选列和上述求和方式，无法在原解析界下靠增大维度解决临界问题。

## 真正需要保留的临界量

原展开中某项用到 $N_i$ 次真实横向偏移
$2\mathrm i j(p_i/q_i-\pi)$，其大小的积中有**精确的**因子
$$
\exp\left(-2\sum_iN_i\log q_i\right)
\prod_i\gamma(p_i,q_i)^{N_i}.
$$
后一个乘积正是单用指数 $2$ 会丢掉的量。若能把完成 $q_i^{-2}$ 的主项抵消后余下的损失压到
$C\sum_iN_i$，并为完全不使用真实横移的项取得独立的严格碰撞余量，便有了用一个固定的极小 $\gamma$ 阈值作反证的具体路径。这里的 $C$ 必须不随所选分母增大。可以把所需结论写得更具体。设 $D_H$ 是这里实际的清分母乘子，$U$ 是某项的全部组合与截断系数的绝对值，并吸收 $|(2j)^{N_i}|$，但剥去上式显示的误差因子。若对每项都能给出
$$
\log D_H+\log U+\log|\det T|-2\sum_iN_i\log q_i
\le C\sum_iN_i-\sigma M_HH,
$$
其中 $C$ 跨分母统一、$\sigma>0$ 对固定几何参数和中心即可，那么选所有中心的 $\gamma_i<e^{-C}$ 后，每一项都至多是 $D_H^{-1}e^{-\sigma M_HH}$。原文的展开项数只贡献 $o(M_HH)$ 的对数，求和便与下界 $D_H^{-1}$ 矛盾。$N_i$ 全为零的项也必须达到式中的严格余量。本笔记尚未建立这种解析控制。

原解析界里需处理的一项是 $100K/w_0$。中心 $j\omega$ 位于虚轴，$\omega=2\pi\mathrm i$，而 $e^{h j\omega}=1$。把行展开成纯多项式算子可确切保留这一信息：
$$
[t^\ell](1+t)^h(j\omega+\log(1+t))^d
= (\partial_h+j\omega)^d\binom h\ell .
$$
右边对 $h$ 求导的是次数 $\ell$ 的二项式多项式。由对 $h$ 求导自然拉出 $\log(1+t)$，式子逐项成立。下一项具体问题是能否用这种作用在 $h,j$ 上的多项式结构估计整个行块的行列式，同时保留重复指标带来的消失。逐元素界很容易利用虚轴上的周期性，全块界仍需证明。

后续构造了原模型中的精确反例，排除了单靠虚轴周期性保留原阶数碰撞、同时去掉线性中心代价的统一上界。按横向指标和切向阶数细分的估计可以严格建立。[反例与细分组估计](PERIOD_JET_COLLISION_COUNTEREXAMPLE.zh.md)
