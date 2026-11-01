# 同一个非零子式怎样少付分母成本

日期：2026-11-01。本文从固定 π017 源码逐式重算，讨论可以合法改写的估计，不认证 π 的主结论，也不宣布坏逼近性结论。

最高分位列费用与精确行均值在[此前仓库判据](../endpoint-joint-bound-2026-10-31/README.md)已有；这里重述以接入新的行列整数化和秩积分选基。有关尾项费用的另一项独立改写见 [tail_elimination.md](tail_elimination.md)，本文第1节仅指整数化本身不会自动消尾。

源固定为 `openai/math@4eaee87548250b2682affda7d05d947c05554bae`。定位：

- [插值定理，interpolation.tex 24–57 行](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/interpolation.tex#L24)。
- [行列定义、非零子式与行均值，determinant.tex 45–81 行](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/determinant.tex#L45)。
- [原整数化，determinant.tex 99–149 行](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/determinant.tex#L99)。
- [原分支比较，determinant.tex 345–431 行](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/determinant.tex#L345)。

## 1. 可直接证明的整数化改写

固定一个使用全部行的非零子式，其行是 `(j,s,β)`，列是 `Y^h X^α`。定义

\[
 c_H=\frac{\sum_{\text{列}}w\alpha}{MH},\qquad
 b_H=\frac{\sum_{\text{行}}w\beta}{MH}.
\]

原稿在 128 行把每个列成本 `wα` 都按最大值 `H` 收费，即用了 `c_H≤1`。这并不等于实际成本。

令 `T_i=ceil(F_0w_i/v_0)`，`L_i=lcm(1,...,T_i−1)`（空集合取 1）。把每个列乘以 `∏(q_iL_i)^{α_i}`，每个行乘以 `∏(q_iL_i)^{−β_i}`。矩阵元素在 `α≥β` 时成为

\[
 \prod_i\binom{\alpha_i}{\beta_i}
 [t^s](1+t)^h
 \prod_i\{L_iq_ijr_i+L_iq_iG_i(t)\}^{\alpha_i-\beta_i}.
\]

这里 `q_i jr_i=2ijp_i∈Z[i]`，`L_iG_i(t)` 的系数为整数，因此每个元素都是高斯整数；在某一坐标 `α_i<β_i` 时元素为零。缩放后的行列式仍然非零，故其模至少为 1。

更关键的是，对每个坐标

\[
 D_i:=\sum_{\text{列}}\alpha_i-\sum_{\text{行}}\beta_i\ge0.
\]

证明：非零行列式的 Leibniz 展开至少有一个非零乘积；对应的匹配里每个非零元素均满足 `α_i≥β_i`，相加即得。这不是假定每个展开项非零。

于是精确算术下界是

\[
 \log|\Delta_H|\ge-\sum_iD_i\log(q_iL_i).
\]

使用 `log q_i≤w_i`，以及 `Λ=4 ln 2` 时

\[
 \log L_i\le\Lambda(T_i-1)
 <\Lambda F_0w_i/v_0,
\]

得到改进的有效下界

\[
 \boxed{\quad
 \frac{\log|\Delta_H|}{MH}
 \ge-\left(1+\frac{\Lambda F_0}{v_0}\right)(c_H-b_H),
 \qquad c_H\ge b_H.
 \quad}
\]

其效果来自**分母按列与行的指数差收费**：行上的微分已经消耗掉的 `β_i` 不必再次付费。它同时消除了原来给每个变量单独预留的 `ΛF_0m/v_0`，并避免了原证明中把 `log q_i` 换成 `w_i` 时的两个额外损失。没有忽略截断多项式的分母；这些分母由 `L_i^{α_i−β_i}` 精确清除。

这条重算只改下界。它没有自动加强非零性定理，也没有消除解析平移中的尾项成本 `ν/F_0`。

## 2. 合法选择更便宜的同一个子式

原插值定理只保证存在满行秩子式。可以在保证非零的子式中选择 `Σ列 wα` 最小者；这是线性代数中的最小费用列基。

令 `r_H(x)` 为保留所有行、仅保留 `wα≤xH` 的列之后的秩，`0≤x≤1`。将列按 `wα` 升序排列，逐列保留增加秩的列，直到取满 `M` 列。其所得子式非零，而且费用最小。严格恒等式为

\[
 \boxed{\quad c_{H,\min}
 =\int_0^1\left(1-\frac{r_H(x)}M\right)\,dx.\quad}
\]

证明：任一基在阈值 `xH` 以上至少需要 `M−r_H(x)` 列；贪心基在每个阈值都恰好满足这个下界。将“超过阈值的列数”积分，恰好恢复列费用总和。这不是把行秩近似成列数。

还有一个能由原插值定理直接推出的渐近界。固定全部参数，再让 `H→∞`。对任意固定 `0<x<1`，投影到度数小于 `xH` 的行，并使用总度数不超过 `xH` 的列，原定理使这个较小矩阵最终满行秩。这些列都包含在 `wα≤xH` 中。因此

\[
 \liminf_{H\to\infty}\frac{r_H(x)}{M_H}\ge x^{m+1}.
\]

对有理 `x` 取适当可整除度数；一般 `x` 由单调性逼近。被积函数有界，Fatou 引理给出

\[
 \boxed{\quad
 \limsup c_{H,\min}\le\frac{m+1}{m+2}.
 \quad}
\]

行单纯形的精确均值则为 `b_H→mθ/(m+2)`，所以可保证

\[
 \limsup(c_{H,\min}-b_H)
 \le\frac{1+m(1-\theta)}{m+2}.
\]

这比给所有列都按 `H` 收费更细，但仍只是上界。若希望进一步逼近必要下界 `c_H≥b_H`，具体缺口已变成：**证明这个特殊对数插值矩阵在每个横向截断处的秩 `r_H(x)` 有多大**。

## 3. 不依赖选基算法的另一个上界

全部列在度数单纯形中，横向费用 `x=wα/H` 的渐近密度为

\[
 f_m(x)=m(m+1)x^{m-1}(1-x),\quad0<x<1.
\]

因此它服从 `Beta(m,2)` 分布。全列平均费用为 `m/(m+2)`，但**不能把这个全列平均直接套给选中的基**：大量便宜列可能线性相关。一般向量配置中这个推论是假的。

令 `ρ=lim M/N=K(w_0/v_0)θ^m∈(0,1)`。任意 `M` 个不同列的平均费用都不超过全部列中最贵 `ρ` 比例的平均。若 `t` 满足

\[
 1-(m+1)t^m+mt^{m+1}=\rho,
\]

那么这个上界为

\[
 U_m(\rho)=\frac{m/(m+2)-mt^{m+1}
       +m(m+1)t^{m+2}/(m+2)}{\rho}.
\]

所以贪心基还满足 `limsup c_min≤min{(m+1)/(m+2),U_m(ρ)}`。在 `ρ=1/2` 时：

| m | 最贵一半列的平均费用 U | 秩积分所给界 |
|---:|---:|---:|
| 1 | 0.5285954792 | 2/3 |
| 2 | 0.6875000000 | 3/4 |
| 10 | 0.9142676242 | 11/12 |
| 50 | 0.9814307181 | 51/52 |

这两个界来源不同，取较小者即可。

## 4. 两类改动为何未自动提供新机制

**只缩小横向单纯形。** 若把列域改成 `w_0h+wα/ζ≤H`，乘以 `ζ` 后就是 `ζw_0h+wα≤ζH`；行条件同时变成 `ζv_0s+wβ/(θ/ζ)<ζH`。当 `θ/ζ<1` 且相应体积条件成立时，这仍然属于原来的参数族。单看成本由 `1` 变成 `ζ` 不能宣布获得新收益，必须同时重算度数、行权重、插值体积和解析增长。上面的最小费用基则是对固定原矩阵做真实优化。

**各方向使用不同 θ_i。** 如果假设已经证明相应各向异性插值定理，行数中出现 `∏θ_i`，行平均成本变成 `Σθ_i/(m+2)`。固定这个平均成本时，AM–GM 给 `∏θ_i≤(Σθ_i/m)^m`。因此把 `θ_i` 分散本身不会在相同平均成本下增加行体积或碰撞类别收益。它可能改变几何秩条件，但那需要一个新的插值定理，不能仅凭维数计数声称成功。此处没有证明各向异性插值定理。

## 5. 精确实验及其边界

运行 `python3 research/pi-endpoint-rewrite-2026-11-01/selected_columns_lab.py`。

程序复用先前已核对的 Taylor 系数生成器，对 `m=K=1,W=(1,1),V=(1/2,4)`、中心 `(1,0)` 的六个有限矩阵计算精确有理数秩；它们均满行秩。按横向费用排序后取 RREF 主元列，独立核对离散秩积分恒等式、`c_H≥b_H`。例如 `H=5` 时有 12 行、21 列，最小基费用 `c_H=7/60`，行均值 `b_H=1/30`，实际差为 `1/12`。这些数值展示了“按最高列费用收费”的浪费确实可以很大。

这六个矩阵是有限模型；它们未满足并验证源定理全部“大权重”门槛，不能替代源定理的渐近证明。程序还用有理二分给出六组 `U_m(1/2)` 的精确包围；显示用小数不是严格区间本身，完整有理上下界在 `selected_columns_results.json`。

现在值得继续证明的对象是秩函数与指数差费用，而不是反证法这种逻辑形式。最终要把改进下界和同一个非零子式的解析上界一起比较，保留全部尾项及取极限顺序。
