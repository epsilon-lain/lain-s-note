# 保留同一个子式的改写估计

2026-11-01。固定 π017 源码版本为 `openai/math@4eaee87548250b2682affda7d05d947c05554bae`。本文是一组独立重推的行列式估计，以上游插值定理能提供所述非零子式为前提；没有独立认证该几何定理或 π 的主结论。

符号与矩阵见 [determinant.tex 45–81](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/determinant.tex#L45)。虚数单位写作 `i`，横向坐标指标为 `i=1,…,m` 时由上下文区分。`H` 最后沿插值定理允许、充分可整除的次数趋于无穷；维数、权重、近似、中心以及下述两个测试半径均先固定。每个 H 可以选择不同的非零列基，但同一个 H 的全部比较共用该列基。

## 1. 同一个子式与精确算术下界

取全部行 `(j,s,β)`，其中 `0≤j<K`、`v₀s+wβ/θ<H`；固定 `M` 个列 `Y^h X^α` 构成一个非零子式 `Δ_H`。令

\[
 D_i=\sum_{\rm columns}\alpha_i-\sum_{\rm rows}\beta_i,
 \quad
 \chi_H=\frac{\sum_{\rm columns}w\alpha}{MH},
 \quad b_H=\frac{\sum_{\rm rows}w\beta}{MH},
 \quad \zeta_H=\frac{w_0\sum_{\rm columns}h}{MH}.
\]

非零元需要 `α≥β`，非零行列式至少有一条非零置换匹配，故每个 `D_i≥0`，且 `χ_H≥b_H`、`χ_H+ζ_H≤1`。

令 `T_i=ceil(F₀w_i/v₀)`、`L_i=lcm(1,…,T_i−1)`，空集取 1。按列乘 `(q_iL_i)^{α_i}`、按行乘 `(q_iL_i)^{−β_i}` 后每个元素都在 `Z[i]` 中。因此

\[
 \boxed{\quad \ln|\Delta_H|\ge-\mathcal A_H,
 \qquad \mathcal A_H=\sum_iD_i\ln(q_iL_i).\quad} \tag{1}
\]

一个便于估计的推论是

\[
 \frac{\mathcal A_H}{MH}
 \le\left(1+\frac{4\ln2\,F_0}{v_0}\right)(\chi_H-b_H). \tag{2}
\]

证明及选取较便宜非零基的方法见 [geometry_rewrite.md](geometry_rewrite.md)。实际推导优先保留 (1)，不必立刻用 (2)。

## 2. 完整行包使截断尾项从解析估计中精确消失

截断矩阵记 `A_G`，把各 `G_i(t)` 换成完整 `ln(1+t)` 后记 `A_log`。若 `v₀T_i>w_i/θ`，则存在按行权重排列的单位下三角矩阵 `U`，满足

\[
 A_G=U A_{\log},\qquad \det U=1.
\]

因此对同一组选中列

\[
 \boxed{\Delta_G=\Delta_{\log}}. \tag{3}
\]

这是有限系数空间中的精确恒等式，不是尾项数值小的近似。完整行包与高权重尾项保证行变换封闭；某些较小但仍封闭的行集也可使用，任意删行后不能直接沿用。

令 `ω=2πi`、`ε_{j,i}=j(2ip_i/q_i−ω)`。分析端直接使用 (3)，便得到

\[
 \mathrm{row}_{j,s,\beta}
 =\sum_{a\ge\beta}\binom a\beta\,
       \epsilon_j^{a-\beta}
       \big([t^s]f_{a,P}(j\omega+\ln(1+t))\big)_P,
 \quad f_{a,P}(z)=\binom\alpha a e^{hz}z^{|\alpha|-|a|}. \tag{4}
\]

这里求和可限制到 `a≤α` 对至少一个选中列成立的指标。其余指标的整行测试函数恒为零。`j=0` 时 `ε_j=0`，只留下 `a=β`。没有尾指标 `d`，也没有原证明的 `ν/F₀` 与尾系数的 `ln2/v₀` 损失。算术端仍用小分母的 `A_G`。详证和真实非平凡尾项模型见 [tail_elimination.md](tail_elimination.md)。

## 3. 按实际列次数收费的碰撞界

固定实数 `0<r<1` 及

\[
 S=2\pi(K-1)-\ln(1-r),\quad R>\max\{1,S\},\quad
 \gamma=S/R,\quad c_R=\tfrac12\ln(R/S)>0.
\]

`r` 是 `t` 系数的 Cauchy 测试半径，`R` 是整个 `z` 函数的 Taylor 展开圆半径，它们承担不同任务。

从每个 (4) 选取一行，横向指标记 `a_ρ`，同组数量记 `n_a`。令 `T` 为删除这些行标量后的测试矩阵。逐列逐组保留

\[
 |f_{a,P}(z)|\le 2^{|\alpha|}e^{hR}R^{|\alpha|-|a|},\qquad |z|\le R.
\]

这些因子可分别从系数行列式的列和行提出；Cauchy 测试再给 `r^{-s}γ^d`。相同 `(a,d)` 产生相同系数行而严格消失。对每组直接求和：

\[
 \sum_{d_1,\ldots,d_n\ {\rm distinct}}\gamma^{\sum d_j}
 =\frac{n!\,\gamma^{n(n-1)/2}}{\prod_{k=1}^n(1-\gamma^k)}.
\]

故得到统一估计

\[
\begin{aligned}
 \ln|\det T|\le{}&-c_R\sum_a n_a^2
 +R\sum_{\rm columns}h
 +\ln(2R)\sum_{\rm columns}|\alpha|\\
 &-\ln R\sum_\rho|a_\rho|
 +\ln(1/r)\sum_{\rm rows}s+o(MH). \tag{5}
\end{aligned}
\]

该 `o(MH)` 对所有行选择和选中列一致：`M!`、各组 `n_a!` 和几何乘积的对数总共为 `O(M ln M)+O(M)`；固定参数时 `M=O(H^{m+1})`，除以 `MH` 后趋零。零行列式自动满足上界。

这里有三个不同于原最坏值包络的收费方式：

- 指数函数增长为 `R Σh/(MH)=Rζ_H/w₀`，并非每一列都付 `R/w₀`。
- 多项式增长按实际 `Σ|α|` 收费；可进一步用 `Σ|α|/(MH)≤χ_H/w_*`，但无需预先丢掉它。
- 因完整行包保持 `s` 不变，测试费按实际 `Σs` 收费，且单纯形均值给
  \[
  \frac{v_0\sum s}{MH}\to\frac1{m+2},\qquad
  b_H\to\frac{m\theta}{m+2}.
  \]
  因而此项极限是 `ln(1/r)/((m+2)v₀)`，原包络 `ln2/v₀` 在 `r=1/2` 时还多付了约 `m+2` 倍。

碰撞求和的证明单独见 [cauchy_rewrite.md](cauchy_rewrite.md)。该文的若干障碍定理保留了旧的统一增长包络和最坏测试费，**不适用于已经采用 (5) 的所有方案**。

## 4. 把两种收益放回同一个优化问题

令 `𝒞_H` 是选中列，目标指标集合为

\[
 \mathcal T_H=\{a\ge0:\ a\le\alpha\text{ 对某个 }(h,\alpha)\in\mathcal C_H\}.
\]

对每个行选择 `a_ρ∈𝒯_H`，要求 `a_ρ≥β_ρ`；`j=0` 时再要求 `a_ρ=β_ρ`。零标量的选择可直接排除。定义精确收益

\[
\begin{aligned}
 \mathcal G_H((a_\rho))={}&c_R\sum_a n_a^2+\ln R\sum_\rho|a_\rho|\\
 &-\sum_\rho\ln\binom{a_\rho}{\beta_\rho}
 -\sum_{\rho,i}(a_{\rho,i}-\beta_{\rho,i})\ln|\epsilon_{j,i}|.
\end{aligned} \tag{6}
\]

指数为零的因子贡献 0；`ε=0` 而指数正的项被排除。若某个非零分母近似恰为 π，需另行处理；π 的无理性可作为已知输入，本轮不重证。

令 `𝒢_H^min` 是这些允许选择上的最小值。允许集合可以包含实际系数行列式为零的选择，这只会降低最小值，使上界更保守。将 (4) 展开、对 (5) 求和，得到

\[
 \ln|\Delta_H|\le \mathcal B_H-\mathcal G_H^{\min}+\mathcal R_H,
\]

其中

\[
 \mathcal B_H=R\sum_{\rm columns}h+
 \ln(2R)\sum_{\rm columns}|\alpha|+
 \ln(1/r)\sum_{\rm rows}s. \tag{7}
\]

可以取完全显式的有限余项

\[
 \mathcal R_H=2\ln(M!)+M\{c_R-\ln(1-\gamma)\}
                +mM\ln(\lfloor H\rfloor+1).
\]

前两项控制系数行列式、各组排列和几何乘积。最后一项来自每行候选个数至多 `(floor H+1)^m`（源权重 `w_i≥1`）。固定参数后 `𝓡_H/(MH)→0`。所求最小值保留了**碰撞收益与误差收益同时发生**，没有做原稿的单阈值二选一。

因此一个明确的、仍待建立的充分条件是

\[
 \boxed{\quad
 \liminf_{H\to\infty}
 \frac{\mathcal G_H^{\min}-\mathcal A_H-\mathcal B_H}{MH}>0.
 \quad} \tag{8}
\]

若已为同一组固定数据证明 (8)，上下界才会冲突。当前没有为 π 的端点假设证明 (8)，也没有证明它永远失败。

## 5. κ 在新模型中有明确位置

设 `κ_i=q_i²|π−p_i/q_i|>0`。对于 `j>0`，

\[
 -\ln|\epsilon_{j,i}|=2\ln q_i+\ln(1/\kappa_i)-\ln(2j).
\]

因此每个方向额外贡献 `(a_i−β_i)ln(1/κ_i)`，本式已不被尾指标消耗。不能将所有方向统一压成最小的 `ln(1/κ_i)/w_i` 后，又声称保留了全部信息。

若 `R≥2`，利用 `binom(a,β)≤2^{|a|}`，(6) 的保守下界是

\[
 c_R\sum_a n_a^2+\ln(R/2)\sum_\rho|a_\rho|
 +\sum_{\rho,i}\lambda_{j,i}(a_{\rho,i}-\beta_{\rho,i}),
 \quad\lambda_{j,i}\le-\ln|\epsilon_{j,i}|.
\]

当这些系数用有理下界给出时，这个保守下界的最小化是一个带凸占用成本的精确最小费用流问题。参见 [refined_transport.md](refined_transport.md)。有限模型求解该保守下界，不等于求解含二项式对数的精确收益，也不等于证明其 `H→∞` 极限或取得 (8)。

## 6. 当前具体缺口

当前真正需要的新结果是：同时控制所选非零基的次数分布与受坐标单调约束的占用优化。选便宜列会改变 `𝒯_H`、`Σh`、`Σ|α|`；减少某项成本可能使其他项更差，必须在同一组选列上比较。

另需保持量词：先由反设取得固定的有限组近似，参数与中心全部固定，再取 `H→∞`。若要让维数随近似质量变化，必须重新证明兼容的选取规则；固定维数的 `o(MH)` 不能直接用作变维数误差。

这把问题从“反证法能不能用”改成了 (8) 这个可精确表达、可被模型否定或支持的构造问题。它仍是研究缺口。
