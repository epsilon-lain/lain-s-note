# 数论与动力系统分支：对 π 逼近项目的定向复查

日期：2026-11-01。上游固定为 `openai/math@4eaee87548250b2682affda7d05d947c05554bae`。本文报告的是源码所陈述的结果及局部证明阅读，不把新预印本的主定理视为已经独立验证的事实。

## 实际覆盖范围

- **目录级筛选：43 个 family**，即数论 `001–031` 共 31 个及动力系统 `143–154` 共 12 个，读了其目录标题和说明。没有逐篇读这些 family 的所有修订稿。
- **定理与证明段落级阅读：3 个 family 各 1 篇**：`005` Catalan 常数、`022` 弱非齐次 Duffin–Schaeffer、`153` Bernoulli 卷积。没有把 `017` 的已有 π 审核重复计入本轮。
- 005：阅读 introduction 的行列式策略、nonvanishing 的命题及模素数分块证明、real-place 的插值与有限维算子证明；foundations、certificate 另有取回，未重算其全部证书。022：阅读 introduction 的主定理及量词说明、direct 的首段比较命题和整数行列式计数设置，没有审完整个直接证明。153：阅读 main 的总述及 classification 中的分类定义、有限阶非消失引理及证明、给定前缀的代数单位实现命题及证明。
- 上游总量是 **721 个 manuscripts、373 个 families**；“772”不是这个固定快照的计数。目录接入、源码取回、定向读证明、独立验证是四个不同的完成度。

## 005：优先迁移“固定证书推出增长矩阵族非零”

来源：`preprints/Catalans-constant-is-irrational-September-24-2026/build/sections/` 下的 [introduction.tex](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/Catalans-constant-is-irrational-September-24-2026/build/sections/introduction.tex)、[nonvanishing.tex](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/Catalans-constant-is-irrational-September-24-2026/build/sections/nonvanishing.tex)、[real-place.tex](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/Catalans-constant-is-irrational-September-24-2026/build/sections/real-place.tex)。

原文构造阶数为 `n=48N` 的混合行列式。先在每个矩阵元中消去 ζ(2)，再在反设 `G∈Q` 下使矩阵元有理；有限素数处给出分母下界，实位给出同一行列式的上界。它特别跟踪 `p²` 和 `p` 两层估值。这里的逻辑范式与 lain 的“同一非零整数/有理行列式同时得到不相容上下界”相近，但常数替换并不自动成立。

### 可精确定位的非零命题

`nonvanishing.tex` 的 `prop:nonvanishing` 陈述：若 `G∈Q`，则对所有充分大的素数 p，`v_p(Δ_p)=−96p`，从而 `Δ_p≠0`。

核心机制是固定有理 `49×48` 矩阵 B，取相邻的两个 `48×48` 块 B₀、B₁。源码的固定证书要求

\[
\det B_0\ne0,\qquad \det(B_0+B_1)\ne0,\qquad \det(B_0-B_1)\ne0.
\]

对奇素数 p，在有限域上用两个回文多项式基得到可逆矩阵 a 与 `b=aΠ`。其中 `Πe₀=0`，`Πeᵢ=e_{p−i}`，非零指标两两成对。因此 Π 的特征值是 0、1、−1，重数分别为 1、`(p−1)/2`、`(p−1)/2`。缩放后的大矩阵满足

\[
L_p=a^T\otimes B_0+b^T\otimes B_1,
\quad
L_p((a^T)^{-1}\otimes I)=I_p\otimes B_0+\Pi^T\otimes B_1.
\]

于是源码的 `eq:nonvanishing-factorization` 给出精确分解

\[
\det L_p=(\det a)^{48}\det B_0
\det(B_0+B_1)^{(p-1)/2}
\det(B_0-B_1)^{(p-1)/2}.
\]

这一步是直接的线性代数：在 p 维因子上对角化 Π 即可，不需要 B₀、B₁ 彼此交换。除去固定矩阵的分母和三个非零行列式涉及的有限个坏素数，大矩阵族全部非零。源码用模 101 证书验证三个固定矩阵；本轮没有重跑该证书，也没有独立审完把原始核降至这个分块式的所有算术估值。

**可迁移工具规格：** 对 lain 的实际行矩阵寻找类似 `A_p⊗B₀+C_p⊗B₁` 的精确分解，明确变化基底与分母，继而将无限增长规模的非零性化为少数固定有理矩阵的模素数证书。这比只枚举有限个 H 强。眼下尚未证明实际 π 矩阵有这种结构；这项结果也不直接给出“最小存活行位移代价”。

### 第二个可用机制：节点合流时避免逆 Vandermonde

`real-place.tex` 的 `lem:real-uniform-interpolation` 要求互异非零实节点 `xᵢ∈(−1,1)`，且

\[
\sum_i\frac{1-x_i^2}{1+x_i^2}\le D_*,\qquad D_*=40N-1,
\]

并为分段值 `cᵢxᵢ^{D_*}`（负节点 `cᵢ=1`、正节点 `cᵢ=−5`）构造闭单位圆附近的全纯插值函数 h，满足 `||h||∞≤10e¹²n`，常数与节点间距无关。

`prop:first-sheet-bound` 把混合求值矩阵写成有限维 Hardy 空间上的 `ELR`，其中 R 是等距反转、L 是 `I` 加压缩乘法算子。求值映射 E 的行列式在比值中**代数消去**，于是剩余因子有 `|det L|≤(1+10e¹²n)^n=exp(O(n log n))`。这是对接近碰撞的节点很有价值的范式：控制算子范数而非逆求值矩阵。

**迁移障碍：** 源码用的是特定实节点、两层求值核、Blaschke 乘积及上述总量条件。lain 的复周期节点、微分行、多变量行列式需另建对应的空间与有界插值函数；不能原样套用这个常数。其贡献是潜在的解析上界工具，不是非零证明，更不能在没有逐项核对的情况下与另一个较便宜子式拼接。

**反设边界：** `G∈Q` 会给出精确代数关系；“π badly approximable”或“存在 κ=q²|π−p/q|→0”只给不等式。后者不能替代前者完成有理性塌缩。必须先用真实的有理近似数据构造可清分母的矩阵，再为同一子式证明非零与节省。

## 022：度量结论不能用于指定 π；可借鉴重合分类

来源：`preprints/The-weak-inhomogeneous-Duffin-Schaeffer-conjecture-September-25-2026/build/sections/` 下的 [introduction.tex](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-weak-inhomogeneous-Duffin-Schaeffer-conjecture-September-25-2026/build/sections/introduction.tex) 与 [direct.tex](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-weak-inhomogeneous-Duffin-Schaeffer-conjecture-September-25-2026/build/sections/direct.tex)。

`thm:main` 的准确量词是：固定任意实数 γ 和有限值函数 `ψ:N→[0,∞)`，若

\[
\sum_{q\ge1}\frac{\varphi(q)}q\psi(q)=\infty,
\]

则对 Lebesgue 几乎处处的 x，`||qx−γ||<ψ(q)` 无穷多次成立。分子不要求与 q 互素，零测例外集可以依赖 γ 和 ψ；这不是定量计数渐近式。

取 `γ=0`、`ψ(q)=ε/q` 可得到几乎处处的非 BA 结论，但**不能把 x 指定为 π**。把 γ 改成 π 也不会使它变成 π 的有理逼近定理。

局部可借鉴的是 `direct.tex` 的 `dir:birth-comparisons`：在其层级、窗口和权重体系内，把精确重合计入一次质量，把非精确重合转化为整数行列式/同余计数并求和。在有限质量块 `W≤1` 中，它给出宽窗口的 `O(1)` 控制与窄窗口的 `O(W)+O_{k₀}(W^{1+c})` 控制。这些估计依赖全文的结构性假设，不能作为任意行矩阵的黑箱。

对 lain 的意义是一个审计规则：将零位移、精确碰撞、非精确近碰撞分别证明，避免让“很多形式项”代替“很多实际不消失项”。若要由此推出固定 π 的逼近，仍缺 π 专属的确定性命中或分布估计；该文没有提供。

## 153：紧性反证给出统一有限阶非消失，但必须控制系数

来源：[classification.tex](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/Arithmetic-classification-and-non-Pisot-singularity-for-Bernoulli-convolutions-October-3-2026/build/classification.tex)，路径为 `preprints/Arithmetic-classification-and-non-Pisot-singularity-for-Bernoulli-convolutions-October-3-2026/build/classification.tex`。

`cls:finite-order`（Uniform finite-order nonvanishing）的完整假设与结论是：固定 `λ∈(0,1)`，存在 `M≥1`、`c₀,x₀>0`，使每个实系数多项式

\[
p(t)=-1+\sum_{k\ge1}a_kt^k,\qquad |a_k|\le1
\]

都满足

\[
\sup_{t\in[\lambda,\lambda+x]}|p(t)|\ge c_0x^M
\quad(0<x<x_0).
\]

证明本身就是简洁的反证法：若所有有限 jet 都能同时趋零，系数紧性给出单位圆内的解析极限；它在 λ 的所有阶导数皆零，却有常数项 −1，违背恒等定理。取得有限 jet 下界后，有限维范数等价和统一 Taylor 余项界给出上述定量形状。本轮读完了这个引理的证明。

固定系数界 `|a_k|≤B` 也可沿同一紧性证明处理，但 M、c₀、x₀ 可能依赖 B 和 λ。故它可启发原生有界数字系数核的非消失证书；对系数随 H 增长的行列式多项式，尚无 H 一致的结论。把多项式缩放到系数≤1又会缩小常数项，必须逐项追踪，不能丢掉。它限制一个受控系数族能有多平坦，并不会自动产生所需高阶消去。

`cls:realization` 另陈述：固定 `λ∈(1/2,1)`、`0<a<(1−λ)⁻¹`，足够大的固定整数 d 及所有充分大的 n 下，若 `w∈{−1,0,1}ⁿ`、`w₀=−1`，并有

\[
\left|\sum_{i<n}w_i\lambda^i\right|\le((1-\lambda)^{-1}-a)\lambda^n,
\]

则有次数 dn 的实代数单位 r，其最小多项式系数在 `{−1,0,1}`、模 2 不可约、前 n 个系数为 w，并满足 `λ<r≤λ+4⁻ⁿ/2`。证明组合了上述有限阶非消失、保前缀的有限域不可约完成与介值定理。

这里近似物是**次数随 n 增长的代数单位**，不是分母 q 受控的有理数；指数距离不能转写成 `q²|π−p/q|`。而 `λ=1/π<1/2` 不在此实现命题的区间。分类主定理本身还含无穷交与 limsup，不能当作任意给定实数的有限判定算法。

## 仅到目录级、暂未升级的候选

| Family | 目录方向 | 与目标缺少的连接 |
|---|---|---|
| 007 | 乘法相关/Chowla | π 的数字或收敛分母并未证明满足相同结构假设 |
| 015 | 算术环面包等分布 | 对一族算术对象平均不等于控制 π 的单个轨道 |
| 025 | Egyptian fractions | 重写一个已知有理数不产生它对 π 的误差保证 |
| 145、154 | 遍历平均与混合 | 需要具体系统假设，几乎处处仍不涵盖指定 π |
| 148 | 自相似测度的熵率与维数 | 测度级维数不等于某个固定常数的丢番图性质 |
| 150 | 无理角三角形混合 | 几何参数无理性不能转成 π 的逼近端点估计 |

## 对下一步的具体建议

先把 **005 的有限模素数块证书** 列为非零性工具候选，把 **005 的算子式求值消去** 列为抗节点合流的解析工具候选，把 **153 的有界系数紧性反证** 列为一致性条件检查。每个候选都必须回到项目的实际行矩阵，并输出清分母、非零、最小存活位移、解析上界这四件事对**同一个子式**的证据。目前三篇均未补上这条连接，因此本轮没有得到 π 为 BA 或非 BA 的证明。
