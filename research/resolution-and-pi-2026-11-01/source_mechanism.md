# π017 的结构性消去与观测分辨率：源码核对

2026-11-01。固定来源 `openai/math@4eaee87548250b2682affda7d05d947c05554bae`，稿件 *The irrationality exponent of π is 2*。本文核对“异常好的内凹是否由分辨率造成”这一解释在实际证明中能对应什么。没有独立认证上游插值几何定理、π 的主结论或坏逼近端点；也没有把论文解释为物理场。

## 可以直接向用户说明的联系

论文先假设有任意大分母的异常好逼近，然后从假设给出的近似中挑选若干个，建立一个辅助检测装置。它并未实际找到一个违反结论的 π 逼近，也不是把空间中的某个位置挖成了真实的凹陷。

这个装置利用周期使某些信息重合。重合不是“精度不够所以看不出差别”：在相应 Taylor 展开项中，重复的系数行完全相同，那个展开项严格等于零。要留下非零项，就必须使用不同的展开次数；许多行同时受到这一限制，便产生很强的衰减。

实际有理中心与精确周期有一点偏离，这些偏离又出现在精确的行平移公式中。偏离越小，相应高次因子越小。但非零性要另外由插值保证，并且有理分母还要清除。论文最终比较的是同一个非零对象的解析上界与算术下界，不是把图上无法分辨的距离认定为零。

因此，“有限分辨率会把平滑运动显示为跳跃，或把不同方向显示为重合”可以单独建模检验；它不能直接解释论文中的严格相等与消去。提高观测精度可能揭示原先看不见的非零差异，而精确算术中的真正相等仍然存在。

## 1. 异常好逼近是反设，不是已经构造出的反例

准确源码定位：

- [`conclusion.tex` L55–60](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/conclusion.tex#L55)：固定 ν>2，反设有任意大的 q 满足 `|π−p/q|≤q^(−ν)`。
- 同文件 [L98–115](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/conclusion.tex#L98)：从这个假设族中依次选择分母，以满足权重分离和误差条件，然后设置中心 `c_ji=jr_i`、`r_i=2ip_i/q_i`。
- [`main.tex` L201–209](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/main.tex#L201)：引言给出同样的逻辑次序。

这里“构造了好位置”只能指构造检测假设数据的辅助多变量配置，不能理解为证明中真的生成了无界的异常好 π 逼近。若把这个反设导成矛盾，结论正是相应假设族不存在。

还应区分辅助中心与轨迹：原中心是 `(1,jr_1,…,jr_m)`；并不是 `(exp(jr),jr,…,jr)`。精确指数周期在解析平移时使用。把近似周期直接放进指数里，会改变原矩阵。

## 2. 有限 H 是精确系数包，不是数值像素

[`main.tex` L168–185](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/main.tex#L168)将 jet 定义为有限 Taylor 系数包。列与行由加权次数截断选出。

[`determinant.tex` L45–60](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/determinant.tex#L45)给出的矩阵元，是代换后 `t^s u^β` 的**精确系数**。不存在先把位置或系数舍入到某个像素大小的步骤。有限 H 限制的是选择了哪些多项式与系数测试，不是这些测试值的小数位数。

这两种“有限”在逻辑上不同：

| 操作 | 保留什么、舍去什么 | 能否把两个不同数直接合并 |
|---|---|---|
| 精确的有限 jet 截断 | 只测试指定有限阶的精确系数 | 未被测试的高阶差异可能不出现，但已保留的系数仍按精确值区分 |
| 数值量化或有限精度显示 | 把数值舍入到网格或误差容限 | 可以把不同的已保留系数显示为相同 |

第一种可以称为有限观测阶数，但不能因此把证明中的精确相等解释为浮点误差。并且证明最后令 H 增大：先固定维数、权重和中心，再令 H 沿允许的充分可整除次数趋于无穷，见 [`conclusion.tex` L133–135](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/conclusion.tex#L133)。

对数截断也不是任意丢掉误差。源码 [`determinant.tex` L28–41](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/determinant.tex#L28)要求尾项有足够高的权重，使其保持有限系数空间的过滤结构与满射。此前本项目还证明：完整行包上相应变换是单位三角的，对同一组选列行列式完全不变，见 [精确消尾记录](../pi-endpoint-rewrite-2026-11-01/tail_elimination.md)。

## 3. Taylor 碰撞是精确复数上的结构性消去

源码 [`determinant.tex` L154–169](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/determinant.tex#L154)定义整函数系数行 `f_(a,P)`；[L196–222](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/determinant.tex#L196)将原行精确展开为周期处测试行的线性组合。

对固定横向指标 a，按 Taylor 次数 d 展开。若一个行列式项中两行选择相同 `(a,d)`，两条系数行完全相同，所以该项为零。这是 [`determinant.tex` L308–316](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/determinant.tex#L308)的具体机制。

同组有 n 行，要使项非零，至少必须使用 n 个不同的非负次数，最小总次数是 `0+1+⋯+(n−1)=n(n−1)/2`。Cauchy 的几何衰减将其变成二次收益。这个等式与交错性在精确的复数运算中成立，提高数字精度不会取消它。

源码 [L394–406](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/determinant.tex#L394)再分两类：小横向指标多，就迫使许多行挤在少数组内；小指标少，就有很多大指标承担小位移的高次幂。它没有仅凭“点很近”就断言整个矩阵秩下降。

## 4. 六列模型中哪些退化真的存在

沿用已经保存的[精确六行模型](../pi-construction-mechanisms-2026-11-01/determinant_mechanism.md)，列固定为

\[
 1,\ Y,\ X_2,\ X_2^2,\ X_2^3,\ X_1,
\]

行固定为两个中心的 0、1、2 阶 t 系数，初始横向指标全为零。它满足一般恒等式

\[
 \det A(r_1,r_2)=-3r_2^2(r_1-r_2).
\]

原因可直接定位：`X_1` 列减 `X_2` 列，只剩一个分量 `r_1−r_2`。这不是计算精度造成的近似相等。

- 若 `r_2≠0` 且 `r_1≠r_2`，这组选列的精确秩为 6。
- 若 `r_1=r_2=Z≠0`，这组选列的精确秩为 5。
- 特别地，代入精确周期 `Z=2πi` 仍然严格为 5。其退化不会在“实数/复数的最好分辨率”下消失。

这只是特定六列子式的结论。它不能推广成“整个插值矩阵在精确周期处都只有五维”，也不能推广成“每个解析展开项都退化”。同一模型的某些横向平移项自身仍为满秩，且不同非零项还可相消。

### 固定坐标下可得一个精确可辨识性界

在标准 Euclidean 向量范数下，令 `δ=r_1−r_2`，选单位列组合 `(e_(X1)−e_(X2))/sqrt(2)`。其像的范数恰为 `|δ|/sqrt(2)`，所以

\[
 \sigma_{\min}(A)\le |\delta|/\sqrt2.
\]

常数列的范数为 `sqrt(2)`，故 `σ_max(A)≥sqrt(2)`。当 A 可逆时，二范数条件数满足

\[
 \operatorname{cond}_2(A)\ge\frac2{|r_1-r_2|}.
\]

因此该固定模型接近对角线时，对某个方向确实越来越难区分；有限分辨率可以更早把它误判为完全重合。但难区分的趋势本身在精确运算中已经存在。这里明确固定了列尺度与坐标范数；不能把任意缩放下的“小行列式”都自动等同于条件数变大。

## 5. 算术下界不是“原行列式小于 1 就矛盾”

源码 [`determinant.tex` L99–126](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/determinant.tex#L99)先缩放行列并清除分母，得到非零高斯整数。模至少 1 的是这个清分母对象；原行列式可以很小，清分母成本必须同时计算。

所以真实证书需要：同一个非零子式的解析衰减严格压过其分母成本，才能产生矛盾。原稿在 [L425–431](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/determinant.tex#L425)进行的正是这项比较。

## 6. 分辨率解释的适用边界

可以把以下两种效果并列研究：

1. **结构性小量：** 精确周期、相同系数行、交错性和小位移幂共同造成。在精确运算中仍存在。
2. **观测伪像：** 舍入、采样或阈值判断把不同点/不同系数合并，或把平滑变化显示成跳跃。它依赖观测规则，改变分辨率后可能消失或移位。

两者可能同时出现，不能只看一张低分辨率图就分清。实际检验应先固定连续对象、采样/量化规则与用于比较的距离或范数，再独立比较精确值和显示值。

论文没有定义天体动力学、真实作用力或速度时间参数。若把接近程度约定为“速度”，那是另建的模型；它不能反过来被用作 π017 已有物理场的证据。当前可以建立的联系是辅助观测的可辨识性及精确相关性，而非已证明的空间内凹机理。

## 本轮核对范围

实际读取 `main.tex` L48–264、`determinant.tex` L1–432、`conclusion.tex` L45–156，以及先前六行模型和精确消尾记录。本文没有重新运行先前完整实验，没有重新认证整个插值依赖闭包或 Lean 内核；六列一般恒等式及相应秩事实沿用已保存的精确模型。新增条件数界由上述固定列组合直接推导，不依赖有限数值采样。
