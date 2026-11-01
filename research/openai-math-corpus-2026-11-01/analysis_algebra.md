# 第二轮语料分析：代数、零项与消失阶

2026-11-01。固定上游版本 `openai/math@4eaee87548250b2682affda7d05d947c05554bae`。

**结论：本轮最可直接迁移的内容来自 family 116 的微分系统／Wronskian 取系数方法。它导出我们实际纯周期行的有限容量和 Taylor 截断深度工具；没有给出 π 的逼近结论，也没有给出所需的渐近最低存活位移。** family 108 的行列式复杂度结论不直接提供本题的算术小量。

## 阅读范围和优先级

本分支对 `overview.tex` 全目录作 interpolation、multiplicity、determinant、valuation、matroid、Wronskian、vanishing、positivity 等关键词筛查。目录命中不算全文审阅。结合条目说明，定向查看 family 108、111、114、116、185、193、194 等候选；其中 111、185、193、194 仅作摘要级排除，114 的主定理和秩权重设置与组合分支交叉核对，详细结论归该分支。

这里的主要原文阅读为：

| family | 阅读到的原文 | 对当前缺口的判断 |
|---|---|---|
| 116 | `introduction.tex`、`multiplicity.tex`、`construction.tex` 的主定理及完整相关证明 | 微分系统的有限零检测深度可以转成精确工具；其主结果是非交换公式，不能直接套到交换行列式恒等式 |
| 108 | `introduction.tex` 主定理与路线；`polar.tex` 的 simple-points、product-count 引理；`transfer.tex` 的 characteristic-transfer 引理；`conductor.tex` 的 conductor-length 命题及其使用条件 | 属于行列式表示复杂度、光滑超曲面与一般参数几何；没有对固定 π 周期矩阵的高度／小量估计 |
| 114 | `main.tex` 主定理、精确 oracle 约定、配对秩缺陷权重与 Lorentzian 设置 | 与候选独立横截集直接相接，但不把正计数混同带符号行列式求和；详见组合分支 |

没有逐篇核验全部语料的证明，没有把上游标题声称的全部结论当作已经独立验证的事实。本分支没有运行这些论文的 Lean 证明。

## 1. Family 116：可迁移的是有限零检测，而非直接的逼近定理

原稿 *One Rational Matrix Hitting Point for Noncommutative Formulas* 的[引言主定理](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/One-Rational-Matrix-Hitting-Point-for-Noncommutative-Formulas-September-24-2026/build/introduction.tex)陈述：对 n 个非交换变量、大小至多 s 的无除法公式，可构造维数至多 2ns² 的有理矩阵元组，检测任意特征零域上该类公式的非零性。

**不能直接搬用的地方。** 我们当前的行列式是交换多项式。将表达式任意提升为非交换公式会改变恒等式：例如 xy−yx 在交换环中为零，在自由非交换代数中非零。因而非交换命中点并不自动认证原交换式非零。当前纯周期符号变量只有一个 Z，直接精确多项式秩已经可用，没必要为了套主定理把对象改掉。

更有用的是[消失阶一节](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/One-Rational-Matrix-Hitting-Point-for-Noncommutative-Formulas-September-24-2026/build/multiplicity.tex)给出的 Moura 型估计。若 F 特征零，q∈F[z] 次数 n、q(0)≠0，P 为 m×m 多项式矩阵且每项次数≤n−1，G∈GLₘ(F[[z]]) 满足 G′=(P/q)G，则对常数行 u、常数列 v，非零 y=uGv 有
\[
\operatorname{ord}_0y\le(n-1)\binom m2+m-1.
\]
原稿用导数行递推、Wronskian 的阶数和最大子式理想不变性给出形式证明。[构造一节](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/One-Rational-Matrix-Hitting-Point-for-Noncommutative-Formulas-September-24-2026/build/construction.tex)再将有限 Taylor 系数变成矩阵计算。

此不等式的方向是**消失阶上界**：它给“检查到多深即可确认零”的停止条件。我们所需的是随 H 增长的**最低非零位移代价下界**。二者不能交换。一般行列式若先转成巨大微分系统，也必须计入系统维数与系数高度的增长。

## 2. 本轮导出的工具 A：按函数类型给展开行设硬容量

固定实际选列集合 \(\mathcal C\)，其中每列 \(P=Y^hX^\alpha\)，h 为非负整数，α 为非负整数多重指标。固定展开位移 a。原整个函数为
\[
f_{a,P}(z)=
\begin{cases}
\binom\alpha a e^{hz}z^{|\alpha|-|a|},&\alpha\ge a,\\
0,&\text{否则}.
\end{cases}
\]
将每个可用列赋予签名 \((h,d)=(h,|\alpha|-|a|)\)，定义
\[
\mathcal S_a=\{(h,|\alpha|-|a|):(h,\alpha)\in\mathcal C,\ \alpha\ge a\},
\qquad r_a=|\mathcal S_a|.
\]

**容量命题。** 任意选择了 nₐ 个同一 a 的纯周期展开行的非零方形行列式，必须满足
\[
\boxed{n_a\le r_a.}
\]

证明：同一 a 的任意行都是对函数 eʰᶻzᵈ 施加线性取值／取导数泛函。每列只是一种签名函数的常数倍。所以该 a 行块分解为“行泛函在 rₐ 个签名函数上的取值矩阵”乘以“签名至列的系数矩阵”，秩至多 rₐ。若 nₐ>rₐ，所选行已有线性相关，整个方阵行列式为零。证毕。

这提供一个比仅计总行数更严格、可直接加入现有整数优化器的约束。它是必要条件；满足所有单组容量仍可能因不同 a 之间的相关关系而行列式为零。原评分的 \(c\sum_a n_a^2\) 是非线性占用成本，本工具没有把它改成逐元素可加权费用，也没有据此声称普通加权拟阵交已经解出整个评分。

在旧四列集合 1、X₂、X₂²、X₁ 中，a=(0,0) 只有 (0,0)、(0,1)、(0,2) 三种签名。因此四个不移动行的行列式必为零。这个容量判断直接恢复上一轮符号秩为 3 的现象。

**增长限制。** 对全部原单项式列，rₐ 可能达到 H² 量级，而每个原 β 的行数只在 H 量级。此时容量可能完全不生效。需要在同一个实际选列集合上核算 rₐ；不能从小例子推断它一定提供正的渐近收益。

## 3. 本轮导出的工具 B：确定有限 Taylor 深度

对每个出现的 h，令 dₕ 为该 a 的最大签名次数，令
\[
D_a=\sum_h(d_h+1).
\]
固定任意中心常数 ζ（也可取形式变量 ζ=jZ），考虑
\[
\mathcal E_a=\operatorname{span}
\{(1+t)^h(\zeta+\log(1+t))^k:0\le k\le d_h\}.
\]
**有限截断命题。** 若 y∈𝓔ₐ 的前 Dₐ 个 Taylor 系数全为零，则 y 恒为零。故这些签名函数的第 0 至 Dₐ−1 阶系数构成一个精确的有限独立性判据。

自包含证明：设 z=log(1+t)，该形式换元首项为 t，保持非零级数的消失阶。换元后 y 成为若干 eʰᶻ(ζ+z)ᵏ 的线性组合，满足常系数微分方程
\[
\prod_h(\partial_z-h)^{d_h+1}y=0.
\]
方程阶数为 Dₐ，最高导数系数为 1。若前 Dₐ 个导数在零点都为零，方程逐次求导给全部更高导数为零，因此 y=0。证毕。这正是原稿微分系统思想在 q(t)=1+t、常矩阵分子下的特例；本证明不依赖其非交换主定理。

完整向下次数闭包的 Dₐ×Dₐ Taylor 矩阵还有非零行列式公式。按 h 递增、同 h 下 k 递增排列，令 mₕ=dₕ+1，则
\[
\det=
\frac{\prod_h\prod_{k=0}^{d_h}k!}
{\prod_{\ell=0}^{D_a-1}\ell!}
\prod_{h<g}(g-h)^{m_hm_g}.
\]
这是正规化的合流 Vandermonde 行列式；ζ 平移和 z=log(1+t) 在各自有限系数包上均为单位三角变换，故不改变该行列式。不同 h 是不同整数，右侧非零。

附 [capacity_tool.py](capacity_tool.py) 计算签名容量、闭包深度，并用精确 SymPy 符号矩阵验证六组实例。所有六组都检验：前 Dₐ 个 t 系数的秩等于签名数；独立使用 z 导数得到相同秩；闭包行列式等于上述公式。[capacity_results.json](capacity_results.json) 是实际运行结果。这个验证不意味着所有实际行包都已包含这些 Dₐ 个导数；深度只是一个有限饱和证书。

## 4. Family 108：哪些几何工具不能混成 π 的算术结论

原稿 *A Cubic Lower Bound for Border Determinantal Complexity of the Permanent* 的[主定理](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/A-cubic-lower-bound-for-border-determinantal-complexity-of-the-permanent-September-24-2026/build/source/sections/introduction.tex)陈述 permanent 的 border determinantal complexity 至少为常数乘矩阵阶数的三次方。这里比较的是“以仿射线性矩阵的行列式表示一个多项式所需的矩阵大小”，允许任意复系数及系数逐点极限。

这一主结论不控制我们的某个整系数／高斯有理行列式在指定超越点的绝对值，不给其分母高度，不给与 q²|π−p/q| 的关系。因此不能把“行列式复杂度大”当成“行列式值小”或“位移必须大”。

可借鉴但尚无现成接口的两段是：

- [polar.tex](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/A-cubic-lower-bound-for-border-determinantal-complexity-of-the-permanent-September-24-2026/build/source/sections/polar.tex)：光滑齐次超曲面的一般极向截面产生指定数量的简单交点，并以多齐次 Bézout 控制简单孤立交点，即便另有多余分量。要使用它，先得为当前矩阵构造适当光滑参数空间和保持算术高度的映射；尚未完成。
- [transfer.tex](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/A-cubic-lower-bound-for-border-determinantal-complexity-of-the-permanent-September-24-2026/build/source/sections/transfer.tex)：利用适当射影族的纤维维数上半连续性，将特征二的一处梯度维数界提升到复参数的非空 Zariski 开集。它给一般参数结论，未保证原问题固定的周期对角或指定有理中心落在该开集；也没有给必要的系数高度控制。

[conductor.tex](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/A-cubic-lower-bound-for-border-determinantal-complexity-of-the-permanent-September-24-2026/build/source/sections/conductor.tex)的长度比较依赖特定有限平坦、泛可分的平面曲线笔。当前缺少将它接到我们的射流矩阵并保持全部分母成本的构造。因此本轮把 family 108 降为备选，优先实现 family 116 诱导的两个有限代数工具。

## 下一步可执行任务

在同一个已认证非零的选列集合上，先计算全部 rₐ，给原占用优化加入 nₐ≤rₐ；仍保留全部非线性平方成本。再与完整符号秩判据比较，检查这些容量在 H 增长时是否持续生效。若只是有限 H 的改善，必须明确保留该限制；本轮没有证明正的归一化收益，更没有完成 P2848。
