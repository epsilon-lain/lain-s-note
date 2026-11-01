# 尾项精确消去之后：保留中心、选列支持与异质逼近收益

本文接续 [完整行 packet 的尾项消去](tail_elimination.md)、[精确 Taylor 碰撞求和](cauchy_rewrite.md) 和 [初版联合预算](collision_rewrite.md)。目的不是把几个改善直接宣布成 \(\pi\) 的结论，而是给这些改善建立一个共同且可核对的有限模型。

原始构造固定为 `openai/math@4eaee87548250b2682affda7d05d947c05554bae` 中 π017 的 `build/sections/determinant.tex`。关键位置：45–64 行是全部行与同一个非零选列子式；154–169 行是测试函数；196–222 行是逐行展开；267–335 行是碰撞；394–406 行是原来的两分支收益丢弃。

## 1. 精确消去后，每行只剩逼近误差平移

固定每个中心的完整加权行 packet。单位下三角行变换给出

\[
A_G=U_H A_{\log},\qquad\det U_H=1.
\]

因此任意固定的满行选列子式都有完全相同的行列式值。对 \(A_G\) 做算术整数化、对 \(A_{\log}\) 做解析估计，仍然是在估计同一个非零数。

令 \(\omega=2\pi i\)、\(r_i=2ip_i/q_i\)、\(\epsilon_{ji}=j(r_i-\omega)\)，测试函数为

\[
f_{a,P}(z)=[u^a]P(e^z,z+u),\qquad
f_{a,Y^hX^\alpha}(z)=\binom\alpha a e^{hz}z^{|\alpha|-|a|}.
\]

这里 \(|a|=\sum_i a_i\)，\(\binom\alpha a=\prod_i\binom{\alpha_i}{a_i}\)。精确展开为

\[
\operatorname{row}_{j,s,\beta}^{\log}
=\sum_{a\ge\beta}\binom a\beta\epsilon_j^{a-\beta}
 \left([t^s]f_{a,P}(j\omega+\ln(1+t))\right)_P.
\]

没有截断尾项的额外多指标，也没有降低的测试阶数。由此得到两个应当保留的零结构：

1. \(j=0\) 时，\(\epsilon_0=0\)，任何非零项都必须有 \(a=\beta\)。此时原来的二项式因子也恰为 1。
2. 设选中列的横向指标集合为 \(\mathcal A\)。若没有任何 \(\alpha\in\mathcal A\) 满足 \(a\le\alpha\)，则整个测试行恒为零。应当只保留
   \[
   \mathcal D(\mathcal A)=\bigcup_{\alpha\in\mathcal A}\{a:0\le a\le\alpha\}.
   \]

这些条件只利用同一个已选子式的信息，不需要换另一个有利子式。

## 2. 把行收益与实际选列增长放在同一条估计里

取测试圆半径 \(0<r<1\)，并令

\[
S=2\pi(K-1)-\ln(1-r),\qquad R>\max\{S,1\},
\qquad c_R=\frac12\ln\frac RS.
\]

在 \(|t|=r\) 上有 \(|j\omega+\ln(1+t)|\le S\)。Taylor 系数的精确碰撞求和保留了完整的几何衰减。对一项选定的逐行横向指标 \(a_\rho\)，令 \(n_a\) 为其占用数，则需要关注的收益是

\[
\begin{aligned}
\mathcal S(a_\rho)
={}&c_R\sum_a n_a^2
 +\ln R\sum_\rho |a_\rho|
 -\sum_\rho\ln\binom{a_\rho}{\beta_\rho}\\
&-\sum_{\rho,i}(a_{\rho,i}-\beta_{\rho,i})\ln|\epsilon_{j_\rho i}|.
\end{aligned}
\]

零误差因子先用精确零规则处理，不能直接数值计算 \(\ln0\)。

选中列的增长应单独保留为

\[
\boxed{G_R=R\sum_{\text{selected columns}}h
       +\ln(2R)\sum_{\text{selected columns}}|\alpha|.}
\]

这是因为每列的测试函数满足

\[
|f_{a,P}(z)|\le
 R^{-|a|}\exp\{Rh+\ln(2R)|\alpha|\},\qquad |z|\le R.
\]

因子 \(R^{-|a|}\) 按行提取，余下增长按选中列提取。这里没有把所有列都替换成最大允许次数。

测试代价是

\[
T_r=\ln(1/r)\sum_\rho s_\rho.
\]

还有 Taylor 求和、行列式置换、有限项展开数量等余项。例如令 \(\gamma=S/R\)，碰撞证明可以用统一余项

\[
B_{\rm coll}\le2M\ln M+c_RM-M\ln(1-\gamma).
\]

若每行最多 \(D=\#\mathcal D(\mathcal A)\) 个候选横向指标，逐项取绝对值再求和还可收取 \(M\ln D\)。这些量必须加入同一个总成本；当前脚本没有假装它们已经为实际 \(\pi\) 参数算好。

## 3. 变成有理凸费用流

设调用者已经严格证明有理下界

\[
0<c_0\le c_R,\qquad
0\le d\le\ln(R/2),\qquad
\lambda_{j,i}\le-\ln|\epsilon_{ji}|\quad(j>0).
\]

利用 \(\binom a\beta\le2^{|a|}\)，得到

\[
\boxed{\mathcal S(a_\rho)\ge
c_0\sum_an_a^2
 +\sum_{\rho,i}\lambda_{j_\rho,i}(a_{\rho,i}-\beta_{\rho,i})
 +d\sum_\rho|a_\rho|.}
\]

这里取 \(R\ge2\)，以保证可使用非负的 \(d\)。\(\lambda\) 可以是负的有效下界，程序也支持；这种输入只意味着所取下界在该坐标不提供正收益。\(j=0\) 必须由前面的精确规则强制 \(a=\beta\)，其 \(\lambda\) 是无须解释为对数的占位输入。

把同一 \((j,\beta)\) 的全部 \(s\) 行合并，供给为

\[
m_{j,\beta}=\#\{s\ge0:v_0s+w\beta/\theta<H\}.
\]

若研究的是完整的原 packet，这个数在各 \(j\) 相同。程序允许一般有限供给，便于独立测试。

供给节点 \((j,\beta)\) 向目标 \(a\) 送一单位的线性成本是

\[
\sum_i\lambda_{j,i}(a_i-\beta_i)+d|a|.
\]

边存在的条件是 \(a\ge\beta\)、\(a\in\mathcal D(\mathcal A)\)，且 \(j=0\) 时 \(a=\beta\)。目标 \(a\) 的第 \(k\) 个单位另收凸边际成本 \(c_0(2k-1)\)。因此所有流的总成本恰好是框内目标。

[refined_transport_lab.py](refined_transport_lab.py) 取所有有理成本分母的最小公倍数，把问题变成整数最小费用流。逐次最短增广路用 Bellman–Ford，因而不需要假设所有 \(\lambda\) 非负。流解的成本还会从最终占用数与实际运输分配重新计算一次，逐项核对。

此松弛仍然可能包含真实系数行列式为零的分配；它保留了两种必要零规则，但没有声称这已经刻画全部非零项。

## 4. 异质 \(\kappa_i\) 应当如何进入

对 \(j>0\)，若

\[
\left|\pi-\frac{p_i}{q_i}\right|=\frac{\kappa_i}{q_i^2}>0,
\]

则精确有

\[
-\ln|\epsilon_{ji}|=2\ln q_i+\ln(1/\kappa_i)-\ln(2j).
\]

不同分母、不同 \(\kappa_i\)、不同中心的收益都不相同。用最小值统一替换是一种合法但可能有损的放松，不能把 \(\kappa_i\to0\) 自动换成固定正的额外幂指数。

实际使用时，应分别为这些对数求严格有理下界，并为 \(G_R,T_r\) 及余项求严格有理上界。程序接受这些有理数；它当前没有生成真实 \(\pi\) 超常逼近，没有把浮点对数当成认证区间。

## 5. 每次只删除一种信息的独立实验

主实验用一维横向指标，三个中心 \(j=0,1,2\) 各供给两行，全部 \(\beta=0\)。这是一个有限供给模型，与完整 packet 在各中心有相同供给数的结构一致；没有声称它来自一个已验证满秩的真实 \(\pi\) 子式。

取候选目标 \(a\in\{0,1,2\}\)，选中列横向指标集合 \(\mathcal A=\{1\}\)，因此真正支持域是 \(\{0,1\}\)。再取

\[
c_0=1/6,\quad d=1/60,\quad
\lambda_1=0,\quad\lambda_2=1/2.
\]

\(\lambda_0=0\) 只是计算占位；在正常模型中 \(j=0\) 被强制留在 \(a=0\)，这个输入不参与逼近误差收益。

| 只改变这一项 | 精确最小收益 | 相对完整信息损失 | 穷举全部分配数 |
|---|---:|---:|---:|
| 保留全部信息 | \(101/30\) | 0 | 16 |
| 删除 \(j=0\Rightarrow a=\beta\) | \(61/20\) | \(19/60\) | 64 |
| 删除实际选列支持域 | \(29/10\) | \(7/15\) | 81 |
| 将异质 \(\lambda_j\) 换成最小值 0 | \(61/20\) | \(19/60\) | 16 |

三个损失在这个模型中都严格为正。删除 \(j=0\) 规则的实验人为允许占位值 \(\lambda_0=0\) 下的不可能分配，度量的是“忘掉精确零结构”会损失多少；没有把这些分配说成真实存在的 Taylor 项。

每一行的最小费用流结果，都用另一段程序独立枚举全部逐行分配核对。另有一个二维、非零 \(\beta\)、含负的 \(\lambda\) 下界的模型，遍历 15 种分配，最优值为 \(3221/4862\)，也精确一致。总共五项穷举交叉核对。

输出见 [refined_transport_results.json](refined_transport_results.json)。运行：

```sh
python3 research/pi-endpoint-rewrite-2026-11-01/refined_transport_lab.py
```

## 6. 与同一非零子式的成本接口

程序的 `compare_budget` 接受五个有理数：最小收益、同一子式的算术预算、实际选列增长、测试代价、其他余项。返回

\[
S_{\min}-C_{\rm ar}-G_R-T_r-B_{\rm residual}.
\]

若要把它变成真正的矛盾，调用者还须证明：所选子式非零；所有输入上下界方向正确；候选域覆盖全部真实非零展开项；有限或渐近余项全部包括；这些数据来自同一个选列子式。

本实验仅用人为有理成本测试接口的精确算术，没有把正的玩具差值解释为实际 \(\pi\) 结论。下一项实际建模工作是从一个经满秩验证的较小原矩阵中提取选列、供给与认证对数界，接到这个模型，再检查真实低成本流是否仍被更深的系数秩关系排除。
