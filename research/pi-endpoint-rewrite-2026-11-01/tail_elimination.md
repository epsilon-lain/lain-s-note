# 完整行 packet 使对数尾项从行列式中精确消失

日期：2026-11-01。这是对 π017 构造的代数改写，不依赖 π 是否坏逼近，也不把反证法当成障碍。

源固定为 `openai/math@4eaee87548250b2682affda7d05d947c05554bae`：[determinant.tex 28–41 行](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/determinant.tex#L28) 已说明尾项替换保持加权过滤与秩；[196–257 行](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/determinant.tex#L196) 则逐行展开尾项并为它收费。下面证明一个更强、直接对选定子式有效的事实。

## 1. 精确恒等式

固定 `H`，保留每个中心的**全部**行

\[
 \mathcal R_H=\{(s,\beta):v_0s+w\beta/\theta<H\}.
\]

设 `A_log` 使用完整形式对数 `ln(1+t)`，`A_G` 使用截断多项式 `G_i(t)`，列集合完全相同。令

\[
 \tau_i(t)=G_i(t)-\ln(1+t),\qquad
 \operatorname{ord}_t\tau_i\ge T_i,\qquad
 v_0T_i>w_i/\theta.
\]

对每个多项式列，形式 Taylor 展开给出

\[
 (A_G)_{(j,s,\beta),P}
 =\sum_{\gamma\ge\beta}\binom\gamma\beta
   \sum_{k=0}^s[t^k]\tau(t)^{\gamma-\beta}
   (A_{\log})_{(j,s-k,\gamma),P}.
\]

右边若系数非零，则 `k≥Σ_i T_i(γ_i−β_i)`。因此相应输入行和输出行的权重差为

\[
 \begin{aligned}
 &v_0(s-k)+w\gamma/\theta
   -(v_0s+w\beta/\theta)\\
 &\quad\le-\sum_i(v_0T_i-w_i/\theta)(\gamma_i-\beta_i).
 \end{aligned}
\]

只要 `γ≠β`，这个差严格小于零；若 `γ=β`，只有 `k=0` 能贡献，系数为 1。由于完整行 packet 包含所有更小权重的行，右边没有跑出已保留的行集合。

按行权重升序排列，每个中心的变换矩阵都是单位下三角矩阵。所有中心拼起来仍有

\[
 A_G=U_H A_{\log},\qquad \det U_H=1.
\]

从而对**任意固定的满行方形选列**都有

\[
 \boxed{\det A_G[:,S]=\det A_{\log}[:,S].}
\]

这不只是两者同时非零或秩相同，而是两个数精确相等。可以对 `A_G` 做有理系数整数化，再对同一个行列式直接用 `A_log` 做解析估计。

## 2. 为什么会省去尾项成本

分析 `A_log` 时，令 `ω=2πi`、`ε_{ji}=j(r_i−ω)`、`z=jω+ln(1+t)`，则

\[
 Y=e^z,\qquad X_i=z+u_i+\epsilon_{ji}.
\]

使用源稿的 `f_{a,P}(z)=[u^a]P(e^z,z+u)`，直接有

\[
 \operatorname{row}^{\log}_{j,s,\beta}
 =\sum_{a\ge\beta}\binom a\beta\epsilon_j^{a-\beta}
    \bigl([t^s]f_{a,P}(j\omega+\ln(1+t))\bigr)_P.
\]

这里没有任何 `τ_i`，没有原展开中的 `d` 和 `k`，测试阶数就是原来的 `s`。若 `|π−p_i/q_i|≤q_i^{−ν}`，`w_i=ceil ln q_i`，则各非零系数满足

\[
 \log\left|\binom a\beta\epsilon_j^{a-\beta}\right|
 \le-\nu w(a-\beta)
 +\frac{H}{w_*}\{\ln2+\ln(2K)+\nu\}.
\]

因此原来的平移误差

\[
 \nu/F_0+\ln2/v_0+
   \{\ln4+\ln(2K)+\nu\}/w_*
\]

可以替换为

\[
 \boxed{E_{\rm tr,new}=
       \{\ln2+\ln(2K)+\nu\}/w_*.}
\]

这里删掉的是**平移阶段**的 `ν/F_0` 和 `ln2/v_0`，不是全部解析误差：碰撞引理测试 `[t^s]` 的 Cauchy 估计仍然有自己的 `ln2/v_0`，整函数增长 `100K/w_0` 也仍在。

原稿要求 `F_0>2/θ`。精确变换仅要求 `v_0T_i>w_i/θ`，若仍令 `T_i=ceil(F_0w_i/v_0)`，则 `F_0>1/θ` 足够。与按列、行指数差整数化相结合，`F_0` 只需控制有限截断的系数分母，不再承担解析尾项的小量预算。

为什么逐行估计会损失，而完整行列式没有损失？尾项把某个输出行改成“该行加上更低权重行的线性组合”；它是一种行剪切，行列式天然不变。先拆成很多项分别取绝对值，就破坏了这层精确消去。这是可定位的估计损失。

## 3. 有真实尾项参与的精确实验

运行 `python3 research/pi-endpoint-rewrite-2026-11-01/tail_elimination_lab.py`。

采用 `m=1,K=2,H=10,w_0=w_1=1,v_0=2,θ=1/4`，所以横向 jet 权重 `w_1/θ=4`；中心为 `0` 和 `44i/7`。取 `F_0=5,T=3,G=t−t²/2`。最大行测试阶数为 4，真实尾项为

\[
 \tau=-t^3/3+t^4/4+O(t^5).
\]

有 18 行、66 列。变换矩阵 `U` 存在 **6 个非零的非对角元素**；例如输出 `(j,3,0)` 接收输入 `(j,0,1)`，系数为 `−1/3`，所以这个例子并非“尾项尚未进入 packet”。程序直接独立生成完整对数矩阵和截断矩阵，验证：

- `A_G=U A_log` 的每个元素精确相等；
- `U` 按权重排序后单位下三角，行列式精确等于 1；
- 截断矩阵精确满行秩；同一组 18 个主元列的两个非零行列式精确相等。

对剪切输入保持封闭的条件也有反向实验。只留下行 `(j=0,s=3,β=0)`，取单列 `P=Y³+X`，截断值为 `1`，完整对数值为 `4/3`，不再相等；被遗漏的 `(0,0,1)` 行正是需要的剪切项。完整加权 packet 自动保证这种封闭；某些更小但同样封闭的行集也可以，任意删行则不能保证。

所有计算用 SymPy 精确有理数/高斯有理数，结果在 `tail_elimination_results.json`。有限例子验证恒等式实现与消去机制，不认证源稿“大权重插值”的全部证明。一般恒等式的证明在第 1 节，独立于这个有限实验。

## 4. 下一步仍需控制什么

这项改写真正省去了原证明中可避免的项；它没有自动给出指数 2 的结论。还需要把改进后的算术成本、完整的逼近误差指数、碰撞收益和非零子式选列一起比较。尤其应保留 `c_H−b_H` 与实际秩函数，避免重新用统一最坏上界把刚省下的结构丢掉。
