# 无穷维延拓必须同时保留什么：极限、非零性与算术成本

2026-11-01。此文核对固定稿件 π017 的维数与极限顺序，并给出一个真正相容的无限维有理算子模型。结论不是“无穷维无用”，而是识别哪些性质需要新的统一估计。未证明 π 的坏逼近或非坏逼近，也没有独立认证上游插值满射或 Lean 内核。

固定来源为 `openai/math@4eaee87548250b2682affda7d05d947c05554bae`，稿件目录 `preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026`。本轮读取 `build/sections/conclusion.tex` 1–156、`determinant.tex` 60–85 与 265–340、`interpolation.tex` 15–85 与 178–214，以及本项目[改写后的估计](../pi-endpoint-rewrite-2026-11-01/rewritten_estimate.md)。

## 1. 原论文不是只算有限小数位

原文的两个增长参数要分清：

- m 是横向变量的个数，维数可以依 ν>2 和误差要求选得很大。
- H 是多项式预算；选定维数、分母、权重和中心后，H 仍然趋于无穷。

[`conclusion.tex` 5–7](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/conclusion.tex#L5)明确规定这个次序；[66–96](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/conclusion.tex#L66)先研究 m 增大如何降低若干费用、提高碰撞收益，然后固定一个足够大的 m；[133–135](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/conclusion.tex#L133)再次说明 H 最后趋于无穷。

因此“拓展无穷维”的实质问题是：建立一组相互兼容的截面，让 m 也增长，同时证明几何、算术及余项控制不失效。它不是简单地把有限位数增加到无限位数。

## 2. 一个所有有限截面都满秩、无限极限却有核的精确模型

令空间为实 Hilbert 空间 `ℓ²`，其元素是平方和有限的实数序列。取

\[
 v_j=\frac{\sqrt3}{2^j},\qquad j=1,2,\ldots.
\]

因为 `3Σ_{j≥1}4^{−j}=1`，所以 `||v||=1`。定义有界算子

\[
 T=I-vv^*,\qquad Tx=x-\langle x,v\rangle v.
\]

它有真实的非零核向量 v：`Tv=0`。令 `P_m` 保留前 m 个坐标，`v_m=P_m v`。前 m 个坐标上的压缩矩阵为

\[
 T_m=I_m-v_m v_m^*,\qquad
 (T_m)_{ij}=\delta_{ij}-\frac3{2^{i+j}},\quad1\le i,j\le m.
\]

**这些矩阵来自同一个无限算子，而不是每个维数重新拼接。** `T_{m+1}` 左上 m×m 块正好等于 `T_m`；所有条目都是有理数。

有限几何级数给

\[
 \|v_m\|^2=1-4^{-m}.
\]

在 `v_m` 方向，`T_m v_m=4^{-m}v_m`；在它的正交补上，`T_m` 等于恒等。故

\[
 \boxed{\det T_m=4^{-m}>0,\qquad \sigma_{\min}(T_m)=4^{-m}.}
\]

每个有限截面都满秩，但独立性越来越不稳定。对 m≥2，最大特征值为 1，所以二范数条件数为 `4^m`；m=1 时条件数是 1。不能把最小特征值的倒数直接写成 m=1 的条件数。

若把有限矩阵在剩余坐标上扩展为恒等，得到 `S_m=I-v_m v_m^*`，则每个 `S_m` 都可逆，且

\[
 \|S_m-T\|\le
 (\|v_m\|+\|v\|)\|v_m-v\|
 \le2\,2^{-m}\longrightarrow0.
\]

因此即使是算子范数收敛，可逆算子仍能趋于有核算子；失败的条件正是逆算子的范数没有统一上界。这里不能将零延拓的 `P_mTP_m` 与上述恒等延拓的 `S_m` 混同：前者本身在无限空间上已有尾部核，且不以算子范数逼近 T。

若另能证明 `||S_mx||≥c||x||`，其中 c>0 不依赖 m，则取极限得到 `||Tx||≥c||x||`，从而 T 没有核。**有限截面的满秩，与统一的独立性下界，是两个不同要求。**

## 3. 同一个模型说明：极限行列式为零不自动给算术矛盾

虽然 `det T_m→0`，但它的最小行列式分母恰为 `4^m`。清掉这个分母后

\[
 \Gamma_m=4^m\det T_m=1.
\]

小量与分母成本精确抵消；没有任何非零整数小于 1。

还可以在矩阵条目层面清分母。令 `D_m=diag(2,2²,…,2^m)`，则

\[
 B_m=D_mT_mD_m,
 \qquad(B_m)_{ij}=4^i\delta_{ij}-3\in\mathbb Z.
\]

于是

\[
 \boxed{\det B_m=4^{m(m-1)/2}\ge1.}
\]

因此，无穷维方向确实可以产生一个真实退化、所有截面却非零的对象；但**只证明某种未归一化体积任意小，仍不等于构造某个给定实数的超好有理逼近**。必须说明与有理 p/q 的对应及其实际高度，并保留清分母代价。

若一列非零高斯整数 Γ_m 真正收敛到 0，当然会在某个有限 m 就出现 `|Γ_m|<1` 而矛盾。困难恰在于，归一化后的实数/复数小量趋零，未必意味着 Γ_m 趋零。上例给出最直接的精确控制。

## 4. 固定一组无限前缀，不能直接照搬原文的分离阈值

这一节只限制原证明的**这组充分阈值**，不证明所有无限维插值都不可能。

固定 `0<θ<1`、正数 `w₀,v₀`，并设想固定一个无限权重序列 `w₁,w₂,…`，让每个前缀都直接满足原文同一形状的分离条件。原文 [`interpolation.tex` 180–203](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/interpolation.tex#L180)要求选择 σ_m>0，满足

\[
 (1+\sigma_m)\theta<1,
 \qquad C_{*,m}=\max_{1\le k\le m}
   2k!\left(\frac{m+2}{\sigma_m}\right)^k.
\]

这里准确的常数是最大值，不能不检查条件就把它替换为 k=m 项的等号。

取分离条件中的单元素集合 `A={1},B={0}`。最大不同正指标为 1，属于 A，所以必须有 `w₁>C_{*,m}v₀`。但由 k=1 项及 σ 的限制，

\[
 \frac{w_1}{v_0}>C_{*,m}
 \ge\frac{2(m+2)}{\sigma_m}
 >\frac{2(m+2)\theta}{1-\theta}.
\]

固定的 `w₁/v₀` 不可能对所有 m 满足这个不等式。这说明不能将“对每个 m，重新选择全部权重”误当成“同一个固定无限对象的全部前缀”。如果允许 θ、v₀ 或坐标规范依 m 变化，或证明更强的插值定理，本节限制便不直接适用；届时需要另外定义截面之间的兼容性。

## 5. 变维数时，显式余项要重新核对

此前[改写估计](../pi-endpoint-rewrite-2026-11-01/rewritten_estimate.md)对每个固定参数组给出显式有限余项

\[
 \mathcal R_{m,H}=2\ln(M!)+M\{c_R-\ln(1-\gamma)\}
             +mM\ln(\lfloor H\rfloor+1),
\]

其中 M 为同一子式的行数、`c_R>0`、`0<γ<1`。因此

\[
 \frac{\mathcal R_{m,H}}{MH}
 \le \frac{2\ln M+c_R-\ln(1-\gamma)
                    +m\ln(\lfloor H\rfloor+1)}{H}.
\]

对固定 m，原稿的格点计数给 `M=O(H^{m+1})`，所以该比值趋零。若 m、权重、半径都随 H 变化，下面是一组可以逐项核对的**充分条件**：

\[
 \frac{\ln M}{H}\to0,\quad
 \frac{m\ln(H+1)}{H}\to0,\quad
 \frac{c_R-\ln(1-\gamma)}H\to0.
\]

必须使用实际 M 或已证明对该数据族统一的计数界。原稿对固定维数的 `M∼C_mH^{m+1}`，不能未经一致误差估计就代入变维数极限。

若坚持要求上面**这一显式余项**除以 MH 趋零，那么由于各项非负，`m ln(H+1)/H→0` 也是它的必要条件。然而，这不是所有无穷维方法的必要条件，更不构成无穷维不可能定理：更精确的项数估计可以改变余项，或者主收益本身可以增长到超过一个不趋零的余项。

此外，原插值定理 [`interpolation.tex` 54–56](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/interpolation.tex#L54)的最终 H 门槛允许依赖中心。变维数方案还需同时保证每个选中 `(m,H)` 真正达到相应门槛，并满足可整除要求。

## 6. 真正有用的无穷维结论应给出有限证书

设对一族允许数据 `(m,H)`，同一非零子式满足

\[
 \ln|\Delta_{m,H}|\ge-\mathcal A_{m,H},\qquad
 \ln|\Delta_{m,H}|\le
 \mathcal B_{m,H}-\mathcal G_{m,H}+\mathcal R_{m,H}.
\]

如果无穷维分析证明某条数据序列上

\[
 \liminf\frac{\mathcal G_{m,H}-\mathcal A_{m,H}
                     -\mathcal B_{m,H}}{MH}>0,
 \qquad \frac{\mathcal R_{m,H}}{MH}\to0,
\]

那么某个充分大的**有限截面**已经满足

\[
 \mathcal G_{m,H}>\mathcal A_{m,H}
                      +\mathcal B_{m,H}+\mathcal R_{m,H},
\]

从而上下界冲突。这是一条充分的有限证书原则，不要求无穷维必须退回有限维才有价值；无穷维结构可以负责发现并证明那个统一正余量。严格正的统一极限余量也不是唯一方式，直接找到一个有限正余量同样足够。

相反，只有“未归一化体积在极限为零”“每个有限截面都有非零子式”或“某个几何收益随维数增大”，还没有闭合反证。上面的有理投影模型同时展示了这些条件为何不能被静默替换。

## 复现与本轮边界

运行 `python3 research/infinite-dimensional-pi-2026-11-01/limit_audit_lab.py`；结果在 [limit_audit_results.json](limit_audit_results.json)。脚本用 Fraction 独立消元计算前 10 个压缩矩阵与整数化矩阵的行列式，核对嵌套块、精确特征向量与归一化恒等式。一般结论来自上文证明，未从有限检查外推。

本轮没有给出 π 的异常逼近实例，没有建立新的无穷维插值定理，没有把无限小行列式直接当成非零整数，也没有将源稿关于固定 ν>2 的结果自动推广到端点。
