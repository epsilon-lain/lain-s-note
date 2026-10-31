# 任意常数改善、Lean 反证链与联合变分判据

日期：2026-10-31。**未证明或证伪 P2848。** 本次实际追读 OpenAI 的 Lean 证明，并尝试在指数 2 处保留常数改善。新内容是第 4 节的联合下界、第 5 节对此具体下界的明确门槛，以及第 6 节指出的被松弛掉的约束。参数端点障碍和连分数反例在旧记录中已有，不作为新发现。

## 1. 要反证的假设及量词

用户提出的假设应写成
\[
\tag{H}
\forall C\ge1\ \exists p\in\mathbb Z,\ q\in\mathbb Z_{>0}:
\quad |\pi-p/q|<\frac1{Cq^2}.
\]
其中 p、q 随 C 变化。这等价于 \(\inf_{q\ge1}q\|q\pi\|=0\)，也等价于 P2848 的坏逼近失败表述。若从 (H) 导出矛盾，得到的是 π 坏逼近，即证伪 P2848；并非证明 P2848。

(H) 允许分母超过任意**事先指定**的门槛 Q：有限个 q≤Q 的 \(q\|q\pi\|\) 有正的最小值，选足够大的 C 即可排除它们。因此可以继续选择很大的、相互分离的分母。

不能由此要求返回的 q 满足 \(C\ge q^\varepsilon\)，其中 ε>0 事先固定。令
\[
\kappa=q^2|\pi-p/q|,\quad \lambda=\log(1/\kappa),
\quad w=\lceil\log q\rceil.
\]
κ→0 仅给出 λ→∞；它不保证 λ/w 有固定正下界。采用 OpenAI 提交的 μ(π)=2 结论时，沿分母趋于无穷、κ≤1 的序列，反而有 λ/w→0。

旧记录第 3 节已有明确对照 \(\alpha=[0;1,2,3,\ldots]\)：其收敛分母 \(q_n\ge n!\)，且
\[
\frac1{n+3}<q_n^2|\alpha-p_n/q_n|<\frac1{n+1}.
\]
故它满足 (H) 的 α 版本，同时 μ(α)=2。这是对“μ=2 自动排除任意常数改善”这一步的反例，**不是 π 的反例**。

## 2. 实际 Lean 构造了什么

核对版本：OpenAI `math` 提交
[`4eaee87548250b2682affda7d05d947c05554bae`](https://github.com/openai/math/commit/4eaee87548250b2682affda7d05d947c05554bae)。

主要调用链为：

1. [`Main.lean`](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/lean/OAI/NumberTheory/PiExponent/Main.lean)：`pi_eventual_lower_bound`、`main`。
2. [`InterpolationConsequence.lean`](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/lean/OAI/NumberTheory/PiExponent/Approximation/InterpolationConsequence.lean)：将插值满射和解析估计交给行列式反证。
3. [`DeterminantContradiction.lean`](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/lean/OAI/NumberTheory/PiExponent/Approximation/DeterminantContradiction.lean)：`no_fixed_data_of_interpolation_and_analytic_aggregate`，以及 `Comparison.determinant_bounds_inconsistent` 的调用。
4. [`Parameters.lean`](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/lean/OAI/NumberTheory/PiExponent/Approximation/Parameters.lean) 与 [`AdmissibleParameters.lean`](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/lean/OAI/NumberTheory/PiExponent/Approximation/AdmissibleParameters.lean)：确保两个界严格交叉的选参。

这里没有构造一个新的“反例数”。反证从固定 ν>2 下分母无界的异常好逼近出发，构造非零满行数子式 Δ。清分母后得到非零高斯整数，模至少为 1；另一方面，将有理近似中心平移到指数函数的精确周期 \(2\pi i j\)，用小位移与重复 Taylor 次数造成的抵消，证明同一个子式过小。这两个界冲突。

原 `Parameters ν` 要求
\[
0<\theta<A<B<1,\quad B<1/C_{\rm geom},\quad
B<C_{\rm geom}\theta,\quad 1-\theta<\nu(A-\theta).
\]
因此 A²<B²<θ。若 ν=2，则所需正间隙满足
\[
2(A-\theta)-(1-\theta)=2A-\theta-1
<-(1-A)^2<0.
\]
所以原来的 `Parameters 2` 不可能存在。注意这里的 \(C_{\rm geom}\) 是原论文的几何参数，不是 (H) 的逼近常数 C。这一障碍在旧记录第 10 节已有；本次补上实际 Lean 对应位置。

本次是源码和相关论文估计的定向核对，**未运行 Lean，也没有逐层审计整个依赖闭包**。下面的新推导尚未形式化。

## 3. 将实际 κ 因子接入原证明

[`MatrixTranslationBounds.lean`](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/lean/OAI/NumberTheory/PiExponent/Approximation/MatrixTranslationBounds.lean) 的 `actual_error_exp` 把指数 ν 的逼近转为中心位移估计。若保留实际 κ_i≤1，则相同初等步骤给出
\[
|\epsilon_i|\le2K\kappa_iq_i^{-2}
\le\exp\{\log(2K)+2-2w_i-\lambda_i\}.
\]
这里先选 q_i≥3，使 \(w_i=\lceil\log q_i\rceil>0\)，并令 \(\lambda_i=\log(1/\kappa_i)\)。

沿用[旧记录第 11–13 节](../OPENAI_PI_EXPONENT_P2848_2026-10-27.md)，行 ρ 的原指标为 βρ，展开中选到 aρ≥βρ、0≤dρ≤aρ−βρ，且 \(w\cdot d_\rho<H/F_0\)。令 n_a 为选到横向指标 a 的行数，M 为总行数，\(c=(\log2)/4\)。保留两种收益的分数为
\[
S_H=\frac{2\sum_\rho w\cdot(a_\rho-\beta_\rho)
+\sum_\rho\lambda\cdot(a_\rho-\beta_\rho-d_\rho)
+c\sum_a n_a^2}{MH}.
\]
原文逐项估计给出
\[
\frac{\log|\Delta_H|}{MH}\le E_{\rm an}(2)+o(1)-\inf S_H,
\qquad
\frac{\log|\Delta_H|}{MH}\ge-\bar c_H+\bar b_H-E_{\rm ar}.
\]
因此，足以矛盾的目标是
\[
\liminf_{H\to\infty}(\inf S_H+\bar b_H-\bar c_H)>E,
\quad E=E_{\rm ar}+E_{\rm an}(2).
\]
以上接口在旧记录已有；下一节对它给出新的显式下界。

## 4. 不再分成两种情况：一个联合变分下界

固定整数 m≥1，所有权重、v₀、F₀ 为正，并取 0<θ<1。设
\[
s_i=\lambda_i/w_i,\quad \tau=\max_i s_i,\quad
v_i=w_i+\lambda_i/2,\quad P=\prod_{i=1}^m(1+s_i/2).
\]
尾项代价满足 \(\sum_\rho\lambda\cdot d_\rho\le\tau MH/F_0\)。因此有
\[
S_H\ge\frac{\sum_a\{c n_a^2+2(v\cdot a)n_a\}
-2\sum_\rho v\cdot\beta_\rho}{MH}-\frac\tau{F_0}.
\]

**有限不等式。** 对任意 t>0、非负整数 n_a（实际上非负实数也可以）、Σn_a=M，配方给出
\[
\tag{4.1}
\sum_a\{c n_a^2+2(v\cdot a)n_a\}
\ge2tHM-\frac1c\sum_{a\in\mathbb Z_{\ge0}^m}(tH-v\cdot a)_+^2.
\]
当 v·a≤tH 时逐项差为 \(c(n_a-(tH-v\cdot a)/c)^2\)；当 v·a>tH 时逐项差非负。实际允许的 a 集合之外令 n_a=0，故对全格点求和仍然有效。

固定全部数据后取 H→∞，Riemann 和给出
\[
\sum_a(tH-v\cdot a)_+^2
\sim\frac{2(tH)^{m+2}}{(m+2)!\prod_i v_i},\qquad
M\sim\frac{K\theta^m H^{m+1}}{(m+1)!v_0\prod_iw_i}.
\]
因此 (4.1) 除以 MH 后的极限下界为
\[
2t-\frac{2v_0}{cK\theta^mP(m+2)}t^{m+2}.
\]
在
\[
T=\left(\frac{cK\theta^mP}{v_0}\right)^{1/(m+1)}
\]
取最大值得 \(2(m+1)T/(m+2)\)。同时
\[
\frac1{MH}\sum_\rho v\cdot\beta_\rho
\longrightarrow\frac\theta{m+2}\sum_i(1+s_i/2),
\qquad\bar b_H\longrightarrow\frac{m\theta}{m+2}.
\]

令 \(\rho=Kw_0\theta^m/v_0\in(0,1)\) 为行数与总列数的极限比例。旧记录第 12 节的通用选列界记为 \(U_m(\rho)\)：它是密度 \(m(m+1)x^{m-1}(1-x)\) 的最高 ρ 比例的条件均值，满足
\[
\limsup\bar c_H\le U_m(\rho),\qquad
U_m(\rho)\ge\frac m{m+2}.
\]
我们得到明确的充分判据
\[
\boxed{\quad
G:=\frac{2(m+1)}{m+2}T
-\frac\theta{m+2}\sum_i(1+s_i)
-U_m(\rho)-\frac\tau{F_0}>E.
\quad}
\tag{4.2}
\]
若数据满足原插值和固定数据下的行列式估计，且 (4.2) 成立，就得到反证。这里不再使用高、低指标两分支，也没有直接丢弃 κ 收益。

## 5. 实际尝试后的结果：这个充分判据仍有严格门槛

**以下是对 (4.2) 这个具体判据的限制，不是对真实 S_H 的上界。** 原论文的误差成本至少为
\[
E\ge\frac{100K}{w_0}+\frac{\log2}{v_0}.
\]
记 \(b_m=m/(m+2)\)、\(\bar s=m^{-1}\sum_i s_i\)，并设
\[
\delta_{m,\rho}=\left(\frac{c\rho}{400(m+2)^2}\right)^{1/m}<1.
\]
则不论选择怎样的 K、θ、w₀、v₀（满足上述比例关系），均有
\[
\boxed{\quad
G-E\le b_m\left\{
\delta_{m,\rho}\left(1+\frac{\tau^2}{4(1+\tau)}\right)-1
\right\}.
\quad}
\tag{5.1}
\]

证明如下。置
\[
B=\frac{100K}{\rho}+\frac{\log2}{K\theta^m}.
\]
因为 \(T^{m+1}=cP\rho/w_0\)，有 \(E\ge BT^{m+1}/(cP)\)。对 T>0 取上确界得到
\[
\frac{2(m+1)}{m+2}T-E
\le 2b_m\left(\frac{2cP}{(m+2)B}\right)^{1/m}.
\]
两项的算术–几何平均不等式给出
\[
B\ge20\sqrt{\frac{\log2}{\rho}}\theta^{-m/2},
\]
利用 \(c=(\log2)/4\)，上式右侧至多
\[
2b_m P^{1/m}\sqrt{\delta_{m,\rho}\theta}.
\]
再扣去 \(b_m(1+\bar s)\theta+U_m(\rho)\)，舍弃非正项 −τ/F₀，关于 √θ 配方，得
\[
G-E\le b_m\delta_{m,\rho}\frac{P^{2/m}}{1+\bar s}-U_m(\rho).
\]
最后由 AM–GM、\(0\le\bar s\le\tau\)，
\[
\frac{P^{2/m}}{1+\bar s}
\le\frac{(1+\bar s/2)^2}{1+\bar s}
=1+\frac{\bar s^2}{4(1+\bar s)}
\le1+\frac{\tau^2}{4(1+\tau)},
\]
再用 \(U_m(\rho)\ge b_m\)，得到 (5.1)。证明中扩大了可选参数范围，因此原插值所需的额外条件不会破坏这个上界。

由此可知：不计 κ 收益时，即 τ=0，(4.2) 对任何固定维数和任何列比例均不可能成立。即使保留 κ，想让这个判据成立也**必须**满足
\[
\boxed{\quad
m\log\left(1+\frac{\tau^2}{4(1+\tau)}\right)
>\log\left(\frac{400(m+2)^2}{c\rho}\right).
\quad}
\tag{5.2}
\]
固定 m、ρ 时，τ 有严格正的门槛。而 μ(π)=2 允许 κ→0，却令足够晚的分母对应 τ→0。因此“C 任意大”本身并不保证这次构造成功。让 m 随 τ 增大也不能被直接宣布可行：维数、插值阈值、分母和实际 τ 的选择依赖需要重新解决。

值得保留的细节是：在这个估计中，一阶收益与一阶代价相消，剩余提升由 \(\tau^2/[4(1+\tau)]\) 控制。第 5 节是对本次新判据的可证失败范围，不是 π 坏逼近与否的证据。

## 6. 下一步有一个具体的、尚未使用的约束

(4.1) 只用了 Σn_a=M。真正的展开还要求每行 aρ≥βρ（逐坐标）。因此对于任意向下闭集 I⊂ℕᵐ，有
\[
\boxed{\sum_{a\in I}n_a
\le\#\{\rho:\beta_\rho\in I\}.}
\tag{6.1}
\]
因为 aρ∈I 且 βρ≤aρ，必有 βρ∈I。它禁止把原本具有高指标的行任意塞进低指标组。自由分布优化会忽略这个限制，可能给出严格过小的收益下界。

例如 s_i=0 时，自由连续极值分布在 x=w·a/H 上是支撑 [0,T] 的缩放 Beta(m,2) 分布；原 β 分布对应支撑 [0,θ] 的同一分布。若 T<θ，自由极值把太多质量放在靠近 0 的位置，违反 (6.1) 的径向下闭集版本。这说明遗漏约束可能实际生效；**尚未证明恢复它后能够完成反证**。

后续有两处明确入口：保留 (6.1) 的单调分配约束，或利用非零子式结构把实际列代价压到通用最坏界 U_m(ρ) 以下。这里没有把任一方向称为已完成的 π 特有估计。

## 7. 检查范围与复现

`verify.py` 用有理数精确检查有限配方不等式，并用独立的最小增量分配算法求有限整数最优值；随后检查格点二次截断和的渐近常数，以及 (5.1) 的数值实例。有限有理检查取任意允许的 c=1/6；第 5 节数值使用 c=log(2)/4。输出见 `results.json`。

这些计算检查公式和常数，不证明 π 上存在使判据闭合的数据。严格论证是上面的配方、格点 Riemann 和与不等式链；没有增加新的 π 有理逼近纪录。

相关前序记录：

- [OpenAI μ(π)=2 与 P2848，尤其第 3、10–13 节](../OPENAI_PI_EXPONENT_P2848_2026-10-27.md)
- [自身进制迁移与旧参数门槛](../self-base-theta-2026-10-31/OPENAI_SELF_BASE_TRANSFER.md)
- [端点行列式约因子与整数成本](../endpoint-reduction-2026-10-31/README.md)
- [原论文行列式估计](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/determinant.tex)
