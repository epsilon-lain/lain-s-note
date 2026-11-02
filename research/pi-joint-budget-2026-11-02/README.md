# 联合参数、非零输运与周期子式的剩余障碍

2026-11-02。接续[完整求和与平均次数集中](../pi-aggregate-support-2026-11-02/README.md)。本轮检验两个具体缺口：让维数与参数一起变化，是否能消除旧成本缺口；满足完整时间约束的目标行，是否仍能实现矩最优的类型分布。

**结果：前者在当前标量估计内仍有障碍；后者可以构造，并且完整输运行列式非零。真正周期子式的非零性、大小，以及完整求和的定量抵消仍未解决。没有证明 π 是坏逼近数或非坏逼近数，没有得到新的有理逼近纪录，也没有独立认证上传原稿的 Theorem 2.1 或 μ(π)=2。**

这里的“碰撞”指不同目标行具有相同横向次数，“时间切割”指选定列空间在每一个 Taylor 阶数上的必要容量条件。输运行列式负责把源行变成目标行；周期行列式负责这些目标行在同一批多项式列上是否真正独立。前一个非零不能代替后一个非零。

## 1. 同一批列的成本不能拆开取最小值

对每个横向类型 α，设选中的原始列 Y^h X^α 有 N_α 个。它们的非负整数 h 互不相同，所以有精确有限不等式

\[
\boxed{\sum_C h\ge\frac12\sum_\alpha N_\alpha(N_\alpha-1).}
\]

这把时间增长成本与类型集中程度直接连接起来：让横向次数集中，会迫使同一类型采用更多不同时间次数。两项不能使用不同列族的最优值。

在半密度、v₀=2Kθ^m w₀、R=14K 的既定归一化下，将这个关系与真实 LCM 成本结合，可以对每个允许的标量列均值，显式找到一个目标均值，使当前上界表达式为正。它覆盖任意联合参数序列 θ_m→1，不再要求 θ、B、C 固定。

这证明的是**标量质量／均值碰撞松弛不能给出统一负余量**。其正值不是实际行列式的下界。任意不均衡列集的逐坐标相容性、全部时间结构、更强范数、减小清分母代价和完整求和相消均未被排除。固定参数的严格正余量，也不等于联合参数下有统一正下界。

详见 [UNIFORM_SCALAR_OBSTRUCTION.md](joint-theta-endpoint/UNIFORM_SCALAR_OBSTRUCTION.md)。

另得到同列联合成本的精确分位数公式和保留 D_i≥0 的对偶公式。空间系数相同时，对所有时间成本系数 a>0 一致有

\[
Q_m(a,1/2)=1+a\tau_m+O(1/m),
\qquad \tau_m\sim\frac{1-\log2}{m}.
\]

即使 a 随维数指数增长，这个以精确 τ_m 写出的绝对误差仍为 O(1/m)。将 τ_m 换成它的首项时，另需计入 O(a/m²)。θ 极快趋近 1 的一个闭式分支违反多中心源包的必要计数条件，已经明确排除为候选构造。

- [JOINT_COLUMN_QUANTILE.md](joint-column-variational/JOINT_COLUMN_QUANTILE.md)
- [SPECIAL_TIME_COST.md](joint-column-variational/SPECIAL_TIME_COST.md)

## 2. 矩最优分布确实能通过所有时间切割

选有理膨胀 λ=L/Q>1，A=λθ<1。在守恒阶数 N=s+|β| 的每个桶内，采用

\[
a_i=\lfloor\lambda\beta_i\rfloor+r_i(N),
\qquad t=N-|a|,
\]

其中相位 r(N) 周期遍历 {0,…,L−1}^m。低阶桶保留原行，保证 t 非负。

在时间预算条件 w₀>[m(L−1)+2](m+1)v₀ 下，对所有充分大的 H，这个构造同时保留不同目标行、逐坐标 a≥β、精确守恒阶数，以及原始列空间的所有时间切割。Newton 插值给出完整 Pascal 输运子式的正乘积公式，因此不只是支撑图上存在匹配，**整个输运行列式确实非零**。

相位循环填补固定膨胀留下的格点空缺。先固定参数取高度 H→∞，再令横向权重充分大，目标类型分布接近先前矩优化的尖锥分布，额外 Pascal 系数成本可以趋零。关键固定指数参数族中，几乎全部行可以离开原源包，碰撞成本也可以趋零。这些结论仍不提供真正周期子式的非零性或大小。

详见 [PHASED_DILATION_COMPATIBILITY.md](time-moment-compatibility/PHASED_DILATION_COMPATIBILITY.md)。有限验证使用完整的 77,360 行二维源包，检查 101 个非平凡时间切割、217 个不同 Pascal 矩阵及同列逐坐标次数余量。它不计算或认证 77,360 阶周期行列式。

## 3. 独立性本身有不可忽略的次数代价

把 K 个中心沿固定方向缩成 jZξ，使用最大空间权重坐标作切向锚点。完整源导数包中的同一时间／法向块含 n 个连续切向导数。该块在 K 个中心合并时至少消失 binom(K,2)n² 阶。因此，对固定且不依赖 Z 的任意合法多项式列，

\[
\boxed{\operatorname{ord}_0 P_C(Z\xi)
\ge V_H:=\binom K2\sum_{s,\nu}n_{s,\nu}^2.}
\]

这比仅逐组作中心有限差分保留了更多导数重叠信息。固定参数后

\[
\frac{V_H}{M_HH}\longrightarrow
\frac{K-1}{(m+2)v_{\max}}.
\]

结合列的加权次数预算，若固定参数下沿无界的 H 存在非零源子式，就得到必要权重分离条件

\[
\boxed{\frac{w_{\max}}{w_{\min}}
\ge\frac{\theta(K-1)}{m(1-\theta)+2}.}
\]

它是必要条件，不是插值充分条件。在既有强分离权重下可能很弱；也不控制剩余因子在真实周期 2πi 处的大小。一维合法 Hermite 阶梯可以达到上述消失阶。

详见 [CENTER_COALESCENCE.md](center-coalescence/CENTER_COALESCENCE.md)。

## 4. 周期多项式的根没有自动的简单半平面保护

完整源包和合法原始列，可以产生

\[
P(Z)=Z^4(Z^2+3),\qquad
P(Z)=\frac{Z^6}{4}(4Z^4-57Z^2-276).
\]

第一个有虚轴根；第二个同时有正负实根和正负虚根。另一个十二行例子的非零根分布于四个开象限。因此，不能仅靠完整行包就假设根全在某个边界经过原点的半平面，或整个虚轴无零点。这些有限例子没有根接近 2πi；也不说明实际周期行列式很小。

对于共同径向扰动，完整聚合的相对值有精确根乘积表示。设 P 的非零根为 ρ_j，n=deg P≤D，ω=2πi，则

\[
\frac{\mathcal F(\zeta)}{\mathcal F(0)}
=(1-\zeta)^{D-n}
\prod_j\left(1+\frac{\rho_j\zeta}{\omega-\rho_j}\right).
\]

这给出一个明确的后续定量对象：根到真实周期的加权接近量

\[
\mathcal S_C=(D-n)+\sum_j\frac{|\rho_j|}{|\omega-\rho_j|}.
\]

当每个扰动因子足够小时，有 log|F(ζ)/F(0)|≥−2|ζ|S_C。尚未证明关键高维列族满足所需根分离，也未取得基准 |P(ω)| 的绝对下界。

详见 [PERIOD_ROOT_GEOMETRY.md](period-root-geometry/PERIOD_ROOT_GEOMETRY.md)。

## 5. 下一步需要什么

本轮把“目标支撑不可能实现”这一简单解释排除了一部分：尖锥分布可以同时满足全部时间切割和非零输运。但无法据此断言它会贡献非零的真实周期项。

有两个仍具体可检验的方向：在相位膨胀目标上，建立同一批列的真实周期秩判据；或在完整源行列式的倒数多项式上，取得保留算术预算的定量根分离／完整相消估计。中心合并的次数成本必须保留，不能只验证支撑匹配或形式非零。

所有高维结论采用分层极限：每个维数先选定权重与相位周期，再取充分大的高度。没有证明可任意交换维数、权重、相位周期和高度极限。以上是本轮推导和有限证书，不声称文献首创。

## 6. 复现与复核

在本目录运行以下命令；全部只使用 Python 标准库，无外部数据下载：

```bash
python joint-theta-endpoint/verify_uniform_scalar_obstruction.py
python joint-column-variational/verify_joint_quantile.py
python joint-column-variational/verify_special_time_cost.py
python time-moment-compatibility/verify_phased_dilation.py
python period-root-geometry/verify_period_root_geometry.py
python center-coalescence/verify_center_coalescence.py
```

每个脚本将精确检查结果写到同目录 JSON。六份独立代理复核记录保留在对应目录；这些复核不是外部同行评审。有限枚举用于检查公式和实例，无限参数结论依据各正文证明。

[CHECKS.json](CHECKS.json) 记录本次实际运行状态、输出摘要和除该清单自身以外的文件 SHA-256。上一轮的[矩碰撞](../pi-aggregate-support-2026-11-02/MEAN_DEGREE_AND_COLLISIONS.md)、[内在时间过滤](../pi-aggregate-support-2026-11-02/FILTERED_TIME_SUPPORT.md)和[倒数多项式恒等式](../pi-aggregate-support-2026-11-02/RECIPROCAL_RESUMMATION.md)是本轮所用前置结果。
