# π 坏逼近性的构造攻击：永久纪录反例与常数增益障碍

2026-11-01。接续当天的无理尾项、无穷维端点和原稿迁移检查。

**没有证明 π 是坏逼近数，也没有构造出否定它的逼近序列。本轮完成的是一个普遍放大方案的严格反例，以及原稿参数方案的新定量障碍。** 没有扩大原先的有限搜索范围，也没有把已知连分数纪录重新包装为新纪录。

## 1. 本轮真正检验的目标

对整数 p、正整数 q，定义标准化误差

$$K_\pi(p,q)=q^2|\pi-p/q|.$$

π 是坏逼近数，意味着存在固定 c>0，使所有 p,q 都满足 K_π(p,q)≥c。反方向则需要构造一列分母趋于无穷的有理数，且 K_π→0。单个极小值不能决定这两个方向。

原稿所声称的 μ(π)=2，允许 K_π 沿某些序列缓慢趋零。这个指数结论即使成立，也没有单独决定坏逼近性。

## 2. 严格反例：极佳的当前逼近，可以永远不再改善

取相邻收敛分数 P/Q、R/S，令 r=S/Q，并把目标的未知完整商记为 α>1。目标写成

$$x=\frac{P\alpha+R}{Q\alpha+S},\qquad |PS-QR|=1.$$

对整数线性组合 (p,q)=(aP+bR,aQ+bS)，有精确式

$$\frac{K_x(p,q)}{K_x(P,Q)}=(a+br)|a-b\alpha|.$$

它把消去效果与最终分母成本同时保留下来。由此得到：若 α=[A;1,1,1,…]、A≥3，则当前 K<1/A，而所有更大分母都不能改善当前 K。取前缀 [0;2] 时，当前项还是一个真正的永久全局纪录，A 可以任意大。

因此，“已经找到很深的内凹，所以一定能由相邻种子自动制造更深内凹”这个一般命题为假。**这是对一般机制的反例，不是对 π 的反例。**

还得到一个有限系数覆盖界：只知道 α∈(A,A+1)，并允许全部原始整数系数 1≤b≤B 时，成功尾部区间的总长度至多 2H_B/A，其中 H_B=Σ_(b=1)^B 1/b。它是确定性的覆盖上界，不能解释为 π 的概率。

- [完整推导](local_amplifier/RESULT.md)
- [精确验证器](local_amplifier/verify.py)与[结果](local_amplifier/verification.json)：561 个恒等式、40 组区间并集检查通过。

## 3. 直接保留质量因子后，原稿的成本仍有门槛

令 s_i=ln q_i、λ_i=ln(1/K_π(p_i,q_i))，并选择共同指数余量 ε>0，使 λ_i≥εs_i。这样才能使用原稿指数 ν=2+ε 的估计。

沿原稿 §4 的参数方案，即使维数和参数随 ε 调整，现有行列式判据若能成立，仍必须有

$$
w_*>
\frac{16(1+\varepsilon)}{\varepsilon^3}
\ln\frac{64(1+\varepsilon)^2}{c\varepsilon^4},
\quad c=\frac{\ln2}{4},\quad w_*=\min_i\lceil\ln q_i\rceil.
$$

核心比较是 ln(B/A)<(ε/2)ln C：碰撞收益增长时，中心数量必须付出更大的增长成本。这个成本再经由误差项 ln K/w_* 进入比较；这里 K 是原稿的中心数，与本页标准化误差 K_π 不同。

设最小分母对应 s=ln q、λ=ln(1/K_π)。如果 λ≤s^(2/3)，该批分母就不可能让上述整套充分条件成立。真实 ln q 与 ceil(ln q) 的差已保留。

非坏逼近只保证某条序列上 λ→∞，并不保证 λ 增长快到满足这个门槛。因此，“让质量因子越来越大，再适配原稿参数”尚不能补齐问题。

**范围：只排除压成共同指数后沿用原 §4 参数和现有误差表达式的方案。** 没有证明实际误差不可改善；也没有排除保留每个坐标质量、改变插值空间或其他估计。

- [完整证明](endpoint/CONSTANT_GAIN_OBSTRUCTION.md)
- [有理区间验证器](endpoint/verify.py)与[结果](endpoint/verification.json)：四组允许参数及对数不等式通过；通用结论依据正文证明。

## 4. 补查旧报告的逐坐标约束缺口

旧报告第 7 节曾放松每行类型 a 必须逐坐标大于等于 β 的条件。本轮构造了满足该条件的真实整数格点分配，证明只把这项限制放回去，仍不足以使现有成本判据在指数 2 处成立。

进一步，对中心 j=0 的行作修补后，甚至可以保证每行选到的平移标量非零；额外成本仍小于现有误差预算。这封住了一个具体的候选修补方案。

**没有证明这些分配对应非零行列式。** 所选列的支持、Taylor 阶数和整个行列式的非零性仍可能排除这些低成本分配；这是本模型尚未捕获的结构。因此，该结果不能推广为“原稿所有改进都不可能”。

[完整格点构造、三个参数区间及非零标量修补](ROW_CONSTRAINT_ENDPOINT.md)。这部分依据解析证明和独立复核；另有 [12 组精确格点二阶矩核验](discrete_assignment_verification.json)，检查分配约束、碰撞计数及中心修补，没有用有限数值试验替代格点极限。

## 5. 有限分辨率的反模型

标准连分数理论允许把任意 π 的有限前缀延续成两类数：两者的无理测度都是 2，一类坏逼近，另一类不是。第二类可令尾部第 2^j 位为 j+1，其余为 1，得到一列

$$K=\Theta(1/\ln\ln q)\longrightarrow0.$$

这展示了慢速改善的完整机制。它改变了目标数，不能作为 π 的构造。

把这个构造接到本项目已经认证的 11747 项 π 前缀后，两种目标均与 π 相差小于 10^(-11994)。在所有 q≤10^5000 的观测中，标准化误差的变化小于 10^(-1994)；不是精确误差值相等。更可证明：q≤10^2498 时，所有最近分子、标准化误差的严格排序，以及四个指定有理阈值的判断，全部相同。

因此有限观测不足以判定坏逼近性；它仍可帮助发现需要无限证明的结构。这一部分使用经典理论，没有宣称新数论定理。

- [完整反模型及有限证书转移](countermodels/finite_resolution_countermodel.md)
- [验证器](countermodels/verify_countermodel.py)与[结果](countermodels/verification.json)
- 依赖此前[π 的有限证书](../irrational-numerator-blast-2026-11-01/results_q_10pow5000.json)；本轮没有重复其完整认证。

## 6. 复现与来源

在本目录运行，均仅依赖 Python 标准库：

```bash
python3 endpoint/verify.py
python3 verify_discrete_assignment.py
python3 local_amplifier/verify.py
python3 countermodels/verify_countermodel.py --pi-certificate ../irrational-numerator-blast-2026-11-01/results_q_10pow5000.json
```

有限检查负责发现实现或代数录入错误；无限结论依据各篇正文证明。

所用原稿是用户上传的《The irrationality exponent of π is 2》，对应 [OpenAI 公开稿件](https://github.com/openai/math/tree/main/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026)。本轮核对所需估计，没有独立认证整篇证明。连分数基础事实参见 [Yann Bugeaud 的原始讲义](https://irma.math.unistra.fr/~bugeaud/travaux/Bamako2.pdf)。
