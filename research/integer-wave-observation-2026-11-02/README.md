# 整数波：从稳定观测到可认证的构造工具

2026-11-02。延续用户早先的“整数波、有限观测、逐层细化”，本轮重新指定对象、观测规则、误差预算和构造接口。研究方向仍是构造更好的有理逼近。

合并前核对仓库最新提交 `7ff9c3d57c67240ccd29007454eb27ba71f39dcb`：并行的 [周期商构造](../pi-period-quotient-design-2026-11-02/README.md) 已实现 Hilbert 范数、Schur 度量和强制公因子子格。本轮的 Gram／固定周期实验与该方向相接，不另称首创；此处增加的主要观察工具是有采样误差证书的拓扑读数、数据驱动的零点计数，以及有限谱观测审计。保存时保留该并行记录与最新方向文件。

**本轮交付是有证明和程序的工具原型，尚未得到 π 的质量趋零子序列，也没有证明物理世界的整数来自坍缩。** 下述拓扑、旋转数、谱理论都是已有数学；本轮工作是推导明确证书、实现并连接本项目，不主张基础定理的新颖性。

## 1. 从什么对象读出整数

旧的整除波核

\[
\Delta_d(x)=d^{-2}\left|\sum_{r=0}^{d-1}e^{2\pi irx/d}\right|^2
\]

能在整数输入上检测整除。本轮使用更一般的对象：一条连续、闭合、不经过原点的复波轨迹 \(\gamma:[0,1]\to\mathbb C\setminus\{0\}\)。沿着轨迹连续展开角度，记它为 \(\theta(t)\)，定义读数

\[
W(\gamma)=\frac{\theta(1)-\theta(0)}{2\pi}\in\mathbb Z.
\]

波形、振幅、局部速度可以连续改变；只要保持闭合且不穿过零，整数读数保持不变。改变读数的事件是越过零或破坏闭合。此处的整数是全局绕行数，不是把振幅四舍五入。

波形相乘使读数相加；取倒数使读数变号；将归一化波形视为圆到圆的映射后，映射复合使读数相乘。由此得到连续对象上的整数运算模型。

模型预先采用连续性、复平面和参考周期，不能用它宣称已从“没有数的世界”推导出数。绕数的标准背景见 Sjamaar §4.3；圆上动量谱的另一种整数读数见 Etingof Examples 8.8–8.9。

## 2. 有限观测什么时候真的足够

设波到零的距离至少为 \(m>0\)，变化速度上界为 \(L\)，每周期取 \(N\) 个等间隔样本，样本误差至多 \(\eta\)。将测得的样本连成闭合折线 \(P\)，则

\[
\|P-\gamma\|_\infty\le\eta+\frac{L}{2N}.
\]

因此

\[
\boxed{\eta+\frac{L}{2N}<m\quad\Longrightarrow\quad W(P)=W(\gamma).}
\]

证明：在每个采样间隔内，精确端点插值的误差至多 \(2Ls(1-s)/N\le L/(2N)\)；再加测量误差。连接真波与折线的直线同伦始终不经过零，绕数保持不变。完整证明及精确有理数实现见 [foundations/NOTE.md](foundations/NOTE.md)。

这是一条能使用的“有限数据→精确整数”证书。若缺少速度约束，\(e^{2\pi i kNt}\) 对所有采样点都读成 1，但绕数为 \(kN\)；采样本身不能替代假设。

**无需预知真实边界距离的版本也已实现。** 记 \(d(P,0)\) 为测量折线到原点的最短距离。只要

\[
d(P,0)^2>\left(\eta+\frac{L}{2N}\right)^2,
\]

便可直接从数据认证真波不经过零且绕数读数正确。折线距离平方由逐边投影用有理数精确计算。把波取为“多项式沿正方形边界的函数值”，再用辐角原理，就得到内部零点个数的证书。6 组无噪声／含噪声实例分别认证 2、0、3 个内部零点；边界零点实例被正确拒绝。这个工具已能接入辅助多项式的零点研究。验证失败只表示当前证书不足，不能据此断言存在边界零点。

## 3. 向有理数、实数、复数扩展

| 对象 | 本模型中的含义 | 必须保留的信息 |
|---|---|---|
| 整数 | 一个闭合观察过程的净绕行数 | 连续相位或经认证的离散样本 |
| 有理数 \(p/q\) | \(q\) 个参考周期内走过 \(p\) 圈 | 最小共同返回周期；需约分 |
| 实数 | 长时间平均相位速度的兼容极限 | 一致的有限观察区间，不能只留末端相位 |
| 复数 | 相位速度加对数振幅增长率 | 相位的展开历史和振幅；属于非单位演化 |

实数层采用递增的圆映射提升 \(F:\mathbb R\to\mathbb R\)，满足 \(F(x+1)=F(x)+1\)。其平移数 \(\rho\) 是迭代位移的长期平均。令 \(m_n=\lfloor F^n(x)-x\rfloor\)，则

\[
J_N=\bigcap_{1\le n\le N}\left[\frac{m_n}{n},\frac{m_n+1}{n}\right]
\]

包含 \(\rho\)，逐层嵌套，宽度至多 \(1/N\)。这提供不以连分数位为输入的一致精化规则。平移数理论及这个观察公式属于经典 Poincaré 理论，见 Milnor Chapter IV §14。

复数有两种不同扩展。对于 \(z=a+ib\)，\(e^{2\pi izt}=e^{-2\pi bt}e^{2\pi iat}\) 记录相位和增长；当 \(b\ne0\) 时没有非零整数时间的精确返回。若要把 \(\mathbb Q+i\mathbb Q\) 解释为有限返回，需要分别记录实部、虚部的两个单位相位。共同整数分母的二维实数逼近，其自然指数是 \(3/2\)；采用 Gaussian 整数分母的复逼近才对应另一种指数 2 问题。两者不得混用。

## 4. 将分母成本写入波的能量

对实目标 \(\alpha\)、整数 \(p\)、正整数 \(q\)、尺度 \(s>0\)，定义模式能量

\[
\mathcal E_{\alpha,s}(p,q)=\frac{q^2}{s}+s(p-\alpha q)^2.
\]

这是二维环面上一个椭圆自伴微分算子的 Fourier 模式能量。第一项记录分母成本，第二项记录失谐。对尺度取下确界，得到

\[
\inf_{s>0}\mathcal E_{\alpha,s}(p,q)=2q|p-\alpha q|
=2q^2\left|\alpha-\frac pq\right|.
\]

这部分是已有格流／Dani 对应的谱表示，本身没有给出新的逼近定理。它提供两个可核验接口：

1. **有限无漏项检测。** 能量不超过 \(L\) 的模式必有 \(q\le\sqrt{sL}\)、\(|p-\alpha q|\le\sqrt{L/s}\)。给定目标的有理包围区间，有限整数枚举即可对所有这些模式同时给出严格上下界。
2. **低谱唯一方向。** 同一尺度下，若两个模式的能量都严格小于 1，则它们必共线；完全约分后是同一个有理数。因为格变换的行列式为 1，整数行列式满足
   \[
   |p_1q_2-p_2q_1|\le\sqrt{\mathcal E_1\mathcal E_2}<1.
   \]
   其值只能为零。阈值必须严格：在能量恰为 1 时可以出现独立方向。

纯正交 Fourier 模式的能量是各模式能量的非负加权平均。仅调整这些模式的相位不会凭空降低最低能量。因此要设计有效相消，需要保留实际积分／周期中的交叉项和整数结构。

Machin 公式的有理区间验证：\(s=4{,}000{,}000\)、\(L=1/10\) 时，穷尽全部 \(q\le632\)，唯一既约低谱方向为 \(355/113\)。未约分时还出现它的二、三次谐波。168 组精确检查通过。这个已知分数只是检测器基准。

## 5. 回接构造：保留交叉项的整数波形选择

本项目已有整数原函数坐标给出

\[
I=\int_0^1\frac{N(x)}{(1+x^2)^c}\,dx
=\frac{R(1)}{2^{c-1}}+\frac{\alpha\pi}{4}.
\]

将 \(-R(1)/(2^{c-3}\alpha)\) 完全约分为 \(p/q\) 后，质量精确为 \(Q=4q^2|I|/|\alpha|\)。候选波形仍在满足端点条件的整数格中选择。

本轮将此前的绝对系数和改为精确 Gram 度量。若自由波形为 \(H(x^2)=\sum_i h_i x^{2i}\)，则

\[
\int_0^1 H(x^2)^2\,dx=h^{\mathsf T}Gh,
\qquad G_{ij}=\frac{1}{2i+2j+1}.
\]

所有交叉项均精确有理，可以进入候选发生器并给出严格误差上界；实际约分后的 \(q\) 同时计入评分。完整构造、证书比较、标准方法归属及局限见 [constructive/NOTE.md](constructive/NOTE.md)。这提供具体的波形改造空间；有限证书改善还不能推出质量趋零。

对旧候选保持同一核界时，\(m=5\) 的 Gram 证书比绝对系数和界收紧约 444.54 倍。这不改变该候选的真实误差；此前的直接有理积分证书仍可更紧。核的 \(16^{-m}\) 界已由并行的旧记录得到，本轮未将它当作新发现。

另一个可复用步骤是固定两个周期读数后，在其余波形自由度上做精确 Schur 补优化。最优波形为有理系数，整体清分母后成为合法整数波形；整体倍数会在最终分数的 gcd 中消掉，所以实际 \(q\) 不增，但构造系数的位数仍必须记录。这一步优化同一分数的波形证书，不能把它计成一个新的逼近。

## 6. 可复现文件与后续接口

- `foundations/`：连续波读数、兼容精化、噪声／采样定理与验证。
- `spectral_certificate.py`、`spectral_checks.json`：精确目标区间的有限共振检测。
- `spectral-audit/`：独立谱验证、低谱唯一方向、Rayleigh 接口和复数指数审计。
- `constructive/`：精确 Gram 波形度量、整数格发生器、同周期波形优化和结果。

在本目录中运行：

```sh
python3 foundations/verify_winding.py
python3 foundations/verify_rational_polygon.py
python3 foundations/verify_polynomial_contour.py
python3 spectral_certificate.py
python3 spectral-audit/verify.py
python3 constructive/gram_wave.py
python3 constructive/fiber_projection.py
python3 constructive/verify_gram_wave.py
```

前五项只用 Python 标准库；后三项使用 SymPy 1.14，并在 `constructive/support/` 附带所需旧模块。显示用浮点数不进入精确证书判定。绕数的浮点样例另外明确标为一致性检查，其精确版本也已运行。

后续可执行的构造问题：在整数周期格内，选择和变形波形，使**实际约分分母的平方乘以经认证的剩余积分**下降；使用绕数证书控制变形中的零点／相位事件。仍需找到一个随规模增长的统一估计，证明这种下降能无限延续。仅指定 \(F(x)=x+\pi\) 或在能量中写入 \(\pi\)，并不提供该估计。

主方向保持构造优先；没有继续推进旧固定族的素数池排除。

## 来源与归属

- Sjamaar, *Manifolds and Differential Forms*, §4.3: https://pi.math.cornell.edu/~sjamaar/manifolds/manifold.pdf
- Milnor, *Dynamics Notes*, Chapter IV §14: https://www.math.stonybrook.edu/~jack/DYNOTES/dn14.pdf
- Etingof, *Mathematical Ideas and Notions of Quantum Field Theory*, Examples 8.8–8.10: https://ocw.mit.edu/courses/18-238-geometry-and-quantum-field-theory-spring-2023/mit18_238_s23_lec_full.pdf
- Gorodnik, *Lectures on Diophantine Approximation and Dynamics*, Lecture I: https://www.math.uzh.ch/gorodnik/tifr/index.html
- Sutherland, *18.783 Elliptic Curves*, Lecture 14 §14.2：辐角原理与零点计数。https://math.mit.edu/classes/18.783/2023/LectureNotes14.pdf
- 本仓库 `research/pi-constructive-redesign-2026-11-02/`、`research/pi-prescribed-cancellation-2026-11-02/`：整数原函数格、既有核界与有理积分证书。

上述来源用于说明已有结构的归属。本轮具体误差预算、代码验证、Gram 发生器与整合过程详见各子目录；没有文献新颖性或形式化证明器验收的声明。
