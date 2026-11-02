# 构造优先：整数原函数坐标、不对称端点与有符号误差证书

2026-11-02。接续 `697aef18b10f969e7bd163c2abe1674db9353279` 后用户明确纠正的研究方向。**本轮不继续主攻旧 Legendre 族的素数池排除问题。** 当前方向另见 [CURRENT_DIRECTION.md](../CURRENT_DIRECTION.md)。

本轮得到了可证明的构造参数选择规则、同预算下的整数格成本改善，以及可直接验证候选质量的有理积分上界程序。**尚未构造出质量趋零的无限序列，也没有新的 π 有理逼近纪录。** 已知 \(355/113\) 等仅作为实现基准。

## 1. 先安排整数系数，再做积分

选整数 $c\ge3$、$\alpha\ne0$ 及奇整数多项式
\[
R(x)=x\sum_{j=0}^{c-1}r_jx^{2j},\qquad r_j\in\mathbb Z.
\]
定义整数分子
\[
N=(1+x^2)R'-2(c-1)xR+\alpha(1+x^2)^{c-1}.
\]
全参数恒等式为
\[
\boxed{I=\int_0^1\frac{N(x)}{(1+x^2)^c}\,dx
=\frac{R(1)}{2^{c-1}}+\frac{\alpha\pi}{4}.}
\]
有理边界分母在构造时已限制为二的幂。设
\[
g=\gcd(R(1),2^{c-3}\alpha),\quad
\frac pq=-\frac{R(1)}{2^{c-3}\alpha},\quad q>0,
\]
其中 (p/q) 完全约分，则
\[
\boxed{q=\frac{2^{c-3}|\alpha|}{g},\qquad
Q:=q^2|\pi-p/q|=\frac{4q^2}{|\alpha|}|I|.}
\]
这改变了构造坐标：自由选择整数原函数数据和端点条件，而非沿用旧 $A_n,B_n) 递推等待大 gcd。奇算术成本仍会影响 (\alpha$ 和所需分子大小；下面对它作精确计算。

## 2. 一个全参数成立的建设性设计规则

要求
\[
N(x)=x^{2a}(1-x^2)^bH(x^2),\quad H\in\mathbb Z[t],
\quad a,b\ge1,\quad c\ge a+b+1.
\]
这些是整数线性条件，允许 $\deg H\le c-a-b$。取其**完整饱和整数解格**，映射到整数对 $ (\alpha,R(1)) $，记像为 $\Lambda_{c,a,b}$。本轮证明对每个奇素数 $p$：
\[
\boxed{\Lambda_{c,a,b}\otimes\mathbb Z_p
=p^{\lfloor\log_p M\rfloor}\mathbb Z_p\times\mathbb Z_p,
\quad M=\max(2a-1,b).}
\]
因此周期像格协体积的奇数部分**恰好**为
\[
\boxed{\operatorname{odd}\operatorname{lcm}(1,\ldots,\max(2a-1,b)).}
\]
这既是必要限制，也证明了在奇素数处没有额外隐藏限制；不是只得到一个粗上界。

证明用两端 Taylor 整性。在零点端，条件等价于 $\alpha/(2j+1)$ 在奇素数处整，$0\le j<a$。在另一端，写
\(1/$2+2u+u^2)=\sum h_j u^j$，则需要 $\alpha h_{k-1}/k$ 为整；当 $k$ 为奇素数幂时 $h_{k-1}$ 是二进单位。这给出 $M$。充分性在 $t=x^2$ 坐标用 $t=0,1$ 的 Hermite/CRT 插值，所得次数至多 $a+b\le c-1$，没有暗中增加次数。完整证明见 [integer-kernel/INTEGER_ANTIDERIVATIVE_FAMILY.md](integer-kernel/INTEGER_ANTIDERIVATIVE_FAMILY.md)。

固定参数预算 $a+b=s$，应优先考察 **$b\approx2a$**，这样 $M$ 从 $a=b$ 时的约 $s$ 降到约 $2s/3$。这里 $a,b$ 是分子中 $x^{2a},(1-x^2)^b$ 的参数；它们对应的实际 $x$ 端点消失阶数是 $2a,b$。

## 3. 同阶数、同端点预算下，确实减少了多少成本

保持 $c=16$、$a+b=14$，精确 Smith/Hermite 计算给出
\[
\Lambda_{16,7,7}=45045\mathbb Z\times1024\mathbb Z,
\qquad
\Lambda_{16,5,9}=315\mathbb Z\times32768\mathbb Z.
\]

| 参数 | 必须承担的奇周期因子 | 边界周期步长 | 周期像格协体积 |
|---|---:|---:|---:|
| $a=b=7$ | 45045 | 1024 | 46126080 |
| $a=5,b=9$ | 315 | 32768 | 10321920 |

奇因子减少 **143 倍**，二进边界成本增加 **32 倍**，完整周期像格协体积净改善 **$143/32$**。这是一项公平预算下的算术设计收益，**不能直接当成 $Q$ 改善了 $143/32$ 倍**：端点形状改变后，积分与 $H$ 的系数成本也要重新计算。

本轮据此选择了具体的增长参数族
\[
\boxed{a=m,\quad b=2m,\quad c=4m+1,\quad\deg H\le m+1.}
\]
它保留随规模增长的自由乘子，而非始终只扫一次乘子。其强制奇周期因子是 $\operatorname{odd}\operatorname{lcm}(1,\ldots,2m)$，可以在选参数时直接计费。

这一预声明参数族已经实际完成 \(m=1,\ldots,6\) 的完整整数核计算，每个规模只用五组预先指定的成本权重，生成候选时不输入 π 或目标分数。下面是各规模在这些有限候选中的最好结果，均另经第 5 节的有理积分证书独立核验：

| \(m\) | \(c\) | 完全约分的候选 | \(Q\) 约值 |
|---:|---:|---:|---:|
| 1 | 5 | \(22/7\) | 0.06195997 |
| 2 | 9 | \(355/113\) | 0.00340631 |
| 3 | 13 | \(355/113\) | 0.00340631 |
| 4 | 17 | \(104348/33215\) | 0.36586381 |
| 5 | 21 | \(1459097/464445\) | 1.62020680 |
| 6 | 25 | \(85563208/27235615\) | 0.63806586 |

规模 6 已产生分母 \(27235615\) 且严格 \(Q<1\) 的候选。这些仍是有限基准，不是新的最佳逼近纪录；质量不随规模单调下降，也没有据此证明无限子序列成功。全部 30 个候选的整数数据见 [integer-kernel/new_family_results.json](integer-kernel/new_family_results.json)，独立积分认证见 [certificates/new_family_integral_checks.json](certificates/new_family_integral_checks.json)。

## 4. 以质量的上界选择构造

令 $h=\sum_j|H_j|$。对 $s=a+b$，直接最大值计算得到
\[
\max_{0\le y\le1}\frac{y^a(1-y)^b}{(1+y)^s}
=\frac{a^ab^b}{2^as^s}.
\]
因此当 $c\ge s$ 时有完全有理、可在计算 π 之前使用的**上界**：
\[
\boxed{Q\le\frac{4q^2h}{|\alpha|}
\frac{a^ab^b}{2^a(a+b)^{a+b}}.}
\]
对上面的增长族，它变成
\[
\boxed{Q_m\le\frac{4q_m^2h_m}{|\alpha_m|}
\left(\frac2{27}\right)^m.}
\]
这是建设性的充分目标：如果实际选出的整数解能满足
\(q_m^2h_m/|\alpha_m|=o$(27/2)^m)$，就能证明 $Q_m\to0$。**这项整数解高度控制尚未证明**，也不把它预设成真。减小强制奇因子只是为这个充分目标改善参数条件。

这个简易上界舍弃了有符号积分的相消，对好候选可能很松。为此本轮还实现了更精确的积分证书。

## 5. 不输入 π 数值的有符号误差证书

程序 [certificates/integral_certificate.py]$certificates/integral_certificate.py) 直接对整数 (N,c$ 验证积分。利用
\[
\frac1{(1+x^2)^c}
=\left(\frac23\right)^c
\sum_{k=0}^K(-1)^k\binom{c+k-1}{k}
\left(\frac{2x^2-1}{3}\right)^k+\text{余项},
\]
其中 $ |(2x^2-1)/3|\le1/3$。所有截断项乘 $N$ 后都可精确积分为有理数。对 $K\ge c$，余项绝对值一致地不超过
\[
\left(\frac23\right)^c
\frac{\binom{c+K}{K+1}}{3^{K+1}}
\left(1-\frac{c+K+1}{3(K+2)}\right)^{-1}.
\]
再乘可计算的 $\sum_j|N_j|/(2j+1)$，便得到 $I$ 的有理上下界，包括符号相消；最后通过质量恒等式给出 $Q$ 上界。

基准整数输入 $c=7,\alpha=-339$、
\[
(r_0,\ldots,r_6)=(339,1921,4473,5455,3560,1200,92)
\]
产生已知 $355/113$，程序在不输入 π 或其连分数的情况下证明
\[
0.00340631193013<Q<0.00340631193015.
\]
**这只是证书实现的已知基准。** 积分本身含有 π 的信息，“未输入 π 的小数”不意味着已经绕过丢番图难题。无限序列仍需上节的整数高度或更强相消估计。

独立验证另用精确原函数恒等式与 Machin 区间核对 $c=1,\ldots,20$ 的 80 组整数输入（含零周期情况），全部通过。证书误差可由参数 $K$ 显式缩小，没有把浮点误差当作逼近。

## 6. 做过哪些新路对照，为什么没有继续扩大它们

- 饱和整数核及自由乘子：完成 $4\le c\le12$ 的 698 个短向量检查。自由乘子相对同 $c$ 的常乘子控制有实质有限改善，但得到的是 $22/7,333/106,355/113$ 等已知基准；没有单调改善序列。
- 移动对称复路径：推导了能自动只保留 π 的整数端点参数族。精确剔除了与旧端点仅反射、缩放等价的情况。近临界 Pell 参数改善解析误差，却在这些试验中增加了更多分母成本；不再为它追求最终失败定理。
- 固定 $\mathbb Q(\sqrt2)$ 及耦合共轭零点：建立了明确的有理迹与范数费用公式，并分别检查两条连续对数分支。有限测试没有找到新的好逼近；作为下一次改变共轭权重时的基准保留。
- 端点耦合 type-I Machin 变体：已做有限精确对照，其中部分只是旧 Padé 表的对偶，没有算术收益。停止扩大该扫描。

这些对照的作用是及时选择下一种构造，不是把“证明旧方法不行”重新变成主线。详细推导、结果范围与代码分别保存在 `integer-kernel/`、`contour-controls/`、`type-I-control/`；简单二次域迹控制另在 `trace-control/`。

## 7. 接下来的明确研究对象

主线固定为第 3 节的增长族及其完整整数核：利用已证明的奇周期成本规则，选择随 $m$ 增长的 $H$，同时控制 $\alpha$、$R(1)$、$H$ 的高度与实际积分相消，用第 4、5 节的上界检验成功。普通的总导数解 $\alpha=0$ 不产生 π 逼近，必须剔除。

不能根据格协体积降低就声称质量下降，也不能根据有限搜索没有成功就宣称此族被排除。当前完成的是一次有具体数学对象、参数设计定理和可复现成功判据的改造；通往无限质量趋零结论的高度估计仍然开放。

## 来源与复现

双极点文献对照核对了 Zeilberger–Zudilin (2020) 的[作者原文](https://sites.math.rutgers.edu/~zeilberg/mamarim/mamarimPDF/pimeas.pdf)及[作者代码页](https://sites.math.rutgers.edu/~zeilberg/mamarim/mamarimhtml/pimeas.html)。该文用于理解如何改造积分与核查已知结构，不把其中的无理指数上界视为本项目的质量趋零结论。

本轮主线证明为整数原函数、Taylor 整性与 Hermite 插值的直接推导。运行以下命令可重现核心验证：

```bash
python integer-kernel/verify_family.py
python integer-kernel/verify_odd_period_image.py
python certificates/integral_certificate.py
python certificates/verify_certificate.py
python integer-kernel/verify_new_family.py
python certificates/certify_new_family.py
```

整数核搜索需要 SymPy；积分证书主程序只需 Python 标准库。SymPy LLL 大有理数向下取整的本地修正、完整饱和整数核与初始子格的差别，以及未完成长跑的范围，均在主证明和代码中标明。
