# 从六维球面构造中的周期系统到有理 π 构造：实际桥接试验

2026-11-02。这轮接续构造主线，完成了从几何工具到有理逼近的具体试验，不恢复旧素数池排除任务。

**联系已落实为可计算的公式链。** 从相关构造中的整数矩阵与椭圆曲线，提取经典二维模周期系统；通过明确的代数变量替换和已有周期恒等式，消去其他常数，得到纯有理系数的 \(1/\pi\) 级数，并实际核算既约分母、余项和一次 Padé 改造。

这条链恢复的是经典椭圆周期／Ramanujan 公式，不是新恒等式或新的好逼近纪录。它只使用原文可独立核算的矩阵与曲线，不依赖整个复流形被确认为六维球面的拓扑结论。

## 1. 从四维整数矩阵抽出二维周期系统

输入为 [Alpöge 构造稿](https://alpo.ge/s6.pdf) 开头的两个整数矩阵。计算找到原始不变旗标
\[
\mathbb Z\gamma\subset
W=\mathbb Z\gamma+\mathbb Zu+\mathbb Zw\subset\mathbb Z^4.
\]
二维子商 \(W/\mathbb Z\gamma\) 上的矩阵为
\[
A_1=\begin{pmatrix}-1&1\\-1&0\end{pmatrix},\quad
A_2=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\quad
A_0=\begin{pmatrix}1&-1\\0&1\end{pmatrix}.
\]
抽取这个子商不需要额外清除有限指数分母。尖点矩阵还满足
\[
N=T_0-I,\quad N^2=0,\quad
\operatorname{Smith}(N)=(1,1,0,0),\quad T_0^k=I+kN.
\]
它给出精确的对数分支差分，绕行系数仅线性增长。

这里得到的是整数局部系统的子商，不是一般复二维环面的椭圆商。还严格检查了：这个二维模块并不直接分裂嵌入原四维模块；直接用常整数列同时消去整个周期列的 \(\mu,\beta\)，只剩常数列。不能跳过解析归一化。详见 [monodromy/INTEGER_SUBQUOTIENT.md](monodromy/INTEGER_SUBQUOTIENT.md)。

## 2. 独立从原曲线导出周期方程

原稿 Remark 3.12 给出
\[
E_t:\ y^2=4x^3-3t^3(t-1)x-t^4(t-1)^2,\qquad j(E_t)=1728t.
\]
对 \(\oint dx/y,\oint x\,dx/y\)，我们直接作了符号微分约化。再换规范
\[
F(t)=t^{3/4}(t-1)^{1/4}\oint dx/y,\qquad z=1/t,
\]
并选解析尖点解的常数归一化，得到
\[
\boxed{z(1-z)f''+(1-\tfrac32z)f'-\tfrac5{144}f=0,\quad
f(z)={}_2F_1(1/12,5/12;1;z),\quad f(0)=1.}
\]
因此 \(f\) 是规范化周期，不是未经缩放的几何积分。独立符号审计见 [audit/verify_gauss_manin.py](audit/verify_gauss_manin.py)。

其系数可以直接归一成整数：
\[
u_n=1728^n[z^n]f\in\mathbb Z,\qquad
(n+1)^2u_{n+1}=12(12n+1)(12n+5)u_n.
\]
这把“存在整数单值矩阵”推进为实际可计费的周期递推。

对数伴随解的归一化也已算清：
\[
g=f\log z+h,\qquad
2\pi i\tau=g/f-\ln1728,\qquad
fg'-f'g=\frac1{z\sqrt{1-z}}.
\]
连接常数 \(\ln1728\) 不可省略。正向 \(z\) 回路给 \(\tau+1\)，原稿所命名的尖点生成元取反向，给 \(\tau-1\)，二者互逆。详见 [periods/MODULAR_SUBSYSTEM.md](periods/MODULAR_SUBSYSTEM.md)。

## 3. 真正接到 π 的显式变换

定义
\[
G(z)=f(z)^2={}_3F_2(1/6,1/2,5/6;1,1;z),
\]
\[
H(v)={}_3F_2(1/2,1/2,1/2;1,1;v)
=\sum_{n\ge0}\binom{2n}{n}^3(v/64)^n.
\]
实际桥梁是
\[
\boxed{G\!\left(\frac{27v^2}{(4-v)^3}\right)
=\sqrt{1-v/4}\,H(v).}
\]
我们核对了完整的有理微分算子变换，比例因子精确为 \(1/(v+8)\)，另核对了展开系数；不是小数拟合。代码见 [audit/verify_pullback.py](audit/verify_pullback.py)。

在 \(v_0=1/64\)，有 \(j=255^3\)，对应 \(\tau=i\sqrt7\) 的特殊值。利用经典模变换和 Legendre 关系，
\[
\boxed{5H(v_0)+42v_0H'(v_0)
=\sum_{n\ge0}\frac{(42n+5)\binom{2n}{n}^3}{4096^n}
=\frac{16}{\pi}.}
\]
其余周期值和额外对数常数在这里被消去。此式是已知的 Ramanujan 公式，证明来源为 [Guillera 公式 (2)](https://arxiv.org/pdf/1807.07394)。

在原 \(G\) 坐标，同一个特殊值是
\[
(432+7182\theta_z)G((4/85)^3)=85\sqrt{255}/\pi.
\]
换到 \(H\) 坐标后，代数平方根被变换吸收，最终级数纯有理。这说明算术坐标的选择具有实际作用。

## 4. 生成有理数并把分母账目算到底

令
\[
S_N=\sum_{n=0}^N\frac{(42n+5)\binom{2n}{n}^3}{4096^n},
\qquad C_N=4096^NS_N,\qquad r_N=16/S_N.
\]
有整数递推
\[
C_0=5,\qquad C_N=4096C_{N-1}+(42N+5)\binom{2N}{N}^3.
\]
相应的整数线性式已经明确：
\[
L_N=C_N\pi-16\cdot4096^N=-\pi\,4096^N T_N.
\]
因此这里确实只留下整数与 π，未夹带其他周期常数。
完全约分公式为
\[
\boxed{\gcd(16\cdot4096^N,C_N)=2^{3s_2(N)},\qquad
q_N=C_N/2^{3s_2(N)}.}
\]
这里 \(s_2(N)\) 是二进制数位和。证明用中央二项式二进赋值及最后一项的唯一最小赋值，无需寻找未知大 gcd。

正项相邻比小于 \(1/32\)，从而
\[
a_{N+1}<T_N:=16/\pi-S_N<\frac{32}{31}a_{N+1},\qquad
r_N-\pi=\frac{16T_N}{S_N(S_N+T_N)}.
\]
因此误差证书全部是有理数，不输入 π 的小数。

## 5. 实际尝试了分母压缩，而不只截断

对整数系数级数
\[
\mathcal H(t)=\sum_{n\ge0}(42n+5)\binom{2n}{n}^3t^n
\]
先做对角 Padé，再于 \(t=1/4096\) 求值和取倒数；比较时固定源级数系数数量。

同样只使用前三个系数：

| 方法 | 完全约分候选 | 既约分母 | \(q^2|\pi-p/q|\) |
|---|---:|---:|---:|
| 直接截断 | \(33554432/10680675\) | 10680675 | 约 \(1.09\times10^9\) |
| \([1/1]\) 周期 Padé | \(3041744/968217\) | 968217 | 约 \(567103.34\) |

改造使既约分母缩小约 **11.03 倍**，质量严格改善 **1892 倍以上**。它依然是差逼近，不能把相对改善称为爆破解。

完成的对照覆盖截断 \(N=0,\ldots,50\) 和 Padé 阶数 \(1,\ldots,10\)。更高的普通 Hankel 方程产生的系数成本，在这次对照中很快压过收益，因此不继续扩大这种无结构扫描。全部数据见 [rational-controls/ramanujan_checks.json](rational-controls/ramanujan_checks.json)。

## 6. 给构造主线留下什么

这次试验已经从相关几何构造走到**整数递推、纯有理 π 构造、完全约分公式和误差证书**。接通的是经典椭圆模周期机制；目前没有从六维球面的整体拓扑中额外获得算术收益。

值得进一步试的是保留二阶周期方程、周期与导数配对以及特殊值关系的分母设计。直接把它们压成高阶泛型 Padé 方程，会丢掉可用于控制整数系数的结构。本记录给出这条支线的精确输入和公平基准，不要求先解决另一个几何问题，也不把差基准的最终排除升级为主攻任务。

并行的整数原函数／规定相消主线继续保留，见 [规定相消研究](../pi-prescribed-cancellation-2026-11-02/README.md) 和 [当前研究方向](../CURRENT_DIRECTION.md)。本记录不覆盖它们。

## 核验与复现

- 原始整数旗标、Smith 不变量、不分裂与周期列消元均已精确检查。
- 独立从曲线作微分约化和代数换规范，完整导出周期方程。
- 40 阶 Frobenius、整数性、Wronskian／Clausen 核查，30 阶变换核查，独立 12 阶模形式展开。
- 完整有理微分算子变换，另有独立 0–17 阶级数比较。
- 51 个截断及 10 个 Padé，所有约分和误差使用精确算术。

运行：

~~~bash
python monodromy/verify_monodromy.py
python audit/verify_gauss_manin.py
python audit/verify_pullback.py
python periods/verify_modular_subsystem.py
python rational-controls/ramanujan_exact.py
~~~

SymPy 用于符号和有理线性代数，Fraction 用于误差区间。通用公式由证明支撑，有限检查不替代全参数论证。

来源：[Alpöge 原稿](https://alpo.ge/s6.pdf)开头、§3.1、Remark 3.12；[Engel 的后续论证](https://arxiv.org/html/2609.38442v1)；[Zagier 的周期与微分方程讲义](https://people.mpim-bonn.mpg.de/zagier/files/tex/ECMHirzebruchLecture/ECMHirzebruchLectureCorrected.pdf)；[Ramanujan 原文](https://ramanujan.sirinudi.org/Volumes/published/ram06.pdf)公式 (29)、(34)；[Guillera 的证明](https://arxiv.org/pdf/1807.07394)公式 (2)。
