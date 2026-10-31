# 自身进制、解析高斯和与 P2848：峰高估计及其限制

2026-10-31。lain 的研究方向，Mirai 协助推导与复现。

**状态：没有证明或证伪 P2848。** 本文给出一个解析模型的显式峰高误差界、峰形极限、连分数参数关系和有限精确证书。它们是本次研究得到的部分结果，不主张文献新颖性；没有运行 Lean 或独立形式化审查。

这里沿用用户澄清后的方向：数值不变，以无理数自身为尺度研究它的有理逼近。没有把十进制数位搬到另一个底数后重新赋值。本文的高斯和使用底数的幂尺度，尚未利用非整数进制的数字可容许性等特有结构。

## 1. 目标与解析对象

仓库中的 P2848 指

\[
\liminf_{n\to\infty}n|\sin n|=0,
\]

等价于 \(\liminf_{q\to\infty}q\|q\pi\|=0\)，也等价于 π 的普通连分数部分商无界。本文的“不曾解决”只描述本次工作的状态，不断言外部最新研究状态。

固定无理数 \(\beta>1\)，定义

\[
F_{p,q}(t)=\frac{\beta^{2t}(p-\beta q)^2+
                       \beta^{-2t}(p+\beta q)^2}{2\beta},
\qquad
\Theta_\beta(t)=\sum_{(p,q)\in\mathbb Z^2}e^{-F_{p,q}(t)}.
\tag{1}
\]

和式包含原点。设

\[
L_\beta(t)=
\left\{\frac1{\sqrt{2\beta}}
\begin{pmatrix}\beta^t(p-\beta q)\\ \beta^{-t}(p+\beta q)\end{pmatrix}:
(p,q)\in\mathbb Z^2\right\}.
\]

对应矩阵的行列式为 1，因此这是面积为 1 的格。\(F\) 是格向量的长度平方。

对任意紧实数区间，\(F_{p,q}(t)\ge c(p^2+q^2)\)，其中 \(c>0\) 只依赖该区间和 β，所以和式及其各阶导数一致收敛。更强地，(1) 在

\[
|\operatorname{Im}t|<\frac{\pi}{4\log\beta}
\tag{2}
\]

内全纯：写 \(t=x+iy\)，有
\(\operatorname{Re}F_{p,q}(t)=\cos(2y\log\beta)F_{p,q}(x)>0\)，并在任意紧子集上一致支配。因此 \(\Theta\) 是实解析函数。

对非零 \((p,q)\)，无理性保证 \(p\pm\beta q\ne0\)。直接极小化得

\[
t_{p,q}=\frac{\log\left|\frac{p+\beta q}{p-\beta q}\right|}{2\log\beta},
\qquad
m_{p,q}=\min_tF_{p,q}(t)=\frac{|p^2-\beta^2q^2|}{\beta},
\tag{3}
\]

以及精确公式

\[
F_{p,q}(t_{p,q}+h)=m_{p,q}\cosh(2h\log\beta).
\tag{4}
\]

当 \(p,q>0\) 时 \(t_{p,q}>0\)。当 \(p/q\to\beta\) 时，

\[
m_{p,q}=q^2\left|\beta-\frac pq\right|
                  \frac{\beta+p/q}{\beta}
\sim2q^2\left|\beta-\frac pq\right|.
\tag{5}
\]

## 2. 本次的定量结果：整个和式的峰由一个方向控制

**命题。** 取任意本原整数向量 \((p,q)\)，即 \(\gcd(p,q)=1\)。令 \(f=F_{p,q}(t)>0\)，则

\[
0\le\Theta_\beta(t)-\sqrt{\frac\pi f}\le B(f),
\tag{6}
\]

其中

\[
B(f)=
\sqrt{\frac\pi f}\frac{2e^{-\pi^2/f}}{1-e^{-3\pi^2/f}}
+2\left(1+\sqrt{\frac\pi f}\right)
\frac{e^{-1/f}}{1-e^{-3/f}}.
\tag{7}
\]

注意 (6) 是整个二维无穷和的上、下界，而不仅是单项或数值截断。它对所有 \(f>0\) 成立，在 \(f\to0\) 时尤其有效。

**证明。** 本原性允许把该向量补成行列式为 1 的格基。旋转坐标后，基可写成

\[
v=(\sqrt f,0),\qquad w=(u\sqrt f,1/\sqrt f),
\]

其中 \(u\in\mathbb R\)。故

\[
\Theta_\beta(t)=\sum_{b\in\mathbb Z}e^{-b^2/f}
                       \sum_{a\in\mathbb Z}e^{-f(a+bu)^2}.
\tag{8}
\]

一维高斯的 Poisson 求和式为

\[
\sum_{a\in\mathbb Z}e^{-f(a+s)^2}
=\sqrt{\frac\pi f}\sum_{k\in\mathbb Z}
                 e^{-\pi^2k^2/f}e^{2\pi i ks}.
\tag{9}
\]

由余弦不超过 1，移位和式不超过未移位和式
\(\vartheta(f)=\sum_a e^{-fa^2}\)。由单调积分比较，
\(\vartheta(f)\le1+\sqrt{\pi/f}\)。另一方面，(9) 在 \(s=0\) 给出
\(\vartheta(f)\ge\sqrt{\pi/f}\)。

在 (8) 中，\(b=0\) 的行就是 \(\vartheta(f)\)，其余行均为正。因 \(j^2\ge1+3(j-1)\) 对 \(j\ge1\) 成立，

\[
\sum_{j\ge1}e^{-c j^2}\le\frac{e^{-c}}{1-e^{-3c}}\quad(c>0).
\]

分别用于 \(b=0\) 行的 Poisson 余项和 \(b\ne0\) 的行，正好得到 (6)—(7)。证毕。

**几何含义。** 当一个本原方向很短时，与它线性无关的方向必须具有至少 \(1/\sqrt f\) 的垂直分量。因此其贡献受到 \(e^{-1/f}\) 的压制。高峰的主要来源是同一个短方向的整数倍。

## 3. 峰形、局部平均与连分数的解析编码

把 (4) 代入 (6)，令 \(m=m_{p,q}\)、\(t_0=t_{p,q}\)，得

\[
\Theta_\beta(t_0+h)=
\sqrt{\frac{\pi}{m\cosh(2h\log\beta)}}+R(h),
\quad 0\le R(h)\le B\bigl(m\cosh(2h\log\beta)\bigr).
\tag{10}
\]

沿任意 \(m\to0\) 的本原向量序列，对有界的 \(h\) 一致成立

\[
\frac{\Theta_\beta(t_0+h)}{\Theta_\beta(t_0)}
\longrightarrow\frac1{\sqrt{\cosh(2h\log\beta)}}.
\tag{11}
\]

\(t_0\) 是这个向量的平衡时刻，本文不声称它精确等于整个和式的局部最大值时刻。对于固定 \(H>0\)，(10) 还给出

\[
\int_{t_0-H}^{t_0+H}\Theta_\beta(t)^2\,dt
\sim\frac{\pi}{m\log\beta}
           \arctan\!\bigl(\sinh(2H\log\beta)\bigr).
\tag{12}
\]

所以足够高的峰在这个时间坐标下不会任意变窄；固定宽度的局部平方平均也能检测它。

再设 \(p_n/q_n\) 是 β 的普通连分数收敛分数，\(n\ge0\)，

\[
\alpha_{n+1}=[a_{n+1};a_{n+2},\ldots],\quad
r_n=\frac{q_{n-1}}{q_n},\quad \delta_n=\frac{p_n}{q_n}-\beta,
\]

约定 \(q_{-1}=0\)。由收敛分数行列式恒等式，

\[
q_n^2|\delta_n|=\frac1{\alpha_{n+1}+r_n},\qquad
\frac1{m_n}=\frac{\alpha_{n+1}+r_n}{2+\delta_n/\beta}.
\tag{13}
\]

结合 (6)，沿 \(a_{n+1}\to\infty\) 的子列得到

\[
\Theta_\beta(t_n)^2
=\frac\pi2\left(a_{n+1}+r_n+\frac1{\alpha_{n+2}}\right)+o(1).
\tag{14}
\]

这里误差是加性 \(o(1)\)：把 \(2+\delta_n/\beta\) 换成 2 产生的误差为 \(O(q_n^{-2})\)，高斯余项产生的误差指数衰减。公式把较大的部分商写成了解析函数的峰高，但尚不能在不知连分数尾部的情况下预测这些峰。

## 4. 与 P2848 的等价性，以及证明两侧的具体目标

对任意无理数 \(\beta>1\)，

\[
\liminf_{q\to\infty}q\|q\beta\|=0
\quad\Longleftrightarrow\quad
\sup_{t\ge0}\Theta_\beta(t)=\infty.
\tag{15}
\]

**证明正向。** 取 \(q^2|\beta-p/q|\to0\)，约分后仍可取这样的序列；分母不能有界。由 (5)，\(m\to0\)。相应 \(t_0>0\)，再由 (6) 的下界，\(\Theta_\beta(t_0)\to\infty\)。

**证明反向。** 记格最短非零长度为 \(\lambda(t)\)。若它在半轴上一致不小于 \(d>0\)，半径 \(d/2\) 的互不相交圆盘给出计数界
\(N_t(R)\le(1+2R/d)^2\)。积分得

\[
\Theta_\beta(t)=\int_0^\infty 2R e^{-R^2}N_t(R)\,dR
\le1+\frac{2\sqrt\pi}{d}+\frac4{d^2}.
\tag{16}
\]

因此高斯和无界必迫使 \(\lambda(t_j)^2\to0\)。有限的时间区间不可能如此，所以可取 \(t_j\to\infty\)。对对应向量有
\(|p_j-\beta q_j|\le\sqrt{2\beta}\lambda(t_j)\to0\)。非零整数性和 β 的无理性迫使 \(|q_j|\to\infty\)；改变整体符号后取 \(q_j>0\)，于是 \(p_j/q_j\to\beta\)。再用 \(m_{p_j,q_j}\le F_{p_j,q_j}(t_j)\to0\) 和 (5) 即得结论。证毕。

两个明确的解析目标是：

- **证明 P2848：** 证明 \(\Theta_\pi\) 在正半轴无界。例如证明 \(\limsup_{T\to\infty}T^{-1}\int_0^T\Theta_\pi(t)^2dt=\infty\) 就足够。后者是充分条件，本文没有证明它，也不声称它与 (15) 等价。
- **证伪 P2848：** 证明存在全局有限上界 \(M\)。此时对所有正整数 \(p,q\)，在平衡时刻有 \(m_{p,q}\ge\pi/M^2\)（非本原向量可先约分）。若 \(|\beta-p/q|<1\)，(5) 给出 \(q^2|\beta-p/q|\ge\beta\pi/((2\beta+1)M^2)\)；其余分数有平凡下界。因此得到有效的坏逼近常数 \(\min\{1,\beta\pi/((2\beta+1)M^2)\}\)。

有限计算只能检验有限个峰；对很多或几乎所有 β 的统计定理，也不能直接用于指定的 β=π。

## 5. 根号数对照与自身底数的限制

令 \(\beta=\sqrt D\)，其中 D 是正非平方整数。若 \(a^2-Db^2=\pm1\)、\(\varepsilon=a+b\sqrt D>1\)，变换

\[
(p,q)\longmapsto(p',q')=(ap+Dbq,bp+aq)
\]

是整数格的双射，且
\(p'+\beta q'=\varepsilon(p+\beta q)\)、
\(p'-\beta q'=(a-b\beta)(p-\beta q)\)，其中
\(|a-b\beta|=\varepsilon^{-1}\)。所以

\[
F_{p',q'}(t+T)=F_{p,q}(t),\qquad
\Theta_\beta(t+T)=\Theta_\beta(t),\qquad
T=\log_\beta\varepsilon.
\tag{17}
\]

| 自身底数 | 采用的单位 \(\varepsilon\) | 一个周期 |
|---|---|---|
| \(\sqrt2\) | \(1+\sqrt2\) | \(\log_{\sqrt2}(1+\sqrt2)\) |
| \(\sqrt3\) | \(2+\sqrt3\) | \(\log_{\sqrt3}(2+\sqrt3)\) |

连续性和周期性给出有界性；\(|p^2-Dq^2|\ge1\) 也给出 \(m\ge1/\sqrt D\)。没有为 π 得到对应的整数范数下界或无穷次更深返回的机制。

**底数的作用需要精确限定。** 令 \(s=t\log\beta\)，(1) 就变成含 \(e^{\pm2s}\) 的同一条格流。若保留目标数 β，只把缩放因子改为另一底数 \(b^{\pm t}\)，得到的函数是 \(\Theta_\beta(t\log b/\log\beta)\)。因此在当前模型中，“自身进制”只是自然的时间参数；解析性或选底数本身并未消除丢番图困难。这个观察仅针对当前高斯和模型，不排除真正利用非整数进制数字结构的其他模型。

**一个计算陷阱。** 固定整数盒子 \(|p|,|q|\le N\) 后，所有非零项都在 \(t\to\infty\) 时趋于零，因为 \(p-\beta q\ne0\)。因此任何固定盒子的截断高斯和都趋于 1，即使完整高斯和实际无界。长时间计算必须移动窗口或做格基约化，不能把这种截断曲线误当成有界性证据。

## 6. 实际验证与有限证书

[verify.py](verify.py) 包含标准库的精确检查，以及可选的 mpmath 数值检查。[results.json](results.json) 是本次实际运行结果。

```bash
python3 verify.py > results.json
```

标准库部分使用 Machin 恒等式 \(\pi=16\arctan(1/5)-4\arctan(1/239)\) 和交错级数余项，得到有理上下界。固定 \(p,q\) 时，\(p^2/x-xq^2\) 严格递减，因此可以用整数和 Fraction 严格包围 \(m=|p^2-\pi^2q^2|/\pi\)。程序验证：

- \(355/113\) 的 \(m<\pi/21^2\)，故由 (6)，\(\Theta_\pi(t_{355,113})>21\)。
- 仓库原有第 430 个收敛分数见证满足 \(m<\pi/180^2\)，故 \(\Theta_\pi(t_{430})>180\)。这里只重新认证该有理数的误差，不重新认证其连分数索引，也没有发现新的逼近纪录。
- 根号 2、根号 3 的两个 Pell 变换，共 7,442 组整数范数恒等式检查通过；普遍有效性由第 5 节代数证明给出。

原有第 430 个见证的 \(m\) 被严格包围在

\[
0.000096261122492578809256234158150580588553989
<m<
0.000096261122492578809256234158150580588553990.
\]

| 分数 | \(m\)，近似 | \(\Theta_\pi(t_0)\)，数值 | (7) 的余项上界，数值 |
|---|---:|---:|---:|
| \(22/7\) | 0.12394488705481634 | 5.03770508079964536 | 0.00378218398688 |
| \(355/113\) | 0.006812624149518618 | 21.4742359919863348 | \(8.02098956557\times10^{-63}\) |

表中小数使用 180 位精度的 mpmath 计算，**不是区间证书**；严格的数值结论是上述 \(>21\)、\(>180\)。程序还在两种根号底数及 π 的四个分数上、各三个偏移时刻检查 (4)、(6)，对这四个分数检查 (13)，并在两种根号底数、各三个时刻检查周期性。数值计算采用 (8) 的格基和 (9) 的 Poisson 形式，避免固定原坐标盒子的陷阱。数值检查不代替本文件的证明。

## 7. 下一步研究应集中在哪里

本次得到的是可检验的定量桥梁：好逼近产生什么高度、什么宽度的解析峰，峰是否可能由其他方向大量叠加造成，以及两侧证明分别需要怎样的全局估计。

若继续推进，应尝试从 π 的某个独立解析表示，直接控制第 4 节的长时间量；若仍先输入完整的 π 连分数尾部再计算峰，则主要是重新编码原问题。本文尚未找到从 π 的积分、级数或函数方程到该长期估计的有效推论。单纯的解析延拓、Poisson 对偶和坐标缩放在根号数对照中同样存在，不能独立承担这一步。

## 来源与范围

- 仓库此前的有限见证：[WAVE_CERTIFICATES_AND_CONTROL.md](../pi-approximation-2026-10-27/WAVE_CERTIFICATES_AND_CONTROL.md) 及该目录中的 `wave_certificate_results.json`。
- NIST DLMF，[§20.7(viii), Transformations of Lattice Parameter](https://dlmf.nist.gov/20.7#viii)：标准 theta 变换背景。本文直接写出需要的高斯 Poisson 公式及后续估计。
- Nicolas Chevallier，[*The natural extension of the Gauss map and the Hermite best approximations*](https://www.numdam.org/item/JTNB_2022__34_2_619_0/)，JTNB 34(2), 2022, 619–636：连分数、格和随参数变化的二次型之间已有成熟联系。本文不把格流或解析编码本身声称为新发现，也不把该文的几乎处处结论套用于 π。

本次未修改仓库原有 P2848 目标，也未声称上述新记录解决了该目标。


## 后续：长期平方平均

[不重叠区间分解与 π 正积分尝试](MEAN_SQUARE_DECOMPOSITION.md)：给出全局线性余项界、单峰充分条件和一次未成功的约分实验。


## 后续：稀疏高峰与平均路线的限制

[构造反例及完整证明](SPARSE_PEAK_OBSTRUCTION.md)：峰高无界但长期平方平均有界；给出不漏掉稀疏高峰的固定长度移动窗口判据。


## 后续：OpenAI 论文结论与证明方法的迁移

[无理指数、峰高增长率及端点尺度限制](OPENAI_SELF_BASE_TRANSFER.md)：证明一般无理指数与高斯峰指数增长率的关系；采用外部论文的 μ(π)=2 后，得到 π 高斯和及固定窗口平方积分的次指数增长。另检查 π 单位归一化中的算术成本，量化原参数比较对缓慢趋零误差的限制。这些结论不判定峰是否无界，也没有解决 P2848。


## 后续：剩余行列式的反例检验

[三个反例与实际中心的误差代价](RESIDUAL_DETERMINANT_COUNTEREXAMPLES.md)：确认二维约尽差值后的非零性；构造整系数指数小量及对数插值子矩阵中的阶乘小量，核算清分母与中心偏移；给出三维不能直接推广的反例。这些有限子矩阵不替代原稿的完整行渐近构造，未解决 P2848。
