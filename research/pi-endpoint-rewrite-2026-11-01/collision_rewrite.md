# 保留碰撞收益与横向误差收益：对 π017 的一次具体改写

这里的障碍来自矩阵的规模、清分母成本和两次估计之间的损失，反证法本身没有造成障碍。本文把这些量放回同一个非零子式中分析。

研究对象固定为 `openai/math` 的快照 `4eaee87548250b2682affda7d05d947c05554bae`，目录 `preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026`。所用实际源码行号如下：

- `build/sections/determinant.tex:45–81`：全列空间、全部行、同一个满行非零子式和 \(\bar b_H\)。
- `determinant.tex:99–149`：清分母及算术下界。
- `determinant.tex:196–257`：逐行精确翻译。标量与选中哪一列无关，所以可以对同一个子式用行列式多线性。
- `determinant.tex:267–335`：同一横向指标的 Taylor 次数碰撞。
- `determinant.tex:345–431`：把碰撞与逼近误差分成两种情况，再分别丢掉一种收益。
- `build/sections/conclusion.tex:45–95`：\(A,B,C,\theta,m\) 的选择及指数规模关系。

[原始 determinant.tex](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/determinant.tex) · [原始 conclusion.tex](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/conclusion.tex)

## 1. 两次具体损失

全部行是

\[
\rho=(j,s,\beta),\quad 0\le j<K,\quad
v_0s+\frac{w\beta}{\theta}<H,
\qquad w\beta=\sum_{i=1}^m w_i\beta_i.
\]

把行数记为 \(M\)，并定义

\[
\bar b_H=\frac{\sum_\rho w\beta_\rho}{MH}.
\]

论文只在最终比较中使用 \(\bar b_H\le\theta\)。但对这个**完整的行单纯形**，实际上有精确极限

\[
\boxed{\bar b_H\longrightarrow b=\frac{m\theta}{m+2}}.
\]

证明：令 \(x_0=v_0s/H\)、\(x_i=w_i\beta_i/(\theta H)\)。行的渐近分布是 \(m+1\) 维单纯形 \(x_i\ge0,\sum x_i<1\) 上的均匀分布，各坐标平均为 \(1/(m+2)\)。而 \(w\beta/H=\theta\sum_{i=1}^m x_i\)。严格边界不影响 Riemann 和极限。这与论文原有计数使用同一个固定参数、\(H\to\infty\) 顺序。

因此，指数 \(\nu=2\) 时，大横向分支的理想间隙应由

\[
2A(1-\eta)-1-\theta
\]

改成

\[
\boxed{g_{\rm mean}=2A(1-\eta)-1-\frac{m\theta}{m+2}}.
\]

这是实质改善，有限维下多出 \(2\theta/(m+2)\)。但它本身尚不能完成比较，见下一节。

第二次损失发生在 `determinant.tex:394–406`：小横向分支把所有逼近误差收益丢掉，大横向分支把所有碰撞收益丢掉。下面保留两者之和。

## 2. 一个严格的障碍：均值改进仍救不了原来的两分支参数方案

这里证明的是**原来的 \(A,B,C\) 指数参数方案**的障碍，不是所有可能构造的障碍。假定

\[
0<\theta<A<B<1,\qquad CB<1,\qquad B<C\theta,
\]

以及论文选择

\[
K=\lfloor C^m\rfloor,\quad w_0=B^{-m},\quad
v_0=2K\theta^m w_0.
\]

碰撞分支中的量仍是

\[
L=\frac{\eta^2(B/A)^m}{2(m+1)},\qquad
c=\frac{\ln2}{4}.
\]

即使把所有误差设为零，并把两条分支的算术成本都改成精确的 \(1-b\)，也不可能同时有

\[
g_{\rm mean}>0,\qquad cL>1-b.
\]

证明如下。\(CB<1\) 与 \(B<C\theta\) 给出 \(B^2<\theta\)。而 \(g_{\rm mean}>0\) 给出

\[
A>\frac{1+b}{2(1-\eta)}.
\]

于是

\[
\frac BA<\frac{2(1-\eta)\sqrt\theta}{1+m\theta/(m+2)}
\le(1-\eta)\frac{m+2}{m+1}.
\]

最后一步成立，因为函数 \(2\sqrt\theta/(1+m\theta/(m+2))\) 在 \(0<\theta<1\) 递增，端点值为 \((m+2)/(m+1)\)。所以

\[
L<\frac{e\eta^2(1-\eta)^m}{2(m+1)}.
\]

对 \(0<\eta<1\)，微分给出

\[
\max\eta^2(1-\eta)^m=\frac{4m^m}{(m+2)^{m+2}}
<\frac4{(m+2)^2}.
\]

同时 \(1-b>2/(m+2)\)、\(e<3\)、\(c<1/4\)。合并即得

\[
\boxed{\frac{cL}{1-b}<\frac{3}{4(m+1)(m+2)}\le\frac18.}
\]

碰撞分支甚至只能获得所需理想成本的八分之一以下。这排除了“只把 \(\bar b\) 算精确，再微调有限维参数”这条修补方式。

另有更粗但直观的统一界：\(\eta^2(1-\eta)^m\le4/27\)，故 \(L<1/9\)、\(cL<1/36\)。不能只用这个统一界和趋零的 \(1-b\) 比较；上面的含 \(m\) 比例界才完成加强后的障碍证明。

这个结论并不排除新的 Cauchy 半径、更强的跨组依赖分析、不同参数方案、选列成本优化或完全不同的插值形状。若更改碰撞常数或参数关系，应重新检查证明中的对应不等式。

## 3. 不再切成两种情况：同一子式的联合预算

对每条原始行 \(\rho\)，在其精确展开中选择一个指标 \(a_\rho\)。所有真正出现的非零项必满足

\[
a_\rho\ge\beta_\rho\quad\text{逐坐标},\qquad wa_\rho\le H.
\]

令 \(n_a=\#\{\rho:a_\rho=a\}\)。在指数 \(\nu=2\) 时，保留论文已有的两种收益得到

\[
\boxed{\mathcal S_H(a_\rho)
=c\sum_a n_a^2+2\sum_\rho w(a_\rho-\beta_\rho).}
\]

第一项来自同一 \((a,d)\) Taylor 系数行的碰撞；第二项来自同一行展开的逼近误差标量。两者可以同时出现。它们不是两个不同子式的估计。

若仅假设 \(|\pi-p_i/q_i|\le q_i^{-2}\)，论文逐行估计及碰撞估计仍适用，只需把其中的 \(\nu\) 换成 \(2\)。若另保留 \(\kappa_i=q_i^2|\pi-p_i/q_i|\)，还有额外收益，但这里没有把它擅自当成固定正指数。

定义有限松弛问题

\[
S_H^{\rm rel}=\min\{\mathcal S_H(a_\rho):
a_\rho\ge\beta_\rho,\ wa_\rho\le H\}.
\]

它允许了某些实际上标量为零、Taylor 行列式为零的选择，因此真正非零项上的最小收益不会小于它。这是一个可以用于**充分下界**的保守松弛。

如果同一个选中非零子式有有效算术下界

\[
\ln|\Delta_H|\ge-C_H^{\rm ar},
\]

而展开的翻译、整函数、求和与有限 \(H\) 误差总量为 \(E_H^{\rm total}\)，那么

\[
\boxed{S_H^{\rm rel}>C_H^{\rm ar}+E_H^{\rm total}}
\]

才给出矛盾。新的整行整列清分母方式、实际选中列的平均权重，都可以通过 \(C_H^{\rm ar}\) 接入。不能把一种子式的低分母成本和另一种子式的高碰撞收益拼接。

这个准则也没有把 \(S_H^{\rm rel}\) 不够大解释成实际 \(\Delta_H\) 不小：松弛问题可能被不可能出现的项拖低。那时需要进一步利用真实项的零结构。

## 4. 这恰好是一个凸最小费用流问题

把具有相同 \(\beta\) 的原始行合并，供给数为

\[
m_\beta=K\#\{s\ge0:v_0s+w\beta/\theta<H\}.
\]

供给节点 \(\beta\) 只允许向 \(a\ge\beta\)、\(wa\le H\) 送一单位流。目标节点 \(a\) 接收 \(n_a\) 单位。其碰撞成本 \(cn_a^2\) 是凸的，第 \(k\) 个单位的边际成本是 \(c(2k-1)\)。向 \(a\) 输送的线性成本是 \(2wa\)，最后减去常数 \(2\sum_\beta m_\beta w\beta\)。

因此无需枚举指数多的逐行选择，可以用精确整数最小费用流得到这个松弛问题的全局最小值。

附带的 [collision_transport_lab.py](collision_transport_lab.py) 使用保守有理常数 \(c_0=1/6<c\)。这可由 \(\ln2=2\operatorname{atanh}(1/3)>2/3\) 精确证明。全部边成本乘以 \(6\)，从而只有整数运算；实际代码也允许有理 \(\nu\)。

它不是求解真正行列式项的工具：当前没有把实际 \(j,s\)、尾项消失和跨横向指标的线性依赖全部加入。其作用是让“碰撞和误差应当合起来估计”成为一个可复核的精确对象。

## 5. 连续模型：可以进一步研究的对象

固定 \(m,K,v_0,\theta,w_i\)，令 \(H\to\infty\)。取缩放坐标

\[
x_i=w_i\beta_i/H,\qquad y_i=w_ia_i/H,
\qquad |x|=\sum_i x_i.
\]

原始供给的形式密度是

\[
r(x)=\frac K{v_0}\left(1-\frac{|x|}{\theta}\right)
\mathbf1_{|x|<\theta}.
\]

若目标占用 \(n_a/H\) 逼近密度 \(f(y)\)，则供给与目标之间必须有坐标单调运输 \(y\ge x\)，目标支持在 \(|y|\le1\)。联合预算的候选连续泛函为

\[
\frac{c\int f(y)^2\,dy+
2\left(\int |y|f(y)\,dy-\int |x|r(x)\,dx\right)}{\int r(x)\,dx}.
\]

这里给出的是从 Riemann 和得到的候选极限模型；本文没有证明离散最优值到连续最优值的收敛，也没有求解一般的坐标单调运输问题。对固定的光滑运输试探，Riemann 和计算可以逐个验证。

不移动指标，即 \(a=\beta\)，对应 \(f=r\)，误差收益为零，碰撞收益的精确形式极限是

\[
\frac{c\int r^2}{\int r}=\frac{2cK}{(m+2)v_0}.
\]

该试探可能对应真实零行列式，所以不能把它直接当成真实最坏项。它指出下一项值得研究的结构：哪些低费用运输根本不可能产生一个非零 Taylor 系数子式？

## 6. 有限核验与适用范围

运行：

```sh
python3 research/pi-endpoint-rewrite-2026-11-01/collision_transport_lab.py
```

结果保存在 [collision_transport_results.json](collision_transport_results.json)。三个极小模型均用独立穷举全部分配交叉核对了全局最优值。

另外用 \(m=2,w=(1,3),\theta=1/2,v_0=1,K=2\) 做了四个预算模型：

| \(H\) | 行数 \(M\) | 联合收益下界的精确最小值（\(c_0=1/6\)） | 理想主算术成本 \(MH-\sum w\beta\) |
|---:|---:|---:|---:|
| 4 | 12 | \(40/3\) | 44 |
| 6 | 24 | \(112/3\) | 128 |
| 8 | 44 | \(248/3\) | 300 |
| 10 | 72 | 160 | 600 |

这四个例子没有产生理想正间隙。权重 \((1,3)\) 没有被声称满足论文所需的渐次分离，插值满秩没有被声称成立，所有实际误差也尚未加入；它们核验的是费用流模型及预算含义。正的有限模型间隙同样不能代替全参数估计或证明 \(\pi\) 的性质。

本轮获得了两项可用结果：一个排除简单均值修补的严格障碍定理；一个保留全部逐行横向指标、允许接入实际清分母成本的联合预算工具。尚待闭合的是实际非零子式的零结构与连续/渐近联合预算，而非反证法的逻辑形式。
