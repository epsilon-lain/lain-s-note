# OpenAI 的 π 无理指数论文：向自身底数高斯模型迁移什么？

2026-10-31。lain 提问，Mirai 推导与核对。

**状态：没有证明或证伪 P2848。** 本文在采用外部论文主定理的前提下，推出 π 高斯和的次指数增长；给出一般无理指数与峰高增长率的关系；检查 π 单位归一化是否改变行列式比较，并量化原比较方案在端点附近的一项必要尺度条件。这里的“新”仅指本仓库新增的推导，不主张文献首创。未运行 Lean，也没有独立审核外部论文的全部证明。

## 1. 先对齐来源和已有工作

这里指 OpenAI 仓库中的 *The irrationality exponent of π is 2*，版本固定为 `4eaee87548250b2682affda7d05d947c05554bae`：

- [论文 PDF](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/paper.pdf)。
- [形式化范围说明](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/lean/docs/017.md)：陈述 \(\mu(\pi)=2\)，即每个固定 \(\nu>2\) 都有最终的 \(|\pi-p/q|\ge q^{-\nu}\) 下界。该说明不把 Flint–Hills 推论包含在所选形式化陈述内。
- [行列式估计源码](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/determinant.tex)。
- [参数选择与最终证明源码](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/conclusion.tex)。

已核对上述范围说明及相关源码；以下对 π 的推论明确依赖该外部主定理。

此前的 [OpenAI 论文迁移记录](../OPENAI_PI_EXPONENT_P2848_2026-10-27.md) 已经做过这些工作，不能重新算作本次突破：第 10 节证明原参数余量至多是 \(O((\nu-2)^2)\)；第 11 节保留实际误差系数 \(\kappa_i\)；第 13 节分析原有分支估计的障碍；第 19 节构造多维差值消去，并发现相同差值因子也进入算术下界。本文接着这些结果推进。

## 2. 可以迁移的主结论：无理指数决定峰高的指数增长率

沿用 [README](README.md) 的定义，固定无理数 \(\beta>1\)：

\[
F_{p,q}(t)=\frac{\beta^{2t}(p-\beta q)^2+\beta^{-2t}(p+\beta q)^2}{2\beta},
\qquad
\Theta_\beta(t)=\sum_{(p,q)\in\mathbb Z^2}e^{-F_{p,q}(t)}.
\]

这里的 \(t\) 是自身底数的尺度时间，不是进制数位的编号。\(\mu(\beta)\) 采用通常的无理指数定义。

**命题。** 对 \(2\le\mu(\beta)\le\infty\)，有

\[
\boxed{\displaystyle
\limsup_{t\to\infty}\frac{\log\Theta_\beta(t)}{t}
=\left(1-\frac2{\mu(\beta)}\right)\log\beta,}
\tag{1}
\]

其中约定 \(2/\infty=0\)。这是本高斯模型的一般关系，不使用 π 的特殊公式。

### 2.1 从最短格向量到高斯和

记 \(\lambda(t)^2=\min_{(p,q)\ne(0,0)}F_{p,q}(t)\)。该格的协体积为 1。README 第 2 节的本原向量估计与第 4 节的圆盘计数给出绝对常数 \(C\)，使

\[
\frac{\sqrt\pi}{\lambda(t)}\le\Theta_\beta(t)
\le C\bigl(1+\lambda(t)^{-1}\bigr).
\tag{2}
\]

具体地，最短向量必本原。若 \(\lambda^2\le1\)，README 的余项函数 \(B(f)\) 在 \(0<f\le1\) 有界，得到上式；若 \(\lambda\ge1\)，直接用圆盘计数得到统一上界。下界来自最短向量的整数倍及一维 Poisson 求和。

### 2.2 上界：有理逼近下界约束所有时刻

先设 \(\mu(\beta)<\infty\)，固定 \(\nu>\mu(\beta)\)。存在 \(c_\nu>0\)，使所有整数 \(p\) 及非零整数 \(q\) 满足

\[
|\beta-p/q|\ge c_\nu |q|^{-\nu}.
\]

最终下界可以这样扩展到所有分母：对有限个较小分母取正的最小误差常数；无理性保证该常数不为零。

令 \(a=\beta^{2t}\ge1\)。由定义的第一项，以及对 \(p\) 配方，分别有

\[
F_{p,q}(t)\ge\frac{c_\nu^2}{2\beta}a|q|^{-2(\nu-1)},
\qquad
F_{p,q}(t)\ge\frac{2\beta}{a+a^{-1}}q^2\ge\beta a^{-1}q^2.
\]

把第一下界乘以第二下界的 \(\nu-1\) 次方，消去 \(q\)，得到

\[
F_{p,q}(t)\ge C_{\beta,\nu}a^{-(\nu-2)/\nu}.
\]

当 \(q=0,p\ne0\) 时，\(F\ge a/(2\beta)\)，也满足相同形式的下界。因此由 (2)，

\[
\limsup_{t\to\infty}\frac{\log\Theta_\beta(t)}t
\le\frac{\nu-2}{\nu}\log\beta.
\]

令 \(\nu\downarrow\mu(\beta)\) 即得所需上界。

若 \(\mu(\beta)=\infty\)，把格矩阵写成
\(\operatorname{diag}(\beta^t,\beta^{-t})A_\beta(0)\)。最小奇异值给出
\(\lambda(t)\ge c_\beta\beta^{-t}\)，所以 (2) 给出增长率至多 \(\log\beta\)。

### 2.3 下界：好逼近产生相应高度的峰

若 \(2<\nu<\mu(\beta)\)，存在无穷多个既约分数 \(p/q\)，\(q\to\infty\)，满足 \(|\beta-p/q|\le q^{-\nu}\)。记

\[
\kappa=q^2|\beta-p/q|,\qquad \ell=\log(1/\kappa).
\]

则 \(\ell\ge(\nu-2)\log q\)。该向量的平衡时刻和最小长度平方为

\[
t_0=\frac{\log|(p+\beta q)/(p-\beta q)|}{2\log\beta}
=\frac{\log q+\ell/2}{\log\beta}+O_\beta(1),
\qquad
m=\frac{|p^2-\beta^2q^2|}{\beta}\sim2\kappa.
\]

沿该本原方向的高斯和至少为 \(\sqrt{\pi/m}\)，故
\(\log\Theta_\beta(t_0)\ge\ell/2+O_\beta(1)\)。于是

\[
\limsup_{t\to\infty}\frac{\log\Theta_\beta(t)}t
\ge\frac{\nu-2}{\nu}\log\beta.
\]

令 \(\nu\uparrow\mu(\beta)\)，或在 \(\mu=\infty\) 时令 \(\nu\to\infty\)。当 \(\mu=2\) 时，下界 0 直接来自 \(\Theta\ge1\)。这就证明 (1)。

### 2.4 对 π 和当前积分路线的含义

采用外部论文的 \(\mu(\pi)=2\)，(1) 给出

\[
\boxed{\lim_{t\to\infty}\frac{\log\Theta_\pi(t)}t=0.}
\tag{3}
\]

等价地，对每个 \(\eta>0\)，存在有限常数 \(C_\eta\)，使全部 \(t\ge0\) 上都有 \(\Theta_\pi(t)\le C_\eta e^{\eta t}\)。本文没有从外部主定理提取可计算的 \(C_\eta\)。

对固定 \(H>0\) 定义 \(J_H(s)=\int_s^{s+H}\Theta_\pi(t)^2\,dt\)。从 \(H\le J_H(s)\le H C_\eta^2e^{2\eta(s+H)}\)，立即得到

\[
\lim_{s\to\infty}\frac{\log J_H(s)}s=0.
\tag{4}
\]

然而，次指数增长与有界性是两个问题。\(\sqrt2\) 的 \(\mu=2\)，对应高斯和周期且有界；[稀疏高峰记录](SPARSE_PEAK_OBSTRUCTION.md) 的显式连分数则有 \(\mu=2\) 且高斯和无界。后者由 \(a_{2^j}=j+1\)、其他 \(a_n=1\) 构造，\(\log a_{n+1}/\log q_n\to0\)，用收敛分数误差公式和 Legendre 判据即可验证其无理指数为 2。

因此 (3)—(4) 不决定 P2848，也不决定长期平方平均是否发散。它们说明：若 π 存在无限增高的峰，这些峰必须在次指数尺度上增高。长期平均有界仍不能证伪 P2848。

## 3. 把周期改成自身单位：哪些步骤真的改变？

原论文使用 \(\omega=2\pi i\)、有理逼近中心 \(r_i=2ip_i/q_i\)，以及 \(e^{j\omega}=1\)。令 \(z'=z/\pi\)，则周期变成 \(2i\)，而指数函数变成 \(e^{\pi z'}\)。写法变简单，但其 Taylor 系数开始携带目标数 π。

可以精确追踪同一个子式的变化。原矩阵的行记为 \(\rho=(j,s,b)\)，列记为 \((h,a)\)，其中 \(a,b\) 为多重指标。原元素形如

\[
M_{\rho,(h,a)}=[t^su^b](1+t)^h\prod_i\bigl(jr_i+G_i(t)+u_i\bigr)^{a_i}.
\]

在每个横向坐标中同时改用 π 单位，取
\(r_i'=r_i/\pi\)、\(G_i'=G_i/\pi\)、\(u_i'=u_i/\pi\)。归一化矩阵满足精确恒等式

\[
M'_{\rho,(h,a)}=\pi^{|b|-|a|}M_{\rho,(h,a)}.
\]

对任何原来选定的方形子式，令
\(d_H=\sum_{\mathrm{columns}}|a|-\sum_{\mathrm{rows}}|b|\)，则

\[
\boxed{\Delta'_H=\pi^{-d_H}\Delta_H.}
\tag{5}
\]

因此原来的解析上界和算术下界同时乘以 \(\pi^{-d_H}\)；比较的余量没有改变。归一化中心 \(2ip_i/(\pi q_i)\) 一般不再是高斯有理数，不能继续声称仅用原有整数分母就能把 \(\Delta'_H\) 变成高斯整数。若使用原清分母因子 \(D\)，正确的关系是 \(D\pi^{d_H}\Delta'_H=D\Delta_H\in\mathbb Z[i]\)。

这不是排除所有 π 进制方法，而是精确排除“只改坐标，就免费获得更强的算术—解析比较”这一步。真正利用非整数进制的可容许数位、不同的插值函数或新算术结构，需要另行证明。

## 4. 原有幂次比较能分辨多慢的小系数？

旧记录第 10 节已经证明，原文的参数余量满足

\[
0<g=\nu\bigl(A(1-\eta)-\theta\bigr)-(1-\theta)
<\frac{(\nu-2)^2}{4(\nu-1)}.
\tag{6}
\]

原文最终误差 \(E=E_{\rm ar}+E_{\rm an}\) 的各项为正，包含 \(\nu/w_*\)，其中 \(w_* =\min_i\lceil\log q_i\rceil\)。原比较要求 \(E<g\)，所以必要条件是

\[
\boxed{w_*>\frac{4\nu(\nu-1)}{(\nu-2)^2}.}
\tag{7}
\]

这是本次从已有二次余量继续得到的尺度限制。它不是充分条件；维数、权重分离及其他误差项还会加上限制。

为了把小系数 \(\kappa_i\) 塞回原来的幂次假设，必须选择共同的 \(\varepsilon=\nu-2>0\)，使每个选中的逼近都满足

\[
\kappa_i\le q_i^{-\varepsilon},\qquad
\ell_i:=\log(1/\kappa_i)\ge\varepsilon\log q_i.
\]

取达到最小权重的指标 \(j\)。由 (7)，

\[
\boxed{\ell_j>
2\sqrt{\nu(\nu-1)}\,\frac{\log q_j}{\sqrt{w_*}}.}
\tag{8}
\]

由于 \(w_*=\lceil\log q_j\rceil\)，在 \(\nu\downarrow2\)、\(q_j\to\infty\) 时，右端约为 \(2\sqrt{2\log q_j}\)。也就是说，若仍保留这套误差估计并把实际误差重新压成共同的幂次条件，则至少需要 \(\log(1/\kappa)\) 达到 \(\sqrt{\log q}\) 的量级。像 \(\kappa=1/\log q\) 这样缓慢趋零的尺度，过不了该必要条件；这只是说明比较方案的分辨能力，不声称 π 的误差实际遵循该例。

**适用范围必须保留：** (7)—(8) 只约束原文这套充分比较条件。它们不排除直接保留 \(\kappa_i\)、改善误差项、选择不同子式或不同函数。也不能把固定 \(\nu\) 的证明直接对角化成 \(\nu(q)\downarrow2\) 的定理；原文先固定维数，再固定逼近尺度，最后才令多项式次数趋于无穷。

## 5. 现在值得继续迁移的部分

原论文最有价值的可迁移部分是：多点插值保证子式不为零；清分母给出离散算术下界；精确周期允许平移；同一横向组内的 Taylor 次数碰撞产生严格相消。单纯的高斯峰图像尚未提供对应的非零算术对象。

下一步需要针对真实目标工作：

\[
\kappa_i\to0,\qquad
\ell_i=\log(1/\kappa_i)\to\infty,
\qquad
\ell_i/\log q_i\to0.
\]

最后一式在采用 \(\mu(\pi)=2\) 后必然成立。它并不提供前两式；前两式就是仍待证明的非坏逼近方向。

旧记录第 11 节已经保留了 \(\kappa_i\) 的精确幂次，所以应继续研究其中的次主项收益，并与全部算术成本比较。第 19.7 节又表明，已找到的差值因子同时整除清分母后的高斯整数；两侧消去共同因子后，真正要估计的是剩余行列式 \(Q_H\)。有效的新进展必须给出未被算术下界同步吸收的收益，或给出从 π 独立表示到固定窗口峰高的新估计。

若要证明 P2848，仍需直接推出 \(\sup_t\Theta_\pi(t)=\infty\)；若要证伪，则需统一上界。本文的增长率关系、坐标缩放恒等式和尺度必要条件均未完成这一步。

## 6. 检查范围

本文件给出了 (1)—(8) 的推导；关键输入是 README 已证明的高斯和估计，以及明确标出的外部论文主定理和误差公式。另用 SymPy 检查了一个四阶插值矩阵的行列式缩放，以及 (6) 所用的配方恒等式。有限符号检查只是防止代数抄写错误，不替代上面的普遍证明，也不属于对外部论文或 P2848 的形式化验证。
