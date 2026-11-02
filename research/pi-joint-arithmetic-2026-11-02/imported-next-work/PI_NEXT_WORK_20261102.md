# π 辅助多项式的近根与算术成本

2026 年 11 月 2 日。核查基准为 [epsilon-lain/lain-s-note 的 befa7cf7 提交](https://github.com/epsilon-lain/lain-s-note/tree/befa7cf7a515a18c3bd00ef17a4c6f53dd9ae6df)。本轮在该构造内补强了近根定位，证明完整反射因子的系数高度主阶，并建立经典 Padé 对照。下一阶段优先研究能够实际输出整数线性形式的构造，同时核算误差与约分后分母。

本轮结论没有证明 π 的坏逼近性或非坏逼近性。下文的新推导是相对该仓库检查点的推进，没有声称其为文献中首次出现的定理。独立复核指其他研究代理逐步检查证明，不是外部同行评审。

## 近根可以加强为唯一简单的纯虚根

沿用原来的双频率一维多项式

\[
P_n(Z)=\det\bigl([z^s]e^{hz}(jZ+z)^a\bigr),\qquad
j,h\in\{0,1\},\quad 0\le s,a<n,
\]

并记 \(\omega=2\pi i\)。充分大的整数 \(n\) 下，最近于 \(\omega\) 的根 \(\zeta_n\) 唯一、简单，并严格位于虚轴，且

\[
\boxed{\zeta_n-\omega
=-\frac{P_n(\omega)}{P_n'(\omega)}\bigl(1+O(n^{-2})\bigr).}
\]

距离的线性项也被确定：

\[
\boxed{\log|\zeta_n-\omega|
=-2n\log n-\bigl[\log(16/\pi^2)-2\bigr]n+O(\log n).}
\]

偶数阶最终从虚轴下方接近，奇数阶最终从上方接近。充分大是渐近条件；本轮没有认证最小起始阶数。

关键新增步骤是对所有端点插入阶数作统一估计。令 \(A_j=\prod_{a=0}^{j-1}a!\)、\(\Delta_j=A_j^4/A_{2j}\)，并约定 \(A_0=\Delta_0=1\)。定义

\[
K_n=\frac{[(2n-2)!]^2(2n-1)!}{[(n-1)!]^4(2\pi)^{2n-1}}.
\]

将行列式写成真实矩矩阵与端点扰动之和，混合 Vandermonde 积分给出第 \(k\) 阶系数的阶乘上界。它们满足几何型控制：首项为 \(K_n\)，其余项的逐阶主比例可用 \(K_{n-1}\) 控制，而

\[
\frac{K_{n-1}}{K_n}
=\frac{(n-1)(2\pi)^2}{8(2n-3)^2(2n-1)}=O(n^{-2}).
\]

所以在半径为固定常数乘 \(K_n^{-1}\) 的圆盘内，高阶余项整体趋于零。Rouché 定理给出恰有一个简单根；实系数与奇偶对称再迫使它位于虚轴。这个论证补上了从大导数到唯一近根所需要的高阶控制。

完整证明见 `root_geometry/UNIQUE_NEAR_PERIOD_ROOT_LATEX.md`，独立复核见同目录 `INDEPENDENT_UNIQUE_ROOT_AUDIT.md`。固定 \(Z=\omega\) 的独立端点变量也有唯一近根，且落在单位圆上；它与上述 \(Z\) 根仍是不同对象。

## 完整反射因子的高度有明确主阶

偶数 \(n=2r\) 时，仓库已经证明

\[
d_nP_n(Z)=(-1)^rZ^nF_+(Z)F_-(Z),\qquad
d_n=B_n^2,\quad B_n=A_n/A_r^2,
\]

其中 \(F_\pm\) 是首一偶整数多项式，次数均为 \(r(r-1)\)。本轮进一步证明其常数项的精确绝对值：

\[
|F_+(0)|=(n-1)B_n,\qquad |F_-(0)|=B_n.
\]

证明通过提取反射行列式的最低次项，再用 Pascal 矩阵共轭化简完成，适用于全部偶数阶。配合显式系数范数界

\[
\|F_\pm\|_1\le B_n n!2^n(3/2)^{n(n-1)/2},
\]

可得

\[
\boxed{\log H(F_\pm)=\frac{n^2}{4}\log n+O(n^2).}
\]

这里 \(H\) 专指最大系数绝对值；这不是将数的绝对对数 Weil 高度与多项式系数高度混用。Mahler 测度的对数也具有相同主阶。对包含最近根的那个完整反射因子，

\[
\frac{-\log|\zeta_n-\omega|}{\log H(F_\pm)}
=\frac{8+o(1)}n\longrightarrow0.
\]

这量化了完整因子的成本，尚未控制可能存在的更小不可约因子。有限模素数证书证明 \(n=4,6,8,10,12\) 的十个反射因子均在 \(\mathbb Q\) 上不可约；一般不可约性仍未证明。

完整证明及证书见 `arithmetic_route/`。有理求值的已知障碍仍保留：对次数为 \(D\) 的首一偶整数因子，在既约 \(ia/b\) 处的约分后分母恰为 \(b^D\)。改写成 \(Z^2\) 的多项式不会再节省这个指数。

## 经典 Padé 给出线性次数的对照

经典指数函数的对角 Padé 公式及 Beta 积分余项见 [Driver 与 Temme，1999](https://ir.cwi.nl/pub/1229/1229D.pdf)。用它们的奇部可以得到下列首一偶整数多项式。对奇数 \(m\ge3\)，定义

\[
C_m(Z)=\sum_{\substack{1\le k\le m\\k\text{ 为奇数}}}
\frac{(2m-k)!}{k!(m-k)!}Z^{k-1}.
\]

它的次数只有 \(m-1\)。本轮从精确余项推导出：其最近根同样满足

\[
\log|\zeta_m-2\pi i|
=-2m\log m+[2+2\log(\pi/2)]m+O(\log m).
\]

同时，

\[
\log|C_m(2\pi i)|
=-m\log m+(1+2\log\pi)m+O(\log m),
\]

\[
H(C_m)=\frac{(2m-1)!}{(m-1)!},\qquad
\operatorname{den}C_m(ia/b)=b^{m-1}.
\]

因此，在直接比较中心小值与有理求值下界的这条方法中，仍需 \(m\) 与 \(b\) 同阶，渐近阈值为 \(m/b=e\pi^2\)；所得下界仍为 \(\exp[-O(b\log b)]\)。这个对照证明：单独把辅助次数从二次降到一次，并不足以跨到临界的 \(q^{-2}\) 尺度。该多项式族是经典对象，本轮没有声称它是原来 \(P_n\) 的因子。

完整推导、文献映射、数值示例和独立复核见 `literature/`。

## 下一轮应完成的算术工作

优先目标是构造并证明整数线性形式

\[
L_j=q_j\pi-p_j,\qquad p_j\in\mathbb Z,\quad q_j\in\mathbb Z_{>0},
\quad \gcd(p_j,q_j)=1,
\]

并联合控制 \(q_j\) 与 \(L_j\)。如果进攻 π 非坏逼近这一方向，需要的是 \(q_j\to\infty\) 且

\[
\boxed{q_j|L_j|
=q_j^2\left|\pi-\frac{p_j}{q_j}\right|\longrightarrow0.}
\]

单个构造的失败不会证明 π 坏逼近；证明坏逼近需要对所有既约有理数建立统一正下界。

下一轮可以按下面的顺序开展，每步都有可检验的产出：

1. **建立整数线性形式的基准。** 重建一个在固定代数参数处求值的对数或反正切 Padé 积分模型，逐步核对系数清分母及整数公因子。以经典 π 无理性测度构造为对照，记录实际的约分后分母和误差，而不只记录解析积分很小。
2. **检验原端点结构能否迁移到这个基准。** 保留精确相消与非零性，把原先产生高次数多项式的步骤替换成对 π 的一次形式。每个候选同时列出整数系数、清分母成本、约分收益和严格误差界。只有联合成本改善才构成算术推进。
3. **继续追踪最近根的最小多项式。** 将当前有限不可约证书推广，或找出确实包含最近根的低成本因子。有限实验用于发现结构；一般结论需要统一证明。进一步的高维构造也应以联合成本改进作为验收条件。

上述第一项的经典参照包括 [Zeilberger 与 Zudilin 的 π 无理性测度论文](https://sites.math.rutgers.edu/~zeilberg/mamarim/mamarimPDF/pimeas.pdf)。它提供的是整数线性形式与算术消因子的范例；本文没有将其中的历史数值当作当前最佳界。

## 文献工具的适用边界

- [The kissing polynomials and their Hankel determinants](https://www.damtp.cam.ac.uk/user/na/NA_papers/NA2015_01.pdf) 与本项目的真实指数权重矩行列式有精确变量替换关系，可借鉴其多重积分与正交多项式结构。它对应 \(G_n(Z,e^Z)\)，不是全局的 \(G_n(Z,1)\)；两者在周期处的值相同，导数需保留端点链式修正。
- [Cullinan 与 Scheel，On the arithmetic of Padé approximants to the exponential function](https://arxiv.org/abs/2007.01329) 提供 Newton 多边形、不可约性与 Galois 群工具。它的不可约性结论不能直接移植到本项目的反射因子或上述奇部。
- 本轮没有独立认证 OpenAI 原稿的全部高维结论。即使某条证明建立无理性指数为 2，也仍需分别研究临界尺度上的统一正下界与趋零子序列。

## 复现与证据

在解压目录中运行：

```bash
python source-audit/verify_nearest_root_scale.py
python source-audit/verify_period_log_derivative.py
python root_geometry/verify_unique_root.py
python arithmetic_route/verify_height.py
python arithmetic_route/verify_modular_irreducibility.py
python literature/verify_pade_benchmark.py
```

原检查点的两组精确验证、本轮混合端点系数恒等式、Laguerre 阶乘行列式、高度恒等式和有限不可约性证书均通过。高精度根计算仅作示例，不承担渐近定理的证明。原检查点和高度验证使用 Python 标准库；新根示例及 Padé 示例使用 mpmath，模素数证书使用 SymPy。

新根定位依赖原检查点已经建立的周期非零下界和真实端点迹渐近。完整证明中的这些输入均显式保留，未用有限测试替代。
