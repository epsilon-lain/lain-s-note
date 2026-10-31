# 两中心插值行列式：显式列约化、分母与系数界

2026-10-31。接续用户提出的“约去共同因子后，显式构造剩余矩阵”。

**状态：得到一般代数恒等式、显式剩余列、整数整除和增长界；没有证明或证伪 P2848。** 本文仅处理原构造的两个中心索引 \(j=0,1\)，即 \(K=2\)。结果不需要采用 \(\mu(\pi)=2\)，也没有运行 Lean。有限精确检查用于核对公式和实现；一般结论由下文证明承担。“新增”指此前仓库没有这份构造，不主张文献首创。

本次核对的仓库基线为 f75e3d37c625c3bee276bbb238a47077106561df。此前 [OpenAI 论文迁移记录第 19 节](../OPENAI_PI_EXPONENT_P2848_2026-10-27.md) 已有特殊列、逐列差值因子及其算术抵消；本文补充列与列之间的相消，并给出实际可构造的剩余矩阵。

上传前另核对了 174bebe9bfef96f1ddf415209569c08f7ad5bd59 新增的[余式反例记录](../self-base-theta-2026-10-31/RESIDUAL_DETERMINANT_COUNTEREXAMPLES.md)。本文的列约化与该记录的小余式反例不同；新记录和旧笔记第 21 节均已保留。

## 1. 结论先行

固定一组旧构造中的 \((h,a)\)，记 \(U=U_h\)、\(B=B_h\)，其中 \(0\le B\le U\)。令 \(\delta=r_2-r_1\) 为两个中心步长的差。该组列可精确写为

\[
C_j(\delta)=\delta^{j+1}\sum_{s=0}^{U}
\binom{U+1+j}{U-s}\delta^sV_s,\qquad 0\le j\le B,
\tag{1}
\]

其中向量 \(V_s\) 与列编号 \(j\) 无关。

定义

\[
c_b=\binom{2b}{b},\qquad
T_{j,b}=\frac{(-1)^{b+j}}{c_b}
\binom{b+j}{j}\binom{U+b+1}{b-j},
\quad 0\le j\le b.
\tag{2}
\]

则 \(T_{b,b}=1\)，且有精确列变换

\[
\boxed{
\sum_{j=0}^{b}T_{j,b}\delta^{b-j}C_j(\delta)
=\delta^{2b+1}R_b(\delta),
}
\tag{3}
\]

其中剩余列有闭式

\[
\boxed{
R_b(\delta)=
\frac{(-1)^b}{c_b}
\sum_{s=b}^{U}
\binom{s}{b}\binom{U+b+1}{U-s}
\delta^{s-b}V_s.
}
\tag{4}
\]

因此这组列的共同因子从逐列得到的
\(\delta^{(B+1)(B+2)/2}\)，加强为
\(\delta^{(B+1)^2}\)。式 (3) 是对角为 1 的三角列变换，变换本身不改变行列式。

对旧记录所有特殊列同时操作，记所得完整剩余矩阵为 \(\mathcal R_H\)，补全用的其他列保持原样，则

\[
\boxed{
\mathcal D_H(r,r+\delta)=\delta^{E_H}\det\mathcal R_H(r,\delta),
\qquad
E_H=\sum_h(A_h+1)(B_h+1)^2.
}
\tag{5}
\]

这里约去的是已经严格保证的因子，**不宣称 \(E_H\) 是整个行列式的最大消失阶**。不同 \((h,a)\) 组的剩余列还可能相关；第 8 节给出继续判定的递推。

## 2. 原矩阵、适用条件与向量 \(V_s\)

沿用旧记录第 19 节：正有理权重满足

\[
0<w_0<v_0,\qquad w_1\le w_2,\qquad 0<\theta<1.
\]

若需要沿用原论文保证完整行列式非零，还必须满足它的插值体积与权重分离条件。本文的列恒等式本身不替代这些条件。

允许行的索引记为 \((\jmath,\ell,\beta_1,\beta_2)\)，其中

\[
\jmath\in\{0,1\},\qquad
v_0\ell+(w_1\beta_1+w_2\beta_2)/\theta<H.
\]

取 \(F_0>2/\theta\)，并定义

\[
T_i=\left\lceil F_0w_i/v_0\right\rceil,\qquad
G_i(t)=\sum_{k=1}^{T_i-1}(-1)^{k+1}t^k/k.
\]

空和取 0。由 \(w_1\le w_2\) 有 \(T_1\le T_2\)，且 \(G_2-G_1\) 的每个非零项的 \(t\)-权重大于 \(w_1/\theta\)。

旧构造选择固定的 \(\tau>0\)，满足

\[
2\tau<1,\qquad
w_0/v_0+\theta\tau(1+2w_2/w_1)<1,
\]

并在

\[
(1-2\tau)H/v_0\le h\le(1-\tau)H/v_0
\]

的整数 \(h\) 上取

\[
U_h=\left\lceil\theta(H-v_0h)/w_1\right\rceil-1,\qquad
A_h=B_h=\left\lfloor\theta(H-v_0h)/(4w_2)\right\rfloor.
\]

这里 \(B_h\le U_h\)。对 \(0\le a\le A_h\)、\(0\le b\le B_h\)，特殊列多项式为

\[
P_{h,a,b}(Y,X_1,X_2)
=(Y-1)^hX_1^a(X_2-X_1)^{U_h+1+b}.
\]

以下取 \(H\) 充分大，并显式要求所有特殊列满足
\(w_0h+w_1a+w_2(U_h+1+b)\le H\)。
上述 \(\tau\) 条件保证这一点最终成立。它确保第 6 节对截断分母的预算适用于每个 \(V_s\)，不把任意小 \(H\) 的实验当成原渐近构造。

固定 \((h,a)\)，令 \(r=r_1\)、\(U=U_h\)，并写

\[
W=u_2-u_1+G_2(t)-G_1(t).
\]

在中心 \(\jmath\) 处，特殊列的系数来自

\[
t^h(\jmath r+G_1(t)+u_1)^a
(\jmath\delta+W)^{U+1+b}.
\tag{6}
\]

每个 \(W\) 因子至少消耗 \(w_1/\theta\) 的行预算，所以可见项中至多有 \(U\) 个 \(W\) 因子。在中心 \(\jmath=0\) 处，(6) 的全部允许行系数因此为零。在中心 \(\jmath=1\) 处，令 \(V_s\) 的该行分量为

\[
[t^\ell u_1^{\beta_1}u_2^{\beta_2}]
t^h(r+G_1(t)+u_1)^aW^{U-s},
\qquad 0\le s\le U,
\tag{7}
\]

并把中心 0 的分量设为零。展开 (6) 就得到 (1)。因此 \(V_s\) 确实独立于列编号。这一步依赖 \(\jmath=0,1\)，不能不加处理地推广到更多中心索引。

## 3. 列相消与闭式变换的证明

先只证明因子数。把 (1) 中的 \(B+1\) 列代入行列式的多线性展开。若两列选择相同的 \(s\)，对应向量都是 \(V_s\) 的倍数，该项为零。非零项必须选择互异的 \(B+1\) 个非负整数，其和至少是 \(B(B+1)/2\)。加上已有的 \(\sum_{j=0}^{B}(j+1)\)，总阶数至少为 \((B+1)^2\)。不同组同时使用这个论证即可得到 (5) 中的整除。

要得到显式剩余矩阵，使用整数恒等式

\[
\sum_{j=0}^{b}(-1)^j
\binom{b+j}{j}\binom{U+b+1}{b-j}
\binom{U+1+j}{U-s}
=
\begin{cases}
0,&s<b,\\
\binom{s}{b}\binom{U+b+1}{U-s},&s\ge b.
\end{cases}
\tag{8}
\]

**证明 (8)。** 提取阶乘后，左侧等于

\[
\frac{(U+b+1)!}{(U-s)!\,(b!)^2}
\sum_{j=0}^{b}(-1)^j\binom bj
\frac{(b+j)!}{(s+j+1)!}.
\]

当 \(s<b\) 时，阶乘比是关于 \(j\) 的次数 \(b-s-1<b\) 的多项式，其 \(b\) 阶有限差分为零。当 \(s\ge b\) 时，

\[
\frac{(b+j)!}{(s+j+1)!}
=\frac1{(s-b)!}\int_0^1x^{b+j}(1-x)^{s-b}\,dx.
\]

求和后出现 \((1-x)^b\)，积分等于
\(b!s!/((s-b)!(b+s+1)!)\)。代回即得 (8)。

将 (1) 代入 (3)，使用 (8)，再除以 \((-1)^bc_b\)，即得到 (4)。所有求和有限，没有交换无穷极限。

旧记录的保证阶数是
\[
E_H^{\mathrm{old}}=\sum_h(A_h+1)(B_h+1)(B_h+2)/2.
\]
固定参数并令 \(H\to\infty\)，新的保证阶数满足
\[
\frac{E_H}{H^4}\longrightarrow
\frac{15\theta^3\tau^4}{256v_0w_2^3},
\qquad
\frac{E_H}{E_H^{\mathrm{old}}}\longrightarrow2.
\tag{9}
\]
这是对保证因子数的加强，不是算术—解析净余量加倍。

## 4. 分母：精确公式与足够界

式 (2) 的所有分母均整除
\[
c_b=\binom{2b}{b}\le4^b.
\tag{10}
\]
设
\[
n_{j,b}=\binom{b+j}{j}\binom{U+b+1}{b-j}.
\]
该变换列所需的**最小共同分母**精确等于
\[
\frac{c_b}{\gcd(c_b,n_{0,b},\ldots,n_{b,b})}.
\tag{11}
\]
式 (4) 的标量系数同样以 \(c_b\) 为足够分母，但不声称其最小分母总与 (11) 相等。

所有特殊列的归一化分母乘积记为
\[
\mathcal C_H=\prod_h\prod_{a=0}^{A_h}\prod_{b=0}^{B_h}c_b.
\]
于是
\[
\log\mathcal C_H
\le(\log2)\sum_h(A_h+1)B_h(B_h+1).
\tag{12}
\]
这个上界通常是 \(H^4\) 量级，不能作为低阶项丢掉。第 7 节说明怎样在算术比较中准确处理它，避免重复计算代价。

## 5. 次数、系数增长与近对角误差

将条目视为关于 \(r,\delta\) 的有理系数多项式。由 (7)，
\(\deg_rV_s\le a\)。式 (4) 因此给出
\[
\boxed{
\deg_\delta R_b\le U-b,\qquad
\deg_rR_b\le a,\qquad
\deg_{r,\delta}R_b\le a+U-b.
}
\tag{13}
\]

令 \(\gamma_1=\sum_{k=1}^{T_1-1}1/k\)，
\(\gamma_{21}=\sum_{k=T_1}^{T_2-1}1/k\)，空和为零。
对一个多项式，\(\|\cdot\|_1\) 表示全部系数绝对值之和；对列向量取其分量中最大的这个范数。系数提取不会增大该范数，故
\[
\|V_s\|_1\le(2+\gamma_1)^a(2+\gamma_{21})^{U-s}.
\tag{14}
\]

令 \(n=U-b\)，把 (4) 中的 \(s\) 写成 \(b+\ell\)，并定义整数
\[
k_\ell=\binom{U+b+1}{2b+\ell+1}\binom{b+\ell}{b},
\qquad 0\le\ell\le n.
\]
它们有精确递推
\[
\boxed{
k_0=\binom{U+b+1}{2b+1},\qquad
k_{\ell+1}
=k_\ell\frac{(U-b-\ell)(b+\ell+1)}
{(\ell+1)(2b+\ell+2)}.
}
\tag{15}
\]
逐项比较得 \(k_\ell\le k_0\binom n\ell\)。结合 (14)，得到无未指定常数的界
\[
\boxed{
\|R_b\|_1\le
\frac{\binom{U+b+1}{2b+1}}{\binom{2b}{b}}
(2+\gamma_1)^a(3+\gamma_{21})^{U-b}.
}
\tag{16}
\]

如果只在给定中心 \(r\) 和差值 \(\delta\) 处估计分量，则同样有
\[
|R_b(r,\delta)_{\mathrm{row}}|
\le\frac{k_0}{c_b}
(|r|+1+\gamma_1)^a
(2+\gamma_{21}+|\delta|)^{U-b}.
\tag{17}
\]
这些是上界，不声称右侧趋零。

对任意固定 \(r\)，也可直接保留向量：
\[
R_b(0)=\frac{(-1)^bk_0}{c_b}V_b,
\]
\[
\left\|R_b(\delta)-\frac{(-1)^bk_0}{c_b}V_b\right\|
\le\frac{k_0}{c_b}
\bigl((1+|\delta|)^{U-b}-1\bigr)
\max_{s>b}\|V_s\|.
\tag{18}
\]
空余项取 0。这清楚标出了“中心合并以后到底剩下什么”。

还有一个有用的算术核对：
\[
\prod_{b=0}^{B}\frac{\binom{U+b+1}{2b+1}}{\binom{2b}{b}}
=\left|\det\left[\binom{U+1+j}{U-s}\right]_{s,j=0}^{B}\right|
\in\mathbb Z_{>0}.
\]
它来自上述三角化后的主元乘积。因此该乘积至少为 1；当 \(U=B\) 时恰为 1。不能单独把 \(c_b\) 置于分母看作净解析收益。

## 6. 完整剩余行列式的可追踪界

设补全后的方阵大小为 \(M=M_H\)，其他列是固定选定的原单项式 \(Y^hX_1^{\alpha_1}X_2^{\alpha_2}\)。记
\[
D_{\mathrm{tr}}=
\prod_{i=1}^{2}\operatorname{lcm}(1,\ldots,T_i-1)^{\lfloor H/w_i\rfloor},
\]
空的最小公倍数取 1。

原条目的多项式系数分母整除 \(D_{\mathrm{tr}}\)；剩余特殊列的分母整除 \(c_bD_{\mathrm{tr}}\)。因此
\[
D_{\mathrm{tr}}^{M}\mathcal C_H\det\mathcal R_H
\in\mathbb Z[r,\delta].
\tag{19}
\]
这是一项足够的系数清分母界，不声称最优。

对每个特殊剩余列，取 \(d_j=a+U-b\)，并以 (16) 的右侧为 \(L_j\)。
对每个补全单项式列，可取
\[
d_j=\alpha_1+\alpha_2,\qquad
L_j=2^h(2+\gamma_1)^{\alpha_1}(3+\gamma_2)^{\alpha_2},
\]
其中 \(\gamma_2=\sum_{k=1}^{T_2-1}1/k\)；第二个括号的 3 来自 \(r+\delta+u_2\)。

于是由行列式展开，
\[
\deg\det\mathcal R_H\le\sum_jd_j,\qquad
\|\det\mathcal R_H\|_1\le M!\prod_jL_j.
\tag{20}
\]
因此 (19) 中整数多项式的系数绝对值之和至多
\[
D_{\mathrm{tr}}^{M}\mathcal C_H\,M!\prod_jL_j.
\tag{21}
\]
固定维数、权重及截断后，\(M=O(H^3)\)，每列次数为 \(O(H)\)；这些显式界给出总次数 \(O(H^4)\)、清分母后对数系数界 \(O(H^4+H^3\log H)\)。这里的粗渐近只总结 (19)—(21)，不是新的负向衰减估计。

## 7. 算术一侧：同一个因子仍然被除掉

现在代入有理中心
\[
r_i=2\mathrm i\,p_i/q_i,\qquad
q_i>0,\qquad
L=\operatorname{lcm}(q_1,q_2),\qquad
d=L\delta\in\mathbb Z[\mathrm i]\setminus\{0\}.
\]

原清分母方案对特殊列使用 \(q_1^aL^{U+1+b}\)，对行使用
\(q_1^{-\beta_1}q_2^{-\beta_2}\)，并对每个条目乘 \(D_{\mathrm{tr}}\)。记完整清分母行列式为 \(\Gamma_H\)。

本节的整除恒等式不需要非零假设；涉及模至少为 1 的结论，则明确假定选定的完整子式满足 \(\Gamma_H\ne0\)。

对 (7) 的向量定义
\[
W_s=
D_{\mathrm{tr}}q_1^aL^{U-s}
\operatorname{diag}(q_1^{-\beta_1}q_2^{-\beta_2})V_s.
\tag{22}
\]
则 \(W_s\in\mathbb Z[\mathrm i]^M\)：\(q_1^a\) 支付第一括号的中心与横向导数分母；\(L^{U-s}\) 支付剩余 \(W^{U-s}\) 的横向导数分母；\(D_{\mathrm{tr}}\) 支付截断对数的系数分母。逐项分配这些幂即可证明整性，与旧记录第 19.7 节的论证一致。

清分母后的特殊列因此精确成为
\[
\widetilde C_j=d^{j+1}\sum_{s=0}^{U}
\binom{U+1+j}{U-s}d^sW_s.
\tag{23}
\]
直接对 (23) 使用第 3 节的多线性论证，得到
\[
\boxed{d^{E_H}\mid\Gamma_H\quad\text{在 }\mathbb Z[\mathrm i]\text{ 中}.}
\tag{24}
\]
这里没有从“乘上一个额外整数后可整除”倒推整除；(24) 是原整数列展开的直接结论，不需要 \(d\) 与任何二项式系数互素。

同一闭式列变换给出算术剩余矩阵 \(\widetilde{\mathcal R}_H\)，满足
\[
\det\widetilde{\mathcal R}_H=\Gamma_H/d^{E_H}
\in\mathbb Z[\mathrm i]\setminus\{0\}.
\tag{25}
\]
其第 \(b\) 个特殊列为
\[
\widetilde R_b=
\frac{(-1)^b}{c_b}\sum_{s=b}^{U}
\binom{s}{b}\binom{U+b+1}{U-s}d^{s-b}W_s.
\]
各条目可能有分母 \(c_b\)，但整个行列式是高斯整数。与解析剩余列的准确关系是
\[
\boxed{
\widetilde R_b
=D_{\mathrm{tr}}q_1^aL^{U-b}
\operatorname{diag}(q_1^{-\beta_1}q_2^{-\beta_2})R_b.
}
\tag{26}
\]
相对于原列的 \(L^{U+1+b}\)，指数下降恰好是 \(2b+1\)，对应已经约掉的差值幂；不能再把这次下降计为独立收益。

若再把每个特殊剩余列乘以 \(c_b\)，得到整数矩阵 \(\widehat{\mathcal R}_H\)，则
\[
\det\widehat{\mathcal R}_H
=\mathcal C_H\,\Gamma_H/d^{E_H}.
\]
因此其正确的已知算术下界是
\[
|\det\widehat{\mathcal R}_H|\ge\mathcal C_H,
\]
而不是只使用较弱的 1。这防止把清掉 \(c_b\) 的成本重复记入比较。

## 8. 若剩余列仍相关：继续剥离的精确递推

式 (4) 只完成同一 \((h,a)\) 组内部的显式约化。不同组的 \(V_b\) 可能相关。下面给出在固定有限矩阵上进一步找出最大差值阶的代数过程；本文没有由它得到新的关于 \(\pi\) 的渐近估计。

把 \(\mathcal R_H(r,\delta)\) 视为系数域 \(\mathbb Q(r)\) 上关于 \(\delta\) 的矩阵。固定选列的非零性保证其行列式不是零多项式。

1. 求在 \(\delta=0\) 处的秩 \(k\)。若等于当前大小 \(n\)，本阶段结束；若为 0，把整个矩阵除以 \(\delta\)，记行列式多了 \(n\) 次差值，再重复。
2. 当 \(0<k<n\) 时，通过行列置换找到 \(A(0)\) 可逆的 \(k\times k\) 块：
   \[
   M(\delta)=\begin{pmatrix}A&B\\C&D\end{pmatrix}.
   \]
   在形式幂级数环 \(\mathbb Q(r)[[\delta]]\) 中定义
   \[
   S=D-CA^{-1}B,\qquad M_{\mathrm{next}}=S/\delta.
   \]
   因为 \(S(0)=0\)，除法合法，并且
   \[
   \det M=\pm\det A\;\delta^{n-k}\det M_{\mathrm{next}}.
   \]
   记录额外的 \(n-k\) 次差值，随后处理更小的矩阵。

逆矩阵的系数无需猜测。若 \(A=\sum_{j\ge0}A_j\delta^j\)、
\(A^{-1}=\sum_{j\ge0}Y_j\delta^j\)，有
\[
Y_0=A_0^{-1},\qquad
Y_m=-A_0^{-1}\sum_{j=1}^{m}A_jY_{m-j}.
\tag{27}
\]

这给出所要求的显式递推。只需有限精度：设
\[
D_\delta=\sum_{\text{special }j}(U-b)
+\sum_{\text{other }j}\alpha_2.
\]
则 \(\deg_\delta\det\mathcal R_H\le D_\delta\)，额外消失阶至多是它。实现时若当前累计阶数为 \(e\)，保留到 \(D_\delta-e\) 阶即足够，不可在每次除以 \(\delta\) 后忘记降低剩余精度。

中间主元运算可能引入关于 \(r\) 的有理函数分母；不能把它们当成新的免费算术收益。若最终额外阶数为 \(e\)，真正的多项式
\[
Q(r,\delta)=\det\mathcal R_H(r,\delta)/\delta^e
\]
仍然继承 (19)—(21) 的分母与系数界，因为除以单项式只移动系数位置，不放大系数。这是对最终对象的控制，不声称中间表达式永不膨胀。

这里计算的是 \(\mathbb Q(r)\) 上的最大差值阶。在特殊的代数中心值代入以后，消失阶可能进一步上升；不能未经检查就把形式变量计算的主元用于所有数值中心。

## 9. 精确例子、代码与核对范围

[verify.py](verify.py) 使用 Python 标准库的有理数与高斯有理数运算。[results.json](results.json) 保存实际运行结果。

[verify_saturation.py](verify_saturation.py) 实现第 8 节的秩与 Schur 补递推；[saturation_results.json](saturation_results.json) 保存八个有限有理多项式矩阵的独立核对和完整步骤。

运行：

    python3 verify.py --output results.json
    python3 verify_saturation.py --output saturation_results.json

最小展示取 \(U=7,\ A=B=1\)。一组二项式系数矩阵及变换为
\[
\begin{pmatrix}8&36\\28&84\end{pmatrix},
\qquad
\begin{pmatrix}1&-9/2\\0&1\end{pmatrix}.
\]
变换后主元为 \(8,-42\)。两组原来的差值次数都是 \(1,2\)，共保证 6 次；约化后的次数是 \(1,3,1,3\)，共 8 次。一个允许的局部四行子式精确等于
\[
112896\,\delta^8.
\]
因此这个局部四列块的最大公差值阶恰为 8。该例只核对局部行条件，不声称它验证了原论文的全局权重分离阈值。

脚本从原二项式矩阵独立执行有理数列消元，再与 (2)—(4) 比较；同时核对小块的精确行列式、分母最小公倍数、清分母后高斯整数整除以及实际截断多项式的系数向量。具体范围和计数以结果文件为准。一般证明是第 2—7 节，不能由有限检查代替。

实际行系数包的检查采用 \(K=2\)、\(v_0=1\)、\(w_1=2\)、\(w_2=4\)、\(\theta=1/2\)、\(H=160,h=128\)，因此 \(U=7,A=B=1\)。两个截断长度分别为 10 和 20，中心为 \(\mathrm i,4\mathrm i/3\)。共有 560 个活动行、四列；直接稀疏展开与分解式逐行比较。保留 \(G_2-G_1\) 改变了 311 个向量系数，因此检查没有把不同截断的影响默认为零。此有限参数不宣称满足原稿全局权重分离条件。

递推脚本独立用逐排列多项式行列式作参考。一个由共享向量公式构造的四阶例子满足
\(\det M_{\mathrm{original}}=\delta^8\det M_{\mathrm{reduced}}\)，但
\(\det M_{\mathrm{reduced}}=-9\delta^2\)；递推正确识别出还需约去两次。这是说明“保证阶数未必最大”的代数模型，不是原稿完整行矩阵的新增证书。

## 10. 本次推进与仍缺的信息

本次新增的是：更强的成组差值整除、保持行列式的闭式列变换、剩余列的显式递推、中心二项式分母界、次数与系数界，以及与原算术清分母的精确对应。

还没有得到的是：约去这些结构后，剩余行列式在精确周期 \(2\pi\mathrm i\) 附近具有足够小的解析上界。当前系数界控制增长，并不证明衰减。增加的差值因子在算术侧也增加，因而没有额外的净逼近收益。

若要沿旧笔记的反证方向完成端点，需要在仅假设
\(\kappa_i=q_i^2|\pi-p_i/q_i|\to0\)、不预设下降速度时，使剩余非零算术对象的解析上界小于其算术下界。那将证明 \(\pi\) 是坏逼近数、证伪 P2848。本文没有完成这一步，也没有证实该方向最终可行。

## 来源

- [仓库旧构造与算术整除，§§19.1—19.7](../OPENAI_PI_EXPONENT_P2848_2026-10-27.md)。
- [自身底数模型与论文迁移范围](../self-base-theta-2026-10-31/OPENAI_SELF_BASE_TRANSFER.md)。
- OpenAI 原论文 [行列式部分](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/determinant.tex) 与 [插值部分](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/interpolation.tex)，固定版本 4eaee87548250b2682affda7d05d947c05554bae。

本文的显式列恒等式及后续推导不依赖外部论文主结论的真实性；使用其完整非零子式时，原插值条件仍作为明确前提。
