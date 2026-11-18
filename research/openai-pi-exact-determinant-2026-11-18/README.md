# OpenAI 的 π 行列式在指数 2 处的精细比较

2026-11-18。延续[端点联合评分](../endpoint-joint-bound-2026-10-31/)与[实际子式选列](../contradiction-rank-2026-11-01/)。目标是在假定存在 $q_i^2|\pi-p_i/q_i|\to0$ 的分数时，沿 OpenAI 的插值证明构造矛盾。得到同一非零子式的更紧算术下界与解析上界，尚未证明两界必然交叉。原论文的正式结论是指数任意大于 2 的下界，其参数选择依赖这个严格不等式。[论文](https://github.com/openai/math/blob/main/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/main.tex) · [原行列式一节](https://github.com/openai/math/blob/main/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/determinant.tex) · [Lean 陈述范围](https://github.com/openai/math/blob/main/lean/docs/017.md)

## 1. 一个子式可以用两种对数展开计算

采用论文的权重，列为 $P=Y^hX^\alpha$，满足 $w_0h+w\cdot\alpha\le H$。行是 $(j,s,\beta)$，其中 $0\le j<K$，$v_0s+w\cdot\beta/\theta<H$。取同一组使全行子式 $\Delta$ 非零的列。记 $r_i=2ip_i/q_i$ 和 $\omega=2\pi i$。

算术侧的行是提取 $P(1+t,jr_i+G_i(t)+u_i)$ 的 $[t^su^\beta]$，可取
$$
 N_i=\left\lfloor\frac{w_i}{\theta v_0}\right\rfloor,qquad
 G_i(t)=\sum_{k=1}^{N_i}(-1)^{k+1}t^k/k.
$$
解析侧可用完整形式级数 $\log(1+t)$。原论文也在证明插值秩不变时使用保持过滤的代换。这里进一步观察到**每个固定选列的全行子式的数值完全相同**。代换 $u_i\mapsto u_i+G_i-\log(1+t)$ 的非零替换项具有严格更高的权重。把各行按权重排序，每一行系数等于原来同一行加若干严格低权重行的线性组合。原矩阵保留了所有这些行，且行变换对角全是 1。因此该结论甚至可在固定高度 $H$ 下成立，且无需比较两个不同的选列。

这样在解析侧不再有任何对数截断余项。位移标量恰为
$$
 \gamma_{\rho,a}=\binom a\beta\prod_i[j(r_i-\omega)]^{a_i-\beta_i},\qquad a\ge\beta.
$$
若 $\kappa_i=q_i^2|\pi-p_i/q_i|$，每次实际使用位移都保留因子 $\kappa_i^{a_i-\beta_i}$，无需引入常数 $F_0$ 或将其降格为一个固定的正额外指数。

## 2. 同一子式的清分母和额外整数因子

令 $M$ 为行数，$S=\sum_{\rm rows}s$，$C_i=\sum_{\rm cols}\alpha_i-\sum_{\rm rows}\beta_i$。后者均非负，因为非零行列式至少有一个非零的排列乘积，每对行列的指标必须满足 $\alpha\ge\beta$。依原论文乘列因子 $\prod_iq_i^{\alpha_i}$ 和行因子 $\prod_i q_i^{-\beta_i}$，得到
$$
 B=\Delta\prod_iq_i^{C_i}.
$$
对普通素数 $\ell$ 定义
$$
 d_{\ell i,k}=(v_\ell(k)-v_\ell(q_i))_+,\qquad
 a_{\ell i}=\max_{1\le k\le N_i}d_{\ell i,k},\qquad
 b_\ell=\max_{i,1\le k\le N_i}d_{\ell i,k}/k ,
$$
空最大值取 0。可取清分母整数
$$
 D=\prod_\ell\ell^{\min\{\sum_i a_{\ell i}C_i,\lfloor b_\ell S\rfloor\}},
 \qquad DB\in\mathbb Z[i].
$$
证明直接展开同一行列式。每个非零排列项在坐标 $i$ 只有 $C_i$ 个因子可贡献对数系数，每个这样的系数为 $\pm q_i/k$，因此最坏赋值成本为 $a_{\ell i}C_i$。另一方面这些因子共占用的 $t$ 次数至多 $S$，而每次 $t$ 次数的成本至多 $b_\ell$。两界取小。普通整数清分母足够处理高斯有理数。即使丢掉第二个界及 $q_i$ 的素因子回扣，也有
$$
 \log D\le(4\log2)\sum_i C_i N_i
 \le\frac{4\log2}{\theta v_0}\sum_i w_i C_i.
$$
这项成本按实际选列与行的横向次数差付，而无需像原粗界按每个坐标的最大可能列次数分别付。

还有一个与上述 $D$ 兼容的整数因子。令中心 $j$ 暂时成为独立变量 $z_j$。$DB$ 仍为高斯整系数多项式。取最大权重的最后一个坐标作主元，且其中心方向分量非零。把该坐标的偏导换成沿向量 $(r_i)_i$ 的导数。对每个中心和固定 $s$，这个变换只会用同权重或更便宜的偏导，故在完整行包上可逆。固定 $s$ 与其余偏导指标 $\beta'$ 后，每个节点具有从 0 开始的 $J=J_{s,\beta'}$ 个连续切向导数。将节点 $z_k$ 的这些行在 $z_j$ 平移展开，若要避开 $z_j$ 已有的 $J$ 个导数行，至少要取 $J$ 个互异且各不小于 $J$ 的导数阶。总共至少付出 $J^2$ 次 $z_k-z_j$。外积逐块相乘，得到
$$
 V_K^{\sum_{s,\beta'}J_{s,\beta'}^2}\mid DB,
 \qquad V_K=\prod_{d=1}^{K-1}d!.
$$
多变量中每一对差是不同的不可约因子，乘积首项系数为单位。因此先在高斯有理多项式环得到的除法，在原高斯整数多项式环仍为整除；最后令 $z_j=j$。行变换的比率为非零常数，不影响对节点差的消失阶。对 $K=1,2$ 这个因子等于 1。

## 3. 解析侧不丢掉一半重复度数收益

固定 $0<r<1$，用半径 $r$ 对局部变量 $t$ 作 Cauchy 估计。令 $J_0=2\pi(K-1)-\log(1-r)$，取 $R>\max(1,J_0)$，$b=J_0/R\in(0,1)$。这使 $|j\omega+\log(1+t)|\le J_0$。对一列 $P=Y^hX^\alpha$，整函数 $f_{a,P}(z)=[u^a]P(e^z,z+u)$ 在 $|z|\le R$ 的上界可取 $e^{Rh}(2R)^{|\alpha|}$。展开每行时仍用论文的观察，即同一横向指标 $a$ 的两条 Taylor 行若选了相同度数，整个项为零。

但无需将几何衰减劈成两半。对有 $n$ 行的一个 $a$ 组，所有互异度数的几何权恰为
$$
 \sum_{d_1,\ldots,d_n\text{ pairwise distinct}} b^{\sum_j d_j}
 =\frac{n!\,b^{n(n-1)/2}}{\prod_{j=1}^n(1-b^j)}.
$$
先排序度数，再用相邻度数减去最小间距，得到独立的几何级数。于是保留整份 $\log(1/b)\binom n2$ 的收益。每组的分母可统一按 $P_b=\prod_{j\ge1}(1-b^j)^{-1}$ 控制，这个积收敛。同一组的 Taylor 系数行列式使用 Leibniz 公式，列的增长成本以实际 $R\sum h$ 计，而不是每列都取 $RH/w_0$。

具体地，选一个非零标量 $\gamma_{\rho,a_\rho}$ 给每行，记 $n_a$ 为该横向指标的行数。设
$$
 \mathcal S=\min_{\text{这些选择}}\left(
    \log(1/b)\sum_a\binom{n_a}{2}
    -\sum_\rho\log|\gamma_{\rho,a_\rho}|\right).
$$
也可先精确剔除其整块行列式为零的选择，再计算该最小值。由上面界、$\prod_a n_a!\le M!$ 及每行至多 $(\lfloor H\rfloor+1)^m$ 个 $a$，得到
$$
 \log|\Delta|\le
 R\sum h+\log(2R)\sum|\alpha|+S\log(1/r)
 +2\log(M!)+M\log P_b+mM\log(\lfloor H\rfloor+1)-\mathcal S .
$$
此式与第 2 节比较的是同一个非零 $\Delta$。对固定参数令 $H$ 增长，最后三个显式组合成本是 $o(MH)$。新的精确比较仍须在指数 2 处证明 $\mathcal S$ 有足够大的统一下界。如果允许 $a_\rho=\beta_\rho$ 的零位移选择进入最小值，$\kappa_i$ 在该项中完全不出现。能不能把足够多这样的项从真正的最低成本中排除，仍需关于实际子式的秩或联合约因子的额外定理。单改清分母和重复度数上界没有给出这个定理。

## 4. 精确复核

运行 `python3 research/openai-pi-exact-determinant-2026-11-18/verify.py`。有限测试取实有理中心，检验的代换、整除和几何级数公式同时适用于复有理中心。每行的截断按真实权重选择。程序寻找非零全行子式，再用同一批列重算完整对数，四个子式的有理值逐一精确相等。旧清分母上界按原论文的合法长截断计算。表中位数表示用于保证整性的**上界**，与实际既约分母分别列出。

| 行数 | 实际的 $B$ 的分母 | 新清分母上界位数 | 原清分母上界位数 |
|---:|---:|---:|---:|
| 26 | 1215 | 37 | 869 |
| 20 | 18225 | 38 | 881 |
| 26 | 5 | 43 | 869 |
| 18 | 1 | 26 | 783 |

最后一例用三个中心，另核验清分母所得整数被 $2^6$ 整除。没有运行 Lean，也没有把这些有限实例当成所有分母的结论。
