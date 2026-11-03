# 文献与项目接口：π 逼近、共同约分和进制复杂度

检索与整理日期：2026-11-03。以作者论文、arXiv 或期刊原始页面为来源；本文区分已发表文章、预印本作者声称、本项目推论和待检验想法。本轮没有逐篇独立审核完整证明。

## 1. 首选：整数系数整体约分

Zeilberger–Zudilin，[The Irrationality Measure of Pi is at most 7.103205334137...](https://arxiv.org/abs/1912.06345v2)，2020，v2；[作者原文 PDF](https://sites.math.rutgers.edu/~zeilberg/mamarim/mamarimPDF/pimeas.pdf)。论文用特定复积分、部分分式和素数集改进 π 的无理性测度上界。这个数字是当时的历史上界，不用于断言今天的最强结论。

实用细节：论文 Lemma 2 针对其完整系数 \(A_j\) 给出整除性，随后 Lemma 3 用素数乘积 \(\Phi_n\) 改善整数化成本。v2 明确修正了旧版本一个更强的错误说法：不能把对合并系数成立的整除性，直接套给每一个求和项。

**与项目的接口是建议，不是现成证明：**对双 arctan 的完整合并分子先做整数化，再研究整体系数和的模素数抵消；保持两角已约分分母及合并后真正 gcd 的账。论文的特殊素数窗口和 0.645... 常数不能直接复制到本项目的不同核。

## 2. 同时 Padé 与共同分母

- Tanguy Rivoal，[Simultaneous Padé approximants to the Euler, exponential and logarithmic functions](https://numdam.org/articles/10.5802/jtnb.914/)，2015：构造共享多项式的显式同时 Padé 近似。其 Theorem 5 中共享多项式用于同一变量的不同函数；它没有直接给出本项目两个不同 arctan 数值的共同既约整数分母。将 arctan 转成复对数后，还要核对分支、两个特化点和真实清分母成本。
- Kawashima–Poëls，[Padé approximation for a class of hypergeometric functions and parametric geometry of numbers](https://arxiv.org/abs/2202.10782v2)，2022 预印、2023 期刊：给出一类超几何函数的显式 Padé 构造、有效 Poincaré–Perron 渐近及同时逼近框架。可能帮助统一控制递推余项与系数增长；本项目参数满足其具体定理全部假设与否尚未核验。
- Bernard de Mathan，[Simultaneous Diophantine approximation with a divisibility condition](https://www.numdam.org/articles/10.5802/jtnb.1303/)，2024 卷、2025 在线：研究带整除条件的二维同时逼近，背景与 Littlewood 型问题相关。它可作为公共整除约束的对照，并不直接证明当前 π 构造的共同分母收益。

## 3. 给现有相消结论一个多项式速率

Matveev，[An explicit lower bound for a homogeneous rational linear form in the logarithms of algebraic numbers II](https://www.mathnet.ru/eng/im314)，2000。对两个固定正代数数的非零整数对数线性形式，可使用其有效下界；本段是把这一经典工具应用于项目构造的推导，不是论文中的 π 结论。

设 \(a=2/5\)、\(b=3/7\)，按[既有推导](../coupled-arctan-2026-11-03/README.md)定义正代数数
\[
\rho_c=\frac{\sqrt{1+c^2}-1}{\sqrt{1+c^2}+1},\qquad
\gamma=\frac{\log\rho_a}{\log\rho_b}.
\]
两个代数数乘法独立：若非零整数幂相等，则相等值属于 \(\mathbb Q(\sqrt{29})\cap\mathbb Q(\sqrt{58})=\mathbb Q\)；其中一个二次域共轭把该数变成倒数，迫使其为 1，矛盾。因此 Matveev 给出有效但本轮未算的常数 \(c>0,C\ge1\)，使
\[
\|q\gamma\|\ge c q^{-C}\quad(q\ge1).
\]

记 \(q_j\le L<q_{j+1}\) 为 \(\gamma\) 的收敛分母。由 \(\|q_j\gamma\|<1/q_{j+1}\) 得 \(q_j\ge(cL)^{1/C}\)。任意连续 \(q_j\) 个旋转点距等间距的有理网至多 \(1/q_{j+1}\)，所以是 \(O(1/q_j)\) 网；平移起点不影响该界。

在此前两个 29 进子窗口的相位兼容区域内，先留出固定正比例边界缓冲，偶数 n 对应的连续索引块长度 \(L\asymp N=29^k\)。相位目标含固定平移、m 为奇数；上述网格仍适用。此前的两项相对误差渐近各有 \(O(1/N)\) 余项，故存在有效的 \(\delta>0\)，使同一列候选满足
\[
\boxed{\frac{|e_n(a)+e_m(b)|}{|e_n(a)|}=O(N^{-\delta}),
\qquad v_{29}(q)\le k-2.}
\]
可取任意不大于 \(\min(1,1/C)\) 的上述指数；隐藏常数与 C 可能极大。本轮未明确计算常数、未用 Lean 形式化，也未将该速率与完整既约分母增长闭合。它只控制相对相消，**不推出 \(q^2|\pi-p/q|\to0\)**。

## 4. 指数 2、π 与进制：量词要分开

OpenAI 的 [2026 预印本](https://github.com/openai/math/blob/main/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/main.tex)主张 \(\mu(\pi)=2\)，[对应 Lean 范围说明](https://github.com/openai/math/blob/main/lean/docs/017.md)说明其选择的形式化结论。这里按作者声明引用，不将本轮检索称为对论文全部证明、依赖或构建的独立审计。即使成立，其 \(\nu>2\) 的下界阈值依赖 \(\nu\)，也不自动判定 \(\inf q^2|\pi-p/q|\) 是否为零。

一个可以自证的量词对照：\(x=[0;1,2,3,\ldots]\) 的下一部分商 \(a_{n+1}=n+1\) 无界，因此沿收敛分数 \(q_n^2|x-p_n/q_n|<1/(n+1)\to0\)。同时 \(q_n\ge n!\)，由标准连分数指数公式得 \(\mu(x)=2\)。黄金比例倒数也有指数 2，但它的部分商全为 1，是坏逼近数。两例都不是 π。

Bugeaud–Kim，[On the b-ary expansion of a real number whose irrationality exponent is close to 2](https://arxiv.org/abs/2510.02059v2)，定理 1.4：对整数底 \(b\ge2\)，若无理数 \(\xi\) 有 \(\mu(\xi)=2\)，令 \(p(n,\xi,b)\) 为其底 b 小数中长度 n 的不同连续数位块个数，则
\[
\liminf_{n\to\infty}\frac{p(n,\xi,b)}n\ge\frac87.
\]
**条件性接口：**若前述 π 指数预印本成立，就可以把该已发表定理用于 π 的每个整数底。这不证明 π 正规，不证明特殊长零段出现，也不直接适用于实数底或 π 自身底，更不解决爆破解。

对用户提出的“刘维尔数在某个底的图线很直”，[Nandakumar–Vangapelli](https://arxiv.org/abs/1204.4104v2) 对任意指定整数底构造了该底正规、同时是刘维尔数的例子。论文另一种有限多底数论构造有广义 Artin 猜想假设；这项对照不应被误述成该文无条件证明全底正规刘维尔数。

## 5. 检索状态

[Chudnovsky 2026 预印本](https://arxiv.org/abs/2610.28079v1) §12 对匹配共轭嵌入的对数比提出更强的作者声称；当前双 arctan 衰减率比的适用假设、定理及形式化依赖尚未逐项核验，故不用于第 3 节结论。

Bai 的 [π 上界预印本](https://arxiv.org/abs/2609.11276v2)是 2026-09-11 的历史结果和局部优化声明，不应与后续 OpenAI 预印本声明混同。[arctan(1/2) 预印本](https://arxiv.org/abs/2608.19782) 的原页面本轮无法打开；仅保留待核验线索，不转述其未经核验的技术细节。

研究的下一关仍是：在相消阶数上控制两角真正既约分母的共有素数与剩余成本。Matveev 的多项式相消不自动胜过指数级分母增长；没有得到新 π 有理逼近纪录。
