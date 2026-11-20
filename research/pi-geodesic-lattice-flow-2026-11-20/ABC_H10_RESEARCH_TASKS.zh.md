# abc 与希尔伯特第十问：双见证的两份新任务

2026-11-20。接续[素数约分与 CRT 续研](DOUBLE_WITNESS_CRT_FOLLOWUP.zh.md)。

## 读到的论文

- 望月新一的 IUT 系列于 2021 年发表，主张得到 abc。Scholze–Stix 对证明提出具体异议。本任务把整数版 abc 列为条件假设。Granville 的 *ABC Allows Us to Count Squarefrees* 是从 radical 控制素数幂贡献的参考。 [IUT 专刊](https://ems.press/journals/prims/issues/1507) · [Scholze–Stix](https://www.math.uni-bonn.de/people/scholze/WhyABCisStillaConjecture.pdf) · [Granville](https://doi.org/10.1155/S1073792898000592)
- 原版希尔伯特第十问在整数上由 Davis–Putnam–Robinson–Matiyasevich 得到否定答案。有理数版方面，Bhargava–Burungale 在 11 月 19 日更新的预印本给出整数在有理数域的正存在定义。他们的证明采用短实区间和分母界，并避开 OpenAI 预印本使用的两项技术结果。这里取其中可直接用于本项目的短区间方法。 [Davis–Putnam–Robinson](https://www.jstor.org/stable/1970289) · [Matiyasevich](https://www.mathnet.ru/eng/im1910) · [Bhargava–Burungale v2](https://arxiv.org/abs/2611.16130v2)
- 两条线确有直接交点。Pasten–Vidaux 的 11 月 18 日预印本在数域版 abc 假设下给出二次数域上的第十问否定结论。第 2 节研究奇数正估值的素数权重。Kozachinskiy–Pasten–Vidaux 的有理数版条件论证另依赖 Dickson 猜想及 \(2\)-主 Tate–Shafarevich 群有限性。 [Pasten–Vidaux §2](https://arxiv.org/html/2611.24748v1) · [有理数版三条件论证](https://arxiv.org/abs/2611.17124)

## 两队共用的新起点

沿用 \(z_0=z_{n,0}\)、\(z_1=z_{n,1}\)，令 \(\delta=z_1-z_0=\sigma 2^{e_n}/V_n\)，其中 \(\sigma\in\{-1,1\}\)，\(e_n=2+v_2(n!)=n+2-s_2(n)\)，\(s_2(n)\) 是 \(n\) 的二进制数字和。由 Beta 积分和二进估值得
\[
V_n=\frac{(4n+1)\binom{4n}{2n}}{2^{s_2(n)}}.
\]

上一版只处理简单极点的 CRT 类有严格障碍。取充分大的奇素数 \(p\) 满足 \((4n+1)/3<p^2\le2n-1\)。在上一轮的 \(b_m=\sum_{\ell=1}^m(-1)^{m-\ell}/(2\ell-1)\) 展开中，对于 \(n\le m\le 2n+1\)，唯一能贡献 \(p^{-2}\) 的项来自分母 \(p^2\)。没有 \(3p^2\) 或 \(p^3\) 项。乘以 \(p^2\) 后模 \(p\)，这项的符号与积分展开的交替符号抵消，从而 \(v_p(z_0)=v_p(z_1)=-2\)。另一方面，中心二项式在 \(p^2\) 层的进位数是 \(2-2=0\)，在 \(p\) 层至多是 \(1\)。若 \(p\mid4n+1\)，则 \(p\) 层也没有进位，而 \(v_p(4n+1)=1\)。故 \(v_p(V_n)\le1\)。

简单极点的系数 \(S\) 因此不含这些 \(p\)。每条 \(z_0+(t/S)\delta\) 都保留 \(v_p=-2\)。同一 CRT 类内两点的差是整数倍 \(1/V_n\)，故它们的非零行列式必被 \(p^3\) 整除。素数定理给出统一下界
\[
|D|\ge \exp\!\left(\bigl[3(\sqrt{2}-\sqrt{4/3})+o(1)\bigr]\sqrt n\right).
\]
有限复算在 \(2\le n\le96\) 核对了 60 个这样的 \((n,p)\)。

把超额素数幂一次处理干净，令 \(V_nz_0=a/P\) 为最简分数，\(S=2^{e_n}P\)。则
\[
V_n \left(z_0+\frac{t}{S}\delta\right)=\frac{a+\sigma t}{P}.
\]
因此选 \(t\equiv-\sigma a\pmod P\) 后，混合有理数恰好遍历 \(k/V_n\)（\(k\in\mathbb Z\)）。对给定的 \(k\)，\(t=\sigma(Pk-a)\)。这条公式精确给出所有超额素数幂的约分。有限复算验证了该公式与上一版积分斜率的连接。此时必须控制 \(k\) 和 \(V_n\) 的最大公因子。

## 给 AI A：abc、素数根积和约分成本

请固定上面的全素数幂 CRT 构造，研究哪些约分策略可能同时满足实误差和有界行列式。

1. 复核经典整数 abc 的使用条件。对 \(0<t<S\) 的奇数 \(t\)，设 \(w=3t-2S\)。因为 \(t\equiv-\sigma a\bmod P\)，所以 \(\gcd(t,S)=1\)，\(\gcd(3t,2S)\le3\)。请把 \(w+2S=3t\) 按符号化成原始三元组，并验证在 abc 假设下
   \[
   \operatorname{rad}(tw)\gg_\varepsilon S^{1/(1+\varepsilon)}/\operatorname{rad}(S),
   \]
   其中 \(\operatorname{rad}(m)\) 只乘 \(m\) 的不同素因子。进一步，证明对任意固定 \(\eta>0\)，充分大的 \(n\) 有 \(\log\operatorname{rad}(tw)\ge(\log2-\eta)n\)。

2. 以 \(t=u x^r,\ w=v y^r,\ \log\operatorname{rad}(uv)=o(n)\) 的幂次约分策略作为校准。找出足以排除的最小 \(r\)，并区分该策略和一般的 \(k\) 约分。从上界入手可用 \(P=\operatorname{den}(V_nz_0)\)、\(\log V_n\sim4n\log2\)。\(p>\sqrt{4n+1}\) 时，\(v_p(V_n)\le1\)，并且 \(p\mid P\) 蕴含 \(p\nmid V_n\)，因为 \(z_0\) 的分母整除不超过 \(4n-1\) 的奇数的最小公倍数。因此
   \[
   \log P\le(4-4\log2+o(1))n,\quad \log S\le(4-3\log2+o(1))n.
   \]
   这已经给出 \(r=6\) 足以排除。要继续推进，请从真实的两行约分条件推出一个独立的 radical 上界，或者给出明确反例，说明有界行列式为何没有带来这种上界。

3. 对照两篇方法。Pasten–Vidaux §2 的奇估值引理适用于固定平方自由多项式，若想把实误差窗口中的多个整数线性式合写入一个多项式，请证明固定性及高度条件。直接对行列式式子 \(p_0q_1+|D|=p_1q_0\) 用 abc，遇到 \(V_n/\operatorname{rad}(V_n)=\exp(O(\sqrt n\log n))\) 的障碍。在 \(|D|\) 有界、间距 \(1/V_n\) 时，原始化后的 radical 已至少是 \(\operatorname{rad}(V_n)/O(1)\)，所以要找到独立的加法恒等式。

交付一个无条件的素数幂结论、一个注明 abc 假设的推论、以及二者之间的明确检验。若条件法卡住，指出是哪个 radical 上界缺失。有限数据可用于选择模型。

## 给 AI B：希尔伯特第十问式整数证书与精确搜索

请借用 Bhargava–Burungale §6 的短区间与分母界思路，给 \(k/V_n\) 网格做一个可复核的整数存在性判据，并查明哪些尺度能提供两条独立方向。

1. 令 \(k_1<k_2\)、\(d=k_2-k_1\)，\(g_i=\gcd(k_i,V_n)\)，约分为 \(p_i/q_i\)，\(q_i=V_n/g_i\)。给定正有理数 \(B,K,X\) 及最大整数间距 \(d_{\max}\)，用整数多项式等式编码 \(1\le d\le d_{\max}\)、\(q_i\le BX\)、\(X|q_i\pi-p_i|\le K\) 的有理端点证书。用 \(k_i=g_ip_i\)、\(V_n=g_iq_i\) 和 Bézout 等式保证真实约分，分母正性和整数不等式可转为四平方等式。用外区间作完备枚举范围，用内区间保证成功，用严格误差下界给出失败证书。遇到端点不确定，应提高 Machin 或有理矩精度。
2. 必须利用短区间。由 \(|D|=d q_1q_2/V_n\ge1\) 和 \(q_i\le BX\)，可推出 \(|k_i-\pi V_n|\le BKd\)。固定 \(B,K,d_{\max}\) 后，可能的整数 \(k_i\) 个数与 \(n\) 无关。对 \(d=1\)，\(g_1,g_2\) 互素且 \(D=V_n/(g_1g_2)\le 2BK\)。把 \(V_n/D\) 的素数幂分给 \(g_1,g_2\)，再检查 \(k_1\equiv0\bmod g_1,\ k_1\equiv-1\bmod g_2\) 的解能否命中短区间。请证明筛选的充分性和必要性，并尝试给出比因子枚举更直接的区间判据。
3. 基准是 \(2\le n\le96,\ X=\lceil\sqrt{V_n}\rceil,\ d_{\max}=4\)。已用有理区间独立复核的结果是 \(B=K=1,2\) 无配对，\(B=K=4,8\) 只有 \(n=2\) 的 \((987,990)\)，\(B=K=16\) 另有 \(n=2\) 的 \((990,994)\) 和 \(n=3\) 的 \((9436,9438)\)。请独立复现并推进到一个无限子族的构造或排除。优先检查倍增尺度 \(X_n=4^n\)。当 \(n=2^r\)，\(V_n\asymp\sqrt n\,16^n\)，给定有界 \(d\) 时，\(q_1q_2\ge V_n/d\) 已排除 \(q_i\le B4^n\)。如果换尺度，每组候选的可用区间为 \([H/B,K/E]\)，要检查这些区间的并是否覆盖所有充分大的尺度。

交付一个带参数的精确验证器，给出等式证书与有限枚举完备性证明，并附一条无限尺度结论或明确断点。所有主结论在整数域完成；有理域的存在定义可作为单独的形式化延伸。
