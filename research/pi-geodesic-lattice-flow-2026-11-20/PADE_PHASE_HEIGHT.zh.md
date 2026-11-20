# 对角 Padé 相位的高斯约分与临界高度

本轮接续两份上传的 Padé 压缩包。两份的 SHA‑256 都是 b409e64fe4ccd265adb3c027fc9698c12b1bf85f11c2940f1702d43e6dd8cb0b。所附脚本重新生成的 CSV 和文本逐字节匹配，共 24 条。本文证明了观察到的 $1+i$ 规律和坐标分母公式，并给出 $n$ 可变时的约分上界。它排除对角 Padé 相位构造的一个截断范围。$\pi$ 的坏逼近性仍未判定。

## 约分的统一公式

取既约 $p/q>0$。写
$$
P_n(z)=\sum_{k=0}^n\binom nk\frac{(2n-k)!}{(2n)!}z^k,\quad
H_n=(2n-1)!!,\quad
C_n=H_n q^n P_n(2ip/q)=a_n+ib_n.
$$
令 $H_0=1$。原报告清分母所得的 $A$ 是 $2^n C_n$。整系数展开和递推为
$$
C_n=\sum_{k=0}^n\binom{2n-k}{k}(2n-2k-1)!!\,q^{n-k}(ip)^k,\qquad
C_n=(2n-1)qC_{n-1}-p^2C_{n-2},
$$
其中 $C_0=1,C_1=q+ip$，$(-1)!!=1$。两式用系数公式直接核对。这个递推也来自 Lambert 正切连分数的分子分母递推。对相邻向量取整数行列式，得到
$$
b_n a_{n-1}-a_n b_{n-1}=p^{2n-1}.
$$
于是 $c_n=\gcd(a_n,b_n)$ 整除 $p^{2n-1}$。这里的 $c_n$ 始终为奇数。若 $p$ 为奇数，由行列式可知。若 $p$ 为偶数，$q$ 必为奇数，而递推给出 $(a_n,b_n)\equiv(1,0)\pmod2$。

令 $g_n=\gcd_{\mathbb Z[i]}(C_n,\overline C_n)$。在除以 $c_n$ 后，余下两个互素整数坐标分别记作 $a',b'$。高斯公因子整除 $2a'$ 和 $2ib'$，故只有 $1+i$ 可能留下。它恰好在两坐标同为奇数时出现，且 $(1+i)^2$ 不再整除。若 $p,q$ 都为奇数，模 2 的递推循环是 $1,1+i,i$。若 $q$ 偶而 $p$ 奇，模 2 在 $1,i$ 之间交替。结合 $p$ 偶的情况，得到
$$
g_n\sim c_n(1+i)^{\varepsilon_n},\qquad
\varepsilon_n=1\ \Longleftrightarrow\ p,q\ \text{都为奇数且}\ n\equiv1\pmod 3.
$$
符号 $\sim$ 表示相差一个高斯单位。

约分后的相位候选
$$
W_n=\frac{P_n(2ip/q)}{P_n(-2ip/q)}-1
    =\frac{(C_n-\overline C_n)/g_n}{\overline C_n/g_n}
$$
的高斯分母记为 $B_n=\overline C_n/g_n$。上面的坐标判据给出 $\gcd_{\mathbb Z[i]}(B_n,\overline B_n)=1$。同时为 $W_n$ 的实部和虚部清分母的最小正整数，是满足 $B_n\mid m$ 的最小正整数 $m$。因 $m=\bar m$，这等价于 $B_n\overline B_n\mid m$，故该最小整数恒为 $|B_n|^2$。

## 奇素数的消因子与统一上界

设 $\ell$ 是 $p$ 的奇素因子，令 $\nu=v_\ell(p)$、$h=\lfloor\log_\ell(2n-1)\rfloor$。取 $C_n/(H_n q^n)$ 的展开。第 $k\ge1$ 项除去 $(ip/q)^k$ 后的系数是
$$
\frac{\binom{2n-k}{k}}{(2n-1)(2n-3)\cdots(2n-2k+1)} .
$$
若 $\nu>h$，各非恒定项的赋值都严格为正，因而 $v_\ell(c_n)=v_\ell(H_n)$。若 $\ell^2>2n-1$，分母中每项至多带一份 $\ell$，而连续 $k$ 个奇数中至多有 $\lceil k/\ell\rceil$ 个 $\ell$ 的倍数。对 $k\ge2$ 这严格小于 $k$。$k=1$ 的比系数恰等于 1，因此也有 $v_\ell(c_n)=v_\ell(H_n)$。

剩余素数满足 $\ell\le\sqrt{2n-1}$ 且 $\nu\le h$。整数行列式给出 $v_\ell(c_n)\le(2n-1)\nu\le(2n-1)h$。结合初等 Chebyshev 素数计数上界，
$$
\log c_n
\le \log H_n+(2n-1)\,\pi_{\mathrm{pr}}(\sqrt{2n-1})\,\log(2n-1)
=O(n^{3/2}),
$$
其中 $\pi_{\mathrm{pr}}(x)$ 计不超过 $x$ 的素数个数。这里的常数与 $p,q$ 无关。也可由 $\sum_{\ell^j\le x}\log\ell=O(x)$ 得到素数计数上界。将素数按 $\ell\le\sqrt{x}$ 和 $\ell>\sqrt{x}$ 分组即可，前一组至多 $\sqrt{x}$ 项，后一组各贡献至少 $(\log x)/2$。

小素数同阶相消确实可以超过双阶乘赋值。例如 $(n,p,q)=(5,3,1)$ 给出 $C_5=-1620+243i$，故 $c_5=81$、$H_5=945$，两者的 $3$ 进赋值分别为 4 和 3。本例排除了把 $v_\ell(H_n)$ 当成普遍上界的做法。

## 相位误差与高度的共同尺度

记 $z=2ip/q$、$V=e^z-1$，假设 $p/q\to\pi$。对固定 $n\ge1$，$W_n$ 的分子对应一个非零有理系数多项式，并在超越点 $z=2i\pi$ 取非零值。分母同样非零。因此任何固定 $n$ 都不能保持 $W_n$ 相对 $V$ 的精度。

对 $n\to\infty$，$P_n(z)$ 在紧集上收敛到 $e^{z/2}$。系数逐项收敛且绝对值被 $2^{-k}/k!$ 控制。在 $z=ix$、$6\le x\le7$ 上，分母绝对值一致趋于 1。另一方面，
$$
|B_n|=\frac{H_n q^n|P_n(-z)|}{|g_n|},
\qquad
\log|B_n|\ge n\log q+\log H_n-O(n^{3/2}).
$$
因此若 $n\to\infty$ 且 $n=o((\log q)^2)$，就有 $|B_n|/q^2\to\infty$。由同一个不等式，在阶数变化的相位近似中若要求 $|B_n|=O(q^2)$，则对充分大的 $q$ 必须有 $n\ge c(\log q)^2$，其中 $c>0$ 是统一常数。

上传报告的 Padé 余项恒等式在偶数阶的符号应改为 $(-1)^n$。积分形式是
$$
e^zP_n(-z)-P_n(z)
=\frac{(-1)^n z^{2n+1}}{(2n)!}
\int_0^1 e^{tz} t^n(1-t)^n\,dt .
$$
逐次分部积分即可得到此式。报告使用绝对值给出的弦误差上界完全保留。将积分按中心对称并归一化为期望。对 $z=ix$，归一化积分模长在
$1-x^2/[8(2n+3)]$ 与 1 之间。下界用 $\cos u\ge1-u^2/2$ 和权 $(1-s^2)^n$ 的二阶矩 $1/(2n+3)$。由此和 Stirling 粗估，
$$
\log|W_n-V|=-2n\log n+O(n)
$$
在 $6\le x\le 7$ 一致成立。对归一化误差不超过 1 的逼近，达到固定相对精度至少需要 $n\log n\ge\log q-O(n)$。用已知的 $\mu(\pi)=2$ 可知 $|V|\gg q^{-5/2}$，从而选 $n=\lceil 2\log q/\log\log q\rceil$ 已足以取得固定相对精度，但实际分母超过 $q^2$。更一般地，在上述 $n=o((\log q)^2)$ 范围，若 $|W_n-V|\le|V|/2$，则约分后分子模至少为 $|B_n||V|/2$，并趋于无穷。[指数 2 的形式化范围](https://github.com/openai/math/blob/main/lean/docs/017.md)

这个障碍针对标准对角 Padé 单相位比较。它把下一条同路线的缺口缩小到 $n\gtrsim(\log q)^2$ 时的小素数约分上界。

## 另一条可验收的整数比较

若对每个充分大的 $q$ 能构造两条整数线性形式 $L_j=A_j\pi-B_j$，满足
$$
A_1B_2-A_2B_1\ne0,\qquad
\max_j|A_j|\le Cq,\qquad
\max_j|L_j|\le \frac1{2q},
$$
则至少一个整数 $qB_j-pA_j$ 非零。对该 $j$，
$$
1\le |qB_j-pA_j|
\le q|L_j|+q|A_j||\pi-p/q|,
$$
即给出 $q^2|\pi-p/q|\ge1/(2C)$。反过来，若 $\pi$ 的所有后续部分商都不超过 $M$，取分母第一次达到 $2q$ 的连分数收敛子及其下一项。两项分母都在 $[2q,2(M+1)^2q]$，相邻行列式为单位，两条线性形式的绝对值都小于 $1/(2q)$。因此独立双形式的尺度条件等价于坏逼近性。可独立测试 Machin 公式组合出的整数线性形式是否有这种共同尺度。这里尚无满足三项条件的显式构造。

[标准对角 Padé 的系数资料](https://citeseerx.ist.psu.edu/document?doi=e300d49488eac386400f5249f10bea35d4f8f7e6&repid=rep1&type=pdf) · [Lambert 递推资料](https://webusers.imj-prg.fr/~michel.waldschmidt/articles/pdf/IMPA2010Cours12.pdf) · [独立复算](verify_pade_phase_height.py)
