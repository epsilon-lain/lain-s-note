# 对角 Padé 的约分峰值与可变阶数

复核上传的累计研究包，其 SHA‑256 为 784d4ab696efee6ec4858093bfc99be83c449aea8af87a7c87718f8d6c6f08a9。24 条基础结果逐字节重现。独立程序对 $n\le200$ 逐奇素数复核可达最大赋值。本轮证明绝对线性界失败。在无限多个阶数上，额外消去至少达到 $(1-o(1))n\log n$。同时，所有输入都满足 $O(n\log n\log\log n)$ 的额外消去上界。这两条都允许 $p,q$ 变化。

## 记号与交接稿勘误

设 $P_n(z)$ 是标准 $[n/n]$ 指数 Padé 的分子，$P_n(0)=1$，分母为 $P_n(-z)$。设 $H_n=(2n-1)!!$，并定义
$$
T_n(t)=H_n P_n(2it)=\sum_{k=0}^n \binom{2n-k}{k}(2n-2k-1)!!(it)^k .
$$
对既约 $p/q>0$，令 $C_n=q^n T_n(p/q)=a_n+ib_n$、$c_n=\gcd(a_n,b_n)$。这里把交接稿的对应量除以了 $2^n$。令 $V_n$ 为 $P_n(2ip/q)/P_n(-2ip/q)-1$ 在 $\mathbb Z[i]$ 中约分后的分母。已有
$$
c_n\mid p^{2n-1},\qquad
|V_n|=\frac{H_n q^n|P_n(-2ip/q)|}{c_n\,2^{\varepsilon_n/2}},\quad \varepsilon_n\in\{0,1\}.
$$
详见[前一轮高度笔记](PADE_PHASE_HEIGHT.zh.md)。

交接稿 §2.4 只用最低赋值层的公式，在该层模 $\ell$ 整体消失后必须加入更高项。精确反例是 $(n,p,q,\ell)=(5,51,16,3)$。两项最低赋值都是 3，最低层的和为 $-3483648000$，赋值为 5。完整 $C_5=-1859999760-62126109i$ 的公因子有七阶。这组分数在目标区间 $3\le p/q\le7/2$ 内。

上传的 $n\le9$ 局部极值可以重现。继续枚举时，$n=61,\ell=11$ 的最大额外赋值已经是 109，因此基于低阶数据的 $C=\log6$ 失效。下面证明任意绝对 $C$ 都失效。交接稿中的固定分子界和固定输入的近似指数结论，只在各自量词下使用。

## 无界的归一化消去率

取奇素数 $\ell\ge5$，取奇整数 $a$ 满足 $1\le a<\ell$ 且 $-2a$ 是模 $\ell$ 的非零平方。令 $2n-1=a\ell^2$。每个这样的奇素数都能用 $a=\ell-2$。在 $\mathbb Z_\ell[i]$ 中考察 $T_n(\ell u)/H_n$。第 $k$ 项除去 $(i\ell u)^k$ 后的系数为
$$
\frac{\binom{2n-k}{k}}{\prod_{j=1}^k(2n-2j+1)}.
$$
对 $k\ge3$，分母赋值至多 $\lceil k/\ell\rceil+\lceil k/\ell^2\rceil<k$，因为所有因子均小于 $\ell^3$，且在奇数序列中每过 $\ell^s$ 项最多出现一个 $\ell^s$ 的倍数。对 $k=1$ 比系数为 1。对 $k=2$ 比系数为 $(n-1)/(2n-1)$。因此整个多项式的系数是 $\ell$ 进整数，实部模 $\ell$ 为
$$
1-\frac{n-1}{a}u^2\equiv1+\frac{u^2}{2a}\pmod\ell.
$$
它有非零单根。Hensel 提升给出 $u_0\in\mathbb Z_\ell^\times$，使 $T_n(\ell u_0)$ 的实部精确为零。

由 Padé 接触条件，
$$
F(t)=e^{2it}T_n(-t)-T_n(t)
$$
的 Taylor 展开从 $t^{2n+1}$ 起。由于 $T_n$ 是高斯整系数多项式，$t^j$ 的系数乘 $j!$ 是高斯整数。令
$$
K_\ell=\left\lceil(2n+1)\left(1-\frac1{\ell-1}\right)\right\rceil.
$$
在 $t=\ell u_0$，由 $v_\ell(j!)\le j/(\ell-1)$ 可知每项都被 $\ell^{K_\ell}$ 整除，且赋值随 $j$ 趋于无穷。实部为零给出 $T_n(-t)=-T_n(t)$，于是
$$
F(t)=-(1+e^{2it})T_n(t).
$$
前因子模 $\ell$ 是 $-2$，为单位，故 $T_n(t)$ 的两个坐标都被 $\ell^{K_\ell}$ 整除。这里在秩 2 的环 $\mathbb Z_\ell[i]$ 中逐坐标理解，$\ell$ 分裂时也成立。

令 $p/q\equiv\ell u_0\bmod\ell^{K_\ell+1}$，并对每个其他奇素数 $s\mid H_n$ 加上 $p/q\equiv0\bmod s^{v_s(H_n)}$。CRT 将条件合成为 $p\equiv r$、$q\equiv1\bmod M$。例如取 $q=(M+1)^j$，选 $p$ 在 CRT 规定的类，且 $p\equiv1\bmod M+1$，再从中取最接近 $\pi q$ 的代表元。这样 $\gcd(p,q)=1$、$p/q\to\pi$，且其余素数至少保留各自的双阶乘指数。因此
$$
\log(c_n/H_n)\ge (K_\ell-v_\ell(H_n))\log\ell
\ge \bigl[(2-3/(\ell-1))n-O(1)\bigr]\log\ell .
$$
最后一步用 $v_\ell(H_n)\le n/(\ell-1)+2$，因为只有前两级素数幂可能贡献。故对于任意绝对 $C$，即使将输入限制在 $3\le p/q\le7/2$，也有 $\log c_n-\log H_n>Cn$ 的反例。

还有无限多个 $\ell$ 可以取 $a=1$。给定有限的奇素数集合，取 $m$ 为它们乘积的两倍，从 $m^2+2$ 取一个奇素因子。它在集合之外，而且 $-2$ 是模它的平方。此时 $n=(\ell^2+1)/2$，上述下界就是 $(1-o(1))n\log n$。

同余条件还能与常数依赖 $n$ 的二次级逼近同时实现。取 $\pi$ 的相邻收敛子向量 $(P_j,Q_j),(P_{j-1},Q_{j-1})$，两列组成整基。模 $M$ 解出系数类 $(u,v)$，使线性组合对应 $(r,1)$，故 $\gcd(u,v,M)=1$。先取 $1\le u\le M$，再用 CRT 选 $v>0$，令其在整除 $u$ 而不整除 $M$ 的素数上同余 1。这样可取互素正系数 $u,v$，大小仅依赖 $M$。由收敛子的余项界，组合后的分数满足 $q^2|\pi-p/q|\le(u+v)^2$。这个常数随同余模数变化，构造也尚未使 Padé 分母达到 $O(q^2)$。

## 短前缀给出统一上界

令 $m_j=2n-2j+1$、$h=\lfloor\log_3(2n-1)\rfloor$、$J=\min(n,2h)$。对每个 $p$ 的奇素因子 $\ell$，写 $\nu=v_\ell(p)$。若 $v_\ell(c_n)>v_\ell(H_n)$，则必须有某个 $k\ge2$ 的项赋值不高于常数项。$k=1$ 比系数是 1，因此不合格。忽略二项式系数的非负赋值，有
$$
k\nu\le \sum_{j=1}^k v_\ell(m_j)
\le \frac{k}{\ell-1}+\lfloor\log_\ell(2n-1)\rfloor .
$$
所以 $k\le J$，且 $\nu\le\sum_{j=1}^J v_\ell(m_j)/j$。设 $p_\ast$ 是由满足上述短前缀存在条件的奇素数构成的 $p$ 的全素幂因子。则
$$
\log p_\ast\le\sum_{j=1}^J\frac{\log m_j}{j}=O(\log n\,\log\log n).
$$
其余素数的赋值至多等于基线。对 $p_\ast$ 中的素数应用 $c_n\mid p^{2n-1}$，得到
$$
\log(c_n/H_n)\le(2n-1)\log p_\ast=O(n\log n\,\log\log n).
$$

当 $p/q\to\pi$ 时，由 $P_n(z)\to e^{z/2}$ 在紧集上一致收敛，可对足够大的 $n$ 取 $|P_n(-2ip/q)|\ge1/2$。固定阶数的候选不能保持相对相位精度，因为对应的分子分母是有理系数非零多项式，在超越点 $2i\pi$ 没有零点。故沿任何同时保持相对相位精度和 $|V_n|\le Cq^2$ 的序列，都有 $n\to\infty$，且
$$
\log p_\ast\ge\frac{(n-2)\log q-\log(2\sqrt2 C)}{2n-1}.
$$
特别地，$\log n\,\log\log n\gg\log q$，于是 $p_\ast\ge(1-o(1))\sqrt q$。另一端，由 $c_n\le p^{2n-1}$ 和 $H_n^{1/n}\sim2n/e$，当 $p/q\in[3,7/2]$ 时，$n\ge20q$ 足以排除 $|V_n|=O(q^2)$。综上，对某个绝对常数 $\delta>0$，残余窗口为
$$
q^{\delta/\log\log q}\le n<20q.
$$

## 同一输入至多一个候选阶

若 $m>n\ge3$，对 $D_j=\operatorname{Im}(C_{n+j}\overline C_n)$ 有 $D_0=0$、$D_1=p^{2n+1}$，且满足与 $C$ 相同的递推。因 $p/q\le7/2$，将第 $j$ 项除以 $p^j$ 后，递推系数大于 2，故后续均为正，且
$$
D_{m-n}\le p^{2n+1}q^{m-n-1}H_m/H_{n+1}.
$$
另一方面，$c_nc_m\mid D_{m-n}$。结合 $|P_n(-2ip/q)|,|P_m(-2ip/q)|\ge1/2$，足够大的 $n$ 给出
$$
|V_n||V_m|\ge\frac{H_nH_{n+1}}{8(p/q)^{2n+1}}.
$$
这个下界与阶数跨度 $m-n$ 无关。符合残余窗口的 $n$ 已大到使右端远超 $q^4$。所以沿 $p/q\to\pi$，对任意固定的高度常数 $C$，充分大的每个输入至多有一个阶数同时保持相对相位精度和 $|V_n|\le Cq^2$。单个孤立候选仍未排除。此处没有判定 $\pi$ 的坏逼近性。

## 下一轮可验收任务

- Kimi 从 $2n-1=\ell^2$ 的单根族出发，研究根在 $p=\ell u,\ q\asymp\ell^2$ 的共振尺度上有没有短有理代表。先选择 $\ell=11,17,19$，交模 $\ell^k$ 的证书和 $v_\ell(c_n)$。将短代表高度和任意 CRT 提升分别列出。
- GLM 独立研究短前缀必要条件对实际 $\pi$ 收敛子分子的限制。对满足 $q^{\delta/\log\log q}\le n<20q$ 且 $p_\ast\ge(1-o(1))\sqrt q$ 的阶数做筛选，并提出可证明的统一排除条件。有限筛选用于确定下一条引理。

[独立复算脚本](verify_pade_gcd_prefix.py)。用参数 --max-order 200 可复跑局部极值。
