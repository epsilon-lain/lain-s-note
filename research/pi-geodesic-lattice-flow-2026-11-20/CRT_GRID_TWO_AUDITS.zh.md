# 两个证明包复核与续研

复核 `abc_reduction_audit_bundle.zip` 和 `pi_double_witness_verifier.zip`。以下沿用 [原 CRT 任务](ABC_H10_RESEARCH_TASKS.zh.md) 的记号。两个续研结论都是无条件的。它们排除了限定构造，π 的坏逼近性仍待判定。

## 1. 双见证的复算与勘误

五组基准在 $2\le n\le96$ 上完整重跑。生成的 JSON 与附件逐字节相同，SHA-256 为 `057b3ea9ff0ec65b912a15d81e840476887151d32915a0ace68d63c2a21d707e`。所有辅助断言也通过。

正文第 6 节取 $X=4^n$ 时，实际应有 $q_1q_2\le B^2 16^n$，印出的指数 $2n$ 多了一倍。改正后仍得到原阈值 $n>2B^4 d_{\max}^2$。第 7 节的有效区间右端是一个保守下界，从它低于 64 无法推出配对在 64 处失败。补跑 `python3 pi_witness_verifier.py single 3 64 16 16 4`，严格枚举并拒绝了全部 8192 个候选，因而该尺度的无解结论成立。

## 2. 网格可以换阶，原分子间距仍必须增大

设 $\alpha$ 为任意实无理数。固定 $B,K>0$，假设每个充分大的正整数尺度 $X$ 都可找到同一网格 $k/V_n$ 上的两条不同既约方向 $p_i/q_i$，满足
$$
 q_i\le BX,\qquad |q_i\alpha-p_i|\le K/X.
$$
网格阶 $n$ 可随 $X$ 变，每一对的原分子间距记作 $d=|k_1-k_2|$。则存在一列整数尺度 $X_m\to\infty$，其上任何合格配对都必须有 $d\gg_{B,K}\log X_m$。特别地，即使针对每个尺度任意挑选网格，$d=o(\log X)$ 的全尺度方案也不成立。

令 $T=\lceil2K\rceil+1$。给定任意既约 $p/q$，在 $X=Tq$ 取与它独立的一条见证。取整数行列式，得到
$$
 1\le BT\,q|q\alpha-p|+K/T .
$$
设 $c=(1-K/T)/(BT)>0$。于是大分母上的加权误差统一下界为 $c$。在 $X\to\infty$ 的合格配对中，两行的分母都趋于无穷。这是因为若分母有界，$\alpha$ 无理便使相应线性形式的最小正误差远离零。因此
$$
 cX/K\le q_i\le BX.
$$
令 $D=|p_1q_2-p_2q_1|$。两行估计给 $1\le D\le2BK$，而网格恒等式是 $D=dq_1q_2/V_n$。令 $A=c^2/(2BK^3)>0$，则 $V_n\ge AX^2$。

$V_n$ 严格递增，且精确递推式为
$$
 V_m/V_{m-1}=\frac{(4m-1)(4m+1)}{(2m-1)\operatorname{odd}(m)} .
$$
取 $m=2^r$，此比值大于 $8m$。令 $X_m=\lfloor\sqrt{V_{m-1}/A}\rfloor+1$，于是任何可用网格必有 $n\ge m$，因为 $AX_m^2>V_{m-1}$。从 $D\ge1$ 及 $q_i\le BX_m$ 得
$$
 d\ge \frac{V_m}{B^2 X_m^2}\gg_{B,K} m .
$$
由中心二项式界知 $\log V_{m-1}\asymp m$，即 $m\asymp\log X_m$。证明中使用的加权误差下界由假设的全尺度覆盖推出，无须假定 $\alpha$ 原先是坏逼近数。

**真正需要改变的是公约数。** 令 $g=\gcd(k_1,k_2,V_n)$，$G=\gcd(q_1,q_2)$。约分后两分数的最小公分母为 $W=V_n/g=\operatorname{lcm}(q_1,q_2)$，并且
$$
 D=(d/g)G .
$$
所以 $d/g$ 仍受 $2BK$ 限制。在上述跳阶尺度，$d\gg\log X$ 只能通过 $g\gg\log X$ 来实现。候选应该按有效网格 $W=V_n/g$ 研究。选到合适大小的 $W$ 后，还要证明 $\pi W$ 附近的两个整数能得到所需的真实分母。

附带有限检查采用另一列容易精确计算的尺度 $\lfloor\sqrt{V_m/(2MB^2)}\rfloor$，其中 $M$ 是 $d$ 的上限。[复算脚本](verify_grid_gap_ladder.py) 枚举每个尺度可能使用的全部 $n$，不限于 $n=m$。在 $m=8,16,32,64,128,256,512$ 和三组 $(B,K,M)=(1,1,1),(4,4,4),(8,8,4)$ 上，共 21 个尺度都没有配对。在 $m=4$，后两组恰有同一对 $n=2,\ (k_1,k_2)=(987,990)$。

## 3. abc 附件的勘误与实际参数

令 $w=3t-2S$，$R=\operatorname{rad}(tw)$。完整乘积的正确公式是
$$
 \operatorname{rad}(w\cdot2S\cdot3t)=\operatorname{rad}(6Stw)\le6 R\operatorname{rad}(S).
$$
例如 $t=1,S=2$ 给原始三元组 $(-1,4,3)$。完整 radical 为 6，而附件写出的上界 $2R\operatorname{rad}(S)$ 只有 4。附件六次幂小节中的 $\operatorname{rad}(uvx^6y^6)=\operatorname{rad}(uv)|xy|$ 应写成不等式。附件给出的通用反例三例均复算通过，其中 $S$ 始终是奇数。

实际模型中，设 $z_0,z_1$ 为原 Beta 积分的两条有理斜率，$V=V_n$，$e=2+v_2(n!)$，$\delta=z_1-z_0=\sigma 2^e/V$，$\sigma\in\{-1,1\}$。取最简 $V z_0=a/P$ 和 $S=2^eP$。候选满足 $0<t<S$、$t$ 为奇数及 $t\equiv-\sigma a\bmod P$，故 $\gcd(t,S)=1$。注意 $P$ 为奇数。原始三元组取公因子 $\gcd(3t,2S)\le3$ 除掉。给定 abc 假设，因 $\operatorname{rad}(S)\le2^{1-e} S$，
$$
 \log R\ge e\log2-\frac{\varepsilon}{1+\varepsilon}\log S-O_\varepsilon(1).
$$
先取充分小的 $\varepsilon$，就有任意 $\eta>0$ 对应的 $\log R\ge(\log2-\eta)n$，充分大的 $n$ 成立。附件通用反例的 $S$ 与实际模型中的偶数 $S$ 不同。

后面的无条件推导只需要两个已有算术界。第一，$P$ 包含区间 $((4n+1)/3,2n-1]$ 的全部素数 $\ell$，当 $n$ 充分大。展开 $z_0$ 的二项式平均时，每一项分母中只有一个奇倍数 $\ell$，所以 $\ell z_0$ 模 $\ell$ 恒为 $4(-1)^{(\ell-1)/2}$，而该区间内 $v_\ell(V_n)=0$。第二，
$$
 \log S\le(4-3\log2+o(1))n.
$$
取 $L$ 为不超过 $4n+1$ 的正奇数的最小公倍数。$z_0$ 的既约分母整除 $L$。对奇素数 $\ell$，$v_\ell(V_n)$ 的阶乘公式在每层 $\ell^j$ 的贡献是 $\lfloor(4n+1)/\ell^j\rfloor-2\lfloor2n/\ell^j\rfloor\in\{0,1\}$，而 $\ell^j>4n+1$ 时贡献为零。因此 $V_n\mid L$，从而 $P=\operatorname{den}(V_nz_0)\mid L/V_n$。由素数定理，$\log L=4n+o(n)$，而 $\log V_n\sim4n\log2$。连同 $e=n+O(\log n)$ 即得估计。

## 4. 实际 CRT 无条件排除有界余因子的所有共同幂次 $r\ge3$

对任意固定整数 $r\ge3$ 和常数 $C$，假设
$$
 t=u x^r,\qquad w=v y^r,\qquad 0<|uv|\le C
$$
在无限多个 $n$ 成立，四个因子都是非零整数。固定一组出现无穷多次的 $(u,v)$。若 $r$ 有奇素因子 $s$，改用 $s$ 次幂。否则改用四次幂。以下把相应幂次记为 $s$。奇数性保证 $u,v,x,y$ 均为奇数。由于 $w\equiv3t\bmod\ell$ 且 $t$ 与 $S$ 互素，$v/(3u)$ 在上述整段素数区间内都必须是模 $\ell$ 的 $s$ 次幂，有限个固定坏素数可忽略。

若 $s$ 是奇素数而 $c_0=v/(3u)$ 不是有理数的 $s$ 次幂，则 $Z^s-c_0$ 在 $\mathbb Q[Z]$ 上不可约。因为若有次数 $1\le j<s$ 的真因子，其根之积给出一个有理数的 $s$ 次幂等于 $c_0^j$，再用 $\gcd(j,s)=1$ 就可推出 $c_0$ 本来是 $s$ 次幂。传递的 Galois 群含有 $s$ 循环，由 Chebotarev 素数计数定理，正比例的素数上没有一次因子。该计数渐近保证每个充分大的固定比例区间都含这样的素数，与上面的模条件矛盾。[Chebotarev 与模素数因子分解](https://docnum.univ-lorraine.fr/public/DDOC_T_2013_0183_BARBULESCU.pdf)

若 $s=4$，先用同样的二次论证推出 $c_0=b_0^2$，取 $b_0>0$。若 $b_0$ 不是有理平方，$\mathbb Q(i,\sqrt{b_0})$ 有固定 $i$ 而反转 $\sqrt{b_0}$ 的元素。于是有一类正比例的素数 $\ell\equiv1\bmod4$ 上，$b_0$ 和 $-b_0$ 都是二次非剩余，$c_0$ 无四次根，矛盾。因此这时 $c_0$ 也是有理四次幂。

写 $c_0=(a_0/b_1)^s$ 为既约形式。$a_0,b_1$ 都为奇数。记 $U=b_1x$、$Y=a_0y$，有精确恒等式
$$
 3u(U^s-Y^s)=2b_1^s S,\qquad v_2(U^s-Y^s)=e+1 .
$$
当 $s$ 为奇数，$v_2(U^s-Y^s)=v_2(U-Y)$。当 $s=4$，$v_2(U^4-Y^4)=v_2(U-Y)+v_2(U+Y)+1$，且其中一个 $U\pm Y$ 恰有估值 1。不论哪种，总有一个非零的 $U\pm Y$ 的绝对值至少是 $2^{e-1}$。但 $|U|+|Y|=O(S^{1/s})$，由高度界知
$$
 \frac13(4-3\log2)<\log2.
$$
因此 $|U|+|Y|=o(2^e)$，矛盾。结论覆盖固定幂次、有界余因子的这一族。放宽到随 $n$ 增长的余因子时，上述素数计数不能按固定数域直接使用。

## 5. 平方情形的精确剩余问题

若 $t=u x^2,w=v y^2$ 且 $|uv|$ 有界，用相同素数区间和二次论证先得到 $v/(3u)$ 是有理平方。由于 $\gcd(t,w)=\gcd(t,2S)=1$，比较 $w/(3t)$ 的素数估值，所有充分大的候选只能是
$$
 (t,w)=(x^2,3y^2)\quad\text{或}\quad(3x^2,y^2).
$$
第一种强制 $3\mid S$，第二种强制 $S\equiv1\bmod3$。分别得到 $(x-y)(x+y)=2S/3$ 和 $(3x-y)(3x+y)=2S$。这两因子的最大公因子恰为 2。一因子只含一个 2，另一因子承载整个 $2^e$。

在 $|k-\pi V_n|=O(1)$ 的网格窗口内，$t/S=\theta_n+O(2^{-e})$，其中 $\theta_n=\mathbb E[1/(1+Z)]$，$Z$ 服从 $\operatorname{Beta}(n+1/2,n+1)$。在 $1/2$ 处二阶展开。由 Beta 矩公式，
$$
 \mathbb E(Z-1/2)=-1/(8n)+O(n^{-2}),\quad
 \mathbb E(Z-1/2)^2=1/(8n)+O(n^{-2}),\quad
 \mathbb E(Z-1/2)^4=O(n^{-2}).
$$
三阶余项期望为 $O(n^{-3/2})$，故 $\theta_n=2/3+5/(54n)+o(1/n)$。依两种平方情形分别令 $(X,Y)=(x,y)$ 或 $(3x,y)$，取正号，则 $Y/X=\sqrt{w/(3t)}\sim\sqrt5/(6\sqrt n)$。于是因子比应满足
$$
 \frac{X+Y}{X-Y}=1+\frac{\sqrt5}{3\sqrt n}+o(n^{-1/2}).
$$
下一项精确任务是检验这一对近等因子能否同时承担规定的素数幂及全量 $2^e$，并使 $t$ 命中原 CRT 同余 $t\equiv-\sigma a\bmod P$。即使平方族全部排除，一般整数约分策略仍需独立分析。


## 附记。上一轮 Niven 任务的双素数模型

沿用 [Niven 递推和 $T_n$](THETA_WAVE_SEED_CONTRADICTION.zh.md)。在 $p$ 是两个大小相差至多一倍的不同奇素数之积、$q\le n\le2q$ 且 $n$ 为偶数时，实际上 $T_n\to\infty$。至少一个素因子 $\ell$ 满足 $\ell^2\nmid 2n-1$，否则 $p^2\mid2n-1$。设 $c_n=(2n)!/n!$。将 $A_n$ 的系数除以 $c_n q^n$ 后，常数项为 1，$k=2$ 的项为 $-(p/q)^2(n-1)/(4(2n-1))$，有正的 $\ell$ 赋值。对偶数 $k\ge4$，分母 $k!(2n)(2n-1)\cdots(2n-k+1)$ 的 $\ell$ 赋值至多 $2+2(k-1)/(\ell-1)<k$，这里 $\ell\ge5$ 且 $\ell^3>2n$ 对充分大的 $q$ 成立。因此余下所有项也有正赋值，得到 $v_\ell(A_n)=v_\ell(c_n)\le n/(\ell-1)+2$。另一素因子 $\ell'$ 由 $g_n\mid p^{2n-1}$ 控制，且 $\ell'\le\sqrt{2p}<\sqrt{8q}$。所以 $g_n\le \exp(o(n))(8q)^n$。又 $S_n\ge q^{2n}$，由已发表的无理性测度上界可用 $|\pi-p/q|\ge q^{-8}$，于是 $S_n|\pi-p/q|/g_n\to\infty$。这关闭了上轮提出的受限模型，任意分子及跨阶线性组合仍需分析。[无理性测度上界](https://doi.org/10.2140/moscow.2020.9.407)
