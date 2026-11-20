# 借鉴 Navier–Stokes 爆破构造的两个 π 临界问题

目标是从既约 $p/q$ 沿某列满足 $\epsilon=q^2|\pi-p/q|\to0$ 的反设出发，得到矛盾。允许 $\epsilon$ 比任何负幂都下降得慢。本文给出两个已落实的子结果及其下一项任务。

[OpenAI 的 Navier–Stokes 论文](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf) 第 3 节先在不同收缩尺度的坐标中构造集中主部，再在环带内用波的协方差实现需要的应力，并用矩条件和逐阶修正控制剩余项。其 Navier–Stokes 构造含光滑外力。这里借鉴核心轮廓、外区匹配和残差预算的组织方法。数论中的相位始终固定为 $2n^2$。

## A. 一个临界尺度内的四次积分

沿用 [theta 波笔记](THETA_STRIP_AND_GAUSSIAN_AUDIT.zh.md)，令 $\alpha=1/\pi$，$w(n)=1+(-1)^n/2$，
$$
W(t)=\sum_{n\in\mathbb Z} w(n)e^{-tn^2+2in^2},
\qquad E(T)=\int_T^{2T}|W(t)|^4\,dt .
$$
设有一列既约 $a/b$ 使 $b^2|\alpha-a/b|\to0$，令 $\delta=\alpha-a/b$、$\eta=2\pi|\delta|$。用模 $2b$ 的高斯 Poisson 展开时，加权 Gauss 均值 $c_0$ 满足 $|c_0|\ge(2b)^{-1/2}$，其余 Fourier 系数为 $O(b^{-1/2})$。在 $t=s\eta$ 且 $s$ 留在固定的正紧区间时，
$$
 W(t)=c_0\sqrt{\frac{\pi}{t-2\pi i\delta}}\bigl(1+O(e^{-c/(b^2\eta)})\bigr).
$$
强度 $\sqrt{t}|W(t)|^2$ 的主部在 $\eta$ 处达峰，除去常数后的轮廓为 $\sqrt{s}/\sqrt{1+s^2}$。积分得到
$$
 E(\eta)=\frac{\pi^2|c_0|^4}{\eta}(\arctan2-\pi/4)(1+o(1))
 \asymp (b^2\eta)^{-1}.
$$
所以 $E(T)=O(1)$ 足以反驳反设。反过来，坏逼近时 $|W(t)|=O(t^{-1/4})$，同一结论成立。

写 $r_w(N)=\sum_{u^2+v^2=N}w(u)w(v)$，从而
$$
 W(t)^2=\sum_{N\ge0}r_w(N)e^{-tN+2iN}.
$$
把 $(N,M)$ 的成对项合并后，四次积分的真正相位为 $\cos(2(N-M))$，对应的权重为
$$
 \frac{e^{-T(N+M)}-e^{-2T(N+M)}}{N+M}
$$
（$N=M=0$ 取 $T$）。$r_w(N)/r_2(N)$ 在 $N\equiv0,1,2\bmod4$ 的系数分别为 $9/4,3/4,1/4$，$3\bmod4$ 没有表示。由 $r_2(2n)=r_2(n)$ 和 [Borwein–Choi，式 (3.16)](https://www.sfu.ca/~schoia/paper/dirichlet.pdf) 的 $\sum_{n\le X}r_2(n)^2=4X\log X+O(X)$，得 $\sum_{N\le X}r_w(N)^2=(25/4)X\log X+O(X)$。理由是 $N$ 为奇数、恰被 2 整除、以及被 4 整除的三部分对原二阶矩的系数分别是 $2,1,1$。因此对角项独自产生
$$
 (25/8)(\log 2)\log(1/T)+O(1).
$$

固定非零移位 $h$ 的贡献是 $O(1)$。还可统一控制短移位窗口。对 $h\ne0$，计数 $u^2+v^2-z^2-y^2=h$，令 $(A,B)=(u-z,v-y)$，则 $A(u+z)+B(v+y)=h$。若原坐标在 $[-R,R]$，把 $d=\gcd(A,B)$ 除去，$d\mid h$。最大坐标为 $M$ 的原始方向有 $O(M)$ 个，对每个方向的方程至多 $O(1+R/M)$ 个盒内解，且 $M=O(R/d)$。求和得到 $O(R^2\sum_{d\mid h}d^{-1})$。因此 $\sum_{N\le X}r_w(N)r_w(N+h)\ll(X+h+1)\sigma_{-1}(h)$，$h>0$。加上指数权重后，每个正移位对 $E(T)$ 的贡献至多 $O((1+Th)e^{-Th}\sigma_{-1}(h))$，其中 $\sigma_{-1}(h)=\sum_{d\mid h}1/d$。由于 $\sum_{1\le h\le H}\sigma_{-1}(h)=O(H)$，任何 $H=o(\log(1/T))$ 的窗口都不足以抵消对角发散。

**任务 A。** 找一个随 $T$ 缩小而扩大的移位窗，在固定相位 $\cos(2h)$ 下证明负向相消，或者给出反例说明所提窗仍不足。可以先核对 [Chamizo 的移位相关定理](https://www.jstage.jst.go.jp/article/jmath1948/51/1/51_1_237/_article/-char/ja) 在带权同余条件下的适用范围，再分别算主项和带 $h$ 的误差。目标是余项足够小到可以和其它区域拼接，最后争取 $E(T)=O(1)$。控制平均角度，或只使用 $\mu(1/\pi)=2$，不能区分所有固定角。对照 $[0;1,2,3,\ldots]$ 的无理性指数是 2，然而它不是坏逼近数。

## B. 指数二行列式证明的参数端点

[OpenAI 的 π 论文参数节](https://github.com/openai/math/blob/main/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/conclusion.tex) 给出一组充分条件。写 $\nu=2+h$，$h>0$。原条件包括
$$
 0<\theta<A<B<1,\qquad CB<1,\qquad B<C\theta,\qquad
 g=\nu(A(1-\eta)-\theta)-(1-\theta)>0 .
$$
前两条含 $C$ 的条件相乘，给 $B^2<\theta$，故 $A<\sqrt\theta$。令 $x=\sqrt\theta$，任何原来的可行参数因此满足
$$
 g<(2+h)x-(1+h)x^2-1\le\frac{h^2}{4(1+h)} .
$$
在 $\nu=2$，第一个上界就是 $-(1-\sqrt{\theta})^2<0$。原证还要求比较误差小于 $g$，而原误差预算至少有 $\Lambda/w_*$，其中 $w_*=\min_i\lceil\log q_i\rceil$。即使逐组适应 $h$，仍需 $h\gg w_*^{-1/2}$。对异常分数 $\delta_i=\epsilon_i/q_i^2$，共同可用的 $h$ 必须不大于所有 $\log(1/\epsilon_i)/\log q_i$ 的最小值。因此分母最小的种子也须满足 $\epsilon_i\le\exp(-c\sqrt{\log q_i})$。这条是沿用原误差预算的必要条件。单单 $\epsilon_i\to0$ 不会提供它。这个结论只排除照搬同一组参数条件，未排除新的行列式方案。

**任务 B。** 追查原来 $CB<1$、$B<C\theta$ 各自由什么计数和算术误差产生。尝试改造一组有限的插值列或行，在端点允许 $A>(1+\theta)/2$，或消掉正的 $1/w_*$ 量级损失并给出新的选择次序。先在明确定义的小块上同时认证三件事：有理中心的秩、消去了哪一批低阶项、所得非零子式的真实清分母费用。理想输出是清分母后非零整数被一个只依赖 $\epsilon_i$ 的趋零量压住。有限秩改进可作为阶段成果，必须标出它覆盖的参数范围。
