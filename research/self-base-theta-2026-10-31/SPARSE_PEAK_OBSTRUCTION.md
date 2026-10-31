# 稀疏高峰反例：长期平方平均有界不能判定坏逼近
2026-10-31。接续 [长期积分分解](MEAN_SQUARE_DECOMPOSITION.md)。
本次得到的是对所提证明路线的限制，不是 π 的证明或反例；不主张文献新颖性。

## 1. 结论
存在明确给定的无理数 β∈(1,2)，同时满足
\[
\liminf_{q\to\infty}q\|q\beta\|=0,\qquad
\sup_{t\ge0}\Theta_\beta(t)=\infty,\qquad
\sup_{T\ge1}\frac1T\int_0^T\Theta_\beta(t)^2\,dt<\infty.
\]
因此上一份笔记的长期平方平均发散目标，严格强于目标数非坏逼近。此前只称它为充分条件是正确的；本次补上不能反向推断的构造证明。

## 2. 构造
取普通连分数
\[
\beta=[1;a_1,a_2,\ldots],\qquad
a_k=\begin{cases}j+1,&k=2^j,\ j\ge1,\\1,&\text{其余 }k.\end{cases}
\]
部分商无界，但累计和增长很慢。若 J=⌊log₂N⌋，则
\[
\sum_{k=1}^N a_k=N+\frac{J(J+1)}2\le2N.
\tag{1}
\]
最后的不等式由 J(J+1)/2≤2^J≤N（J≥1）得出；N=1单独成立。

令 p_n/q_n 为收敛分数，q_(-1)=0、q_0=1。
q_n=a_nq_(n−1)+q_(n−2)，所以
\[
q_n\ge2^{\lfloor n/2\rfloor}.
\tag{2}
\]
前文恒等式给出
\[
q_n^2|\beta-p_n/q_n|
=\frac1{\alpha_{n+1}+q_{n-1}/q_n},
\quad \alpha_{n+1}=[a_{n+1};a_{n+2},\ldots].
\]
沿 n=2^j−1，右侧小于1/(j+1)，趋零。故 β 非坏逼近；由已有峰下界，Θ_β无界。此处的 β 是人为构造的数，不是 π。

## 3. 证明长期平均仍然有界
采用前文分解，固定 c=1/4：
\[
\int_0^T\Theta_\beta^2\le K_cT+\pi S_c(T).
\tag{3}
\]
每个本原方向在 S_c 中最多计算一次，其完整区间贡献不超过
\[
\frac{\pi}{2m\log\beta}.
\tag{4}
\]

先识别所有可能贡献的方向。因为 1<β<2、t≥0：
- q=0的本原方向满足 F≥1/β>1/2，不进入F<c；
- p=0的方向满足F≥β>1；
- pq<0时，β^(2t)(p−βq)²/(2β)>1/4，也不进入。
因此变更整体符号后，可以取p,q>0。

若某时刻F_(p,q)(t)<c，则m_(p,q)<c，并有
\[
q^2|\beta-p/q|=\frac{m_{p,q}}{1+p/(\beta q)}<c<\frac12.
\]
由Legendre判据，该既约分数必为一个收敛分数。

若相应区间与[0,T]相交，则
\[
c>F_{p,q}(t)\ge
\frac{\beta^{-2t}(p+\beta q)^2}{2\beta}
>\frac{\beta q^2}{2}\beta^{-2T},
\]
所以q<√(2c/β)β^T<β^T。
由(2)，所有有贡献的收敛分数编号均不超过
\[
N_T=\left\lfloor\frac{2\log\beta}{\log2}T+1\right\rfloor .
\tag{5}
\]

对任意收敛分数，
\[
\frac1{m_n}
=\frac{\alpha_{n+1}+q_{n-1}/q_n}{1+p_n/(\beta q_n)}
<a_{n+1}+2.
\]
于是由(1)、(4)，即使把所有相交区间都扩大为完整区间，也有
\[
S_c(T)
\le\frac{\pi}{2\log\beta}\sum_{n=0}^{N_T}(a_{n+1}+2)
\le\frac{2\pi}{\log\beta}(N_T+1).
\]
结合(3)、(5)得到显式线性界
\[
\int_0^T\Theta_\beta(t)^2dt
\le
\left(K_c+\frac{4\pi^2}{\log2}\right)T
+\frac{4\pi^2}{\log\beta}.
\tag{6}
\]
因此T≥1的平均一致有界。这里没有交换无限求和与不受控极限；紧时间段中有贡献的方向只有有限个。

## 4. 应调整的解析目标
对任意固定窗口长度H>0，保留移动窗口的最大值可避免长等待的稀释：
\[
\liminf_{q\to\infty}q\|q\beta\|=0
\iff
\sup_{s\ge0}\int_s^{s+H}\Theta_\beta(t)^2dt=\infty.
\tag{7}
\]
证明：左侧成立时，取m_j→0的平衡时刻t_j→∞。在[t_j−H/2,t_j+H/2]上，
\[
\Theta_\beta(t)^2
\ge\frac{\pi}{m_j\cosh(H\log\beta)},
\]
积分趋无穷。反向，若Θ有界，则所有长度H的窗口积分一致有界；再用前文Θ无界与非坏逼近的等价性即可。

(7)仍是等价重述，不是新的π算术输入，但它不会因稀疏高峰而漏判。
对π，长期平方平均若能证明发散仍然有用；若只能证明有界，不能因此证伪P2848。要证伪仍需统一峰高上界，或等价的全局丢番图下界。

## 5. 验证、来源与未完成部分
[sparse_peak_check.py](sparse_peak_check.py) 用标准库整数运算检查前4096项的累计和公式、2N上界和分母增长下界。有限检查只检验实现；无限结论由第2—3节证明承担。
Legendre判据参考 [Stanford continued fraction notes](https://crypto.stanford.edu/pbc/notes/contfrac/converge.html)。
高斯峰估计和线性余项常数K_c来自本目录前两份笔记，未作Lean形式化。

此前正积分族的gcd下界仍未闭合。本次没有把有限样本推广为该积分族永远失败的结论，也没有获得π的无界峰或有界峰结论。
