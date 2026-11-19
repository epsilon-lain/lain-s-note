# 二次振荡和的精确尺度判据

对 $t>0$ 定义两个绝对收敛的和
$$
 S_+(t)=\sum_{n\in\mathbb Z}e^{-tn^2+2\mathrm i n^2},\qquad
 S_-(t)=\sum_{n\in\mathbb Z}(-1)^n e^{-tn^2+2\mathrm i n^2}.
$$
相位 $2n^2$ 是以弧度计的整数角。**$\pi$ 坏逼近，等价于两个和都满足 $|S_\pm(t)|=O(t^{-1/4})$，其中常数须对 $0<t\le1$ 统一。** 下文推导精确恒等式、统一估计的等价性以及单个共振的形状。

## 正项恒等式

令
$$
 Z(t)=\sqrt{\frac{t}{2\pi}}\bigl(|S_+(t)|^2+|S_-(t)|^2\bigr).
$$
在平方和的双重求和里，同奇偶的 $(n,m)$ 权重为 $2$，异奇偶权重为零。写 $n=u+v,\ m=u-v$，得到
$$
 |S_+(t)|^2+|S_-(t)|^2
 =2\sum_{u,v\in\mathbb Z}e^{-2t(u^2+v^2)+8\mathrm iuv}.
$$
采用 $\widehat f(\xi)=\int_\mathbb R f(v)e^{-2\pi \mathrm iv\xi}\,dv$。对 $f(v)=e^{-2tv^2+8\mathrm iuv}$ 做泊松求和，得
$$
 \sum_{v\in\mathbb Z} e^{-2t v^2+8\mathrm iuv}
 =\sqrt{\frac{\pi}{2t}}\sum_{k\in\mathbb Z}
       \exp\!\left(-\frac{(4u-\pi k)^2}{2t}\right).
$$
所以
$$
 \boxed{\displaystyle Z(t)=
 \sum_{u,k\in\mathbb Z}\exp\!\left(-2tu^2-\frac{(4u-\pi k)^2}{2t}\right).}
$$
原求和绝对收敛。变换后的正项和收敛，因为指数是固定 $t>0$ 上由可逆线性变换得到的正定二次型。[泊松求和公式](https://dlmf.nist.gov/1.8.iv)

## 格点及最短向量

右式为平面格
$$
 \Lambda_t=\left\{\left(\sqrt{2t}\,u,\frac{4u-\pi k}{\sqrt{2t}}\right):u,k\in\mathbb Z\right\}
$$
上的高斯和。其基本平行四边形面积为 $\pi$。若最短非零向量的欧氏长度记为 $\lambda(t)$，则半径 $R$ 的圆内至多有 $(1+2R/\lambda(t))^2$ 个格点，因为以这些点为圆心画半径 $\lambda(t)/2$ 的开圆，圆内部互不重叠。对整数半径分层可知，只要 $\lambda(t)\ge c_0>0$，$Z(t)$ 就有统一上界。反过来，最短向量的全部整数倍均出现在正项和中，故当 $\lambda(t)\le1$ 时
$$
 Z(t)\ge \left(2\lfloor1/\lambda(t)\rfloor+1\right)e^{-1}.
$$
因此统一有界性等价于 $\inf_{0<t\le1}\lambda(t)>0$。

置 $\beta=4/\pi$。对于 $u\ne0$，向量长度的平方至少为
$$
 2tu^2+\frac{(4u-\pi k)^2}{2t}\ \ge\
 2|u||4u-\pi k|=2\pi u^2|\beta-k/u|.
$$
对于 $u=0$ 的非零向量，平方长度至少为 $\pi^2/(2t)$。所以 $\beta$ 坏逼近时，最短长度有统一正下界。若 $\beta$ 非坏逼近，取 $u_j^2|\beta-k_j/u_j|\to0$，其中 $u_j\to\infty$。令 $t_j=|4u_j-\pi k_j|/(2u_j)$，则 $t_j\to0$，并且对应平方长度恰为 $2\pi u_j^2|\beta-k_j/u_j|\to0$。

最后，正实数 $x$ 与 $4/x$ 的坏逼近性相同。只须比较足够接近 $4/x$ 的正既约分数 $p/q$，此时 $p/q$ 远离零且有界，而
$$
 q^2|4/x-p/q|=\frac{q}{xp}\,p^2|x-4q/p|.
$$
离 $4/x$ 较远的分数自动满足一个统一正下界。这个变换自身为逆，也保持通常的无理性指数，因为近似分数经过变换后的既约分母介于 $p/4$ 和 $p$，与 $q$ 同阶。取 $x=\pi$ 完成判据证明。

## 一次共振的形状与误差上界

取一对互素整数 $(u,k)$ 且 $u\ne0$，令 $\eta=u^2|\beta-k/u|$，$t_*={|4u-\pi k|}/{(2|u|)}$，并用 $s=\log(t/t_*)$ 观察该对。对应向量长度 $\ell$ 满足
$$
 \ell^2=2\pi\eta\cosh s.
$$
该向量是格中的本原向量，可以扩充为一组格基。把格点按平行于它的直线分层。编号 $b\in\mathbb Z$ 的直线与原点的垂直距离为 $|b|\pi/\ell$，每条线上的格点间距为 $\ell$。再用一次一维泊松求和，任意相位平移 $\xi$ 上的纵向高斯和满足
$$
 \sum_{a\in\mathbb Z}e^{-(a\ell+\xi)^2}
 =\frac{\sqrt{\pi}}{\ell}
 \left(1+2\sum_{j\ge1}e^{-\pi^2j^2/\ell^2}\cos(2\pi j\xi/\ell)\right).
$$
令 $\rho=e^{-\pi^2/\ell^2}$。$b=0$ 的直线含原点，其纵向平移为零，因此单这一条线已贡献至少 $\sqrt\pi/\ell$。对所有直线，余弦的绝对值至多为 1。将垂直方向的高斯和乘进去，得到
$$
 1\le\frac{\ell Z(t)}{\sqrt\pi}
 \le \left(1+\frac{2\rho}{1-\rho}\right)^2.
$$
特别地，当 $\ell\le1$ 时，
$$
 Z(t)=\frac{1+\varepsilon(t)}{\sqrt{2\eta\cosh s}},
 \qquad 0\le\varepsilon(t)\le5\exp\!\left(-\frac{\pi}{2\eta\cosh s}\right).
$$
此时该本原向量就是唯一最短方向。与其独立的格向量离它的直线至少 $\pi/\ell>1$。若用固定门槛 $0<L\le1$ 看最短向量长度何时小于 $L$，这一段行程的对数间距为 $2\operatorname{arcosh}(L^2/(2\pi\eta))$，条件是 $2\pi\eta<L^2$。不同本原方向的这些开区间互不相交。归一化平方振幅的共振形状因此有精确的单向主项。

还可以直接读出它与穿越折线的对应。若 $k/u=p_n/q_n$ 是 $\beta$ 的第 $n$ 个收敛分数，记 $r_n=[0;a_{n+2},a_{n+3},\ldots]$、$s_n=q_{n-1}/q_n$，其中 $a_j$ 是部分商。连分数余项给 $\eta^{-1}=a_{n+1}+r_n+s_n$。于是当 $a_{n+1}+r_n+s_n\ge2\pi$ 时，在该尺度中心有
$$
 Z(t_*)=\sqrt{\frac{a_{n+1}+r_n+s_n}{2}}\,(1+\varepsilon(t_*)),
 \quad
 0\le\varepsilon(t_*)\le5e^{-\pi(a_{n+1}+r_n+s_n)/2}.
$$
前一份折线研究中的穿越位置比 $R_n=1+a_{n+1}/(1+r_n)$，可以在这里读作归一化能量尖峰平方的尺度。

## 无理性指数给出的增长率

写 $\mu=\mu(\beta)=\mu(\pi)$。有精确关系
$$
 \limsup_{t\downarrow0}\frac{\log Z(t)}{\log(1/t)}
   =\frac{\mu-2}{2\mu},
$$
其中 $\mu=+\infty$ 时约定右端为 $1/2$。证明如下。对每个 $\gamma>\mu-2$，存在 $c_\gamma>0$ 使 $|\beta-k/u|\ge c_\gamma|u|^{-2-\gamma}$。在 $|u|\ge t^{-1/(2+\gamma)}$ 时用向量平方长度的第一项，其余非零 $u$ 用第二项，于是 $\lambda(t)^2\gg_\gamma t^{\gamma/(2+\gamma)}$。最短向量本原，若 $\lambda(t)\le1$，用上一节的高斯估计得 $Z(t)\ll1/\lambda(t)$。否则用计数上界。这给出上极限至多 $\gamma/(2(2+\gamma))$。当 $\mu=+\infty$，只用第一项和 $u=0$ 的情形便有 $\lambda(t)^2\ge2t$ 对充分小的 $t$ 成立。

对每个 $0<\gamma<\mu-2$，取无穷多个 $u$ 满足 $|\beta-k/u|<|u|^{-2-\gamma}$。在其共振时间 $t_*$ 有 $T=\log(1/t_*)\ge(2+\gamma)\log|u|+O(1)$，而 $\ell^2=4u^2 t_*$。上一节下界给 $\log Z(t_*)\ge T/2-\log|u|+O(1)$，所以相应上极限至少为 $\gamma/(2(2+\gamma))$。在 $\mu=2$ 时已有 $Z(t)\ge1$，于是上下界一致。

因而无理性指数等于 $2$ 只推出 $Z(t)\ll_\epsilon t^{-\epsilon}$ 对每个 $\epsilon>0$ 成立。指数仍为零的无界尖峰完全符合这个条件。要完成坏逼近性证明，需要在整数角 $2n^2$ 的两个振荡和上证明统一的平方根抵消。上面的正项恒等式和共振估计准确给出了应排除的尖峰及其时间尺度。
