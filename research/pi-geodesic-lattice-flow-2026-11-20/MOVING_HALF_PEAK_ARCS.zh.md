# 半高峰允许移动多远

记 $S_Q$ 为圆周上分母在 $[Q,2Q)$ 的既约分数，$H,Q$ 为正整数，$L=HQ^2$。令 $K_L(y)=|\sum_{j=0}^{L-1}e^{2\pi i j y}/L|^2$，$F_\beta=\sum_{x\in S_Q}K_L(\beta-x)$，$E_{Q,H}=\{\beta:F_\beta\ge1/2\}$。用 $\|\cdot\|$ 表示圆周距离，$B(x,d)$ 为圆周闭球。取 $c_*$ 为 $(0,1)$ 中 $\operatorname{sinc}^2(\pi c_*)=1/2$ 的唯一根，$\operatorname{sinc}(z)=\sin z/z$ 连续延拓到零。

**半高集合夹逼。** 对所有 $H\ge20$ 有
$$
\bigcup_{x\in S_Q} B(x,c_*/L)\ \subseteq\ E_{Q,H}\ \subseteq\
\bigcup_{x\in S_Q} B\big(x,(c_*+31/H^2)/L\big).
$$
左侧包含对所有 $H\ge1$ 成立。证明时，令 $t=L\operatorname{dist}(\beta,S_Q)$。最近核值是 $\operatorname{sinc}^2(\pi t)/\operatorname{sinc}^2(\pi t/L)$。当 $t\le c_*$，分母不超过 1，且其余核都非负。对右侧，由[此前的 Farey 间距估计](FIXED_ANGLE_WINDOWS_AND_SPARSE_PEAKS.zh.md)，其余核的贡献小于 $27/H^2$，并且 $F_\beta\ge1/2$ 蕴含 $t<1$。用 $\sin u/u\ge1-u^2/6$，得 $\operatorname{sinc}^{-2}(\pi t/L)\le (1-5/(3L^2))^{-2}\le1+4/L^2$。若 $1/2\le t<1$，第一瓣的单调性给出 $F_\beta\le(4/\pi^2)(1+4/L^2)+27/H^2<1/2$。故 $t<1/2$，并有 $\operatorname{sinc}^2(\pi t)\ge1/2-31/H^2$。

当 $c_*\le t\le1/2$，置 $u=\pi t$。由 $\sin v\ge 2v/\pi$ 在 $[0,\pi/2]$ 上成立，故 $\sin u-u\cos u=\int_0^u v\sin v\,dv\ge 2u^3/(3\pi)$，又 $\operatorname{sinc}(u)\ge2/\pi$。微分得到 $-\frac d{dt}\operatorname{sinc}^2(\pi t)\ge8t/3\ge1$，末步用 $c_*>3/8$，可用 $\operatorname{sinc}(3\pi/8)\ge1-(3\pi/8)^2/6>1/\sqrt2$ 核对。积分便得 $t\le c_*+31/H^2$。

**某个中心的最小移动。** 给定任意中心 $\alpha$，设 $D_Q=Q^2\operatorname{dist}(\alpha,S_Q)$，$u_+=\max(u,0)$。球并的距离直接给出
$$
(HD_Q-c_*-31/H^2)_+\le L\operatorname{dist}(\alpha,E_{Q,H})\le(HD_Q-c_*)_+.
$$
给定固定 $C\ge0$，允许移动 $C/L$ 后仍可见半高峰的中心为 $E_{Q,H}$ 的这个半径的闭邻域。固定 $C$ 且 $H$ 足够大时，两侧夹逼中扩张的球各自互不相交，因为既约节点相距超过 $1/(4Q^2)$。因此随 $H,Q\to\infty$，该中心集合的 Lebesgue 长度乘 $H$ 趋于 $18(C+c_*)/\pi^2$，这里用了 $|S_Q|/Q^2\to9/\pi^2$。该长度计算同时考虑全部相位。

**巨型纪录。** 取 $\alpha=1/\pi$，令 $r$ 为巨型纪录分数倒数的分母，锚定窗口 $Q=r$，$A=845188380288$。按[纪录报告](PI_RECORD_PEAK_AND_GAP.zh.md)，$1/(A+2)<D_r<1/A$。选 $C=1/50$，这个数是专门为给出可无歧义核验的移动档位选的。令 $h=391276977722$，则
$$
h/A < C+c_*,\qquad (h+1)/(A+2)> C+c_*+31/(h+1)^2 .
$$
[有理数验证器](verify_moving_record_threshold.py)把 $c_*$ 严格夹在 $0.44294647068945233$ 和 $0.44294647068945235$ 之间，直接验证这两条不等式。因此对于每个整数 $1\le H\le h$，都存在 $\|\beta-1/\pi\|\le1/(50Hr^2)$ 使 $F_\beta\ge1/2$。对于每个 $H>h$，该邻域内的每个观察点都低于半高。这个有限档位以公开的 $A$ 为输入，完全不需要事先重建数百亿位的分母。原先 $C=0$ 的未决档位仍依赖邻项。

在巨峰后的某个二进分母窗口已知 $D_Q>A/18$，本结论还给出 $L\operatorname{dist}(1/\pi,E_{Q,H})\ge HA/18-c_*-31/H^2$。也就是说，要从中心移到半高峰，在那个窗口至少须走将近 $HA/18$ 个核分辨率单位 $1/L$。要判决坏逼近性仍须控制无穷多个独立分母尺度。
