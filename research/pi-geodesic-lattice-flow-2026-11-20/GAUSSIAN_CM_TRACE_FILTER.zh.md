# Gaussian-CM trace filter: a three-endomorphism order-11 quasi-period approximant

Date: 2026-11-20  
Status: **proved formal expansion, exact rationality, minimal norm-4 support, and exact elliptic experiment**.

## 1. Use the full CM endomorphism ring, not only integer multiplication

Keep the square-CM curve
\[
E:\quad y^2=x^3-36x,\qquad P=(12,-36),
\]
with formal parameter \(z=\xi=-x/y\). Its CM ring is \(\mathbb Z[i]\). The order-four automorphism is
\[
\tau(x,y)=(-x,iy),
\qquad \tau^*z=iz.
\]
For a Gaussian integer \(\alpha=a+bi\), define
\[
[\alpha]Q=[a]Q+[b]\tau(Q).
\]
For rational \(Q\), the point \([\alpha]Q\) is defined over \(\mathbb Q(i)\), and complex conjugation sends it to \([\bar\alpha]Q\).

Define the normalized trace multiplier
\[
\boxed{
T_{a,b}(Q)=
\frac{\xi([a+bi]Q)+\xi([a-bi]Q)}
{2a\,\xi(Q)}
\in\mathbb Q,\qquad a\ne0.
}
\tag{1}
\]
The value is rational because the two summands are Galois conjugates. This is the key extra freedom: the arithmetic filter can use Gaussian endomorphisms while its output remains rational.

## 2. The two nonzero local moments

The formal logarithm is
\[
L(z)=\int_0^z\frac{dt}{\sqrt{1+144t^4}}
=z-\frac{72}{5}z^5+864z^9+O(z^{13}),
\]
and its inverse begins
\[
\exp_F(T)=T+\frac{72}{5}T^5+\frac{864}{5}T^9+O(T^{13}).
\]
Since
\[
[\alpha]_F(z)=\exp_F(\alpha L(z)),
\]
composition gives
\[
[\alpha]_F(z)
=\alpha z+\frac{72}{5}(\alpha^5-\alpha)z^5
+\frac{864}{5}(\alpha^9-6\alpha^5+5\alpha)z^9
+O(z^{13}).
\]
Taking the trace and dividing by \(a z\), we obtain
\[
\boxed{
T_{a,b}(Q)
=1+\frac{72}{5}r_1(a,b)z^4
+\frac{864}{5}r_2(a,b)z^8+O(z^{12}),
}
\tag{2}
\]
where
\[
r_1(a,b)=a^4-10a^2b^2+5b^4-1,
\]
\[
\begin{aligned}
r_2(a,b)={}&a^8-36a^6b^2+126a^4b^4-84a^2b^6+9b^8\\
&-6a^4+60a^2b^2-30b^4+5.
\end{aligned}
\]

The target local quasi-period series is
\[
\boxed{
z\zeta(L(z))
=1+12z^4-\frac{3888}{7}z^8+O(z^{12}).
}
\tag{3}
\]
Thus it suffices to find rational weights \(\lambda_\nu\) and Gaussian multipliers \((a_\nu,b_\nu)\) satisfying
\[
\sum_\nu\lambda_\nu=1,\qquad
\sum_\nu\lambda_\nu r_1(a_\nu,b_\nu)=\frac56,\qquad
\sum_\nu\lambda_\nu r_2(a_\nu,b_\nu)=-\frac{45}{14}.
\tag{4}
\]

## 3. A minimal-norm three-term solution

Use the three Gaussian integers
\[
1,\qquad 1+i,\qquad 2.
\]
Their norms are \(1,2,4\), and their moment pairs are
\[
\begin{array}{c|rr}
\alpha & r_1 & r_2\\ \hline
1 & 0&0\\
1+i&-5&45\\
2&15&165
\end{array}
\]
The unique solution of (4) on these three nodes is
\[
\boxed{
(\lambda_0,\lambda_1,\lambda_2)
=\left(\frac{233}{210},-\frac{13}{105},\frac1{70}\right).
}
\tag{5}
\]
Consequently, with \(T_{1,0}=1\),
\[
\boxed{
F_{\rm CM}(Q)
=\frac{233}{210}
-\frac{13}{105}T_{1,1}(Q)
+\frac1{70}T_{2,0}(Q)
=1+12z^4-\frac{3888}{7}z^8+O(z^{12}).
}
\tag{6}
\]
Therefore
\[
\boxed{
Z_{\rm CM}(Q):=\frac{F_{\rm CM}(Q)}{\xi(Q)}
=\zeta(u)+O(\xi(Q)^{11}).
}
\tag{7}
\]
Every value \(Z_{\rm CM}(Q)\) at a rational point \(Q\) is rational.

This support is optimal in maximum Gaussian norm. Up to units and conjugation, the Gaussian integers of norm at most \(4\) with nonzero real part are exactly \(1,1+i,2\). Their moment points are \((0,0),(-5,45),(15,165)\). The target \((5/6,-45/14)\) lies on none of the three lines through pairs of these points (the corresponding nonzero determinants are \(-150/7,-1300/7,-11650/7\)). Hence one or two trace multipliers cannot work with maximum norm at most \(4\), while (5) uses all three. No Gaussian trace filter of maximum norm below \(4\) can match both local moments.

## 4. Insert the filter into the exact Legendre bridge

Use the same consecutive convergents \(2/7,3/11\) of the elliptic logarithm ratio for \(P=(12,-36)\). Put
\[
z_7=\xi(7P),\qquad z_{11}=\xi(11P),
\]
and let \(R_n=\zeta(nu)-n\zeta(u)\) be the rational addition cocycle. The exact identities are
\[
\Omega=7L(z_{11})-11L(z_7),
\]
\[
H=11R_7-7R_{11}-11\zeta(L(z_7))+7\zeta(L(z_{11})),
\qquad \pi=\Omega H.
\]
Truncate only the formal logarithm at
\[
L_1(z)=z-\frac{72}{5}z^5
\]
and replace each local zeta value by the rational filter (7):
\[
\widehat\Omega=7L_1(z_{11})-11L_1(z_7),
\]
\[
\widehat H=11R_7-7R_{11}
-11Z_{\rm CM}(7P)+7Z_{\rm CM}(11P).
\]
Then
\[
\boxed{\widehat\pi_{\rm CM}=\widehat\Omega\,\widehat H\in\mathbb Q.}
\tag{8}
\]

## 5. Exact rational experiment

All elliptic group operations, Gaussian endomorphism evaluations, traces, filters, and the rational number (8) are computed exactly. The error against \(\pi\) is independently enclosed using Machin's alternating arctangent formula.

For the two local zeta values:

| Input | Numerator bits | Denominator bits | \(\log_{10}|Z_{\rm CM}-\zeta|\) |
|---|---:|---:|---:|
| \(7P\) | 446 | 442 | \(-8.04147457\) |
| \(11P\) | 1092 | 1088 | \(-8.88890783\) |

The resulting rational approximation is:

| Construction | Reduced denominator bits | \(\log_{10}|\widehat\pi-\pi|\) |
|---|---:|---:|
| Previous fixed integer three-point filter | 7270 | \(-3.25519747\) |
| **Gaussian-CM trace filter (5)** | **2500** | **\(-6.38634348\)** |

Thus the CM-ring filter cuts the denominator bit length by about \(65.6\%\), while improving the absolute error by a factor of about \(1.35\times10^3\) in this same finite Legendre experiment.

The improvement is not obtained by increasing the number of local moments: it comes from using the endomorphism \([1+i]\) and taking its Galois trace. The three support norms are only \(1,2,4\).

## 6. What to optimize next

The new design space is the Gaussian integer lattice:
\[
(a,b)\in\mathbb Z^2,\qquad N(a+bi)=a^2+b^2.
\]
For each \((a,b)\), the two local moments \(r_1,r_2\) are explicit integer polynomials. The next optimization can therefore search supports in increasing Gaussian norm, solve the exact rational moment system, and score the resulting filters by their **actual reduced rational height** on the chosen elliptic return points.

This is a more faithful optimization than restricting to positive integer multiplication indices. It also suggests testing other CM orders, where the endomorphism ring may offer additional low-norm directions.

## References and reproducibility

- Square-CM Legendre bridge and exact period/cocycle identities: CM_LEGENDRE_PERIOD_PI_BRIDGE.md.
- Exact Gaussian trace filter and Machin-interval verification: verify_gaussian_cm_trace_filter.py.
- Exact output summary: gaussian_cm_trace_filter_results.json.
