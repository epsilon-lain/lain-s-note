# CM periods give an exact bridge to pi: Legendre relation and rational elliptic certificates

Date: 2026-11-20  
Status: **exact period identity, exact two-convergent elimination formulas, and exact rational approximant construction**. This is the first direct bridge from a CM elliptic curve to rational approximations of \(\pi\).

## 1. Choose the square-CM curve because its period product is exactly pi

Take
\[
E:\quad y^2=x^3-36x,\qquad P=(12,-36).
\]
The Weierstrass invariants are \(g_2=144\), \(g_3=0\), so \(j(E)=1728\). The lattice is square. Let \(\Omega\) be the positive full real period and \(H\) its Weierstrass zeta quasi-period:
\[
\wp(u+\Omega)=\wp(u),\qquad
\zeta(u+\Omega)=\zeta(u)+H.
\]
The order-four CM automorphism
\[
\tau(x,y)=(-x,iy)
\]
acts on the uniformizing parameter by \(u\mapsto iu\). Hence the second full period and quasi-period are \(i\Omega\) and \(-iH\).

The full-period Legendre relation is
\[
H(i\Omega)-(-iH)\Omega=2\pi i.
\]
Therefore
\[
\boxed{\Omega H=\pi.}
\tag{1}
\]
This is the exact bridge: \(\pi\) is the symplectic area of the square-CM period/quasi-period lattice.

The point \(P\) is non-torsion. The good reductions have
\[
\#E(\mathbb F_5)=8,\quad
\#E(\mathbb F_7)=8,\quad
\#E(\mathbb F_{11})=12,\quad
\#E(\mathbb F_{13})=20.
\]
Prime-to-\(p\) injectivity of reduction bounds the rational torsion order by the gcd \(4\). But
\[
2P=(25/4,35/8)
\]
is not a 2-torsion point, so \(P\) is not torsion.

## 2. Rational elliptic points encode the small period errors

Let \(u\in\mathbb R\) be the elliptic logarithm of \(P\), chosen so that
\[
\wp(u)=12,\qquad \wp'(u)=-72.
\]
Numerically,
\[
\Omega=1.0704505140376155993\ldots,
\]
\[
u=0.2967744354148354546\ldots,\qquad
\theta=\frac{u}{\Omega}=0.2772425549084343844\ldots.
\]
Two consecutive convergents are
\[
\frac{m_1}{n_1}=\frac27,\qquad
\frac{m_2}{n_2}=\frac3{11},
\]
with determinant
\[
D=n_2m_1-n_1m_2=11\cdot2-7\cdot3=1.
\]
Their period errors are
\[
\delta_7=7u-2\Omega=-0.06347998017138\ldots,
\]
\[
\delta_{11}=11u-3\Omega=0.05316724745034\ldots.
\]
For every \(n\), the local parameter
\[
z_n=\xi(nP)=-\frac{x(nP)}{y(nP)}
\]
is rational. Near the identity, the formal logarithm \(L\) satisfies
\[
\delta_n=L(z_n).
\]

## 3. The zeta addition cocycle is rational

Define
\[
R_n=\zeta(nu)-n\zeta(u).
\]
The Weierstrass zeta addition formula gives \(R_n\in\mathbb Q\), with exact recurrence
\[
R_1=0,\qquad
R_2=\frac{3x(P)^2-36}{2y(P)}=-\frac{11}{2},
\]
and, for \(n\ge2\),
\[
\boxed{
R_{n+1}=R_n+
\frac{y(nP)-y(P)}{x(nP)-x(P)}.
}
\tag{2}
\]
Exact rational group arithmetic gives
\[
R_7=-\frac{14566275648}{441041329},
\]
\[
R_{11}=
-\frac{253434733433445222234144}
{28941184902620566636633}.
\]

Since \(nu=m\Omega+\delta_n\), quasi-periodicity gives
\[
n\zeta(u)+R_n=\zeta(\delta_n)+mH.
\]
Eliminate \(\zeta(u)\) between the two convergents. Together with the corresponding elimination for \(\delta_n=nu-m\Omega\), this gives the exact identities
\[
\boxed{\Omega=7\delta_{11}-11\delta_7,}
\tag{3}
\]
\[
\boxed{
H=11R_7-7R_{11}
-11\zeta(\delta_7)+7\zeta(\delta_{11}).
}
\tag{4}
\]
Combining (1), (3), and (4) expresses \(\pi\) entirely in terms of two small elliptic logarithms, two rational cocycles, and two local zeta values.

## 4. The CM symmetry makes both local series sparse

Put \(z=-x/y\). For \(E:y^2=x^3-Ax\), \(A=36\), and \(w=-1/y\), the curve equation gives
\[
w=z^3-Azw^2.
\]
The invariant differential simplifies exactly to
\[
\frac{du}{dz}=\frac1{\sqrt{1+4Az^4}}.
\]
Thus
\[
\boxed{
L(z)=\int_0^z\frac{dt}{\sqrt{1+144t^4}}
=\sum_{r\ge0}
(-1)^r\binom{2r}{r}\frac{36^r}{4r+1}z^{4r+1}.
}
\tag{5}
\]
The first terms are
\[
L(z)=z-\frac{72}{5}z^5+864z^9
-\frac{933120}{13}z^{13}
+\frac{117573120}{17}z^{17}+\cdots.
\]

The Weierstrass zeta expansion is also sparse:
\[
\zeta(u)=\frac1u-\frac{12}{5}u^3
-\frac{432}{175}u^7
-\frac{31104}{17875}u^{11}-\cdots.
\]
After substituting \(u=L(z)\), the normalized local zeta function becomes
\[
\boxed{
z\zeta(L(z))
=1+12z^4-\frac{3888}{7}z^8
+\frac{466560}{11}z^{12}+O(z^{16}).
}
\tag{6}
\]

## 5. A CM-sparse arithmetic filter directly approximates the zeta value

The multiplication series is
\[
\frac{\xi([j]V)}{\xi(V)}
=j+\frac{72}{5}(j^5-j)z^4
+\frac{864}{5}(j^9-6j^5+5j)z^8+O(z^{12}).
\tag{7}
\]
Matching (6) with integer combinations of (7) gives two explicit filters.

### Two-point filter

The weights
\[
(\lambda_2,\lambda_3)=(285,-34),\qquad C=468
\]
satisfy
\[
\sum\lambda_jj=C,\qquad
\sum\lambda_j(j^5-j)=\frac{5C}{6}.
\]
Therefore
\[
\boxed{
Z_7(V):=
\frac{285\,\xi(2V)-34\,\xi(3V)}
{468\,\xi(V)^2}
=\zeta(u)+O(z^7).
}
\tag{8}
\]

### Four-point filter

The weights
\[
(\lambda_2,\lambda_3,\lambda_4,\lambda_5)
=(-15064,17375,-5291,539),\qquad C=3528
\]
satisfy all three exact integer moments
\[
\sum\lambda_jj=3528,
\]
\[
\sum\lambda_j(j^5-j)=2940,
\]
\[
\sum\lambda_j(j^9-6j^5+5j)=-11340.
\]
Consequently
\[
\boxed{
Z_{11}(V):=
\frac{-15064\,\xi(2V)+17375\,\xi(3V)
-5291\,\xi(4V)+539\,\xi(5V)}
{3528\,\xi(V)^2}
=\zeta(u)+O(z^{11}).
}
\tag{9}
\]
For rational \(V\), both filters are rational numbers. This is the desired arithmetic transfer: CM symmetry has turned the local quasi-period function into a rational combination of a few elliptic multiples.

## 6. An explicit rational sequence converging to pi

Truncate (5) after \(r=K\), writing \(L_K\). Truncate the zeta Laurent series after \(u^{4K-1}\), writing \(Z_K\). For the convergent pair \(2/7,3/11\), define
\[
d_{7,K}=L_K(z_7),\qquad d_{11,K}=L_K(z_{11}),
\]
\[
\Omega_K=7d_{11,K}-11d_{7,K},
\]
\[
H_K=11R_7-7R_{11}
-11Z_K(d_{7,K})+7Z_K(d_{11,K}),
\]
and
\[
\boxed{\Pi_K=\Omega_KH_K\in\mathbb Q.}
\tag{10}
\]
The analytic expansions and exact identities above imply \(\Pi_K\to\pi\).

The experiment uses exact Fraction arithmetic for every elliptic point, cocycle, and \(\Pi_K\). Only the error against \(\pi\) is displayed numerically.

| \(K\) | \(\log_{10}|\Pi_K-\pi|\) | Denominator bits |
|---:|---:|---:|
| 0 | -1.3093 | 325 |
| 1 | -4.3037 | 5,322 |
| 2 | -7.2105 | 17,526 |
| 3 | -10.0344 | 36,672 |
| 4 | -12.8122 | 62,836 |
| 5 | -15.5612 | 95,869 |
| 6 | -18.2909 | 135,919 |

A second adjacent pair, \(5/18,28/101\), was also tested. At \(K=3\), its rational approximant has error about \(1.19\times10^{-21}\), with denominator 2,273,673 bits. These are calibration data for the height/error optimization.

## 7. The exact uniform-coverage criterion

Here is the final bridge from a rational approximant sequence to **all** rational approximations.

Let \(r_k=a_k/b_k\) be reduced rationals with strictly increasing denominators \(b_k\to\infty\), and let
\[
\varepsilon_k=|\pi-r_k|.
\]
Suppose there are constants \(c_0>0\) and \(c_1<1/2\) such that
\[
\boxed{
c_0\le b_k^2\varepsilon_k\le c_1,
\qquad
b_{k+1}^2\varepsilon_k\le\frac12
}
\tag{11}
\]
for every \(k\). Then
\[
\boxed{
q^2\left|\pi-\frac pq\right|
\ge \min(c_0,1/2)
}
\tag{12}
\]
for every reduced rational \(p/q\).

**Proof.** Choose \(k\) with \(b_k\le q<b_{k+1}\). If \(p/q=r_k\), then \(q=b_k\) and (11) gives the lower bound \(c_0\). Otherwise distinct reduced rationals satisfy
\[
\left|\frac pq-r_k\right|\ge\frac1{qb_k}.
\]
Hence
\[
q^2\left|\pi-\frac pq\right|
\ge \frac q{b_k}-q^2\varepsilon_k
\ge 1-\frac12=\frac12.
\]
This proves (12). \(\square\)

This criterion gives a concrete target for the CM construction: produce rational period approximants whose normalized errors \(b_k^2\varepsilon_k\) stay in a fixed interval below \(1/2\), while the next denominator is covered by the current error window.

## 8. Next optimization target

The square-CM Legendre identity supplies the exact occurrence of \(\pi\); the two-convergent elimination turns it into rational arithmetic; and the sparse zeta filters supply low-order rational approximations to the quasi-period. The next optimization is to improve the rational height/error ratio of these filters until the covering inequalities (11) become available.

This is the concrete bridge from CM symmetry to a uniform constraint on rational approximations to \(\pi\).
