# Shifted six-point construction along convergents of the tail-slope

Date: 2026-11-19  
Status: **proved construction lemma for the two-arctangent tail model**. It produces six integer points with all six tail ratios tending to 1 and a nonzero bounded integer determinant. It does not prove bad approximability of \(\pi\).

## 1. Exact tail asymptotics

For \(0<x<1\), define
\[
r_m(x)=\int_0^x\frac{t^{2m}}{1+t^2}\,dt.
\]
Integration by parts gives
\[
r_m(x)=\frac{x^{2m+1}}{(2m+1)(1+x^2)}
+\frac{2}{2m+1}\int_0^x\frac{t^{2m+2}}{(1+t^2)^2}\,dt.
\]
Consequently
\[
\frac{x^{2m+1}}{(2m+1)(1+x^2)}
\le r_m(x)
\le
\frac{x^{2m+1}}{(2m+1)(1+x^2)}
\left(1+\frac{2x^2(1+x^2)}{(2m+3)}\right).
\tag{1}
\]
The upper bound follows by replacing \((1+t^2)^{-2}\) by \(1\) in the remainder integral and dividing by the leading term. In particular, uniformly for \(x\) in a fixed compact subinterval of \((0,1)\),
\[
r_m(x)=\frac{x^{2m+1}}{(2m+1)(1+x^2)}
\left(1+O(m^{-1})\right).
\tag{2}
\]

Put
\[
a=\frac25,\qquad b=\frac37,\qquad
\lambda=\frac{\log a}{\log b}.
\]
The number \(\lambda\) is irrational: if \(\lambda=u/v\in\mathbb Q\) with positive integers \(u,v\), then \(a^v=b^u\), impossible by unique prime factorization.

Using (2), if \(m\to\infty\), \(n/m\to\lambda\), and \(n-\lambda m\to d\), then
\[
\frac{r_m(a)}{r_n(b)}
\longrightarrow
C\,b^{-2d},
\qquad
C=\lambda\frac ab\frac{1+b^2}{1+a^2}.
\tag{3}
\]
Define the unique real shift
\[
d_*=\frac{\log C}{2\log b}.
\tag{4}
\]
Then the limiting ratio in (3) is exactly 1 when \(d=d_*\).

For the chosen \(a,b\),
\[
\lambda\approx1.08142693934742458,\quad
C\approx1.0299304184261184,\quad
d_*\approx-0.01740311548764688.
\]
These decimals are orientation only; the definitions (3)–(4) are exact.

## 2. Six integer points from convergents

Let \(p_j/q_j\) be the continued-fraction convergents of \(\lambda\). Then
\[
\left|\frac{p_j}{q_j}-\lambda\right|<\frac1{q_j^2},
\qquad q_j\to\infty,\qquad \gcd(p_j,q_j)=1.
\tag{5}
\]

Choose an integer \(c_j\) nearest to \(q_jd_*\), and set
\[
c_{j,-}=c_j-1,\quad c_{j,0}=c_j,\quad c_{j,+}=c_j+1.
\]
For each \(c\) in this triple, choose \(u_c\in\{0,1,\ldots,q_j-1\}\) satisfying
\[
p_j u_c\equiv-c\pmod {q_j}.
\tag{6}
\]
This is possible and unique because \(p_j\) is invertible modulo \(q_j\). Let
\[
x_c=u_c+q_j.
\]
Thus \(q_j\le x_c\le2q_j-1\). Define
\[
\kappa_j=\frac{2x_{c,0}-x_{c,-}-x_{c,+}}{q_j}\in\mathbb Z.
\tag{7}
\]
It is an integer because (6) implies
\(x_c\equiv-p_j^{-1}c\pmod {q_j}\), so the numerator in (7) is divisible by \(q_j\).

If \(\kappa_j=0\), replace \(x_{c,+}\) by \(x_{c,+}+q_j\); this changes \(\kappa_j\) to \(-1\). Hence we may and do arrange
\[
1\le|\kappa_j|\le 4.
\tag{8}
\]
(The upper bound follows from \(q_j\le x_c\le3q_j-1\) after the possible adjustment.)

For each \(c\in\{c_{j,-},c_{j,0},c_{j,+}\}\), take the two points
\[
(m,n)=(x_c,n_c),\qquad (x_c+q_j,n_c+p_j),
\quad
n_c=\frac{p_jx_c+c}{q_j}.
\tag{9}
\]
Equation (6) makes every \(n_c\) an integer. Since \(x_c\asymp q_j\), \(p_j/q_j\to\lambda>1\), and \(c_j/q_j\to d_*\), all six coordinates are positive for sufficiently large \(j\), and all \(m,n=O(q_j)\).

For every point in (9),
\[
n-\lambda m
=\left(\frac{p_j}{q_j}-\lambda\right)m+\frac c{q_j}
=d_*+O(q_j^{-1}),
\tag{10}
\]
uniformly over the six points, by (5), \(m=O(q_j)\), and
\(c_j/q_j=d_*+O(q_j^{-1})\). Also \(m,n\to\infty\) and \(n/m\to\lambda\). Therefore (3) gives
\[
\boxed{
\max_{\text{six points}}
\left|\frac{r_m(2/5)}{r_n(3/7)}-1\right|\longrightarrow0.
}
\tag{11}
\]
This is a genuine infinite family of increasingly accurate paired-tail cancellations, not merely a finite numerical fit.

## 3. Exact nonzero determinant

For each point use the quadratic evaluation row
\[
(1,m,n,m^2,mn,n^2),
\]
and let \(A_j\in M_6(\mathbb Z)\) be the resulting matrix, with the rows ordered by the three \(c\)-levels and the two points on each level.

Change the coordinate from \(n\) to \(c=q_jn-p_jm\). The induced change on the six-dimensional quadratic feature basis has determinant \(q_j^4\). In the new coordinates the rows are
\[
(1,m,c,m^2,mc,c^2).
\]
A direct determinant expansion for the three levels \(c_-=c_0-1,c_0,c_0+1\), with paired \(m\)-coordinates \(x_c,x_c+q_j\), gives
\[
\det A_j
=\pm4\,\frac{2x_{c,0}-x_{c,-}-x_{c,+}}{q_j}
=\pm4\kappa_j.
\tag{12}
\]
Thus (8) yields the exact certificate
\[
\boxed{
\det A_j\in\{\pm4,\pm8,\pm12,\pm16\}
\quad\text{for every sufficiently large }j.
}
\tag{13}
\]
In particular the determinant is nonzero, and it is nonzero modulo any prime not dividing \(4\kappa_j\). The determinant does not grow with the point height.

## 4. What this construction does and does not achieve

It improves the finite six-point sample in a specific way:

- all six paired arctangent tails have ratio tending to 1;
- the six points have integer coordinates of size \(O(q_j)\);
- the quadratic evaluation determinant is exactly nonzero and uniformly bounded.

But it does **not** create the required integer contradiction. The line containing the points is only approximated to width \(O(q_j^{-1})\), while the point height is \(O(q_j)\). The generic quadratic six-point determinant bound is therefore not forced below 1. More importantly, (12) shows the exact determinant is already a small nonzero integer: the tail-cancellation information has not yet entered the determinant as an additional small factor.

The missing ingredient is a new column or row functional that depends on the exact tail (or Padé endpoint) and vanishes on the limiting shifted line, while retaining an independent modular nonvanishing certificate. Merely adding more polynomial monomials in \(m,n\) will not automatically solve this: the exact lattice determinant must be recalculated, and the denominator cost must be included.

**Conclusion:** (11)–(13) are a proved construction lemma for the auxiliary two-arctangent tail geometry. They do not imply any approximation statement about \(\pi\) beyond the identity \(\arctan(2/5)+\arctan(3/7)=\pi/4\), and they do not prove that \(\pi\) is badly approximable.
