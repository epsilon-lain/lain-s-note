# Quantitative audit of the shifted six-point tail construction

Date: 2026-11-19  
Status: **proved quantitative auxiliary estimate + explicit scope limitation**.

This note sharpens the shifted-convergent construction. It gives a uniform \(O(q_j^{-1})\) rate for the six tail ratios, and isolates why that fact alone is not yet a Diophantine approximation theorem for \(\pi\).

## 1. A uniform endpoint estimate

For \(0<x<1\), let
\[
r_m(x)=\int_0^x\frac{t^{2m}}{1+t^2}\,dt.
\]
The integration-by-parts identity gives
\[
r_m(x)=\frac{x^{2m+1}}{(2m+1)(1+x^2)}
\left(1+\rho_m(x)\right),
\]
where
\[
0\le \rho_m(x)\le
\frac{2x^2(1+x^2)}{2m+3}.
\tag{1}
\]
Write \(a=2/5,\ b=3/7\), and
\[
\lambda=\frac{\log a}{\log b},\quad
C=\lambda\frac ab\frac{1+b^2}{1+a^2},\quad
d_*=\frac{\log C}{2\log b}.
\]
For \(d=n-\lambda m\), substitution in (1) gives the exact factorization
\[
\frac{r_m(a)}{r_n(b)}
=C b^{-2d}
\left(1+\frac{2d+1-\lambda}{\lambda(2m+1)}\right)
\frac{1+\rho_m(a)}{1+\rho_n(b)}.
\tag{2}
\]
Because \(Cb^{-2d_*}=1\), this becomes
\[
\frac{r_m(a)}{r_n(b)}
=b^{-2(d-d_*)}
\left(1+\frac{2d+1-\lambda}{\lambda(2m+1)}\right)
\frac{1+\rho_m(a)}{1+\rho_n(b)}.
\tag{3}
\]

## 2. Rate along the six-point family

Let \(p_j/q_j\) be convergents to \(\lambda\). Choose \(c_j\) nearest to \(q_jd_*\), and use the three offsets \(c_j-1,c_j,c_j+1\). Solve
\[
p_j u_c\equiv-c\pmod{q_j},\qquad 0\le u_c<q_j,
\]
and put \(x_c=u_c+q_j\). If necessary replace the \(+\) point's \(x_{c_j+1}\) by \(x_{c_j+1}+q_j\) to make
\[
\kappa_j=\frac{2x_{c_j}-x_{c_j-1}-x_{c_j+1}}{q_j}
\]
nonzero. For each offset use the paired points
\[
(m,n)=(x_c,n_c),\quad(x_c+q_j,n_c+p_j),
\quad n_c=(p_jx_c+c)/q_j.
\]

For all six points, for sufficiently large \(j\),
\[
q_j\le m\le4q_j,\qquad
|c|\le |d_*|q_j+2,
\]
and
\[
d-d_*=
\left(\frac{p_j}{q_j}-\lambda\right)m
+\left(\frac c{q_j}-d_*\right).
\]
The convergent estimate \(|p_j/q_j-\lambda|<q_j^{-2}\), together with
\(|c_j/q_j-d_*|\le(2q_j)^{-1}\), yields
\[
|d-d_*|\le\frac{4+3/2}{q_j}=\frac{11}{2q_j}.
\tag{4}
\]
Also \(d\) remains in a fixed bounded interval and \(m,n\asymp q_j\). Equations (1), (3), and (4) therefore imply a constant \(K>0\), independent of \(j\) and of the six selected points, such that
\[
\boxed{
\max_{1\le i\le6}
\left|\frac{r_{m_i}(2/5)}{r_{n_i}(3/7)}-1\right|
\le \frac K{q_j}
}
\tag{5}
\]
for all sufficiently large \(j\). This is stronger than mere convergence, although the constant \(K\) is not optimized here.

## 3. The integer determinant remains bounded

The six quadratic evaluation rows are
\[
(1,m,n,m^2,mn,n^2).
\]
Using \(c=q_jn-p_jm\), the coordinate change on quadratic features has determinant \(q_j^4\). The three offsets are consecutive and each has two points separated by \((q_j,p_j)\). Direct expansion gives
\[
\det A_j=\pm4\kappa_j.
\tag{6}
\]
The residue relation
\[
x_c\equiv-p_j^{-1}c\pmod{q_j}
\]
shows that \(\kappa_j\in\mathbb Z\). After the optional one-step adjustment, the six \(x\)-coordinates lie in \([q_j,3q_j-1]\); hence
\[
0<|\kappa_j|<4.
\]
Consequently
\[
\boxed{|\det A_j|\in\{4,8,12\}.}
\tag{7}
\]
In particular, this exact integer determinant does not tend to zero. The \(O(q_j^{-1})\) tail cancellation has not yet been inserted as a determinant factor.

## 4. Why this does not yet approximate \(\pi\) well enough

The rational number
\[
R_{m,n}=4\bigl(T_m(2/5)+T_n(3/7)\bigr)
\]
satisfies
\[
\frac\pi4-\frac{R_{m,n}}4
=(-1)^m\bigl(r_m(2/5)-r_n(3/7)\bigr)
\]
when \(m,n\) have opposite parity. The estimate (5) controls the *relative difference* of two tails, but not the denominator cost of \(R_{m,n}\) and not a lower bound for that difference.

A denominator upper bound alone is not enough to show \(q^2|\pi-p/q|\to0\); one needs a lower bound for the reduced denominator together with a sufficiently strong upper bound for the actual linear-form error. Conversely, the relative estimate (5) is only \(O(q_j^{-1})\), and the common tail size is exponentially small in \(m\), while the rational truncation denominator contains large prime-power contributions from both \(5\) and \(7\). The present construction has not shown that the error beats the square of the *actual reduced denominator*.

Thus the valid output is (5)–(7), not a new rational approximation record and not a proof of bad approximability.

## 5. Next mathematical target

The next useful construction must make the determinant itself sensitive to the signed tail difference
\[
r_m(2/5)-r_n(3/7),
\]
rather than merely selecting points for which the ratio of the two tails is close to 1. A candidate row/column operation must satisfy both:

1. its determinant contribution is an explicitly controlled small factor from the exact tail difference or a Padé endpoint;
2. after clearing denominators, the resulting integer determinant has an independent nonzero certificate (for example, a valid modular rank computation).

Without both, a small analytic row is not an integer contradiction.
