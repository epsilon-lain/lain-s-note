# Exact tail-column determinant: the six-point bridge and its denominator budget

Date: 2026-11-19  
Status: **exact determinant identity proved; current denominator and tail bounds do not certify the required inequality.**

This note makes the signed two-arctangent tail difference an actual determinant column. It also audits the denominator needed to turn that real determinant into an integer. This is the missing logical bridge between “small tail cancellation” and “integer determinant”; the audit shows why the raw Taylor truncations do not yet close it.

## 1. The exact tail column

Let
\[
a=\frac25,\qquad b=\frac37,\qquad
T_m(x)=\sum_{j=0}^{m-1}(-1)^j\frac{x^{2j+1}}{2j+1},
\]
and
\[
f_{m,n}:=\frac\pi4-T_m(a)-T_n(b).
\]
Since \(\arctan a+\arctan b=\pi/4\),
\[
f_{m,n}=R_m(a)+R_n(b)
=(-1)^m r_m(a)+(-1)^n r_n(b),
\]
where \(r_m(x)=\int_0^x t^{2m}/(1+t^2)\,dt\). If \(m,n\) have opposite parity, then
\[
\boxed{|f_{m,n}|=|r_m(a)-r_n(b)|.}
\tag{1}
\]
Thus the small signed tail difference is exactly a linear form \(\pi/4-\text{rational}\), not merely an approximate ratio.

Choose six integer pairs \((m_i,n_i)\), and define the six-by-six real matrix
\[
M_{i,\cdot}=
(1,m_i,n_i,m_i^2,m_in_i,f_{m_i,n_i}).
\tag{2}
\]
Let \(R_i=T_{m_i}(a)+T_{n_i}(b)\), and let \(L\) be any positive common denominator of the six rationals \(R_i\). Since
\[
f_{m_i,n_i}=\pi/4-R_i,
\]
and the first column of \(M\) is the all-ones column, multilinearity gives
\[
\det M=-\det
\begin{pmatrix}
1&m_1&n_1&m_1^2&m_1n_1&R_1\\
\vdots&&&&&\vdots\\
1&m_6&n_6&m_6^2&m_6n_6&R_6
\end{pmatrix}.
\tag{3}
\]
The right side is rational. Consequently
\[
\boxed{L\det M\in\mathbb Z.}
\tag{4}
\]
This is an exact integer determinant whose last real column is the signed tail error. It is the desired algebraic interface; no approximation has been used in (3)–(4).

## 2. A rigorous smallness bound

Suppose
\[
1\le m_i,n_i\le H,\qquad |f_{m_i,n_i}|\le\delta\le1.
\]
The six column Euclidean norms in (2) are bounded respectively by
\[
\sqrt6,\quad \sqrt6H,\quad\sqrt6H,\quad
\sqrt6H^2,\quad\sqrt6H^2,\quad\sqrt6\delta.
\]
Hadamard's inequality therefore gives
\[
\boxed{|\det M|\le216H^6\delta.}
\tag{5}
\]
(The exponent is \(1+1+2+2=6\) from the four nonconstant polynomial columns \(m,n,m^2,mn\); the constant column contributes no power of \(H\).)

Combining (4) and (5), if
\[
216\,L H^6\delta<1,
\tag{6}
\]
then \(L\det M=0\), hence \(\det M=0\). To get a contradiction one must independently prove \(L\det M\ne0\), for example by reducing the rational determinant in (3) modulo a prime not dividing the chosen clearing denominator. Thus the exact closure conditions are
\[
L\det M\in\mathbb Z\setminus\{0\},
\qquad 216LH^6\delta<1.
\tag{7}
\]

## 3. Denominator cost for the raw Taylor truncations

For \(x=s/t\) in lowest terms,
\[
T_m(s/t)=
\sum_{j=0}^{m-1}(-1)^j
\frac{s^{2j+1}}{t^{2j+1}(2j+1)}.
\]
A valid common denominator is
\[
t^{2m-1}\Lambda_m,\qquad
\Lambda_m=\operatorname{lcm}(1,3,5,\ldots,2m-1).
\tag{8}
\]
Let \(H=\max_i(m_i,n_i)\). A common denominator for all six \(R_i\) is therefore bounded by
\[
L_H=5^{2H-1}7^{2H-1}\Lambda_H.
\tag{9}
\]
This is a deliberately simple upper bound; actual cancellation can reduce it, but cannot be assumed without an exact gcd proof.

The least common multiple satisfies
\[
\log\Lambda_H\le\log\operatorname{lcm}(1,2,\ldots,2H)
\le 2H\log4,
\]
where the last inequality follows from the elementary bound
\(\operatorname{lcm}(1,\ldots,N)\le4^N\). Hence
\[
\log L_H
\le(2H-1)(\log5+\log7)+2H\log4.
\tag{10}
\]
The bound is exponential in \(H\), with exponent
\(2\log35+2\log4\approx 9.128\).

By contrast, for the shifted six-point family, \(m,n\asymp q_j\) and the matched tail size is
\[
r_m(2/5)\asymp \frac{(2/5)^{2m}}m,
\]
whose exponential rate in \(m\) is \(2\log(5/2)\approx1.833\). The relative mismatch established in the shifted construction is \(O(q_j^{-1})\), so the absolute signed error is at most a polynomial factor times this tail. Even this optimistic error scale is exponentially much larger than \(1/L_H\), because the crude but valid clearing-denominator cost (10) has a much larger exponential rate.

Therefore the present estimates do **not certify** the sufficient condition (6) along the current shifted family: the available upper bound for the clearing denominator is far too coarse relative to the available tail estimate. This is not a proof that the route fails. In particular, because (9) is only an upper bound, the actual reduced common denominator may be much smaller; and the actual signed tail difference may be smaller than the current upper estimate. A valid impossibility claim would require a lower bound for the product \(L\delta\), which has not been proved.

## 4. A correction to the shifted quadratic determinant

For the six-point family with three consecutive offsets \(c_-=c_0-1,c_0,c_0+1\), choose representatives \(x_-,x_0,x_+\in[q,2q-1]\) satisfying
\[
p x_c+c\equiv0\pmod q.
\]
Then
\[
\kappa=\frac{2x_0-x_--x_+}{q}\in\mathbb Z.
\]
Since
\[
|2x_0-x_--x_+|\le2(q-1)<2q,
\]
we have \(\kappa\in\{-1,0,1\}\). If \(\kappa=0\), replace \(x_+\) by \(x_++q\); this changes \(\kappa\) to \(-1\). Hence the construction can always be normalized so that
\[
\boxed{|\kappa|=1.}
\tag{11}
\]
For the quadratic feature rows \((1,m,n,m^2,mn,n^2)\), with paired points \((x_c,n_c)\) and \((x_c+q,n_c+p)\), the exact determinant identity is
\[
\boxed{\det A=-4\kappa}
\tag{12}
\]
for the stated row order (reordering rows changes only the sign). Thus the determinant is exactly \(\pm4\), not merely one of \(\pm4,\pm8,\pm12,\pm16\). It is nonzero modulo every odd prime.

## 5. What the next construction must improve

The determinant (3) is the correct way to make the signed analytic cancellation appear in a determinant. The remaining difficulty is no longer conceptual; it is quantitative and arithmetic:

1. Replace the raw Taylor truncations by rational/Padé approximants whose **actual reduced common denominator** is much smaller than (9), with a proved bound.
2. Or construct a multi-row determinant in which the large denominator factors cancel algebraically before taking absolute values.
3. Prove a modular nonzero certificate for that exact reduced integer determinant.
4. Prove a tail-difference bound strong enough that the determinant is less than the reciprocal of the true clearing denominator.

The present note proves (1)–(5), identifies the exact sufficient condition (7), and records why the first naive denominator choice is too expensive. It does not prove \(\pi\) badly approximable, nor does it prove that all Padé or joint-gcd improvements are impossible.

## 6. Exact finite replay

The accompanying [verify_exact_tail_column.py](verify_exact_tail_column.py) computes the six rational truncations with `Fraction`, forms their exact common denominator, builds the integer matrix obtained by clearing only the final column, and checks its determinant by Bareiss elimination. For the six-point sample in `TWO_ARCTAN_SIX_POINT_CANDIDATE.md`, it produces a nonzero integer determinant and verifies it modulo 5. This is a finite certificate only; it does not establish a uniform lower bound or an infinite-family contradiction.
