# Prime blocks, digit residues and a prime-square lifting criterion

2026-11-02. Work from repository checkpoint
`dff048f73da83e63978e1c01576d38ea71af1546`. All new statements concern the
specific integer recurrence below, not arbitrary rational approximations.

Let
\[
A_0=A_1=1,\quad B_0=0,\quad B_1=1,\qquad
X_{n+1}=(2n+1)X_n+n^2X_{n-1},
\]
and \(g_n=\gcd(A_n,4B_n)\). Set \(a_n=A_n/n!\), \(b_n=B_n/n!\).
For odd primes the factor 4 does not affect valuations.

## 1. Exact generating functions and a Lucas property

Formal power series satisfy
\[
f(z):=\sum_{n\ge0}a_nz^n=(1-2z-z^2)^{-1/2},
\qquad
\sum_{n\ge0}b_nz^n=f(z)\int_0^z f(s)\,ds.
\tag{1}
\]
Indeed the recurrence gives
\((1-2z-z^2)f'=(1+z)f\); the second series satisfies the same differential
equation with an extra constant 1 on its right side. Variation of constants,
using \((1-2z-z^2)f^2=1\), proves the second identity. Consequently
\[
\boxed{b_n=\sum_{j=1}^n\frac{a_{n-j}a_{j-1}}j.}
\tag{2}
\]

Writing \(\operatorname{CT}\) for the coefficient of \(x^0\), one has
\[
a_n=\operatorname{CT}\left(1+x+\frac1{2x}\right)^n.
\tag{3}
\]
For example, summing the geometric series in the Laurent polynomial gives
\(1/\sqrt{(1-z)^2-2z^2}\), which is (1). Thus \(a_n\in\mathbb Z[1/2]\).
For every odd prime \(p\), every \(m\ge0\), and \(0\le r<p\), Frobenius and
the exponent support \([-r,r]\) of the last factor give the Lucas identity
\[
\boxed{a_{mp+r}\equiv a_ma_r\pmod p.}
\tag{4}
\]
Only the zero exponent of that last factor is a multiple of \(p\); hence
no other exponents contribute to the constant term. Iterating (4) expresses
\(a_n\bmod p\) as the product of its base-\(p\) digit coefficients.

For
\(s_p=(-1)^{(p-1)/2}\), the coefficients below degree \(p\) in \(f\bmod p\)
equal those in \((1-2z-z^2)^{(p-1)/2}\). Therefore
\[
a_{p-1}\equiv s_p\pmod p,
\qquad
a_{p-1-r}\equiv s_p(-1)^ra_r\pmod p.
\tag{5}
\]
The second identity follows by reversing this polynomial and using
\(z^2(1+2/z-1/z^2)=-(1-2z-z^2)\).

## 2. The exact leading denominator residue for every odd prime

Let \(n\ge p\), \(e=\lfloor\log_p n\rfloor\), and write
\(n=Np^e+r\), with \(1\le N<p\) and \(0\le r<p^e\). Formula (2) shows
\(v_p(b_n)\ge-e\). Multiplying it by \(p^e\), only the terms
\(j=\ell p^e\) survive modulo \(p\). Applying (4) repeatedly yields
\[
\boxed{p^e b_n\equiv s_p^e a_r b_N\pmod p.}
\tag{6}
\]
Here \(b_N\) is \(p\)-integral because \(N<p\), so the right side is an
ordinary finite-field residue.

In particular, if \(a_rb_N\not\equiv0\pmod p\), then
\[
\boxed{v_p(g_n)=v_p(n!)-\lfloor\log_p n\rfloor.}
\tag{7}
\]
If that residue vanishes, one can only conclude
\[
v_p(g_n)\ge v_p(n!)-\lfloor\log_p n\rfloor+1.
\tag{8}
\]
Higher lifting remains necessary. Thus the denominator deficit at its
largest possible order is detected by a fully explicit digit residue.
The factor \(a_r\) itself factors into the lower digit values by (4).

For every \(p>n\), \(p\nmid g_n\): the integer Wronskian
\[
A_nB_{n-1}-A_{n-1}B_n=(-1)^n((n-1)!)^2
\tag{9}
\]
would otherwise vanish modulo \(p\). Also \(p\nmid g_p\), by (6).

## 3. An exact transfer block modulo an odd prime

Let \(M_k=\begin{pmatrix}2k+1&k^2\\1&0\end{pmatrix}\). The product order
is fixed by \((X_{k+1},X_k)^T=M_k(X_k,X_{k-1})^T\). Equations (4)--(6)
and Wilson's theorem give
\[
A_p\equiv0,\qquad B_p\equiv A_{p-1}\equiv-s_p\pmod p.
\]
Also \(B_{p-1}\equiv0\pmod p\). To see the last assertion, pair in (2)
the term indexed \(j\) with the term indexed \(p-j\); their numerators
coincide and their inverse denominators sum to zero modulo \(p\).
The fundamental matrices therefore give
\[
\boxed{M_{p-1}\cdots M_1
\equiv s_p\begin{pmatrix}-1&1\\0&-1\end{pmatrix}\pmod p,}
\tag{10}
\]
and, after multiplying by \(M_0\),
\[
\boxed{M_{p-1}\cdots M_0
\equiv\begin{pmatrix}0&0\\-s_p&0\end{pmatrix}\pmod p.}
\tag{11}
\]
The full block is rank-one nilpotent. The scalar in (10) is governed by
\((-1/p)\), not by \((2/p)\). Since the recurrence coefficients are
periodic modulo \(p\), applying (11) from index \(p\) to index \(2p\)
shows \(X_{2p}\equiv0\). The following singular step gives
\(X_{2p+1}\equiv X_{2p}\equiv0\); the recurrence then propagates zero.
Thus \(p\mid g_n\) for every \(n\ge2p\). This argument starts at
index \(p\), so it does not assume the false initial relation \(B_1=B_0\).

For the incomplete first block, (6) and \(b_1=1\) show
\[
\boxed{p\mid g_{p+r}\quad\Longleftrightarrow\quad
a_r\equiv0\pmod p\quad(0\le r<p).}
\tag{12}
\]
Because \(r!\) is a unit, one may replace \(a_r\) by \(A_r\) here.

## 4. A complete prime-square criterion in the top-half range

Define the integer polynomials \(C_r(u)\) by
\[
C_0(u)=1,\quad C_1(u)=2u+1,
\quad C_{r+1}(u)=(2u+2r+1)C_r(u)+(u+r)^2C_{r-1}(u).
\]
Then \(C_r(0)=A_r\). Put \(D_r=C_r'(0)\); these integers satisfy
\[
D_0=0,\quad D_1=2,
\quad
D_{r+1}=(2r+1)D_r+r^2D_{r-1}+2A_r+2rA_{r-1}.
\tag{13}
\]
The singular coefficient \(p^2\) at index \(p\) vanishes modulo \(p^2\),
so for either solution \(X\),
\[
X_{p+r}\equiv C_r(p)X_p\pmod{p^2}.
\tag{14}
\]
Since \(B_p\) is a unit and
\(C_r(p)\equiv A_r+pD_r\pmod{p^2}\), this controls the second valuation.
For \(0<r<p\), (9) bounds \(v_p(g_{p+r})\le2\). Hence
\[
\boxed{
v_p(g_{p+r})=
\begin{cases}
0,&A_r\not\equiv0\pmod p,\\
1,&A_r\equiv0\pmod p,\quad A_r/p+D_r\not\equiv0\pmod p,\\
2,&A_r\equiv0\pmod p,\quad A_r/p+D_r\equiv0\pmod p.
\end{cases}}
\tag{15}
\]
This range is exactly \(p<n<2p\), where \(v_p(n!)=1\). Therefore the last
case gives a genuine extra prime factor in \(g_n/n!\).

For example, \(p=73,r=34,n=107\) gives
\[
A_{34}/73\equiv37,\qquad D_{34}\equiv36\pmod{73}.
\]
Their sum is zero, so \(v_{73}(g_{107})=2\), whereas
\(v_{73}(107!)=1\). By contrast \(p=17,r=4\) gives residues 12 and 6;
the sum is 1 modulo 17, so \(v_{17}(g_{21})=1\).

## 5. Consequences and remaining barrier

These results identify both the usual loss and a concrete recovery mechanism.
Outside the digit-residue exceptional set, the gcd loses exactly the full
\(p^{\lfloor\log_p n\rfloor}\) factor compared with \(n!\). Exceptional
digits recover at least one power; the shifted continuant can recover an
additional power through a second congruence. Thus an upper bound based on
the false assumption \(g_n\mid n!\) cannot work.

No bound for the aggregate exceptional prime powers is proved here. In
particular these local criteria neither prove nor disprove that
\(g_n\sqrt n/n!\) is unbounded. A global proof still needs control over
the simultaneous exceptional residues and their higher lifts as \(n\)
varies, rather than a count at finitely many indices.
