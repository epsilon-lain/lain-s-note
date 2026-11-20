# Explicit future-lookahead bound from two partial quotients

Date: 2026-11-20  
Status: **proved explicit bound** for the future-veto stopping algorithm.

## 1. Setup

Let \(\alpha\) be irrational, with convergents \(p_n/q_n\), errors
\[
E_n=|q_n\alpha-p_n|,\qquad 0<E_{n+1}<E_n.
\]
For a fixed candidate \(n\), its future crossing with candidate \(j>n\) is
\[
X_{n,j}=\frac{q_j^2-q_n^2}{E_n^2-E_j^2}.
\]
Let
\[
F_n=\min_{j>n}X_{n,j}
\]
be its exact future-veto threshold. After examining \(j=n+1,\ldots,m\), let
\[
U_{n,m}=\min_{n<j\le m}X_{n,j}.
\]
The finite stopping certificate says that if
\[
\frac{q_{m+1}^2-q_n^2}{E_n^2}\ge U_{n,m},
\tag{1}
\]
then \(F_n=U_{n,m}\), so the future has been completely resolved.

## 2. A bound using only the first future crossing

Since \(U_{n,m}\le X_{n,n+1}\), it is enough to ensure
\[
q_{m+1}^2\ge q_n^2+E_n^2X_{n,n+1}.
\]
Put
\[
z_n=\frac{E_{n+1}}{E_n},\qquad r_n=\frac{q_n}{q_{n+1}}.
\]
Then
\[
q_n^2+E_n^2X_{n,n+1}
=q_n^2+\frac{q_{n+1}^2-q_n^2}{1-z_n^2}.
\]
Dividing by \(q_{n+1}^2\) gives
\[
r_n^2+\frac{1-r_n^2}{1-z_n^2}
=\frac{1-r_n^2z_n^2}{1-z_n^2}
\le\frac1{1-z_n^2}.
\]
Thus the sufficient condition is
\[
\boxed{q_{m+1}\ge\frac{q_{n+1}}{\sqrt{1-z_n^2}}.}
\tag{2}
\]

## 3. Fibonacci growth turns this into an explicit number of future terms

The denominator recurrence is
\[
q_{k+1}=a_{k+1}q_k+q_{k-1}\ge q_k+q_{k-1},
\]
because every regular partial quotient \(a_{k+1}\ge1\).

Let \(F_1=F_2=1\) be the Fibonacci numbers. Induction gives
\[
q_{n+k}\ge F_kq_{n+1}+F_{k-1}q_n
\ge F_kq_{n+1}.
\]
Therefore (2) is guaranteed as soon as
\[
\boxed{F_k\ge\frac1{\sqrt{1-z_n^2}}.}
\tag{3}
\]
If \(k\) is the least integer at least 2 satisfying (3), then examining the \(k-1\) future convergents
\[
n+1,n+2,\ldots,n+k-1
\]
is sufficient to determine \(F_n\) exactly.

This is an explicit finite-lookahead bound, not merely eventual termination.

## 4. Replace the error ratio by two partial quotients

The complete-quotient identity gives
\[
z_n=\frac{E_{n+1}}{E_n}=\frac1{\alpha_{n+2}},
\qquad
\alpha_{n+2}=[a_{n+2};a_{n+3},\ldots].
\]

### Case A: \(a_{n+2}\ge2\)

Then \(\alpha_{n+2}\ge2\), so \(z_n\le1/2\). Hence
\[
\frac1{\sqrt{1-z_n^2}}\le\frac2{\sqrt3}<2=F_3.
\]
Take \(k=3\). Therefore:

\[
\boxed{a_{n+2}\ge2\quad\Longrightarrow\quad
\text{at most 2 future convergents need be examined}.}
\tag{4}
\]

### Case B: \(a_{n+2}=1\)

Write \(A=a_{n+3}\ge1\). Then
\[
\alpha_{n+2}=1+\frac1{\alpha_{n+3}},
\qquad
A<\alpha_{n+3}<A+1.
\]
It follows that
\[
z_n=\frac{\alpha_{n+3}}{\alpha_{n+3}+1}
<\frac{A+1}{A+2}.
\]
Consequently
\[
\frac1{\sqrt{1-z_n^2}}
<
\frac{A+2}{\sqrt{2A+3}}.
\]
Define
\[
\boxed{
k(A)=\min\left\{k\ge2:
F_k\ge\frac{A+2}{\sqrt{2A+3}}\right\}.
}
\tag{5}
\]
Then:

\[
\boxed{a_{n+2}=1,\ a_{n+3}=A
\quad\Longrightarrow\quad
\text{at most }k(A)-1\text{ future convergents need be examined}.}
\tag{6}
\]

Since \(F_k\ge\varphi^{k-2}\) for \(k\ge2\), where
\(\varphi=(1+\sqrt5)/2\), a simpler closed upper bound is
\[
\boxed{
\text{lookahead}\le
\left\lceil
1+\log_\varphi\frac{A+2}{\sqrt{2A+3}}
\right\rceil.
}
\tag{7}
\]
In particular the required lookahead is \(O(\log A)\), with an absolute implied constant.

## 5. What is genuinely new here

The earlier stopping theorem said that every fixed candidate can be certified after finitely many future convergents, but gave no explicit bound on their number.

Equations (4) and (6) are a data-dependent lookahead theorem:
- if the second future partial quotient is at least 2, two future terms suffice;
- if it equals 1, the next partial quotient gives a logarithmic bound.

The bound is conservative because it uses only the first future crossing as an upper bound for the true minimum \(F_n\). In practice the algorithm often stops much earlier.

It is not a uniform constant-memory theorem: partial quotients \(a_{n+3}\) can be arbitrarily large for general irrational numbers, so this bound can grow.

## 6. Pi experiment

For each of the first several convergents of \(\pi\), the companion script:
1. reads \(a_{n+2}\) and, when necessary, \(a_{n+3}\);
2. computes the Fibonacci lookahead bound;
3. runs the exact finite-prefix stopping test at high precision;
4. records the actual number of future terms examined.

The theorem is exact for every irrational \(\alpha\). The decimal experiment for \(\pi\) is only a finite numerical illustration; it is not an interval-certified computation of every crossing.

## 7. Scope

This gives an explicit future-lookahead bound in terms of two partial quotients. It does not bound the partial quotients of \(\pi\), does not produce a uniform lookahead independent of \(n\), and does not prove that \(\pi\) is badly approximable. It does turn the two-sided selection rule into a terminating, quantitatively controlled procedure for each fixed candidate.
