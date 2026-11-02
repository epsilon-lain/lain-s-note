# Retaining deficits across many prime bands

2026-11-02. The fixed recurrence and reduced approximants are unchanged. This note gives an exact divisibility certificate for ordinary deficits across many prime intervals, followed by a global necessary prime-concentration condition. It does not prove the concentration cannot occur, and does not settle the proposed eventual S(n)<=n/6 bound.

## 1. Objects and established local input

Let A_0=A_1=1, B_0=0, B_1=1, with

\[
X_{n+1}=(2n+1)X_n+n^2X_{n-1}.
\]

Write ell_n=A_n/n!, w_n=B_n/n!, g_n=gcd(A_n,4B_n), and

\[
e_p(n)=v_p(g_n)-v_p(n!),\quad
E_n=\prod_{p\text{ odd}}p^{\max(e_p(n),0)},\quad
D_n=\prod_{p\text{ odd}}p^{\max(-e_p(n),0)}.
\]

The integer U_n=2^n ell_n has positive odd part C_n. Thus p divides C_n precisely when ell_n vanishes modulo the odd prime p. The established all-index cap and leading-pole congruence are

\[
e_p(n)\le\lfloor\log_p(n-1)\rfloor-v_p(n),
\]
\[
\ell_{bp+r}\equiv\ell_b\ell_r\pmod p,
\qquad pw_{bp+r}\equiv(-1)^{(p-1)/2}\ell_r w_b\pmod p
\tag{1}
\]

when p>sqrt(n), n=bp+r, 1<=b<p, 0<=r<p. Here all denominators used to reduce ell_b,ell_r,w_b are prime to p. A nonzero residue on the right of the second congruence proves e_p(n)=-1.

## 2. A finite deficit certificate

Fix integers n>=2 and B>=1, and define

\[
t=\max\{\sqrt n,n/(B+1)\},\quad
H_{n,B}=\prod_{t<p\le n\atop p\text{ odd prime}}p,
\]
\[
P_{n,B}=\gcd(C_n,H_{n,B}),\qquad
F_{n,B}=\gcd\left(H_{n,B},\prod_{b=1}^{B}B_b\right).
\]

These are finite positive integers; H,P,F are squarefree. The letter B without a subscript is a band count, whereas B_b is the second recurrence sequence. Empty prime products are one.

Put M_{n,B}=H_{n,B}/gcd(H_{n,B},P_{n,B}F_{n,B}). Then

\[
\boxed{M_{n,B}\mid D_n,
\qquad H_{n,B}\mid D_nP_{n,B}F_{n,B}.}
\tag{2}
\]

Proof: take p dividing H but not P F. Write n=bp+r. Since p>n/(B+1), b<=B; since p>sqrt(n), b<p. The condition p not dividing P says ell_n is nonzero modulo p. The first congruence in (1) then forces ell_r to be nonzero. The condition p not dividing F gives B_b nonzero modulo p, hence w_b nonzero because b! is a p-unit. The second congruence proves e_p(n)=-1, so p divides D_n. Distinct such primes give (2).

This proof also covers r=0 and p=n. No upper-band endpoint is silently dropped. It does not require factoring the huge A_n or B_n to justify the theorem.

Let J_n=ln(E_n/D_n). The finite consequence is

\[
\boxed{J_n\le\ln E_n+\ln P_{n,B}+\ln F_{n,B}-\ln H_{n,B}.}
\tag{3}
\]

It explicitly retains genuine positive excess E_n, rather than replacing every first-layer zero by a positive excess.

There is also a finite improvement which does not need an asymptotic estimate for P. Write d_n=v_2(g_n/n!). From the established strict analytic bound Q_n>2R_n^{-2}/(2n+1), where R_n=g_n/n!=2^{d_n}E_n/D_n, equation (2) gives

\[
Q_n>\frac{2^{1-2d_n}H_{n,B}^2}
{(2n+1)E_n^2P_{n,B}^2F_{n,B}^2}.
\tag{3a}
\]

The squarefree overlap can be retained exactly: since D_n>=M_{n,B}, any proved multiple G_n of E_n yields

\[
\boxed{Q_n>\frac{2^{1-2d_n}M_{n,B}^2}{(2n+1)G_n^2}.}
\tag{3b}
\]

For example, G_n=gcd(C_n,odd(lcm(1,...,n-1))) is valid. The sharper gain-only divisor Ghat_n from the preceding local-pole sieve is valid too. With s_n=v_2(U_n), the elementary bound d_n<=s_n-n gives the directly computable certificate

\[
\boxed{Q_n>\frac{2\,2^{2(n-s_n)}M_{n,B}^2}
{(2n+1)\widehat G_n^2}.}
\tag{3c}
\]

This multiplies the existing gain-only sieve floor by M_{n,B}^2>=1, so it is an unconditional finite improvement whenever M_{n,B}>1. It requires U_n and only the prefix B_1,...,B_B, rather than the full B_n. The asymptotic inequality (5) below should not be confused with this finite, always-improving certificate.

## 3. The troublesome small leading coefficients have sublinear total cost

Positivity of the recurrence gives 1<=B_b<=A_b for b>=1. With lambda=1+sqrt(2), the established coefficient bound ell_b<=lambda^b implies

\[
\ln F_{n,B}\le\sum_{b=1}^{B}\ln B_b
\le\frac{B(B+1)}2(\ln B+\ln\lambda).
\tag{4}
\]

For B=floor(n^(1/3)), the right side is O(n^(2/3)ln n)=o(n). Thus primes for which a head coefficient w_b vanishes can be removed at a proved sublinear aggregate cost. No independence assumption about those primes is used.

For this choice B, the threshold t is n^(2/3)(1+o(1)). By the classical prime number theorem,

\[
\ln H_{n,B}=n+o(n).
\]

Define the full high-prime first-layer weight

\[
W(n)=\ln P_{n,B}
=\sum_{t<p\le n\atop p\mid C_n}\ln p,
\qquad B=\lfloor n^{1/3}\rfloor.
\]

Equations (3)--(4) therefore give the all-index asymptotic estimate

\[
\boxed{J_n\le W(n)+\ln E_n-n+o(n).}
\tag{5}
\]

This differs from the upper-half S(n): W spans all prime bands down to approximately n^(2/3).

## 4. Eliminating E at the cost of a first-layer upper bound

For every selected prime p>sqrt(n), positive excess e_p(n) is at most one and implies p divides C_n. Therefore the selected-prime part of E_n divides P_{n,B}. All other positive primes are at most t; there are none above n. Each of their contributions e_p(n)ln p is at most ln(n-1). Consequently

\[
\ln E_n\le W(n)+\lfloor t\rfloor\ln(n-1)
=W(n)+o(n).
\tag{6}
\]

Combining (5) and (6),

\[
\boxed{J_n\le2W(n)-n+o(n).}
\tag{7}
\]

For a fixed finite band count B, the same argument with the usual lcm bound for the lower primes yields

\[
J_n\le2\ln P_{n,B}-\frac{B-1}{B+1}n+o(n).
\tag{8}
\]

For fixed B, the head-prime product is fixed. The lower-prime positive budget is at most theta(n/(B+1))+O(sqrt(n)ln n)=n/(B+1)+o(n); the error accounts for the prime powers beyond their first powers. Thus (8) recovers the upper-half baseline at B=1 and retains an additional negative main term when more bands are included. It is not an automatic improvement for a specified n, because P_{n,B} also includes more exceptional primes.

## 5. What would suffice, and what is still missing

The established quadratic quality relation and binary saturation give

\[
\ln Q_n=n\ln2-2J_n+O(\ln n),
\qquad Q_n=q_n^2|\pi-p_n/q_n|.
\]

Consequently,

\[
\boxed{\ln Q_n\ge(2+\ln2)n-4W(n)+o(n).}
\tag{9}
\]

A sufficient global exclusion would therefore be

\[
\limsup_{n\to\infty}\frac{W(n)}n
<\frac12+\frac{\ln2}{4}=0.6732867951\ldots.
\tag{10}
\]

Equivalently, any unbounded sequence of indices with Q_n bounded above must have

\[
W(n)\ge\left(\frac12+\frac{\ln2}{4}\right)n-o(n).
\tag{11}
\]

In prime-product language, at least this exponential weight of the high-prime primorial must divide the same C_n. Since each selected prime contributes at most ln n, at least (1/2+ln2/4-o(1))n/ln n distinct high primes must be first-layer zeros along such successful indices. These zeros need not all be positive excess primes.

Compared with ln C_n<=ln(2+sqrt(2)) n, the necessary selected-prime weight is at least the fraction

\[
\frac{1/2+\ln2/4}{\ln(2+\sqrt2)}=0.5483\ldots
\]

of this entire available logarithmic height. This ratio measures a necessary concentration requirement; it is not an estimate asserting typical factorization behavior.

This is a necessary structural condition, not a proof that it is impossible. The available height C_n<=(2+sqrt(2))^n gives ln C_n<=1.227947...n, which does not imply (10). Even the bound P<=H only gives W<=n+o(n). Likewise (5) does not automatically beat J_n<=ln E_n: its usefulness requires actual information about W or the complementary nonzero prime product.

Let V=H/P. Equation (2) also gives V|D F. In the growing-band choice, (11) is equivalent to the necessary condition

\[
\ln V\le\left(\frac12-\frac{\ln2}{4}\right)n+o(n)
=0.3267132049\ldots n+o(n)
\]

along any sequence with bounded Q. A uniform lower bound for this complementary nonzero prime product exceeding that constant would therefore exclude the family. This is still a cross-modulus arithmetic problem. Neither the requested upper-half S bound nor the new sufficient W bound is proved here.

## Verification and dependencies

The companion standard-library script passed 5,715 finite divisibility certificates for 2<=n<=1000 and band counts drawn from 1,2,3,8,16,floor(n^(1/3)). It checked 355,049 generic selected-prime valuations directly against the actual recurrence gcds and 42 selected-prime positive-excess support events. Of these certificates, 5,672 have M>1. Every computed new quality floor is at least the preceding gain-only floor and at most the exact analytic floor based on the actual gcd.

For B=floor(n^(1/3)), the previous gain-only sieve left the indices 2,3,4,5,6,7,11,15 below its threshold one. The new finite certificate (3c) leaves only 2,3; thus it proves Q_n>1 for every 4<=n<=1000. This is a finite certified improvement, not an extrapolation to all n. The verification does not prove (10).

The local formulas are from the preceding odd-prime Wronskian and leading-pole proofs; the coefficient height, binary saturation and quality normalization are from the existing fixed-family checkpoints. The only asymptotic prime-distribution input is the classical PNT; no zero-distribution theorem is assumed.
