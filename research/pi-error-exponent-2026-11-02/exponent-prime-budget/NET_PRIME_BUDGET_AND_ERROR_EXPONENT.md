# Exact net prime budget for fixed error prefactors and exponents

2026-11-02. This note translates the actual reduced gcd into error bounds of the form C/q^c, with arbitrary fixed C>0 and c>0. Positive odd-prime excess, negative odd-prime deficits, and the exact binary contribution are kept separately. The main additional consequence is a necessary simultaneous-lift density: when the desired exponent exceeds the binary-loss threshold, a successful index needs order n/ln n positive odd-prime contributions, predominantly from primes larger than sqrt(n). This is a necessary condition within this specified Padé family, not a statement about all rational approximations to pi.

## 1. Exact objects and the signed budget

Use the integer recurrence

\[
A_0=A_1=1,\quad B_0=0,\ B_1=1,\qquad
X_{n+1}=(2n+1)X_n+n^2X_{n-1}.
\]

For n>=1 put F_n=n!, ell_n=A_n/F_n, g_n=gcd(A_n,4B_n), and let

\[
\frac{p_n}{q_n}=\frac{4B_n}{A_n},\quad
p_n=4B_n/g_n,\quad q_n=A_n/g_n>0,
\qquad \varepsilon_n=\left|\pi-\frac{p_n}{q_n}\right|.
\]

All fractions are fully reduced. For each odd prime p define

\[
e_p(n)=v_p(g_n)-v_p(F_n),
\quad E_n=\prod_{p\text{ odd}}p^{\max(e_p(n),0)},
\quad D_n=\prod_{p\text{ odd}}p^{\max(-e_p(n),0)}.
\]

The finite products E_n,D_n are coprime odd integers. Define the exact binary quantities

\[
d_n^A=v_2(A_n)-v_2(F_n),\quad
t_n=v_2(A_n)-v_2(g_n)\ge0,\quad
\beta_n=d_n^A+n/2.
\]

Unique factorization gives the exact identity

\[
\boxed{\frac{g_n}{F_n}=2^{-n/2+\beta_n-t_n}\frac{E_n}{D_n}.}
\tag{1}
\]

Write

\[
J_n=\ln E_n-\ln D_n-t_n\ln2.
\tag{2}
\]

All ln symbols denote the real natural logarithm of a positive real number. Equation (1) keeps the negative odd-prime budget: knowing E_n alone never identifies g_n/F_n.

The proved exact valuation formula for A_n makes beta explicit. With s_2(n) the number of ones in the binary representation,

\[
\beta_n=
\begin{cases}
s_2(n),&n\equiv0\pmod4,\\
s_2(n)-1/2,&n\equiv1\pmod4,\\
s_2(n)+1,&n\equiv2\pmod4,\\
s_2(n)+v_2(n+1)-1/2,&n\equiv3\pmod4.
\end{cases}
\tag{3}
\]

In particular beta_n<=2 log_2(n+1)+1. The separate saturation theorem proves t_n=0 for every n except n=3, where t_3=1. Thus for all n>=4,

\[
\boxed{J_n=\ln(E_n/D_n).}
\tag{4}
\]

Equation (1) and the subsequent bounds remain correct even before applying saturation, by retaining t_n explicitly.

There is also a useful unconditional upper budget for the deficits. The exact convolution w_n=sum_(j=1)^n ell_(j−1)ell_(n−j)/j and ell_k in Z[1/2] show v_p(w_n)>=−floor(log_p n) for odd p. Since ell_n is p-integral, this yields

\[
\boxed{D_n\mid\operatorname{odd}\bigl(\operatorname{lcm}(1,\ldots,n)\bigr).}
\tag{4a}
\]

It supplies a denominator upper bound as well as the lower-budget obstructions studied here. The companion exponent-map note develops its consequences for weak attainable exponents; no additional asymptotic statement about the lcm is needed in this note.

## 2. Two uniform bounds for ell_n

Set lambda=1+sqrt(2), a=ln lambda and b=(ln2)/2. For every n>=0,

\[
\boxed{\frac{\lambda^n}{2n+1}\le\ell_n\le\lambda^n.}
\tag{5}
\]

Here is an elementary coefficient proof. Coefficient extraction from (1+2x+2x²)^n gives its central coefficient U_n=2^n ell_n. Substitute x=y/sqrt(2). The central coefficient of (1+sqrt(2)y+y²)^n is then (sqrt(2))^n ell_n, while the sum of all its 2n+1 coefficients is (2+sqrt(2))^n.

These coefficients are symmetric and attain their maximum at the center. To check the latter without invoking a general convolution theorem, let an already symmetric unimodal sequence be b_k, indexed about zero and zero outside its finite support. Its convolution c_k=b_(k−1)+sqrt(2)b_k+b_(k+1) satisfies, for k>=0,

\[
c_k-c_{k+1}
=b_{k-1}-b_{k+2}+(\sqrt2-1)(b_k-b_{k+1})\ge0.
\]

For k=0 use b_(−1)=b_1. Starting with the single central coefficient proves the assertion for every power. The central coefficient therefore lies between the average and the total. Dividing by (sqrt(2))^n gives (5).

## 3. All-index necessary and sufficient budget tests

The exact positive Legendre error identity gives, for every n>=1,

\[
\frac{2}{(2n+1)\ell_n^2}<\varepsilon_n
<\frac{4}{(2n+1)\ell_n^2}.
\tag{6}
\]

Let c>0 and use x_+=max(x,0). Define the desired linear net-budget rate

\[
\Theta_c=a+b-\frac{2a}{c}.
\]

Since q_n=ell_n F_n/g_n, equations (1), (5), and (6) imply the strict bounds

\[
\boxed{\begin{aligned}
\ln(\varepsilon_n q_n^c)
&>c\Theta_c n-cJ_n-c\beta_n\ln2
-[1+(c-2)_+]\ln(2n+1)+\ln2,\\
\ln(\varepsilon_n q_n^c)
&<c\Theta_c n-cJ_n-c\beta_n\ln2
+[(2-c)_+-1]\ln(2n+1)+\ln4.
\end{aligned}}
\tag{7}
\]

For c<2 the direction of the ell-power bound reverses; the positive-part terms in (7) account for that reversal. No integer c assumption is used.

For a fixed C>0, the target epsilon_n<=C q_n^(−c) consequently has this necessary condition:

\[
\boxed{J_n>
\Theta_c n-\beta_n\ln2
-\frac{1+(c-2)_+}{c}\ln(2n+1)
+\frac{\ln(2/C)}c.}
\tag{8}
\]

The following condition is sufficient:

\[
\boxed{J_n\ge
\Theta_c n-\beta_n\ln2
+\frac{(2-c)_+-1}{c}\ln(2n+1)
+\frac{\ln(4/C)}c.}
\tag{9}
\]

They are uniform finite-index tests, with an explicit logarithmic gap. For c=2 the polynomial terms match, and the remaining gap is only (ln2)/2. Taking C=e changes only the constant terms ln(2/C) and ln(4/C), since ln e=1; it does not change the linear threshold Theta_c.

## 4. Sharp asymptotic normalization and its boundary

The direct-integral and Legendre asymptotics establish

\[
\ell_n=L\lambda^n n^{-1/2}(1+O(n^{-1})),
\quad \varepsilon_n=K\lambda^{-2n}(1+O(n^{-1})),
\]

where

\[
L=\frac{\lambda^{1/2}}{\sqrt{2\pi\sqrt2}},
\qquad K=\frac{4\pi}{\lambda}.
\]

For fixed c>0 define

\[
\mathcal B_{n,c}=J_n+\beta_n\ln2+\frac12\ln n-\Theta_c n.
\]

Substitution into the exact denominator formula yields

\[
\boxed{\varepsilon_n q_n^c
=K L^c\exp(-c\mathcal B_{n,c})(1+O(n^{-1})).}
\tag{10}
\]

Thus epsilon_n q_n^c tends to zero on a chosen subsequence if and only if B_(n,c) tends to +infinity there. For a fixed target C, its asymptotic boundary is

\[
\mathcal B_{n,c}=\ln L+\frac{\ln(K/C)}c+O(n^{-1}).
\]

This boundary statement does not decide an inequality when the difference itself is O(n^(−1)); the uniform tests (8)--(9) or a more precise error certificate are then needed.

At c=2, K L²=sqrt(2), so the boundary becomes

\[
\boxed{\mathcal B_{n,2}
=J_n+\beta_n\ln2+\tfrac12\ln n-bn,
\qquad
\varepsilon_n q_n^2\sim\sqrt2 e^{-2\mathcal B_{n,2}}.}
\tag{11}
\]

For C=e the fixed-quality boundary is B_(n,2)=(ln2)/4−1/2+O(n^(−1)). For quality tending to zero along a specified subsequence one needs B_(n,2) to tend to +infinity there, not merely a fixed positive excess above this constant. Being unbounded above suffices only to extract a further subsequence tending to +infinity.

## 5. What subexponential positive excess would exclude

Discarding D_n>=1 and t_n>=0, and using beta_n<=2 log_2(n+1)+1 in the lower bound (7), gives the all-index inequality

\[
\boxed{\varepsilon_n q_n^c>
\frac{2^{1-c}\exp(c\Theta_c n)}
{(2n+1)^{1+(c-2)_+}(n+1)^{2c}E_n^c}.}
\tag{12}
\]

If ln E_n=o(n), this proves epsilon_n q_n^c tends to infinity for every

\[
\boxed{c>c_*=\frac{2\ln(1+\sqrt2)}{\ln(1+\sqrt2)+(\ln2)/2}.}
\tag{13}
\]

The constant c_* is approximately 1.4355. More generally, an upper rate limsup (ln E_n)/n<=sigma excludes every c for which Theta_c>sigma. These statements also apply along any specified subsequence on which the assumed upper rate holds. No such subexponential excess bound is asserted here for the actual full family.

At equality Theta_c=sigma, the exponential comparison is inconclusive. The explicit beta_n, the ln n term, odd deficits D_n, and the prefactor C determine the remaining scale. Conversely, c<c_* does not by itself prove success: the discarded deficits can still dominate. These thresholds are properties of this candidate family and are not a global irrationality exponent for pi.

## 6. Necessary simultaneous-lift density

Let N_+(n) be the number of distinct odd primes for which e_p(n)>0. The proved prime cap

\[
e_p(n)\le\lfloor\log_p(n-1)\rfloor-v_p(n)
\]

implies that each positive contribution obeys

\[
e_p(n)\ln p\le\ln(n-1).
\]

For n>=4, equations (4) and (8), together with J_n<=ln E_n, therefore force

\[
\boxed{N_+(n)\ln(n-1)>
\Theta_c n-\beta_n\ln2
-\frac{1+(c-2)_+}{c}\ln(2n+1)
+\frac{\ln(2/C)}c}
\tag{14}
\]

at every successful index. In particular, when Theta_c>0, any unbounded successful sequence must have

\[
\boxed{N_+(n)\ge(\Theta_c+o(1))\frac n{\ln n}.}
\tag{15}
\]

There is a stronger localization of that budget. Primes p<=sqrt(n−1) contribute at most sqrt(n−1) ln(n−1) in total, since there are at most sqrt(n−1) such primes and each contributes at most ln(n−1). This is o(n), with no prime-number theorem needed. Every positive prime p>sqrt(n−1) has e_p(n)=1. Consequently the indispensable linear budget must come from many distinct large primes:

\[
\sum_{\substack{p>\sqrt{n-1}\\e_p(n)>0}}\ln p
>\Theta_c n-O(\sqrt n\ln n).
\tag{16}
\]

The number of those large primes is at least

\[
(\Theta_c+o(1))n/\ln n.
\]

For c=2, Theta_2=(ln2)/2. Thus an isolated extra factor, or even o(n/ln n) simultaneously lifted odd primes, cannot reach any fixed target C/q² at arbitrarily large indices within this family. The finite observations of single-prime excess do not establish such sparsity for all n.

## 7. A finite integer-only illustration at n=50000

Take C=e,c=2. Since e<3, the exact lower error bound and E_n<=(n−1)^(N_+(n)) imply the necessary integer inequality

\[
3(2n+1)(n-1)^{2N_+(n)}>2^{1-2d_n^A}.
\tag{17}
\]

This remains necessary even if binary saturation is not invoked, because v_2(g_n)<=v_2(A_n). At n=50000 the exact binary formula gives s_2(n)=6 and d_n^A=−24994. Therefore

\[
3\cdot100001\cdot49999^{2N_+(50000)}>2^{49989}.
\]

Exact integer comparison proves that N=1600 fails this inequality and N=1601 passes. Hence a fraction from this family at index 50000 with epsilon<=e/q² would require **at least 1601 distinct positive-excess odd primes at that same index**. Passing this coarse necessary bound does not ensure existence.

The small-prime budget can be removed more accurately. There are 47 odd primes p<=floor(sqrt(49999))=223. Form the integer

\[
M=\prod_{\substack{p\le223\\p\text{ odd prime}}}
p^{\max(0,\min\{\lfloor\log_p49999\rfloor-v_p50000,
\lfloor\log_p50001\rfloor-v_p50001\})}.
\]

The product M has 188 decimal digits. Every possible small-prime part of E_50000 divides M. If N_large counts the positive-excess primes above 223, success requires

\[
3\cdot100001\cdot M^2\cdot49999^{2N_{\rm large}}>2^{49989}.
\tag{18}
\]

Exact comparison shows failure at N_large=1560 and passage at 1561. Thus **at least 1561 odd primes greater than 223 would have to lift simultaneously**, each contributing its one permitted excess factor. These are unconditional necessary counts derived from the error and valuation bounds, not a search for actual A_50000 or B_50000.

## Verification and dependencies

The companion standard-library verifier checks the exact factorization (1), beta residues, binary saturation, uniform ell bounds in Q(sqrt(2)), and finite C=e targets at rational exponents. It also checks the excess caps and the combinatorial contribution bounds. All infinite-index claims follow from the displayed proofs.

Dependencies are the [positive Legendre error identity](../../pi-joint-arithmetic-2026-11-02/joint-log-pade/LEGENDRE_ARCTAN_JOINT_COST.md), [direct-integral asymptotics](../../pi-joint-arithmetic-2026-11-02/joint-arctan-integral/DIRECT_INTEGRAL_GCD_THRESHOLD.md), the [supplied 2-adic valuation theorem](../../pi-gcd-prime-adic-2026-11-02/gcd_arithmetic/TWO_ADIC_VALUATION.md), the [binary-saturation proof](../exponent-audit/TWO_ADIC_GCD_SATURATION.md), and the [odd-prime Wronskian cap](../../pi-gcd-blast-2026-11-02/gcd-prime-blocks/ODD_PRIME_BLOCKS_AND_EXCESS.md). The [companion exponent map](../exponent-error-map/REDUCED_DENOMINATOR_EXPONENT_MAP.md) develops the sharp normalized-gcd version. No unproved equality v_2(g_n)=v_2(A_n) is silently assumed; it is the explicitly cited saturation theorem, with n=3 retained as an exception.
