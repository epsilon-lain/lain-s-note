# Independent review of the joint integral and reduced-denominator argument

2026-11-02. Reviewed independently:

- `../joint-log-pade/LEGENDRE_ARCTAN_JOINT_COST.md`
- `../joint-arctan-integral/DIRECT_INTEGRAL_GCD_THRESHOLD.md`

**Verdict:** the exact identities, constants, asymptotic equivalence, reduced-denominator obstruction, and neighboring-gcd argument are sound. The results concern a specified family. They do not decide whether pi is badly approximable, or whether this family has a subsequence whose normalized error tends to zero.

## 1. Independent rational-coordinate derivation of the integral identity

Write

\[
r(x)=\frac{x(1-x)}{1+x^2},\qquad
v(x)=\frac{x^2+2x-1}{1+x^2},\qquad
I_n=\int_0^1\frac{r(x)^n}{1+x^2}\,dx.
\]

Direct rational differentiation gives, for n >= 1,

\[
\frac d{dx}\bigl(v(x)r(x)^n\bigr)
=\frac{r(x)^{n-1}}{1+x^2}
\left[4(n+1)r(x)^2+2(2n+1)r(x)-n\right].
\tag{1}
\]

Since r vanishes at both endpoints, integration gives

\[
4(n+1)I_{n+1}+2(2n+1)I_n-nI_{n-1}=0.
\tag{2}
\]

Also v'=(4r+2)/(1+x^2), v(1)-v(0)=2, and I_0=pi/4. Thus I_1=(1-I_0)/2. The sequence D_n=(-2)^n I_n therefore obeys

\[
(n+1)D_{n+1}=(2n+1)D_n+nD_{n-1},
\quad D_0=\pi/4,\quad D_1=\pi/4-1.
\]

This proves, independently of a contour deformation or special-function asymptotic,

\[
(-2)^n I_n=\ell_n\pi/4-w_n,
\qquad
\left|\pi-\frac{4w_n}{\ell_n}\right|
=\frac{4\,2^n I_n}{\ell_n}.
\tag{3}
\]

The positivity of I_n and ell_n makes the sign and the nonzero error exact. Formula (1) was also reduced to two polynomial identities, independent of n, in the companion audit script.

## 2. Independent Laplace constants

Let s=sqrt(2)-1, so lambda=1/s and rho=s/2. The rational kernel r has its unique maximum in (0,1) at x=s. Direct differentiation gives

\[
r(s)=s/2,\qquad
(\log r)''(s)=-\frac1{\sqrt2\,s^3},\qquad
\frac1{1+s^2}=\frac1{2\sqrt2\,s}.
\]

The interior Gaussian contribution is consequently

\[
I_n\sim\left(\frac s2\right)^n
\sqrt{\frac{\pi s}{2\sqrt2\,n}}
=\rho^n\sqrt{\frac{\pi\rho}{\sqrt2\,n}}.
\tag{4}
\]

For ell_n, the real integral

\[
\ell_n=\frac1\pi\int_0^\pi(1+\sqrt2\cos\theta)^n\,d\theta
\]

has the same initial values and recurrence: integrate the derivative of sin(theta)(1+sqrt(2)cos(theta))^n. Its dominant endpoint is theta=0, where the base is lambda and the second derivative of its logarithm is -sqrt(2)/lambda. The half-Gaussian coefficient is

\[
\ell_n\sim\lambda^n\sqrt{\frac{\lambda}{2\pi\sqrt2\,n}}.
\tag{5}
\]

The negative-base part causes no competing dominant contribution: its absolute base is at most s, whereas lambda>1. If a positive integrand is desired, pair theta and pi-theta on [0,pi/2]. The sum of the two nth powers is positive, including for odd n, because the positive base has greater magnitude.

Multiplying (4) and (5), with 2*rho*lambda=1, gives the constant without a normalization ambiguity:

\[
4\,2^n\ell_n I_n\sim\frac{\sqrt2}{n}.
\tag{6}
\]

Taking their ratio instead gives the author's equivalent estimate

\[
\left|\pi-4w_n/\ell_n\right|
\sim4\pi\lambda^{-2n-1}.
\]

The stated relative O(1/n) errors are justified by the smooth nondegenerate maxima. In the author's angular coordinates both expansions are even about their respective maxima. Splitting at a fixed neighborhood and then at distance n^(-2/5) yields an exponentially small complement, Gaussian tails of size exp(-c n^(1/5)), and an integrated quartic remainder O(1/n). The finite checks do not prove these asymptotic estimates; this local analysis does.

## 3. Integer normalization and the exact arithmetic threshold

Multiplying the ell,w recurrence by n! gives integer sequences

\[
A_n=n!\ell_n,\quad B_n=n!w_n,\qquad
X_{n+1}=(2n+1)X_n+n^2X_{n-1},
\]

with A_0=A_1=1 and B_0=0,B_1=1. For

\[
g_n=\gcd(A_n,4B_n),\quad
p_n=4B_n/g_n,\quad q_n=A_n/g_n,
\]

these are the genuinely reduced fractions, and q_n/ell_n=n!/g_n exactly. Equation (6) proves

\[
\boxed{Q_n:=q_n^2|\pi-p_n/q_n|
=\frac{\sqrt2}{n}\left(\frac{n!}{g_n}\right)^2
\bigl(1+O(n^{-1})\bigr).}
\tag{7}
\]

The error term in (7) is a relative analytic error, independent of the size of g_n. Thus for every diverging index subsequence,

\[
Q_{n_j}\to0
\quad\Longleftrightarrow\quad
\frac{n_j!}{g_{n_j}\sqrt{n_j}}\to0.
\tag{8}
\]

No statement that the right-hand condition occurs follows from the equivalence itself. It is a concrete arithmetic target, not a completed construction of a good subsequence.

## 4. Neighboring gcds and the stronger asymptotic product

The integer Wronskian obeys the same first-order recursion as the determinant of the transfer matrix. Its initial value is -1, so

\[
A_nB_{n-1}-A_{n-1}B_n=(-1)^n((n-1)!)^2.
\]

Both terms after multiplying by four are divisible by g_n*g_(n-1). Therefore

\[
g_ng_{n-1}\mid4((n-1)!)^2.
\tag{9}
\]

For h_n=g_n sqrt(n)/n!, n >= 1, this implies for n >= 2

\[
h_nh_{n-1}\le4\sqrt{\frac{n-1}{n}}.
\tag{10}
\]

Since (7) says Q_n=(sqrt(2)/h_n^2)(1+O(1/n)),

\[
\boxed{\liminf_{n\to\infty}Q_nQ_{n-1}\ge\frac18.}
\tag{11}
\]

The passage to a limit inferior is correct even if gcds fluctuate: the relative analytic factors approach one uniformly in that arithmetic choice. Equation (11) is not a claimed finite-n inequality with constant 1/8. It excludes simultaneous quality-zero limits at neighboring diverging indices, while permitting isolated good indices in principle.

## 5. Independent audit of the finite-n Legendre bound

For rational a,c>0, let alpha=c arctan(1/a), ell_n=i^(-n)L_n(ia), w_n=i^(1-n)W_n(ia), and let c*w_n/ell_n=p_n/q_n in lowest terms. Orthogonality gives

\[
(-1)^n\left(\arctan(1/a)-w_n/\ell_n\right)
=\frac{a}{2\ell_n^2}
\int_{-1}^{1}\frac{L_n(t)^2}{a^2+t^2}\,dt.
\]

The phase is correct: C(ia)=-i arctan(1/a), W_n(ia)/L_n(ia)=-i w_n/ell_n, and L_n(ia)^2=(-1)^n ell_n^2. Using the exact norm 2/(2n+1) gives the strict lower error

\[
\left|\alpha-cw_n/\ell_n\right|
>\frac{ca}{(a^2+1)(2n+1)\ell_n^2}.
\]

Strictness holds also at n=0: the weighted denominator attains its extrema only on sets of measure zero and the polynomial square is positive almost everywhere.

The rational Wronskian is ell_n*w_(n-1)-ell_(n-1)*w_n=(-1)^n/n. Adjacent reduced fractions therefore satisfy

\[
q_nq_{n-1}\ge\frac{n\ell_n\ell_{n-1}}c.
\]

This follows from the nonzero **integer** cross-difference after reduction, so it does not assume that n!*ell_n is the reduced denominator. Multiplication gives

\[
Q_nQ_{n-1}>
\frac{a^2n^2}{(a^2+1)^2(4n^2-1)}.
\tag{12}
\]

The factors c cancel. For a=1,c=4 this is >1/16 for every n>=1. The theorem applies separately to an arctangent family and cannot be transferred automatically to a same-index Machin combination. The author's note observes this distinction.

The added same-parity identity is also correct. Adding the two signed adjacent differences and using
\(n\ell_n-(n-1)\ell_{n-2}=(2n-1)a\ell_{n-1}\) gives, for n>=2,

\[
|r_n-r_{n-2}|=\frac{c(2n-1)a}{n(n-1)\ell_n\ell_{n-2}}.
\]

The same reduced-denominator argument and the two strict error bounds then give

\[
Q_nQ_{n-2}>
\frac{n^2(n-1)^2}
{(2n-1)^2(a^2+1)^2(2n+1)(2n-3)}.
\tag{13}
\]

Its limit is 1/[16(a^2+1)^2], equal to 1/64 for the pi family. This prevents either complete parity subsequence from having quality tending to zero, while leaving sparse subsequences possible.

## 6. Every fixed index gap has a positive limiting obstruction

For the pi family put e_n=|pi-r_n|. Equation (3) gives the exact sign pi-r_n=(-1)^n e_n, while (4)--(5) give

\[
\frac{e_n}{e_{n-k}}\longrightarrow\lambda^{-2k}<1
\]

for every fixed integer k>=1. The approximants at indices n and n-k are therefore distinct for all sufficiently large n. Their integer cross-difference after reduction implies

\[
Q_nQ_{n-k}\ge
\frac{e_ne_{n-k}}{|r_n-r_{n-k}|^2}
=\frac{e_n/e_{n-k}}
{\left[1-(-1)^k e_n/e_{n-k}\right]^2}.
\]

The denominator tends to a nonzero number, so

\[
\boxed{\liminf_{n\to\infty}Q_nQ_{n-k}
\ge\frac{\lambda^{-2k}}
{\left[1-(-1)^k\lambda^{-2k}\right]^2}>0.}
\tag{14}
\]

For k=1 this is 1/(lambda+lambda^(-1))^2=1/8. For k=2 it is 1/(lambda^2-lambda^(-2))^2=1/32, improving the limiting consequence of (13).

Suppose n_j is a strictly increasing sequence and Q_(n_j) tends to zero. If n_j-n_(j-1) failed to tend to infinity, some fixed positive gap k would occur infinitely often, by the finite pigeonhole principle. The products at those pairs would tend to zero, contradicting (14). Therefore every quality-zero subsequence in this family must have consecutive index gaps tending to infinity. This necessary sparsity condition does not rule out such a subsequence.

## 7. Validation and remaining limits

`independent_joint_arithmetic_checks.py` checks the rational-coordinate telescoping identity as two integer polynomial identities; independently builds the rational and integer recurrences; and checks strict error bounds, actual reduced denominators, Wronskians, neighboring gcd divisibility, and both finite-n product bounds against rational arctangent/pi intervals. It includes a fractional a<1 case to exercise the complementary-angle branch. Both author verifiers were also rerun and passed.

The direct-integral representation and the Legendre construction are the same family. The arithmetic equivalence (8) and the obstruction (12) are genuine new statements within this investigation. They establish neither the existence nor the impossibility of the sparse gcd gains required by (8), and supply no claim about all rational approximations to pi.
