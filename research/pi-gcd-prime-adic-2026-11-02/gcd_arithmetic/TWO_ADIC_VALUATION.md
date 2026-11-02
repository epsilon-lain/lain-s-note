# Exact 2-adic valuation and the remaining odd-prime surplus

2026-11-02. This concerns precisely the integer recurrence A_0=A_1=1 and A_{n+1}=(2n+1)A_n+n²A_{n-1}, together with B_0=0, B_1=1 obeying the same recurrence and g_n=gcd(A_n,4B_n). It does not assume any Diophantine property of pi.

## Theorem

For every k>=0,

v_2(A_{4k})=2k,
v_2(A_{4k+1})=2k,
v_2(A_{4k+2})=2k+2,
v_2(A_{4k+3})=2k+3+v_2(k+1).

Here v_2 is the additive 2-adic valuation and A_0=1 has valuation zero. Consequently

v_2(g_n)-v_2(n!) <= -n/2+2 log_2(n+1)+1.

If E_n=product over odd primes p of p^{max(0,v_p(g_n)-v_p(n!))}, then

g_n sqrt(n)/n! <= 2(n+1)^2 sqrt(n) 2^{-n/2} E_n.

Thus any subsequence making g_n sqrt(n)/n! unbounded must overcome an exponentially large, universally present 2-adic loss through odd-prime surplus. In particular E_n=exp(o(n)) would exclude success for this entire specified family. This note does not prove that missing odd-prime bound.

## 1. Binomial reversal

Put a_n=A_n/n! and c_m=binom(2m,m)/2^m. The established generating function, or direct coefficient extraction, gives

a_n=sum_{j=0}^{floor(n/2)} binom(n,2j)binom(2j,j)/2^j.

Reverse j=m-l. This yields

a_{2m}=c_m S_m,
S_m=sum_{l=0}^m l! binom(m,l)^2/(2l-1)!!,

a_{2m+1}=(2m+1)c_m T_m,
T_m=sum_{l=0}^m l! binom(m,l)^2/(2l+1)!!.

The convention is (-1)!!=1. Every double factorial denominator here is odd, so both sums are 2-adic integers. Legendre's factorial formula gives v_2(c_m)=s_2(m)-m, where s_2(m) is the number of ones in the binary expansion of m.

## 2. The even-index sum

Modulo 2, S_m=1+m². Thus S_m is odd for even m. For odd m, determine it modulo 4: all l>=4 terms vanish, since v_2(l!)>=3. The l=0,1 contribution is 1+m²=2 mod4. The l=2,3 terms are respectively 2 binom(m,2)^2/3 and 2 binom(m,3)^2/5. Since binom(m,2)=binom(m,3) mod2 for odd m (by the first two binary digits, or Pascal's identity), these two contributions sum to zero modulo4. Therefore

v_2(S_m)=0 if m is even, and 1 if m is odd.

For even m, T_m is odd as well: again only l=0,1 survive modulo2, giving 1+m²=1.

## 3. A convergent 2-adic identity at m=-1

For a 2-adic integer x, define

T(x)=sum_{l>=0} l! binom(x,l)^2/(2l+1)!!.

The series converges uniformly on the 2-adic integers because binom(x,l) is a 2-adic integer and v_2(l!)->infinity. At a nonnegative integer x=m it terminates and equals T_m. At x=-1, binom(-1,l)=(-1)^l, so

T(-1)=sum_{l>=0} l!/(2l+1)!!=0 in Q_2.

Here is an independent proof of the last identity. Formally integrate the geometric expansion

1/[1-2x(1-x)] = sum_{l>=0} 2^l x^l(1-x)^l

from 0 to1, using the algebraic integral of polynomials. Its l-th integral is 2^l(l!)²/(2l+1)!=l!/(2l+1)!!. This operation commutes with summation 2-adically: polynomial coefficients of the l-th term have valuation at least l, and integration reduces it by at most log_2(2l+1).

In Q_2(i), factor the denominator as (1-(1+i)x)(1-(1-i)x). A primitive of its reciprocal is

[log(1-(1-i)x)-log(1-(1+i)x)]/(2i).

The logarithm series converge at x=1 since v_2(1+i)=v_2(1-i)=1/2>0. The endpoint values are log(i) and log(-i), both zero: the logarithm homomorphism on 1+the maximal ideal kills roots of unity. Hence the integral and T(-1) equal zero. Equivalently, this final fact follows by applying log(u^4)=4log(u) to u=±i and using characteristic zero. This identity is entirely 2-adic; it makes no assertion that a real sum is zero.

## 4. Exact valuation for the odd-index sum when m is odd

Let m be a positive odd integer and r=v_2(m+1). The generalized binomial identity implies

v_2(binom(m,l)-binom(-1,l)) >= r-floor(log_2 l), for l>=1.

Indeed, expand the difference using m=(-1)+(m+1) and use binom(m+1,j)=(m+1)binom(m,j-1)/j. All other binomial coefficients are integers, including those with upper argument -1.

For P_l(x)=l!binom(x,l)^2/(2l+1)!! this gives

v_2((P_l(m)-P_l(-1))/(m+1)) >= v_2(l!)-floor(log_2 l).

For every l>=6, the right side is at least2. (Check l=6,7 directly; thereafter Legendre's formula or elementary monotonic intervals gives the inequality.) Thus, using T(-1)=0,

T_m/(m+1) = R(m) mod4,
R(m)=sum_{l=1}^5 [P_l(m)-P_l(-1)]/(m+1).

This is a finite rational polynomial identity. Set m=2u-1; exact expansion gives

R(2u-1)-2 = (4/155925)*(16u^9-240u^8+1800u^7-8200u^6+25853u^5-57235u^4+95560u^3-112150u²+106091u-133890).

The denominator is odd, so R(m)=2 mod4 for every odd m. It follows that

v_2(T_m)=1+v_2(m+1).

This supplies the potentially delicate cancellation case, with no extrapolation from finite data.

## 5. Translate back to A_n and to g_n

Use v_2(n!)=n-s_2(n), s_2(2m)=s_2(m), and s_2(2m+1)=s_2(m)+1. Combining Sections1--4 gives the four exact formulas in the theorem. Since g_n divides A_n, its 2-adic valuation is at most that of A_n. The piecewise formula gives v_2(A_n)<=n/2+log_2(n+1)+1, and s_2(n)<=log_2(n+1); subtracting v_2(n!) proves the displayed deficiency bound.

For odd primes retain only positive valuation surplus and discard all deficits. Multiplying valuations yields the stated bound involving E_n. Any stronger claim requires separate control of those odd primes.
