# Positive endpoint expansion: beta cancellation before denominator clearing

2026-11-02. This is a constructive alternative kernel truncation for the same a=m,b=2m,c=4m+1 family. It gives automatic nonzero periods, a common-denominator envelope without a growing power of three, and exact signed moment-tail bounds. The bounded controls below are poor in actual quadratic quality; the aim of recording this rule is its proved arithmetic mechanism and explicit all-scale error bill, not a record claim.

## 1. Use the endpoint Bernstein powers

For m>=1 put c=4m+1, w(x)=x^(2m)(1-x^2)^(2m), f(t)=(1+t)^(-c). Expand at t=1:

\[
S_K^+(t)=2^{-c}\sum_{k=0}^K
\binom{c+k-1}{k}\left(\frac{1-t}{2}\right)^k,
\qquad K\ge0.
\]

These are the endpoint powers of a Bernstein basis. Every summand is nonnegative on [0,1], with strictly positive constant term. Set

\[
\mu_j=\int_0^1w(x)x^{2j}S_K^+(x^2)\,dx,
\quad j=0,1,\qquad H_K(t)=\mu_1-\mu_0t.
\]

As in the full-truncation construction, the complete truncated integral vanishes exactly: integral w H_K S_K^+ = 0. The coefficients are explicit rational finite sums; no pi value, target fraction or LLL choice enters them.

## 2. The beta moment cancels the apparent dyadic growth

Let T_(j,k) be the kth summand in mu_j. The exact beta moment gives

\[
\boxed{T_{j,k}=2^{-2m-1}\binom{4m+k}{k}
\frac{(2m+k)!}{\prod_{r=0}^{2m+k}(2m+2j+2r+1)}.}
\]

In particular the apparent 2^(-k) in the polynomial coefficient cancels the 2^k in the integrated beta moment. Each T_(j,k), and each mu_j, has dyadic denominator exponent at most 2m+1, independently of K.

The odd denominator follows also from expanding the polynomial before integrating: it divides the odd lcm of the odd monomial denominators through degree 3m+K+1. Consequently

\[
\boxed{L_K^+:=\operatorname{lcm}(\operatorname{den}\mu_0,
\operatorname{den}\mu_1)
\mid2^{2m+1}\operatorname{odd}\operatorname{lcm}
(1,\ldots,6m+2K+3).}
\]

The previous centered expansion had an additional universal factor 3^(4m+K+1). This is an improvement in the common-denominator *envelope*, not a claim that the actual reduced denominators improve by the ratio of those envelopes.

One efficient exact recurrence is

\[
T_{0,0}=\frac{(2m)!}{2^{2m+1}\prod_{r=0}^{2m}(2m+2r+1)},
\quad
\frac{T_{j,k+1}}{T_{j,k}}=
\frac{(4m+k+1)(2m+k+1)}{(k+1)(6m+2j+2k+3)},
\]
\[
T_{1,k}=T_{0,k}\frac{2m+1}{6m+2k+3}.
\]

Thus both moments can be accumulated with rational recurrence rather than expanding a high-degree polynomial.

## 3. Positivity gives a nonzero pi coefficient for every K

Use the already proved closed base periods

\[
J_j=\int_0^1w(x)x^{2j}f(x^2)dx
=\beta_j+\alpha_j\pi/4,
\]
\[
\alpha_0=\binom{2m}{m}/64^m>0,\quad
\alpha_1=r_m\alpha_0,\quad r_m=(6m+1)/(2m-1)>3,
\quad\beta_0=0,\quad\beta_1=-\frac1{(2m-1)\binom{4m}{2m}}.
\]

Because S_K^+ is positive, 0<mu_1<mu_0. Therefore

\[
\boxed{\alpha_H=\alpha_0(\mu_1-r_m\mu_0)<0}
\]

for every m>=1 and every K>=0. There is no need to demand a small truncation error merely to establish a nonzero period. More generally this observation applies to any positive polynomial approximation S, including every even-degree centered Taylor truncation.

Let L=L_K^+, u_j=L mu_j, A=binom(2m,m), C=binom(4m,2m). The generated reduced fraction is exactly the gcd reduction of

\[
\boxed{\frac pq=\frac{2^{6m+2}u_0}
{AC[(6m+1)u_0-(2m-1)u_1]}.}
\]

Every denominator in this displayed ratio is positive before reduction.

## 4. One-sided cancellation and sharp rational tails

Put z=(1-t)/2. Direct differentiation gives

\[
\frac{d}{dz}\left[(1-z)^c\sum_{k=0}^K
\binom{c+k-1}{k}z^k\right]
=-c\binom{c+K}{K}z^K(1-z)^{c-1}<0.
\]

Hence S_K^+(t)/f(t) is strictly increasing in t. Positive covariance with t under the measure w(x)f(x^2)dx gives

\[
I_H=\mu_1J_0-\mu_0J_1>0.
\]

For a direct proof, symmetrize the double integral: its integrand is the positive weight times (t-s)(S_K^+(t)/f(t)-S_K^+(s)/f(s)). Since alpha_H<0, every generated fraction is above pi.

Let delta_j=J_j-mu_j>0. Then

\[
0<I_H=\mu_1\delta_0-\mu_0\delta_1<\mu_1\delta_0.
\]

Write z=k+1 and d=6m+2j+1. The term ratio is

\[
R_{j,k}=\frac12+
\frac{6m-2j-1}{2(2z+d)}+
\frac{8m^2}{z(2z+d)}.
\]

For j=0,1 it is strictly decreasing in k. At k=3m it is less than one: denominator minus numerator in the unreduced expression is m^2+(6j+9)m+2j+2. Thus for K>=3m,

\[
\delta_j\le\frac{T_{j,K+1}}{1-R_{j,K+1}},
\qquad
\boxed{Q\le\frac{4q^2\mu_1}{|\alpha_H|}
\frac{T_{0,K+1}}{1-R_{0,K+1}}.}
\]

The q in this bound is the actual gcd-reduced denominator, not an unreduced coefficient estimate.

## 5. A fixed all-scale schedule with ordinary convergence

Take K=4m. The exact ratio satisfies R_(0,4m)<=6/7, whence delta_0<=6T_(0,4m). Also

\[
T_{0,4m}=2^{-8m-1}\binom{8m}{4m}
\int_0^1x^{2m}(1-x^2)^{6m}dx
\le\tfrac12\left(\frac{6^6}{7^7}\right)^m.
\]

Since |alpha_H|>2 alpha_0 mu_0 and alpha_0>=16^(-m)/(2m+1),

\[
\boxed{0<\frac{p_m}{q_m}-\pi
\le6(2m+1)\left(\frac{746496}{823543}\right)^m\to0.}
\]

Reduced denominators tend to infinity because pi is irrational. This rule cuts the fixed truncation order to 4m and improves its denominator envelope, but these facts alone do not give quadratic-quality decay.

## 6. Bounded design control and scope

`endpoint_truncation.py` evaluated the 24 prescribed pairs m=1,...,6 and K in {m,2m,3m,4m}. It generates solely from rational beta recurrence, verifies the common-denominator divisor, and then uses a rigorous Machin interval to certify the output quality. The endpoint candidates are poor: the smallest actual quality in this bounded control is about7698.5 at m=1,K=2. The sharp tail bounds are also recorded when available.

Accordingly this note does not recommend a larger scan of this plain endpoint truncation. The useful construction mechanisms are the beta-integrated denominator cancellation, automatic nonzero period for positive approximants, rational moment-tail bounds, and explicit cost formula. A further redesign can reuse those ingredients while modifying or telescoping the moment pair. No conclusion about pi bad/nonbad approximability is asserted.
