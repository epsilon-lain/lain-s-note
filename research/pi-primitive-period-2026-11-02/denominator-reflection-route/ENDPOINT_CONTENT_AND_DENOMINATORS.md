# Endpoint splitting, coefficient content, and the least denominator

2026-11-02. This note concerns the same sparse one-dimensional two-frequency raw-column family as the preceding Hankel analysis. It proves an exact least-denominator formula in even dimensions and bounds the odd-dimensional least denominator within a polynomial factor. In particular, the least-coefficient-denominator-cleared actual-period value tends to zero. This is different from the earlier valid but larger row-factorial clearing, whose actual-period value grows. No rational approximation or bad-approximability conclusion follows.

## 1. Objects and the result

For n>=2, use the complete two-center packet and the 2n raw columns from the preceding note:

\[
P_n(Z)=\det\left[[z^s]e^{hz}(jZ+z)^a\right]_{
 j=0,1;\ 0\le s<n\ ;\ h=0,1;\ 0\le a<n}.
\]

Rows are center then time; columns are frequency then spatial degree. For n>=3 these are legal with H=2n,w0=1,w1=2,v0=2,v1=2n,theta=1/n. With this parameter choice n=2 is an algebraic boundary case: K*theta=1, whereas the original admissibility condition is strict. All polynomial identities below remain valid at n=2. The family has fixed transverse dimension and ambient density tending to zero, rather than the critical high-dimensional half-density parameters.

Define A_0=1 and

\[
A_k=\prod_{a=0}^{k-1}a!\quad(k\ge1),\qquad
r=\lfloor n/2\rfloor,\quad s=\lceil n/2\rceil,
\]

\[
G_n=A_r^2A_s^2,\qquad D'_n=\frac{A_n^2}{G_n}.
\]

Let d_n be the least positive integer for which every coefficient of d_nP_n is an integer. The main assertions are

\[
\boxed{D'_nP_n\in\mathbb Z[Z],\qquad d_n\mid D'_n.}
\tag{1}
\]

For even n=2r,

\[
\boxed{\deg P_{2r}=2r^2,\qquad
\operatorname{LC}(P_{2r})=\frac{(-1)^r}{D'_{2r}},\qquad
d_{2r}=D'_{2r}.}
\tag{2}
\]

For odd n=2r+1>=3,

\[
\boxed{\deg P_{2r+1}=2r^2+2r-1,\qquad
\operatorname{LC}(P_{2r+1})
=\frac{(-1)^{r+1}4r(r+1)}{D'_{2r+1}}.}
\tag{3}
\]

Consequently, with c_r=4r(r+1),

\[
\boxed{\frac{D'_{2r+1}}{\gcd(D'_{2r+1},c_r)}
\ \mid\ d_{2r+1}\ \mid\ D'_{2r+1}.}
\tag{4}
\]

The exact odd least-denominator formula is not asserted. For r>=2, c_r divides D'_(2r+1), so (4) also yields the simpler lower bound D'/c_r<=d<=D'. The n=3 case is P_3=2Z^3 and d_3=1.

## 2. An exact signed endpoint functional

Write t=1/Z. The preceding formal Hankel identity says

\[
A_n^2P_n(Z)=Z^{n^2}H_n(1/Z),\qquad
H_n(t)=\det[\widehat\mu_{i+j}(t)]_{i,j=0}^{n-1},
\]

where muhat_0=0 and muhat_k=t-kt*muhat_(k-1). For every polynomial f define

\[
\mathcal L_t f
=t\int_0^\infty e^{-v}\{f(1-tv)-f(-tv)\}\,dv.
\tag{5}
\]

These integrals are finite sums of factorial moments and therefore polynomials in t. Taylor expansion gives

\[
\mathcal L_t f
=t\sum_{j\ge0}(-t)^j\{f^{(j)}(1)-f^{(j)}(0)\}.
\]

For f(u)=u^k this is exactly muhat_k. Expanding the Hankel determinant by the usual alternant identity, and grouping variables assigned to endpoint 1 or endpoint 0, gives an exact endpoint decomposition.

Let V=(V_1,...,V_k) have probability density proportional to

\[
e^{-\sum V_i}\prod_{i<j}(V_j-V_i)^2
\quad\text{on }[0,\infty)^k,
\]

and independently let W have the corresponding l-variable density. Empty groups have expectation one. Put

\[
R_{k,l}(t)=\mathbb E\prod_{i=1}^k\prod_{j=1}^l
\bigl(1-t(V_i-W_j)\bigr)^2.
\]

The Laguerre normalization is

\[
\int_{[0,\infty)^k}e^{-\sum V_i}\Delta(V)^2\,dV
=k!\det[(i+j)!]_{i,j=0}^{k-1}=k!A_k^2.
\]

For example, the factorial determinant follows by factoring i! from row i and j! from column j, leaving the determinant-one binomial matrix binom(i+j,i). Thus

\[
\boxed{H_n(t)=\sum_{k+l=n}(-1)^l A_k^2A_l^2
t^{k^2+l^2}R_{k,l}(t).}
\tag{6}
\]

The exponent includes one t per functional, the within-endpoint Vandermonde powers, and no leading t from the cross-endpoint distances. Formula (6) is polynomial algebra with finite factorial integrals; it is not an asymptotic expansion of a signed probability measure.

## 3. Why the normalized cross coefficients are integers

Every coefficient of the product defining R_(k,l) is an integer polynomial, symmetric separately in V and W. It is enough to prove that the expectation of every symmetric integer polynomial in a Laguerre group is an integer.

The Schur polynomials form an integral basis of symmetric integer polynomials. One direct justification uses their alternant quotients: division by the Vandermonde gives integer coefficients, and the expansion in monomial symmetric functions is triangular with diagonal one. Thus integer change of basis uses no division.

For a partition lambda with at most k parts, apply the alternant formula and integrate using integral x^a exp(-x) dx=a!. If a_i=lambda_i+k-i, 1<=i<=k, one obtains

\[
\begin{aligned}
\mathbb E s_\lambda(V)
&=\frac{\det[(a_i+j)!]_{1\le i\le k,\ 0\le j<k}}
{\det[(k-i+j)!]_{1\le i\le k,\ 0\le j<k}}\\
&=s_\lambda(1,\ldots,1)
\prod_{i=1}^k\frac{(\lambda_i+k-i)!}{(k-i)!}\in\mathbb Z.
\end{aligned}
\tag{7}
\]

To see the second equality, factor a_i! from each numerator row; its remaining j-th entry is the monic degree-j polynomial (a_i+1)...(a_i+j). The resulting determinant is the Vandermonde in the a_i. The denominator is the same computation for lambda=0. Their Vandermonde ratio is the alternant limit s_lambda(1^k), an integer. Every remaining factorial ratio is an integer as well.

Independence of V and W now proves

\[
\boxed{R_{k,l}(t)\in\mathbb Z[t],\qquad R_{k,l}(0)=1.}
\tag{8}
\]

## 4. Common content and the extreme coefficients

For k+l=n, the factorial product A_kA_l is divisible by A_rA_s. Indeed, after interchanging k,l assume k<=r and l>=s. Then

\[
\frac{A_kA_l}{A_rA_s}
=\prod_{j=0}^{r-k-1}\frac{(s+j)!}{(k+j)!}\in\mathbb Z.
\]

Equations (6) and (8) imply that every coefficient of H_n is divisible by G_n. Since its degree is at most n², reversing powers in the identity A_n²P_n=Z^(n²)H_n(1/Z) proves (1). In particular D'_n is an integer.

For even n=2r, k²+l² has the unique minimum 2r² at k=l=r. Its coefficient in (6) is (-1)^r A_r^4. This proves (2), including minimality of the common denominator because its highest coefficient is exactly a signed reciprocal of D'.

For odd n=2r+1, the two balanced splits (r,r+1) and (r+1,r) have the same exponent 2r²+2r+1 and opposite leading signs, so those coefficients cancel. The next unbalanced splits start four powers later.

The Laguerre trace mean is k²: scaling the weight exp(-a*sum V_i) makes its normalization a^(-k²), and differentiating at a=1 gives this mean. Therefore

\[
R'_{k,l}(0)=-2(lk^2-kl^2)=2kl(l-k).
\]

At the next power of t the two balanced contributions sum to

\[
(-1)^{r+1}4r(r+1)A_r^2A_{r+1}^2\ne0.
\]

This proves (3). Reducing that rational leading coefficient gives (4). For r>=2, write explicitly

\[
B=\frac{A_{2r+1}}{A_rA_{r+1}}
=r!(r+1)!\prod_{j=2}^{r}\frac{(r+j)!}{j!}.
\]

All displayed ratios are integers, and the first two factorials already show 2r(r+1) divides B. Hence 4r(r+1) divides B²=D', as used in the simpler bound.

## 5. The least-denominator-cleared value tends to zero

The preceding actual-period theorem gives

\[
\frac16U_n<|P_n(2\pi i)|\le U_n,
\qquad U_n=(2\pi)^{n^2}\frac{A_n^2}{A_{2n}}.
\]

Using log A_n=(n²/2)log n-(3n²/4)+O(n log n), equations (2)--(4) give

\[
\boxed{\log d_n
=\frac12n^2\log n+
\left(\frac12\log2-\frac34\right)n^2+O(n\log n).}
\tag{9}
\]

For odd n the gap between the upper and lower logarithmic bounds is only O(log n). Combining (9) with the two-sided actual-period estimate proves

\[
\boxed{\log\bigl(d_n|P_n(2\pi i)|\bigr)
=-\frac12n^2\log n+
\left[\frac34+\log\left(\frac\pi{\sqrt2}\right)\right]n^2
+O(n\log n).}
\tag{10}
\]

In particular d_n|P_n(2pi i)| tends to zero. Even the explicit valid factor D'_n has the same asymptotic behavior, so computing the exact odd content is unnecessary for this conclusion.

There is no conflict with the earlier finite table: all n=4,...,10 values exceeded one, but a finite initial segment does not determine the asymptotic behavior. There is also no conflict with growth after the larger A_n² clearing. The ratio A_n²/D'_n=G_n has logarithm of order (n²/2)log n and is not negligible.

The number d_nP_n(2pi i) is not an integer. The polynomials d_nP_n are integer polynomials of increasing degree, and small values at a transcendental point do not alone supply exceptionally good rational approximants. At a rational imaginary center 2ip/q, an additional q power up to the actual polynomial degree must be paid, and the change from that center to 2pi i must be bounded. Neither the needed comparison nor a high-dimensional critical construction is proved here. Thus (10) removes this particular coefficient-denominator obstruction, but does not prove or disprove pi's bad approximability.

## Verification

The standard-library companion `verify_endpoint_denominators.py` passed all checks: 25 endpoint splits for n=2,...,6; 188 Laguerre Schur expectations by exact factorial determinants; complete Hankel coefficient and least-denominator checks for n=2,...,10, including 27 evaluations against the original 2n-by-2n raw determinant; and 1,888 balanced-factorial divisibility checks through n=60. Results are recorded in `endpoint_denominator_verification.json`. A separate independent audit and verifier in [independent audit](../minimal-denominator-general/structure/INDEPENDENT_CONTENT_AUDIT.md) also passed. The all-n statements follow from (5)--(10), not from fitting finite coefficients.
