# Exact moment cancellation: periods, residual freedom, and a joint bound

2026-11-02. This note concerns the user's constructive redesign, not the paused prime-pool route. It gives explicit generators and joint arithmetic/analytic upper certificates. It does not prove that the quality tends to zero.

## 1. Prescribed moments give an exact parameter-transfer rule

Fix m>=1, c=4m+1, d=m+1, and

\[
w(t)=t^{m-1/2}(1-t)^{2m}.
\]

The original integral is

\[
I_H=\frac12\int_0^1\frac{w(t)H(t)}{(1+t)^c}\,dt,
\qquad \deg H\le d.
\]

For an integer 0<=L<=d impose the L rational linear conditions

\[
\int_0^1w(t)H(t)t^j\,dt=0\qquad(0\le j<L).
\tag{1}
\]

They are equivalent to cancelling the first L terms of any polynomial kernel expansion based on an affine coordinate of t, including (2t-1)/3. An explicit parametrization of their **entire rational solution space** is

\[
\boxed{H(t)=w(t)^{-1}\frac{d^L}{dt^L}
\left[t^{m+L-1/2}(1-t)^{2m+L}G(t)\right],
\quad G\in\mathbb Q[t],\quad\deg G\le d-L.}
\tag{2}
\]

The right side is a polynomial with rational coefficients. Endpoint derivatives through order L-1 vanish, so integration by parts proves (1). The moment functionals are independent: their Gram matrix on degree<L polynomials is positive definite because w>0 in (0,1). Thus the solution space has dimension d+1-L. The map in (2) is injective: for a nonzero highest term of G of degree g, its image has degree L+g with nonzero leading coefficient (the relevant falling factorial has half-integer nonzero factors). The domain has dimension d+1-L, proving surjectivity onto (1).

Integrating by parts L more times gives the exact identity

\[
\boxed{I_H=(c)_L\int_0^1
\frac{x^{2A}(1-x^2)^BG(x^2)}{(1+x^2)^C}\,dx,}
\tag{3}
\]

where (c)_L=c(c+1)...(c+L-1) and

\[
A=m+L,\quad B=2m+L,\quad C=4m+1+L,
\quad C-A-B=d-L.
\]

Rational overall scaling of (2) is allowed and scales both sides of (3). Clearing all rational denominators yields an integer antiderivative input. Such scaling leaves the resulting reduced p/q and quality unchanged. The rational and pi period coefficients transfer by exactly the same factor in (3). One can see this from unique Hermite reduction plus vanishing integration-by-parts boundaries; alternatively equality of the rational-plus-pi integrals and irrationality of pi imply equality of both rational coefficients.

This is a repeatable cancellation ladder. It also records its cost: increasing L transfers endpoint/pole orders and reduces the remaining multiplier degree. Moment cancellation alone is not a free new improvement.

## 2. A maximal-cancellation generator valid at every scale

Take L=m+1 and G=1. Multiplying (2) by 2^L gives integer coefficients; divide by their gcd to make H primitive. An explicit Bernstein expression before that content division is

\[
\widetilde H(t)=\sum_{j=0}^L(-1)^j\binom Lj2^j
(2m+L)_{\underline j}
\left[\prod_{k=0}^{L-j-1}(2m+2L-1-2k)\right]
t^j(1-t)^{L-j}.
\]

Here empty products are one. This prescribes all coefficients without pi, target fractions, LLL, or a search. If h is the content of this integer polynomial, the factor K=2^L/h gives H=K times the Rodrigues expression. Its integral equals K(c)_L times the constant-multiplier kernel with

\[
\boxed{(A,B,C)=(2m+1,3m+1,5m+2),\qquad C=A+B.}
\tag{4}
\]

Thus its reduced rational approximant is exactly the one from that shifted constant kernel. This is an identity for the complete period pair, not just a comparison of error sizes.

## 3. Exact period extraction without a matrix inverse

For an even numerator N(x)=n(x^2), deg n<=c, first expand

\[
n(t)=\sum_{j=0}^c\nu_j(1+t)^j,\qquad
\nu_j=\sum_{k=j}^cn_k\binom kj(-1)^{k-j}.
\]

Use R(x)=x sum_{j=0}^{c-1} rho_j(1+x^2)^j. Put rho_(-1)=0 and compute

\[
\rho_j=\frac{\nu_j-(2j-2c+1)\rho_{j-1}}{2(c-j-1)}
\quad(0\le j\le c-2),
\]
\[
\boxed{\alpha=\nu_{c-1}+\rho_{c-2},\qquad
\rho_{c-1}=\nu_c,\qquad R(1)=\sum_{j=0}^{c-1}2^j\rho_j.}
\tag{5}
\]

These formulas are exact rational arithmetic. Their derivation is the triangular identity

\[
L_c[x(1+x^2)^j]=2(c-j-1)(1+x^2)^j
+(2j-2c+3)(1+x^2)^{j+1}.
\]

Clear the lcm of the rho and alpha denominators, then divide the full (R,alpha) coefficient vector by its content. This gives primitive integer coordinates in the user's family. The reduced approximation is

\[
p/q=-R(1)/(2^{c-3}\alpha),\quad q>0,
\]

and can be computed before any numerical value of pi is introduced.

## 4. A universal exponential arithmetic certificate

Let now n(t)=t^A(1-t)^B and C=A+B, with A,B>=1. Let alpha and beta be the rational coefficients in

\[
J=\int_0^1\frac{x^{2A}(1-x^2)^B}{(1+x^2)^C}\,dx
=\beta+\alpha\pi/4.
\]

Set D0=2^(2C-2), L_C=lcm(1,...,C-1). Then

\[
\boxed{\alpha=a_0/D0,\quad
\beta=b_0/(D0L_C),\quad a_0,b_0\in\mathbb Z.}
\tag{6}
\]

Proof: subtract n_C(1+t)^C, n_C=(-1)^B, from n(t). The remainder has degree<C. With x=tan theta, each remaining monomial becomes sin^(2j)(theta) cos^(2C-2j-2)(theta). Its finite Fourier expansion has denominator D0 and only cos(2k theta), k<=C-1. Integrating from 0 to pi/4 gives a constant pi term and rational terms with denominators k; sin(k pi/2) is 0 or +-1. The subtracted term integrates to n_C. This proves (6). In fact D0 L_C also clears the full odd rational primitive R: sin(2k arctan x) has an integer-polynomial numerator over (1+x^2)^k.

Consequently, if alpha!=0, the fully reduced denominator satisfies

\[
\boxed{q\le L_C|a_0|=D0L_C|\alpha|.}
\tag{7}
\]

There is no degeneracy at any scale in (4). To see this, put s=sin^2 theta. The finite trigonometric polynomial after subtraction is

\[
P(s)=\frac{s^A(1-2s)^B-(-1)^B}{1-s}.
\]

Its mean is alpha:

\[
\alpha=\frac2\pi\int_0^{\pi/2}P(\sin^2\theta)\,d\theta.
\]

For 0<s<1, |s^A(1-2s)^B|<1. Thus P has strict sign (-1)^(B+1), proving alpha!=0. Also |P|<=A+2B by the mean-value theorem applied to g(s)=s^A(1-2s)^B, since |g'|<=A+2B on [0,1]. Hence

\[
\boxed{0<|\alpha|\le A+2B.}
\tag{8}
\]

## 5. Joint quality upper bound and its limitation

For C=A+B the positive kernel has the rational maximum

\[
J\le\frac{A^AB^B}{2^AC^C}=:V_{A,B}.
\]

Thus the maximal-cancellation generator has the exact computable upper certificate

\[
\boxed{Q\le\frac{4q^2}{|\alpha|}V_{A,B},}
\tag{9}
\]

and the completely explicit all-scale arithmetic certificate

\[
\boxed{Q\le 2^{4C-2}L_C^2(A+2B)V_{A,B}.}
\tag{10}
\]

The latter includes the primitive period/denominator cost and does not assume an unproved coefficient estimate. It is exponentially growing in the tested parameter direction rather than a proof that Q tends to zero. The sharper exact-q bound (9) is also large for the first eight scales. A weak or increasing upper bound is not by itself a lower bound on actual quality, and is not being treated as a global failure theorem.

For m=1,...,8 the exact script generates the primitive integer data, checks (5), verifies both period coefficients in (3), proves equivalence with (4), and checks the arithmetic statement (6). The first two outputs are 11264/3585 and 5214208/1659735; their simple upper bounds are respectively about 1.44e4 and 2.02e7. These are verification outputs, not improved pi approximation candidates.

The constructive lesson is precise: prescribing polynomial moments supplies an explicit generator and a positive remainder, but maximal cancellation spends all multiplier freedom and returns a shifted constant kernel. To seek net quality improvement one should retain G or prescribe conditions that exploit the full rational kernel without merely applying this Rodrigues transfer. The joint bound must still be charged in either case.

Companion files: explicit_periods.py and explicit_periods.json. Standard-library exact arithmetic only; floating-point displayed upper values are descriptive, while the certificates are stored as rational numbers.
