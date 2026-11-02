# Exact moment cancellation, its parameter shift, and a sharper kernel norm

2026-11-02. The construction goal remains a sequence of reduced approximants with quality tending to zero. This note proves what exact polynomial-moment cancellation does in the proposed family and gives a sharper unconditional analytic certificate. It does not prove a quality-zero sequence.

## 1. Full parametrization of the moment conditions

For integers a,b>=1, c>=a+b, and d=c-a-b, write

\[
I_{a,b,c}(H)=\frac12\int_0^1
\frac{t^{a-1/2}(1-t)^bH(t)}{(1+t)^c}\,dt,
\qquad H\in\mathbb Q[t],\quad\deg H\le d.
\]

For 0<=k<=d, impose the exact rational conditions

\[
\int_0^1t^{a+j-1/2}(1-t)^bH(t)\,dt=0,
\qquad 0\le j<k.                         \tag{1}
\]

These are also the conditions cancelling the first k polynomial terms in a Taylor expansion in 2t-1, since those terms span the same space as 1,t,...,t^(k-1).

Every solution of (1), with no exception, has the form

\[
\boxed{H(t)=\frac1{t^{a-1/2}(1-t)^b}
\frac{d^k}{dt^k}\left[t^{a+k-1/2}(1-t)^{b+k}P(t)\right],
\quad P\in\mathbb Q[t],\quad\deg P\le d-k.}       \tag{2}
\]

The expression is a polynomial: applying the Leibniz rule to a monomial P=t^ell gives

\[
\sum_{s=0}^k\binom ks
(a+k+\ell-\tfrac12)_{\downarrow k-s}
(-1)^s(b+k)_{\downarrow s}
t^{\ell+s}(1-t)^{k-s}.                         \tag{3}
\]

Here z_(downarrow h)=z(z-1)...(z-h+1), with empty product one. Thus 2^kH has integer coefficients whenever P does.

To prove completeness, integration by parts k times shows (2) satisfies (1); all boundary terms vanish. The map P->H is injective: for a leading term t^ell, its resulting leading coefficient is a nonzero falling factorial, because a+b+2k+ell-1/2 is not an integer. Finally the k moment conditions have rank k, since their Gram matrix on polynomials of degree<k is positive definite for the positive weight t^(a-1/2)(1-t)^b. Both spaces have dimension d-k+1, proving the asserted bijection.

## 2. The resulting integral and reduced fraction are exactly a shifted construction

Integration by parts gives the exact identity

\[
\boxed{I_{a,b,c}(H)=(c)_k I_{a+k,b+k,c+k}(P),}
\qquad(c)_k=c(c+1)\cdots(c+k-1).               \tag{4}
\]

In the original rational primitive coordinates write I=beta+alpha*pi/4. Equation (4) implies alpha_old=(c)_k alpha_new and beta_old=(c)_k beta_new, by uniqueness of a rational-plus-rational-pi representation. In particular, when alpha is nonzero,

\[
\boxed{-4\beta_{\rm old}/\alpha_{\rm old}
=-4\beta_{\rm new}/\alpha_{\rm new}.}           \tag{5}
\]

These are the same reduced p/q, with the same true quality. Clearing the rational coefficients to reach the integer antiderivative lattice multiplies alpha,beta,I by one common integer and does not change (5). The rising factorial in (4) is a common scale; it must not itself be charged as an actual denominator increase after reduction.

For the proposed a=m,b=2m,c=4m+1,d=m+1, maximal cancellation k=m+1 leaves constant P and parameters

\[
(a',b',c')=(2m+1,3m+1,5m+2).
\]

In this endpoint case c'=a'+b'. The earlier odd-period image theorem assuming c'>=a'+b'+1 cannot be applied to it. Partial cancellation k=m leaves residual linear P and (a',b',c')=(2m,3m,5m+1).

Thus this particular exact-moment proposal is fully explicit and can be implemented without searching, but it reparametrizes higher endpoint/pole orders. It does not supply a new cancellation advantage for free. This is a statement about conditions (1); it does not cover every possible signed multiplier construction.

## 3. A sharper constructive bound in the original family

The fixed proposed family has a particularly simple angular kernel. With x=tan(theta),

\[
\boxed{I_m(H)=16^{-m}\int_0^{\pi/4}
\sin^{2m}(4\theta)H(\tan^2\theta)\,d\theta.}     \tag{6}
\]

For H=1 the integral is exactly

\[
I_m(1)=\frac{\pi}{4\,64^m}\binom{2m}{m}.
\tag{7}
\]

Its rational boundary is zero and its pi period coefficient is binom(2m,m)/64^m, which is nonzero. Constant H in these original parameters therefore gives the trivial approximant 0/1, with quality pi; this is distinct from the alpha=0 degeneracy.

For arbitrary H, let h=sum|H_j|. On 0<=tan²theta<=1, |H(tan²theta)|<=h, so the exact positive kernel mass proves

\[
\boxed{Q_m\le
\pi\frac{q_m^2h_m}{|\alpha_m|}
\frac{\binom{2m}{m}}{64^m}
<\frac{22}{7}\frac{q_m^2h_m}{|\alpha_m|}
\frac{\binom{2m}{m}}{64^m}.}                  \tag{8}
\]

The last certificate is entirely rational. By the central binomial asymptotic, its analytic factor is O(16^(-m)/sqrt(m)), sharper than the previous O((2/27)^m) factor. A sufficient joint height target becomes

\[
q_m^2h_m/|\alpha_m|=o(16^m\sqrt m).
\tag{9}
\]

Condition (9) is not proved for a sequence of chosen integer multipliers. The sharper norm does not retain additional signed cancellation beyond (6).

## 4. Exact Fourier period coordinates

This supplies a repeatable arithmetic extraction rule without pi input. Let N(t) have integer coefficients and degree<=c, and put n=c-1,D=4^n,L=lcm(1,...,n). If N_c is its top coefficient, subtract N_c(1+t)^c to obtain tildeN of degree<=n. Its angular integrand is a trigonometric polynomial:

\[
D\sum_{j=0}^{n}\widetilde N_j
\sin^{2j}\theta\cos^{2(n-j)}\theta
=A_0+\sum_{k=1}^{n}A_k\cos(2k\theta),
\quad A_k\in\mathbb Z.
\]

Integrality follows by using sin²theta=(2-z-z^(-1))/4 and cos²theta=(2+z+z^(-1))/4 with z=exp(2i theta). The symmetric Laurent coefficients also show that each A_k for k>=1 is even. Define

\[
B=2DLN_c+\sum_{k=1}^{n}\frac Lk A_k\sin(k\pi/2)\in\mathbb Z.
\]

Only the values 0,1,-1 of the displayed sine are needed. Then

\[
I=\frac{A_0\pi}{4D}+\frac B{2DL},\qquad
\boxed{\frac pq=-\frac{2B}{LA_0},\quad
q=\frac{|LA_0|}{\gcd(2B,LA_0)}}.              \tag{10}
\]

The formula requires A_0!=0. It makes the actual reduction explicit; q is not inferred from a norm or from the unreduced period coefficient. This is equivalent to rational antiderivative extraction, not a claim of new approximation quality.

## Verification

The companion standard-library verifier checks (1)--(5) on exact rational polynomial data, including all residual monomial bases in the stated bounded range, verifies (7) against independent rational antiderivative reduction, and compares the Fourier periods (10) against that reduction. Its finite computations audit the identities; the universal statements follow from the proofs. No pi numerical value is supplied to multiplier generation or extraction.
