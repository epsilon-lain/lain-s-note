# Independent review: profiled remainders, beta denominators, and parity-free primitives

2026-11-02. This is an internal independent proof and exact-arithmetic audit, not external peer review. The three constructive mechanisms pass within the scopes stated below. None proves quadratic quality tending to zero.

## 1. The exact differential residual and its pointwise bound

For the centered polynomial S_K, direct coefficient cancellation gives

\[
(1+t)S_K'+cS_K=\gamma_K(2t-1)^K,
\qquad
\gamma_K=\frac{c2^c(-1)^K\binom{c+K}{K}}{3^{c+K}}.
\]

Multiplying the error equation by (1+t)^(c-1), with zero error at t=1/2, proves the author's exact integral representation. On the left half of the interval, the relevant antiderivative comparison is

\[
\frac{d}{du}\{u^{K+1}(3-u)^{c-1}\}
\ge\frac{2K-c+3}{2}u^K(3-u)^{c-1}.
\]

On the right half, replacing 3-u by 3+u gives a still stronger bound. Therefore, whenever 2K-c+3>0,

\[
|S_K(t)-(1+t)^{-c}|\le
\frac{|\gamma_K|}{2K-c+3}\frac{|2t-1|^{K+1}}{1+t}.
\]

The sufficient condition K>=c used initially is not necessary. In the user's c=4m+1 family, K=2m already meets the required strict positivity condition.

The profile integral certificate with Z_j is correct, including its coefficient 4q^2/|alpha| and the use of the actual reduced q. The exact ratio to the preceding uniform epsilon is

\[
\frac{|\gamma_K|}{(2K-c+3)\varepsilon_K}
=1+\frac{c-1}{(K+2)(2K-c+3)}.
\]

The author's K=13m bounds Z_j<=64^(-m) and ordinary-error coefficient 256/3 are correct. These improve the certificate for unchanged fractions and introduce no new multiplier-denominator cost.

## 2. A sharper independent K=2m convergence certificate

The following strengthens the coarse maximum bound in the author note. For 0<=t<=1,

\[
\boxed{t(1-t)^2(2t-1)^2\le\frac7{128}.}
\tag{1}
\]

For t<=1/2 put u=1-2t. The left side is h(u)/8, where h(u)=u^2(1-u^2)(1+u). Its derivative is

\[
h'(u)=u(1+u)(2+u-5u^2).
\]

The last factor is negative and decreasing for u>=3/4, so the global maximum of h occurs at u<=3/4. There, u^2(1-u^2)<=1/4 and 1+u<=7/4, giving h<=7/16. For t>=1/2 put u=2t-1; the left side is u^2(1-u^2)(1-u)/8<=1/32. This proves (1) without a numerical maximization.

Now choose the centered schedule K=2m. This is even, so S_K is positive on [0,1]: for a negative centered variable every Taylor term is positive; for a positive centered variable the even Taylor truncation lies above the positive target function. Therefore 0<mu_1<mu_0 and

\[
|\alpha_H|>2\alpha_0\mu_0,
\qquad\alpha_0=\binom{2m}{m}/64^m.
\]

Here E=K=2m, so (1) gives Z_0<=(7/128)^m and Z_1<=Z_0. The residual coefficient simplifies to

\[
\frac{|\gamma_{2m}|}{2}
=\frac{(6m+1)2^{4m}\binom{6m}{2m}}{3^{6m+1}}
\le\frac{6m+1}{3}.
\]

The inequality follows from the nonnegative binomial probability term
binom(6m,2m)(1/3)^(2m)(2/3)^(4m)<=1. Cancelling mu_0 between the integral bound and the period lower bound proves

\[
\boxed{\left|\pi-\frac{p_m}{q_m}\right|
\le\frac43(6m+1)(2m+1)\left(\frac78\right)^m\longrightarrow0.}
\tag{2}
\]

Thus this shorter centered schedule has nonzero period at every scale, ordinary error tending to zero, and reduced denominators tending to infinity. The denominator still must be charged in any Q=q^2 times error estimate. Equation (2) does not prove quadratic-quality decay or monotonic denominator growth.

## 3. Telescoping must retain the correct boundary

For an arbitrary rational primitive R, integration by parts gives the exact boundary term

\[
2R(1)S_K(1)-R(0)S_K(0).
\]

The originally inherited odd-R family has R(0)=0. Once even coefficients are allowed, either impose that normalization explicitly or retain the second term. The corrected author note states this condition. The residual-moment term and alpha term in its telescoping identity have the correct signs and factors.

For the specific zero-pi-period numerator w[(6m+1)-(2m-1)x^2], the closed base periods give rational integral 1/binom(4m,2m), hence R(1)=2^(c-1)/binom(4m,2m) under R(0)=0. The denominator-combination identity follows exactly. It is an arithmetic handle, not a proof that an extra common factor exists.

The centered W and positive increment recurrences also pass: they follow by integrating the Pearson identity of the transformed beta weight. The coefficient signs imply (-1)^kW_(j,k)>0. All increment coefficients are positive for j=0,1,m>=1. No monotonicity of the ratio mu_1/mu_0 is inferred from that fact.

## 4. Positive endpoint truncation

The endpoint beta formula cancels its apparent growing dyadic denominator exactly. Term by term the 2-adic denominator exponent is at most 2m+1; odd denominator control comes independently from polynomial expansion followed by rational monomial integration. Combining the two proves

\[
L_K^+\mid2^{2m+1}\operatorname{odd}\operatorname{lcm}
(1,\ldots,6m+2K+3).
\]

The positive polynomial gives 0<mu_1<mu_0 and hence a nonzero negative pi coefficient for every K>=0. The derivative of S_K^+/f has the stated sign, so its covariance with t is positive; the resulting generated fraction lies above pi. The one-sided tail certificate therefore has the stated strict sign and upper factor.

The term ratio decomposition, monotone decrease, geometric tail bound for K>=3m, and R_(0,4m)<=6/7 all pass. They prove the author's ordinary convergence bound

\[
0<p_m/q_m-\pi\le6(2m+1)
(746496/823543)^m\to0.
\]

The denominator-envelope improvement is genuine but cannot be converted directly into a factor improvement in the actual reduced denominator. The bounded candidates being poor is not an asymptotic failure theorem.

## 5. Full parity freedom does not introduce a hidden logarithm

For deg N<=2c-2, write its odd part as xP(x^2), deg P<=c-2. The substitution t=x^2 reduces it to one half of P(t)/(1+t)^c, whose expansion in powers of 1+t has exponents at most -2. Its antiderivative is rational and contains no logarithm. The even part integrates to a rational number plus a rational multiple of pi. Thus the all-parity integral still has only these two periods.

In the enlarged degree<=2c space the precise log-free condition is [x^(2c-1)]N=0. The odd primitive columns with alpha span all even numerators; even R columns with R(0)=0 span all permitted odd numerators by a triangular matrix with diagonal 2,4,...,2c-2. This proves the stated complete image, rather than merely checking examples.

For N=x^A(1-x)^B H, the highest odd coefficient condition is H_(D-1)=B H_D, D=2c-A-B. The new fair parameters give D=4m+2 and H_(D-1)=2m H_D. The direct lower-degree subfamily and dimension 4m+2 are correct.

The odd-prime period-image proof also passes. At zero the needed cutoff is the largest odd integer <=A; at one it is B. The corresponding integral jets over Z_p are simultaneously realizable by CRT with degree<=A+B+1<=2c-1. The normalization R(0)=0 is included in those jets. Thus the claimed exact odd period bill for A=B=2m is unchanged, with no uncharged degree enlargement or inferred two-adic formula.

The new norm inequality follows from
2x^4+15x^2-11x+2>0. Outside [1/3,2/5] the quadratic part is nonnegative; inside, its minimum is -1/60 while 2x^4>=2/81>1/60. Consequently the factor 4/121 and the candidate-wise Q upper bound are correct. Re-embedding an old numerator can increase the multiplier norm by 4^m, so this is not a free gain for an unchanged old candidate.

The saved c=5 benchmark was independently recomputed: R(1)=-4260, alpha=339, H=(0,-8,83,-286,376,-118,-59), and the reduced fraction is 355/113. It satisfies the derivative, endpoint and log-free conditions exactly. This is a known benchmark, not a new record.

## 6. Independent bounded checks

`verify_denominator_independent.py` uses only the Python standard library and imports no author verifier. It passes:

- 12 centered-profile cases, including K=2m,4m,13m for m=1,...,4;
- 252 rational pointwise remainder checks;
- 44 endpoint polynomial/beta moment and denominator cases, including K=0;
- 4 full-parity telescoping checks with R(0) nonzero, retaining its boundary term;
- the exact c=5 parity-free 355/113 benchmark.

These bounded checks audit arithmetic implementation. The all-scale claims above follow from the written identities and inequalities. No Q-to-zero theorem, new approximation record, or assertion about bad approximability is established.
