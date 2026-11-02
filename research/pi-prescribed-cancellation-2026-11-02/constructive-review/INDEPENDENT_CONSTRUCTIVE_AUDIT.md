# Independent audit of prescribed cancellation and its arithmetic bill

2026-11-02. This is an internal independent mathematical and exact-arithmetic audit, not external peer review. The fixed family is a=m,b=2m,c=4m+1, deg H<=m+1. The reviewed work remains constructive and does not resume the old prime-pool exclusion problem.

## 1. Individual polynomial moments have an exact reparameterization cost

Let w(t)=t^(a-1/2)(1-t)^b. A degree-k Rodrigues multiplier can be chosen as

\[
H(t)=\frac{4^k}{k!w(t)}\frac{d^k}{dt^k}
\{t^{a+k-1/2}(1-t)^{b+k}\}.
\]

It annihilates the first k polynomial moments. The generalized half-binomial denominator rule makes this H integral. Integrating by parts k times gives the exact identity

\[
\int_0^1\frac{x^{2a}(1-x^2)^bH(x^2)}{(1+x^2)^c}\,dx
=4^k\binom{c+k-1}{k}
\int_0^1\frac{x^{2(a+k)}(1-x^2)^{b+k}}{(1+x^2)^{c+k}}\,dx.
\]

Consequently the rational and pi periods scale by the same nonzero integer. Whenever the pi period is nonzero, the fully reduced fraction is exactly unchanged by this representation. At maximal k=m+1 the effective constant-multiplier parameters are (2m+1,3m+1,5m+2). This is not a new degree of arithmetic freedom; it does not show that every signed multiplier or all cancellation conditions fail.

The accompanying `audit_rodrigues.py` independently verifies exact rational period pairs and reduced fractions in 52 cases through m=8, including the maximal case. No pi approximation is used in those identity checks.

## 2. Cancelling the complete rational truncation is a different explicit rule

The complete-truncation construction takes K=13m and defines mu_j by integrating the parent's rational kernel truncation against x^(2m+2j)(1-x^2)^(2m), j=0,1. It sets H=mu_1-mu_0 t. This imposes one rational equation involving the whole truncation; it does not impose K separate zero moments. Thus the degree-one multiplier is not the maximal Rodrigues specialization above.

The determinant proof in `COMPLETE_TRUNCATION_GENERATOR.md` is valid: a nonzero linear multiplier with both periods zero would force R=x^(2a+1)(1-x^2)^(b+1)S(x^2), and its differential numerator has multiplier degree at least two. The contradiction proves independence of the two base period vectors without presuming pi irrationality or numerical cancellation.

The closed base periods in the arithmetic note were also independently checked:

\[
\alpha_0=\binom{2m}{m}/64^m,\quad\beta_0=0,
\quad\alpha_1=\frac{6m+1}{2m-1}\alpha_0,
\quad\beta_1=-\frac1{(2m-1)\binom{4m}{2m}}.
\]

The beta substitutions give these identities directly. The K=13m inequality certifies the pi coefficient is nonzero for every m; it is not just a finite test. The positivity argument for S_K, mu_0 and the displayed fraction denominator is valid too.

## 3. Costs and the exact status of convergence

The universal moment denominator divisor, the Fourier common-period denominator divisor, and the exact gcd formula for q are all valid. Integer realization by the adjugate of M_c is legitimate. Its common scaling cancels from the period ratio exactly; no extra denominator advantage is inferred from that scaling.

Using the closed periods gives q_m <= (3/2)binom(2m,m)L_(13m). The moment denominator L_(13m) is retained. Both the actual-q signed-remainder certificate and the simplified joint bill with coefficient 12 are correct.

The new proof establishes

\[
|\pi-p_m/q_m|\le\frac{128}{3}(4m^2-1)
\left(\frac{2^{33}}{3^{23}}\right)^m\to0.
\]

Since pi is irrational, this implies q_m tends to infinity. It does not imply monotonic denominator growth, and it does not establish q_m^2|pi-p_m/q_m| tending to zero. The universal denominator envelope yields a poor quality upper bill; this must not be described as a theorem that the actual denominators or actual quality fail asymptotically.

## 4. Independent bounded audit

`audit_complete_truncation.py` imports only exact rational utilities from the companion `audit_rodrigues.py`. For m=1,2,3 it independently constructs S_K, integrates its polynomial moments, verifies the whole-truncation cancellation, checks the closed periods, tests the proven nonzero inequality, recovers the rational primitive coefficient by coefficient, and verifies that T_cL_K clears every primitive coefficient to an integer.

It also checks the actual reduced fraction against the closed binomial formula, verifies the moment-denominator divisor, the simplified q bound, and the signed-integral upper bounds.

| m | K | reduced-q bit length | moment-denominator bit length | actual-q remainder certificate below one? |
|---:|---:|---:|---:|---|
| 1 | 13 | 56 | 71 | No |
| 2 | 26 | 122 | 141 | No |
| 3 | 39 | 177 | 210 | No |

These last entries describe the available certificates, not the actual quality. No assertion of a new pi record, a quality-zero sequence, or a proof of bad/nonbad approximability is made.
