# Independent review of the complete normal-block extension

2026-11-02. Independently reviewed `COMPLETE_NORMAL_BLOCK_PRODUCT.md`, equations (1)--(14). The complete source, legal raw columns, exact block factorization, least coefficient denominator, rational specialization cost, and the two different limits all check. No mathematical blocker was found. This review does not certify the original manuscript's critical high-dimensional interpolation theorem.

## Completeness and admissibility

For N>=2 and d>=1, the displayed parameters give v_i=w_i/theta in every spatial direction and v_0=2>w_0=1. They also give 2theta=1/N<1 and 2(w_0/v_0)theta=1/(2N)<1. The complete strict source inequality is exactly

\[
2s+4|\nu|+4N\beta_m<4N.
\]

All indices are nonnegative integers, so beta_m=0 and s<2(N-|nu|); there are no omitted source rows. The chosen columns satisfy

\[
h+2a+(2/N)|\kappa|\le4N-1-(4-2/N)|\kappa|<4N.
\]

Thus they are genuine legal raw monomials. Each normal multiindex of degree t contributes 4(N-t) rows and columns. The count follows from the generating function

\[
\sum_{N\ge1}\sum_{t=0}^{N-1}(N-t)
\binom{t+d-1}{d-1}x^N
=\frac{x}{(1-x)^{d+2}},
\]

which gives M=4 binom(N+d,d+1). The shallow P_2 blocks cause no parameter exception: they are diagonal blocks inside the full legal N>=2 packet, rather than standalone packets using theta=1/2.

## The raw factorization has no hidden mixed-column assumption

Normal coefficient extraction in a raw column is zero unless kappa>=nu coordinatewise. If equality holds in total degree, it must hold in every coordinate. Therefore ordering the normal blocks by total degree really gives a block upper triangular matrix, with square diagonal block P_{2(N-|nu|)}. The explicit within-block row and column ordering fixes the sign exactly as in (6).

The optional U_i=X_i-T substitution is also valid. Every term in T^a(X-T)^kappa has T exponent a+|kappa-lambda| and normal exponent lambda<=kappa; its combined index is a+|kappa|+|lambda|<=a+2|kappa|<2N. Its transformation has diagonal entries one and preserves the selected finite column space. This confirms the block diagonal interpretation, but the raw triangular proof already suffices.

## Arithmetic and degree

Each even block becomes a primitive integer polynomial with unit leading coefficient after multiplication by its exact denominator. Consequently their product is an integer polynomial with unit leading coefficient. The leading coefficient of the original rational polynomial is the signed reciprocal of the product of the denominators. Any integer clearing factor must therefore be divisible by that product. This proves exact minimality, not just an upper bound.

The degree sum is checked by

\[
\sum_{N\ge1}\sum_{k=1}^{N}k^2
\binom{N-k+d-1}{d-1}x^N
=\frac{x(1+x)}{(1-x)^{d+3}}.
\]

Multiplication by two gives precisely

\[
D=4\binom{N+d}{d+2}+2\binom{N+d}{d+1}.
\]

All powers of the polynomial are even. At ia/b with gcd(a,b)=1, multiplying its value by b^D produces an integer. Modulo each prime dividing b only the unit leading term remains, and it is nonzero. Thus the reduced value denominator is exactly b^D. No extra cancellation appears merely by adding normal directions.

## Both asymptotic regimes and the root set

For fixed d and N tending to infinity, the leading logarithm is minus D log N. Replacing log k by log N in the sum of -2k^2 log k costs only O_d(N^{d+2}); this follows from the integrability of x^2 log(1/x) on [0,1]. The lower-order block error sums to O_d(N^{d+1} log N). Hence the coefficient -4/(d+2)! in (11) is correct.

For fixed N>=2 and d tending to infinity, the k=1 multiplicity has degree N-1 in d, while every k>=2 multiplicity has degree at most N-2. The first block has exact value |R_2(2pi i)|=4pi^2>1. Its contribution therefore gives

\[
\log|\mathcal Q_{N,d}(2\pi i)|
=\frac{\log(4\pi^2)}{(N-1)!}d^{N-1}+O_N(d^{N-2}),
\]

with the positive leading sign stated in the note. The degree has leading term 2d^{N-1}/(N-1)!, so the normalized limit is log(2pi). The two asymptotic regimes cannot be interchanged; their error constants have different fixed parameters.

Every block multiplicity is positive for d>=1. The zero set is consequently exactly the union of the finitely many one-dimensional block zero sets and is independent of d at fixed N. Only multiplicities change. This rules out a new closer root from dimension growth within this particular extension, and does not rule out genuinely coupled high-dimensional constructions.

## Independent exact checks

The companion `verify_normal_block_review.py` uses only the standard library. It independently enumerates the packet and columns, checks the displayed count and degree identities, and compares complete raw matrix determinants with independent formal Hankel block determinants modulo two primes at several centers. Its finite checks supplement the all-parameter proofs above; the two infinite limits are not inferred from those computations.
