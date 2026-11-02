# Independent review of the finite coupled raw-column certificate

2026-11-02. Independently reviewed `COUPLED_RAW_CERTIFICATE.md`, its complete packet, determinant identity certificate, new-root comparison, irreducibility, and rational-value denominators. No mathematical blocker was found. This is an independent agent audit, not external peer review or an assertion of a new irrationality result for pi.

## Complete source and legal columns

The strict inequality s+4 beta_1+5 beta_2<6 gives exactly the nine stated rows per center. In particular both nonzero spatial derivative directions are retained, and no other nonnegative multiindex is admissible. The chosen eighteen raw columns have maximum weight 22/5, less than H=6, so the minimum slack 8/5 is correct. The two numerical parameter inequalities equal 1/50<1. These checks do not establish the separated-weight threshold hypotheses of the original manuscript, and the note correctly does not claim them.

For the X_1 normal grading, the row counts are 14 at index zero and four at index one. The column counts are six and twelve. Normal-index-one rows cannot use index-zero columns, so every nonzero determinant permutation uses four of the twelve index-one columns there, leaving eight for index-zero rows. The surplus is exactly eight. This is a genuinely imbalanced profile in the chosen grading; it is not a proof of irreducibility under all possible coordinate changes.

## Exact polynomial identity and irreducibility

The ordinary Taylor coefficient extraction gives the displayed rational polynomial. Independent direct Fraction arithmetic, with the note's row and column orders, agrees at Z=1,-2,1/2,7/3. Its full identity certificate is stronger than these supplementary checks: the degree bound 44 is the column spatial-degree sum 50 minus the source spatial-derivative sum six, and checking 45 distinct centers proves the polynomial equality exactly. The integer row scaling is nonzero and independent of Z, so it preserves this inference.

The least coefficient denominator is 864, because all displayed coefficients clear with it and 13/864 is reduced. The cleared polynomial has coefficient gcd one. It is not monic; retaining the leading coefficient 13 is important for rational evaluation.

The cubic f(x)=13x^3+720x^2+6480x-63936 remains degree three modulo 11 and has no root in that field. Its primitive integer coefficients therefore give irreducibility over Q. For a cubic root alpha, its field norm is 63936/13. The numerator is not divisible by 13, so this norm has 13-adic valuation -1 and cannot be a rational square. If alpha were a square in its cubic field, its norm would be a rational square. Thus adjoining sqrt(alpha) has degree two over that field and degree six over Q. This proves the sextic f(Z^2) irreducible, rather than inferring sextic irreducibility merely from irreducibility of the cubic.

The comparison with products of whole old P_n is valid. After removing powers of Z, the residual degrees of the old family are 0,0,4,6, and at least 12 for n=2,3,4,5, and n>=6 respectively. The only way to obtain a total residual degree six from such whole factors would involve the n=5 residual, which is reducible. It cannot equal this irreducible sextic. No analogous claim about all lacunary one-dimensional constructions is justified or made.

## Root location and approximation quality

The exact signs of S(it) at 6302063/10^6 and 6302064/10^6 are positive and negative. On that interval, t^2 lies between 39 and 40. The quadratic 13y^2-480y+2160 is positive there because its value at 39 and its derivative throughout the interval are positive. Thus the derivative of S(it) is negative, proving a unique simple imaginary root in the stated interval.

The consecutive-exponent comparator has its relevant positive imaginary ordinate between 6.2954 and 6.2956, as certified by the quartic signs. Both endpoints exceed 44/7, and its upper endpoint is below the new root's lower endpoint. Since 2pi<44/7, the claimed strict comparison of distances to the actual first period follows without a decimal approximation of pi. The new certified root is worse than that comparator root. The note properly separates a new algebraic zero from an improvement in approximation.

## Exact rational-value denominators

For coprime a,b with b>=2,

\[
b^6S(ia/b)=-13a^6+720a^4b^2-6480a^2b^4-63936b^6.
\]

At a prime dividing b other than 13, the first term is a unit modulo that prime and the remaining terms vanish. At 13 dividing b, the first term has valuation exactly one and every other term has valuation at least two. Therefore the gcd with b^6 is exactly gcd(b,13), giving reduced denominator b^6/gcd(b,13). Multiplication by (ia/b)^21 and division by i multiplies the real numerator by a^21 and increases the denominator power to 27. Since a is coprime to b, the same cancellation factor remains. Both nonvanishing and the exact denominator formulas follow, including powers of 13 in b.

The fourteen-row control has the complete stated packet, legal columns, and normal surplus two. Independent direct determinants agree with -2Z^13 at four rational centers. Thus positive surplus by itself does not imply a nonzero root or a useful approximation.

## Independent exact checks

`verify_coupled_independent_review.py` uses only the standard library and independently reconstructs ordinary Taylor entries. It passed eight direct raw determinant checks, 980 rational-value denominator checks, the exact root/comparator signs, the modulo-11 cubic test, and the nonsquare norm check. Results are in `coupled_independent_review_verification.json`. These finite checks supplement the main full polynomial identity certificate and the general structural arguments; no asymptotic quality claim is inferred from them.
