# Independent audit of the endpoint-content proof

2026-11-02. Reviewed [endpoint-content proof](../../denominator-reflection-route/ENDPOINT_CONTENT_AND_DENOMINATORS.md), equations (1)–(10), against the integer Hankel representation and independent exact computations. The all-dimension denominator bounds, exact even denominator, odd leading coefficient, and asymptotic decay all check out. Exact odd content is not used or proved.

## The integrality step

The critical point is stronger than rationality of the endpoint integrals: normalized Laguerre expectation sends every symmetric polynomial with integer coefficients to an integer.

For clarity, the Schur-basis argument can be made entirely algebraic. The alternant defining a Schur polynomial vanishes when two variables coincide. Each difference divides the alternant in the integer polynomial ring, because that difference is monic in one variable. The differences are pairwise nonassociated primes in that unique-factorization ring, so their product divides the alternant there. The quotient is therefore a symmetric polynomial with integer coefficients. Its lexicographically leading monomial has coefficient one. Ordered by partition, the Schur-to-monomial-symmetric change of basis is triangular with diagonal one, giving an integer basis in each homogeneous degree.

For a partition of length at most k, use increasing exponents a_i=i+lambda_(k-i), with 0<=i<k. Then

\[
\frac{\det[(a_i+j)!]_{i,j=0}^{k-1}}{A_k^2}
=\left(\prod_{i=0}^{k-1}\frac{a_i!}{i!}\right)
\frac{\prod_{i<j}(a_j-a_i)}{A_k}.
\]

Every factorial ratio is an integer since a_i>=i. The Vandermonde ratio is s_lambda(1,...,1), hence an integer by the preceding polynomial construction. This validates each separately symmetric cross coefficient after independent integration of the two groups. It avoids assuming that normalization of an arbitrary integer-valued integral preserves integrality.

## Content and the least denominator

The endpoint weights A_k^2 A_(n-k)^2 are all divisible by the balanced weight G_n=A_floor(n/2)^2 A_ceil(n/2)^2. Indeed, with k<=r=floor(n/2) and s=ceil(n/2),

\[
\frac{A_kA_{n-k}}{A_rA_s}
=\prod_{j=0}^{r-k-1}\frac{(s+j)!}{(k+j)!}\in\mathbb Z.
\]

Thus D'_n=A_n^2/G_n is a valid integer clearing factor. For n=2r the unique balanced endpoint term has coefficient (-1)^rG_n, proving exact content G_n and exact least denominator d_n=D'_n. For n=2r+1 the balanced constant terms cancel; their first-order terms add to (-1)^(r+1)4r(r+1)G_n. The next endpoint split begins four powers later, so no omitted split can change this leading coefficient. The resulting lower bound comes only from reducing that coefficient; no unproved content equality is needed.

The draft's ancillary odd divisibility is correct. An explicit product that makes it immediate is

\[
\frac{A_{2r+1}}{A_rA_{r+1}}
=r!(r+1)!\prod_{j=2}^{r}\frac{(r+j)!}{j!},\qquad r\ge2.
\]

The first factor is divisible by 2r(r+1). Its square is therefore divisible by 4r(r+1). The exceptional n=3 case is correctly handled separately.

## Asymptotics and scope

The odd least denominator differs from D'_n by at most the factor 4r(r+1); hence its logarithm differs by O(log n). Substituting the factorial-product asymptotics gives

\[
\log(d_n|P_n(2\pi i)|)
=-\tfrac12 n^2\log n+
\left[\tfrac34+\log(\pi/\sqrt2)\right]n^2+O(n\log n).
\]

The sign of the leading term is negative. This overturns an inference of universal growth from the n=4,...,10 table while remaining consistent with growth under the larger row-factorial factor. The value at 2pi i is not an integer, and no rational-approximation or bad-approximability theorem follows from the decay alone.

## Independent checks

`verify_content_structure.py` uses only the Python standard library. It checks complete integer Hankel polynomials for n=2,...,12, including coefficient gcds and both signed leading-coefficient formulas; 183 normalized Schur expectations by independent exact factorial determinants; and balanced-factorial divisibility through n=40. All checks passed. Exact odd content beyond the proved leading-coefficient divisibility remains unproved; no experimental formula is used here.
