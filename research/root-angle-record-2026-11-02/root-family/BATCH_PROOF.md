# # Exact batch treatment of the root family

The search uses only a certified rational interval `L < pi < U`, the rational
threshold `tau = 1/250`, and the cap `N = 10^100`. Continued-fraction data are
read only after the independent search and the saved tree replay finish.

## Why a whole denominator block can be rejected

For a primitive fraction p/q, the target inequality is

    q^2 |pi - p/q| < tau,

equivalently `|q*pi-p| < tau/q`. If `A <= q <= B`, this implies

    q*L - tau/A < p < q*U + tau/A.

At each q, the number of integer numerators in this strict open interval is

    ceil(q*U + tau/A) - floor(q*L - tau/A) - 1.

The interval has positive width, so this number is nonnegative. Consequently
the sum of these counts is zero if and only if the whole widened strip has no
integer lattice point. A zero sum excludes every good root at every denominator
in the block at once.

The implementation takes initial disjoint dyadic blocks
`[1,1], [2,3], [4,7], ...`, truncated at N. A positive count triggers bisection.
At a singleton q, its interval has the exact threshold `tau/q`; each possible p
is scored using the full rational pi interval. Coprimality and reduced scores
are also recorded. The target coverage concerns primitive labels, so different
unreduced representations never count as new rational approximations.

All possible strip numerators are inside `3q < p < 4q` in this run: the stronger
input inequalities `L-tau > 3` and `U+tau < 4` hold. Thus the search does not
leave the user's root family.

## Why floor_sum is a genuinely batched computation

For integers n>=0, m>0, a,b, define

    S(n,m,a,b) = sum_{i=0}^{n-1} floor((a*i+b)/m).

Integer quotient/remainder reductions of a and b contribute explicit polynomial
terms and leave `0<=a,b<m`. Put `n' = floor((a*n+b)/m)` and
`b' = (a*n+b) mod m`. If n'=0, the remaining sum is zero. Otherwise, lattice-point
counting gives

    S(n,m,a,b) = S(n',a,m,b').

Indeed, for each j=1,...,n', the count of i for which `j*m <= a*i+b` is
`n-ceil((j*m-b)/a)`. With k=n'-j this is `floor((k*m+b')/a)`.
This proves the swap step used in the Euclidean loop. After quotient reduction
the positive modulus decreases, so the computation terminates without visiting
each i. Negative a,b are handled by Python's integer divmod before this loop.

Rational affine sums are reduced to this integer form using a common denominator.
The block count uses two such sums because `ceil(t) = -floor(-t)`.

## Actual finite result

The Machin enclosure uses 200 terms for atan(1/5) and 70 terms for atan(1/239),
with alternating next-term bounds. Its width is <10^-281. The condition

    N*(U-L) + 2*tau < 1

is checked; it ensures each leaf uncertainty strip has at most one numerator.

Exact search statistics:

- 333 initial dyadic blocks;
- 345 total tree nodes, all saved and replayed;
- 338 zero-count rejected blocks;
- 6 bisected blocks;
- 690 floor_sum calls and 64294 Euclidean loop iterations during search;
- a single leaf: p=355, q=113;
- its full quality interval is strictly below 1/250;
- no discarded leaf candidates and no unresolved candidates.

Therefore among all primitive root fractions with `q<=10^100`, the *only* one
with `q^2 |pi-p/q| < 1/250` is 355/113. Since its quality is itself below this
threshold, it is also the unique global minimum of that quality over this finite
range. This argument does not require Legendre's theorem or a continued-fraction
search. The separately constructed 194-convergent certificate agrees afterwards.

The floor_sum implementation was first checked against direct summation for
31768 small integer cases, including negative coefficients; 448 small rational
strip cases; all 14 actual-pi dyadic strip counts through q=10000; and a direct
per-denominator candidate scan through q=10000. No comparison decision uses a
floating-point pi value.

This is a finite algorithmic result. The small node count in this particular run
is not asserted as a uniform complexity theorem. The enclosure precision and
large-integer arithmetic costs are part of the method. No infinite quality-zero
subsequence or new pi approximation record has been proved.
