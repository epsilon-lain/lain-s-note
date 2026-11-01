# Independent check of the boundary-rank filtration

2026-11-01. Reviewed `BOUNDARY_RANK_FILTRATION.md` independently of its finite checker. No mathematical error found in equations (1)--(7). This review concerns the stated formal transport claims, not a quantitative bound for the period determinant.

## Cap threshold

Write c=v_1-v_0>0 and A=floor(H/(theta v_1)), with theta<=1/2. If N>A, integrality gives

\[
N>2H/v_1.
\]

Indeed A>=floor(2H/v_1), so an integer strictly above A is strictly above 2H/v_1. Hence the source length satisfies

\[
b_N\le\left\lceil\frac{H-v_0N}{c}\right\rceil
 <\frac{H(v_1-2v_0)}{v_1c}+1
 =\frac H{v_1}-\frac{Hv_0}{v_1c}+1.
\]

At the stated sufficient threshold H>=v_1c/v_0 this is at most H/v_1. Moreover A+1>2H/v_1, so 2b_N<=A+1. Thus the capped outside interval has A-b_N+1>=b_N columns, and the capped and uncapped ranks both equal b_N. When N<=A, the cap is inactive directly. If the source interval is empty, no block is counted. The threshold is sufficient; sharpness is not claimed.

## Rank of each boundary block

For rows beta in {0,...,b-1} and outside columns a>=b, its scalar coefficient apart from the common nonzero error powers is

\[
(-1)^{b-\beta-1}\binom a\beta\binom{a-\beta-1}{b-\beta-1}
=\frac{(-1)^{b-\beta-1}}{\beta!(b-\beta-1)!}
\frac{a!}{(a-b)!}\frac1{a-\beta}.
\]

The first two factors are nonzero row and column scalings. The remaining matrix is Cauchy, with distinct row nodes beta, distinct column nodes a, and a-beta>0. Every square submatrix is nonsingular. Multiplication by a nonzero error to the power a-beta is another row and column scaling. This proves the stated minimum-of-dimensions rank, including the polynomial type cap.

## Integral and normalization

With x=N/H, the limiting source count in a block is

\[
b(x)=\min\left(x,\frac{1-v_0x}{c}\right).
\]

The boundary-rank density is zero before x=1/v_1. Afterwards it is the minimum of (v_1x-1)/c and (1-v_0x)/c. These lines meet at x=2/(v_1+v_0), with height 1/(v_1+v_0). The region is a triangle with base

\[
\frac1{v_0}-\frac1{v_1}=\frac{c}{v_0v_1}
\]

and area c/(2v_0v_1(v_1+v_0)). Dividing by the source-area constant 1/(2v_0v_1) gives c/(v_1+v_0), as claimed. The order of limits is correctly fixed weights first, and then a separate large weight-ratio limit.

## Formal degree filtration

The support implication c dot gamma>=g and c dot gamma<=c_max |gamma| gives the row order ceil(g/c_max). Under normalized source counting, slack x=g/H has density (m+1)(1-x)^m. Its lowest rho mass ends at a=1-(1-rho)^(1/(m+1)); direct integration gives

\[
\int_0^a x(m+1)(1-x)^m\,dx
=\rho-\frac{m+1}{m+2}
\left(1-(1-\rho)^{(m+2)/(m+1)}\right).
\]

The ceilings contribute O(M), which disappears after division by MH. This establishes the guaranteed row-order sum; it does not identify the exact determinant valuation. The formal determinant reduction is the determinant lemma after deleting terms of degree greater than D, with the nonsingular base matrix fixed in the error variables. The inverse matrix remains in the reduced expression, so no conditioning estimate is obtained.

## Exact checks rerun

The companion `verify_rank_filtration.py` passed: 6 truncated determinant reductions, 348 boundary-support checks, 864 exact rational boundary-block ranks, and 15 finite complete-packet counts. Three small-height cases change under the cap, so retaining the stated large-H qualification is necessary for the simpler formula. These checks complement the proof and do not establish any asymptotic claim by extrapolation.

The positive-fraction transport verifier also passed independently: 175 exhaustive target-set cost checks, 12 nonzero fiber determinant checks, and 12 capped complete-packet constructions. The only clarification added to that note during mutual review was that the error values are fixed before the H limit in its transport-norm inequality.
