# Independent audit of the all-size root-distance argument

2026-11-02. Reviewed `../near-root-general/GENERAL_NEAR_PERIOD_ROOTS.md`, equations (1)--(10). No mathematical blocker found. This audit is independent of the note's finite root calculations.

The degree bound is valid: in each reflected determinant factor, two columns of the same parity selecting the same power j select exactly the same coefficient vector. A nonzero determinant term therefore uses distinct nonnegative j values inside each parity group. Their minimum sum is the stated sum of two binomial coefficients. Possible further dependencies can only lower the degree. The resulting bound agrees with the exact leading coefficients proved independently by endpoint splitting.

The low coefficient normalization and sign are correct. The most negative Laurent coefficient of moment j>=1 is (-1)^(j-1)j!. Row and column sign extraction leaves an overall factor (-1)^n. The factorial Hankel determinant is A_n^2 and its (0,0) cofactor is n A_n^2. Zeroing that entry thus gives (1-n)A_n^2, proving [Z^n]P_n=(-1)^(n+1)(n-1). No lower-degree term can occur because every selected Laurent exponent is bounded below by minus its moment index.

The root-product lemma uses the degree inequality in the correct direction: for 0<=x<1 and d<=D, x^d>=x^D. Its geometric inequality |rho|<=|z|+|z-rho| is sufficient and uses no location assumptions on roots. Normalizing R_n(0)=1 excludes zero roots in its product. If the actual period were itself a root, the distance conclusion would still hold with distance zero; transcendence is used only in the optional positivity observation.

The absolute Hankel upper bound applies at every integer exponential period, even where the special cosine lower bound at the first period is unavailable. Writing t=2*pi*|k| gives

\[
\log v_{n,k}
=-n^2\log n+\left[\frac32+\log(t/4)\right]n^2+O_k(n\log n).
\]

Division by D_n=n^2/2-n+O(1) therefore yields

\[
r_{n,k}=\frac{e^3t^2}{16n^2}(1+o_k(1)),
\qquad
\frac{t r_{n,k}}{1-r_{n,k}}
=\frac{e^3t^3}{16n^2}(1+o_k(1)).
\]

This verifies the constant in the stated limsup bound. Replacing t by 44|k|/7 and the D_n-th root by a larger rational q preserves all inequalities. The conclusion concerns algebraic roots of increasing allowed degree, not rational approximations, and the note correctly retains that distinction.

## Audit of the uniform corollary (11)

The added uniform corollary also checks. Before division by D_n, the exact expression is

\[
\log v_{n,k}
=(n^2-n)\log t+2\log A_n-\log A_{2n}-\log(n-1).
\]

For \(2\pi\le t\le n\), the extra term \(-n\log t\) is bounded by O(n log n) with a constant independent of k. The factorial remainder is itself independent of k. Dividing by \(D_n=n^2/2-n+O(1)\) consequently gives the asserted uniform error O(log n/n) in log r. For all sufficiently large n this error has exponential at most 2, so

\[
r_{n,k}\le\frac{e^3t^2}{8n^2}.
\]

The chosen range \(1\le|k|\le n/(8\pi)\) means \(t\le n/4\), giving \(r_{n,k}\le e^3/128<1/2\). Thus the root-product bound is at most

\[
2t r_{n,k}\le\frac{e^3t^3}{4n^2},
\]

exactly as in (11). If \(K_n=o(n^{2/3})\), then this range includes every \(1\le|k|\le K_n\) for all sufficiently large n, and the maximum radius is \(O(K_n^3/n^2)=o(1)\). Distinct period centers have separation at least \(2\pi\), so these disks are eventually disjoint and supply at least \(2K_n\) distinct nonzero roots. This is a uniform extension of the existing argument, not an interchange of a fixed-k asymptotic with an unbounded k-limit.
