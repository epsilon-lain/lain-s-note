# Exact elliptic test of height-aware weights: the quadratic proxy fails

Date: 2026-11-20  
Status: **exact rational comparison on one elliptic curve** plus a general minimal-support theorem. The tested alternative weights satisfy every cancellation moment, but they lose to the binomial weights in both actual rational height and real error.

## 1. The candidates

On
\[
E:y^2=x^3-2,\qquad P=(3,5),\qquad \xi=-x/y,
\]
we previously found alternative integer weights for orders \(m=4,5,6\) that reduce
\[
C_2(\lambda)=\sum_j|\lambda_j|j^2
\]
by 12.1%, 16.8%, and 32.3%, respectively.

The alternative weights are:

\[
\begin{aligned}
m=4:\quad &(4,-5,0,5,-4,1)
&&\text{at }j=2,\ldots,7,\\
m=5:\quad &(4,-5,1,1,1,1,-5,4,-1)
&&\text{at }j=2,\ldots,10,\\
m=6:\quad &(4,-4,-3,6,0,0,-6,3,4,-4,1)
&&\text{at }j=2,\ldots,12.
\end{aligned}
\]
Each satisfies the exact moments
\[
\sum_j\lambda_jj^r=1,\qquad 0\le r\le m.
\]
Thus the formal-group cancellation order is unchanged.

## 2. Exact rational experiment

For \(s\in\{15,26,41\}\), compute every elliptic multiple \([js]P\) by the exact rational group law, and form
\[
B_{m,\lambda}([s]P)
=\sum_j\lambda_j\frac{\xi([js]P)}{\xi([s]P)}.
\]
All group-law operations and all resulting rational values are exact. Only the displayed logarithms of heights and errors are decimal summaries.

For the most interesting recurrence \(s=41\):

| Order | Weights | \(\log_{10}|B-1|\) | Log height \(h(B)\) | Numerator/denominator bits |
|---:|---|---:|---:|---:|
| 4 | binomial | -23.7271 | 285,847 | 412,390 |
| 4 | alternative | -23.3069 | 394,742 | 569,492 |
| 5 | binomial | -23.9267 | 449,189 | 648,043 |
| 5 | alternative | -23.0102 | 1,238,677 | 1,787,033 |
| 6 | binomial | -24.5288 | 666,978 | 962,247 |
| 6 | alternative | -23.4496 | 1,810,374 | 2,611,818 |

The same direction occurs for \(s=15\) and \(s=26\): in all nine \((s,m)\) comparisons, the alternative has both **larger actual logarithmic height** and **larger error**.

At \(s=41\), the alternative errors are approximately 2.63, 8.25, and 12 times larger for orders 4, 5, and 6. The alternative logarithmic heights are about 1.38, 2.76, and 2.71 times the binomial heights.

So the lower \(C_2\) proxy did not predict the actual rational output. The alternatives use larger multipliers, and the denominators of \(\xi([j]V)/\xi(V)\) grow rapidly with \(j\). The common-denominator cost dominates the coefficient-norm saving in this experiment.

## 3. A structural theorem: binomial weights uniquely minimize the largest multiplier

There is a general obstruction to improving the binomial scheme while keeping all multipliers small.

**Theorem.** Suppose finitely supported integer weights \(\lambda_j\), supported on integers \(j\ge2\), satisfy
\[
\sum_j\lambda_j f(j)=f(1)
\]
for every polynomial \(f\) of degree at most \(m\). Then:

1. at least \(m+1\) distinct indices \(j\) have nonzero weight;
2. if every index satisfies \(j\le m+2\), then the support must be exactly \(\{2,3,\ldots,m+2\}\), and the weights are uniquely
   \[
   \boxed{\lambda_j=(-1)^j\binom{m+1}{j-1}.}
   \]

**Proof.** If the support \(S\) had at most \(m\) elements, take
\[
f(x)=\prod_{j\in S}(x-j).
\]
Then \(\deg f\le m\), the weighted sum is zero, but
\[
f(1)=\prod_{j\in S}(1-j)\ne0,
\]
a contradiction. Thus there are at least \(m+1\) support indices. There are only \(m+1\) integers in \(\{2,\ldots,m+2\}\), so a solution supported there must use all of them. The evaluation functionals at these \(m+1\) distinct nodes form an invertible Vandermonde matrix, hence the weights are unique. Lagrange interpolation at \(1\) gives the displayed binomial coefficients. \(\square\)

This says the binomial construction is already optimal for minimizing the maximum multiplication index among all integer (indeed rational) moment solutions. Any genuinely different solution must introduce at least one multiplier \(j\ge m+3\).

## 4. Revised research direction

The experiment rejects the simple strategy “minimize \(\sum|\lambda_j|j^2\) and expect the rational height to improve.” That proxy ignores the denominator structure of the elliptic multiples.

The next useful approaches are therefore:

1. **Denominator-aware optimization:** optimize the actual reduced rational height \(h(B_{m,\lambda}(V))\), not only a canonical-height proxy.
2. **Divisibility-aware supports:** exploit elliptic divisibility sequences and common factors among denominators of \(\xi([j]V)\). A larger multiplier may be acceptable if its denominator is already absorbed by previous terms.
3. **Different algebraic constructions:** seek a single rational function on \(E\) with high-order contact at the identity, rather than a linear combination of many multiplication maps. Any such function must still preserve the required \(p\)-adic congruences.

The present negative experiment is useful: it prevents us from treating a lower coefficient/height proxy as a real improvement.

## 5. Scope

This is a finite experiment on one curve, one generator, three recurrence indices, and three cancellation orders. It does not prove that binomial weights minimize actual rational height among all possible constructions. It does not produce a new approximation to \(\pi\), and it does not prove or disprove that \(\pi\) is badly approximable.
