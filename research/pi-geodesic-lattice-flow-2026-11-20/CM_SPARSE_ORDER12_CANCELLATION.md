# CM-sparse order-12 cancellation: a three-point elliptic construction

Date: 2026-11-20  
Status: **proved formal expansion and minimum-multiplier construction; exact rational experiment at one recurrence index**. This is a curve-specific improvement over the generic binomial finite-difference filter.

## 1. The order-six automorphism removes five-sixths of the formal coefficients

Take
\[
E:\ y^2=x^3-2,\qquad P=(3,5),\qquad z=\xi=-x/y.
\]
Let \(\zeta\) be a primitive sixth root of unity. The map
\[
\tau(x,y)=(\zeta^2x,\zeta^3y)
\]
is an automorphism of \(E\), and
\[
\tau^*z=\zeta^{-1}z.
\]
Since multiplication by \(j\) commutes with \(\tau\), its formal series
\[
[j]_F(z)=\sum_{r\ge1}c_r(j)z^r
\]
satisfies
\[
[j]_F(\zeta^{-1}z)=\zeta^{-1}[j]_F(z).
\]
Thus
\[
\boxed{c_r(j)=0\quad\text{unless }r\equiv1\pmod6.}
\tag{1}
\]
Only the powers \(z,z^7,z^{13},z^{19},\ldots\) occur.

This is the key point: for this CM curve, cancelling every polynomial moment up to degree \(m\) is wasteful. To reach order \(z^{12}\), it is enough to cancel the linear and degree-seven coefficients.

## 2. Compute the first two non-linear coefficients exactly

Set
\[
z=-x/y,\qquad w=-1/y.
\]
The curve equation becomes
\[
w=z^3-2w^3.
\]
Solving recursively,
\[
w=z^3-2z^9+12z^{15}-96z^{21}+O(z^{27}).
\]
The normalized invariant differential is
\[
\omega=\frac{dx}{2y}
=\frac{zw'(z)-w(z)}{2w(z)}\,dz
=\left(1-6z^6+60z^{12}-672z^{18}+O(z^{24})\right)dz.
\]
Therefore the formal logarithm and its inverse are
\[
\log_F(z)=z-\frac67z^7+\frac{60}{13}z^{13}
-\frac{672}{19}z^{19}+O(z^{25}),
\]
\[
\exp_F(T)=T+\frac67T^7+\frac{48}{91}T^{13}+O(T^{19}).
\]
Since \([j]_F(z)=\exp_F(j\log_F(z))\), direct composition gives
\[
\boxed{
[j]_F(z)
=jz+\frac67(j^7-j)z^7
+\frac{12}{91}(4j^{13}-39j^7+35j)z^{13}
+O(z^{19}).
}
\tag{2}
\]

## 3. A three-point integer filter

Choose the integer weights
\[
\boxed{
(\lambda_2,\lambda_3,\lambda_4)=(-364,291,-36).
}
\tag{3}
\]
They satisfy the two sparse moment identities
\[
2\lambda_2+3\lambda_3+4\lambda_4=1,
\]
\[
(2^7-2)\lambda_2+(3^7-3)\lambda_3+(4^7-4)\lambda_4=0.
\]
Equivalently,
\[
\sum_{j=2}^4\lambda_jj=1,
\qquad
\sum_{j=2}^4\lambda_jj^7=1.
\tag{4}
\]
Define
\[
\mathcal C(V)=
-364\frac{\xi([2]V)}{\xi(V)}
+291\frac{\xi([3]V)}{\xi(V)}
-36\frac{\xi([4]V)}{\xi(V)}.
\]
Substituting (2), the coefficient of \(z^6\) vanishes by (4). The next coefficient is
\[
\frac{12}{91}\sum_{j=2}^4\lambda_j
(4j^{13}-39j^7+35j)
=-1\,031\,184\,000.
\]
Hence
\[
\boxed{
\mathcal C(V)
=1-1\,031\,184\,000\,\xi(V)^{12}
+O(\xi(V)^{18}).
}
\tag{5}
\]
The error is order twelve, although only three multiplication maps are used.

## 4. Minimum possible largest multiplier for this sparse cancellation

The two required constraints are
\[
\sum_j\lambda_jj=1,\qquad
\sum_j\lambda_j(j^7-j)=0.
\tag{6}
\]
No solution to (6) can use only two nonzero weights.

Indeed, suppose the support were \(\{a,b\}\), with \(2\le a<b\). Solving the two equations gives
\[
\lambda_b
=-\frac{a^6-1}{b(b^6-a^6)}.
\]
But
\[
b(b^6-a^6)
\ge(a+1)((a+1)^6-a^6)
>6a^6>a^6-1,
\]
so \(0<|\lambda_b|<1\), contradicting the fact that \(\lambda_b\) is a nonzero integer.

Thus at least three distinct multipliers are needed. Since all multipliers are at least 2, every such scheme has largest multiplier \(J\ge4\). The weights (3) attain \(J=4\). Therefore:

\[
\boxed{\text{The construction (3) has the minimum possible largest multiplier }J=4.}
\tag{7}
\]

This is an optimality statement for the two sparse moment constraints (6), not for every possible rational function on the elliptic curve.

## 5. Consequence for the EDS denominator barrier

The previous primitive-EDS argument applies to every fixed finite integer combination. Here the largest multiplier is \(J=4\), and its coefficient is \(\lambda_4=-36\). For all sufficiently large \(s\), the primitive part of \(B_{8s}\) survives in the reduced denominator of \(\mathcal C([s]P)\).

Streng's Proposition 1.3 gives
\[
\log\|D_N^{\mathbb Z}\|
=\widehat h(P)\left(\prod_{p\mid N}(1-p^{-2})\right)N^2
+O_\epsilon(N^\epsilon).
\]
Consequently,
\[
\boxed{
\liminf_{s\to\infty}
\frac{h(\mathcal C([s]P))}{s^2}
\ge\frac{64\,\widehat h(P)}{\zeta(2)}.
}
\tag{8}
\]
For comparison, the binomial order-seven filter has largest multiplier \(J=9\), giving the corresponding primitive-denominator coefficient \(324\,\widehat h(P)/\zeta(2)\). The sparse construction lowers this unavoidable denominator floor to \(16/81\) of that value.

This is a lower bound on the denominator height, not an upper bound and not a proof about \(\pi\).

## 6. Exact rational experiment at \(V=[41]P\)

All elliptic group operations and all values of \(\mathcal C(V)\) below were computed using exact rational arithmetic on \(E:y^2=x^3-2\). Only the displayed logarithms are decimal summaries.

For \(s=41\),
\[
\log_{10}|\xi([41]P)|=-4.694041041408119\ldots.
\]

| Filter | Largest multiplier \(J\) | Leading formal order | \(\log_{10}|\mathcal C([41]P)-1|\) | Logarithmic rational height |
|---|---:|---:|---:|---:|
| Sparse three-point (3) | 4 | 12 | -47.315156 | 95,282.82 |
| Binomial order 4 | 6 | 6 | -23.727130 | 285,846.86 |
| Binomial order 6 | 8 | 6 | -24.528763 | 666,978.16 |
| Binomial order 7 | 9 | 12 | -45.290221 | 911,991.13 |

The exact rational numerator and denominator of the sparse output each have 137,465 bits. For the binomial order-seven output they each have 1,315,726 bits.

The leading term in (5) predicts
\[
\log_{10}|\mathcal C([41]P)-1|
\approx -47.315156331073\ldots,
\]
matching the exact rational computation to the displayed precision. This is a strong finite consistency check of the formal expansion.

Compared with the binomial order-seven filter, the sparse filter has:
- about \(10.45\%\) of its logarithmic rational height,
- about \(0.944\%\) of its absolute error,
- a smaller maximum multiplier (4 instead of 9).

Compared with the binomial order-four filter, the sparse filter has about one third of the logarithmic height and an error smaller by a factor of about \(3.87\times10^{23}\).

The follow-up counts all surviving poles. The exact coefficient of \(\widehat h(P)s^2\) in the log denominator is 84 for the sparse filter and 804 for the binomial order-seven filter. The coefficient 84 is optimal among fixed finite integer filters satisfying (6). [Exact coefficients and optimality proof](EDS_EXACT_POLE_HEIGHT.zh.md).

## 7. What this changes in the research plan

The generic binomial interpolation filter was designed to cancel every polynomial coefficient up to a chosen degree. That ignores the order-six automorphism of this curve.

The improved strategy is:

1. determine the exponents permitted by the CM automorphism;
2. compute only the corresponding coefficient polynomials \(c_{6k+1}(j)\);
3. solve the resulting sparse integer moment problem;
4. minimize the largest multiplier first, then optimize the actual rational height;
5. test the resulting filter with exact elliptic arithmetic.

The three-point filter (3) is the first concrete success of this strategy. It improves both the error and the actual rational height in the finite experiment, while attaining the minimum possible largest multiplier for the required sparse constraints.

It still does **not** prove that \(\pi\) is badly approximable. The remaining bridge from this elliptic construction to a uniform lower bound for \(q\|q\pi\|\) is absent.

## References

- M. Streng, *Divisibility sequences for elliptic curves with complex multiplication*, Algebra & Number Theory 2 (2008), Proposition 1.3, https://msp.org/ant/2008/2-2/ant-v2-n2-p04-p.pdf.
- The EDS primitive-denominator argument for finite integer combinations is recorded in EDS_PRIMITIVE_DENOMINATOR_BARRIER.md.
