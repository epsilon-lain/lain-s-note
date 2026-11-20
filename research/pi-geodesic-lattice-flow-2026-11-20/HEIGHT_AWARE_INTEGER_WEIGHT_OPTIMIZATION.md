# Height-aware integer weights for high-order formal-group cancellation

Date: 2026-11-20  
Status: **exact integer identities and improved finite candidates**. No global optimality claim is made.

## 1. Why the objective must include the multiplication index

For a rational elliptic point \(V\), the canonical height satisfies
\[
\widehat h([j]V)=j^2\widehat h(V).
\]
A crude height bound for a weighted rational combination therefore naturally involves
\[
\mathcal C_2(\lambda)=\sum_j|\lambda_j|j^2,
\]
not merely the coefficient norm
\[
\mathcal C_0(\lambda)=\sum_j|\lambda_j|.
\]
Minimizing \(\mathcal C_0\) alone can push support to large \(j\) and make the height proxy worse. The more relevant first optimization target is \(\mathcal C_2\).

## 2. Polynomial encoding of all moment constraints

For integer weights \(\lambda_j\), put
\[
L(x)=\sum_{j\ge2}\lambda_jx^{j-2}.
\]
The identity
\[
\sum_{j\ge2}\lambda_j f(j)=f(1)
\quad(\deg f\le m)
\tag{1}
\]
is equivalent to
\[
\frac{L^{(r)}(1)}{r!}=(-1)^r
\quad(0\le r\le m),
\]
and hence to
\[
\boxed{xL(x)-1=(x-1)^{m+1}Q(x)}
\tag{2}
\]
for some \(Q(x)\in\mathbb Z[x]\) with \(Q(0)=(-1)^m\).

Conversely, any integer polynomial \(Q\) with this constant term makes
\[
L(x)=\frac{1+(x-1)^{m+1}Q(x)}x
\]
an integer polynomial satisfying all moment identities (1). Its coefficients are the required integer weights.

This turns weight design into a polynomial construction problem.

## 3. Explicit candidates with lower quadratic height proxy

The original binomial weights are
\[
\lambda^{\rm bin}_{m,j}=(-1)^j\binom{m+1}{j-1},
\qquad 2\le j\le m+2.
\]
The following alternative choices of \(Q\) give exact integer weights:

### Order \(m=4\)

Take
\[
Q_4(x)=x+1.
\]
Then
\[
L_4(x)=x^5-4x^4+5x^3-5x+4,
\]
so the nonzero weights are
\[
\boxed{
\lambda_2=4,\quad
\lambda_3=-5,\quad
\lambda_5=5,\quad
\lambda_6=-4,\quad
\lambda_7=1.
}
\]
All other \(\lambda_j\) vanish.

### Order \(m=5\)

Take
\[
Q_5(x)=-(x+1)(x^2+x+1).
\]
Then
\[
L_5(x)=-x^8+4x^7-5x^6+x^5+x^4+x^3+x^2-5x+4,
\]
so
\[
\boxed{
(\lambda_2,\ldots,\lambda_{10})
=(4,-5,1,1,1,1,-5,4,-1).
}
\]

### Order \(m=6\)

Take
\[
Q_6(x)=(x+1)^2(x^2+x+1).
\]
Then
\[
L_6(x)=x^{10}-4x^9+4x^8+3x^7-6x^6+6x^3-3x^2-4x+4,
\]
so
\[
\boxed{
(\lambda_2,\ldots,\lambda_{12})
=(4,-4,-3,6,0,0,-6,3,4,-4,1).
}
\]

Direct exact substitution verifies
\[
\sum_j\lambda_jj^r=1,\qquad 0\le r\le m
\]
in each case. The verifier checks these moment identities and the polynomial divisibility (2).

## 4. Cost comparison

For the three explicit choices of \(Q\) above, we obtain:

| Order \(m\) | Binomial \(\mathcal C_2\) | Alternative \(\mathcal C_2\) | Reduction |
|---:|---:|---:|---:|
| 4 | 431 | 379 | 12.1% |
| 5 | 1119 | 931 | 16.8% |
| 6 | 2815 | 1905 | 32.3% |

These are exact integer comparisons. They are improvements over the binomial construction, not proofs of the globally minimum possible cost.

## 5. Why the moment identity preserves formal-group cancellation

For a one-dimensional commutative formal group,
\[
[j]_F(T)=jT+\sum_{r\ge2}c_r(j)T^r,
\]
where \(c_r(j)\) is a polynomial in \(j\) of degree at most \(r\), and \(c_r(1)=0\). If (1) holds through degree \(m\), then
\[
\sum_j\lambda_j[j]_F(T)=T+O(T^{m+1}),
\]
and hence
\[
\boxed{
\sum_j\lambda_j\frac{[j]_F(T)}T
=1+O(T^m).
}
\]
Because every \(\lambda_j\) is an integer, the same \(p\)-adic integrality/congruence argument as before remains valid.

## 6. Important limitation

The quadratic proxy \(\mathcal C_2\) is motivated by canonical-height growth, but it is not the exact logarithmic height of the reduced rational combination. Denominator cancellation, common factors, and addition-height inequalities can change the true output height. The candidates must still be tested on an actual elliptic curve before claiming an improvement in the real height/error tradeoff.

No claim about \(\pi\)'s bad approximability follows from this weight optimization.


## 7. Actual elliptic test: the proxy improvement does not transfer

The explicit candidates above were tested on \(E:y^2=x^3-2\), \(P=(3,5)\), with exact rational group arithmetic, for \(s=15,26,41\). In all nine comparisons (orders \(m=4,5,6\)), the alternative weights produced both a larger reduced rational logarithmic height and a larger \(|B_m-1|\) than the binomial weights.

For \(s=41\), the alternative logarithmic heights were approximately \(1.38\), \(2.76\), and \(2.71\) times the binomial heights at orders 4, 5, and 6, respectively. The alternative errors were approximately \(2.63\), \(8.25\), and \(12\) times larger. See `ELLIPTIC_HEIGHT_AWARE_WEIGHTS_ACTUAL_TEST.md` and its exact-arithmetic reproduction script.

**Revised conclusion:** \(\sum|\lambda_j|j^2\) is only a coarse proxy. It does not predict the denominator lcm or the actual reduced height of the rational combination. The binomial weights are also uniquely forced if one requires all multipliers to be at most \(m+2\), so any alternative must introduce larger multiplication indices.


## 8. Structural explanation from primitive EDS divisors

The later note `EDS_PRIMITIVE_DENOMINATOR_BARRIER.md` proves that if \(J\) is the largest multiplier with nonzero weight, then the primitive part of the EDS denominator \(B_{2Js}\) survives in the reduced denominator of the weighted combination for all sufficiently large \(s\). For our CM curve this gives
\[
\liminf_{s\to\infty}\frac{h(B_{m,\lambda}([s]P))}{s^2}
\ge \frac{4J^2\widehat h(P)}{\zeta(2)}.
\]
Thus the actual denominator cost has a rigorous quadratic lower-bound component controlled by the *largest multiplier*. Since any non-binomial integer moment solution must use \(J\ge m+3\), the binomial scheme is optimal for this unavoidable primitive-divisor contribution. This does not assert it minimizes every part of the full rational height.
