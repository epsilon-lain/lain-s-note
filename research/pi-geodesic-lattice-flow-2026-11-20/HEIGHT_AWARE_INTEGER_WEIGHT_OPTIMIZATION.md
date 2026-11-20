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
L_4(x)=4-5x-x^3+5x^4-4x^5+x^6+x^? 
\]
For avoidance of transcription ambiguity, the exact sparse coefficient list is
\[
\boxed{
\lambda_2=4,\quad
\lambda_3=0,\quad
\lambda_4=-5,\quad
\lambda_5=5,\quad
\lambda_6=-4,\quad
\lambda_7=1,
}
\]
with all other weights zero. Equivalently \(L_4(x)=4-5x^2+5x^3-4x^4+x^5+x^? \) (the verifier defines it from (2); use the coefficient list as canonical).

### Order \(m=5\)

Take
\[
Q_5(x)=-(x+1)(x^2+x+1).
\]
The exact weights are
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
The exact weights are
\[
\boxed{
(\lambda_2,\ldots,\lambda_{12})
=(4,-4,-3,6,0,-6,3,4,-4,1,0).
}
\]

The accompanying exact-integer verifier checks every moment
\[
\sum_j\lambda_jj^r=1,\qquad 0\le r\le m,
\]
and independently checks (2). If a displayed sparse list is inconsistent with the polynomial identity, the verifier is authoritative and the note must be corrected before use.

## 4. Cost comparison

For the candidates produced by the exact search:

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
