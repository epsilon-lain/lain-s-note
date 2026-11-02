# Primitive reflection factors and the remaining quadratic degree cost

2026-11-02. This note continues the established sparse one-dimensional raw-column family. It turns the two rational reflection determinants into explicitly normalized monic integer factors. This reduces the rational-value degree cost by about a factor of two, but the resulting degree is still quadratic in the size parameter. It does not prove a rational-approximation result at exponent two, or a low-degree factorization beyond these two factors.

## 1. The even-size family and its factors

Let n=2r>=4, and define

\[
A_j=\prod_{a=0}^{j-1}a!,\qquad
B_n=A_n/A_r^2,\qquad d_n=B_n^2,
\qquad \Delta_n=r(r-1)=n^2/4-n/2.
\]

The established endpoint-content theorem proves that d_n is the least integer coefficient denominator of

\[
P_n(Z)=\det\bigl([z^s]e^{hz}(jZ+z)^a\bigr),
\quad j,h\in\{0,1\},\quad 0\le s,a<n.
\]

Its leading coefficient is (-1)^r/d_n, its degree is n^2/2, and its first nonzero coefficient is [Z^n]P_n=-(n-1). Rows and columns have the orders used in the previous notes.

Write

\[
L_\pm(Z)=\det\left([z^s](Z/2+z)^a
\bigl(e^{z/2}\pm(-1)^a e^{-z/2}\bigr)\right)_{s,a=0}^{n-1}.
\]

The exact reflection identity gives P_n=L_+L_-. We prove below that

\[
\operatorname{ord}_0L_\pm=r,\quad \deg L_\pm=r^2,
\quad \operatorname{LC}(L_+)=\frac{(-1)^{r(r-1)/2}}{B_n},
\quad \operatorname{LC}(L_-)=\frac{(-1)^{r(r+1)/2}}{B_n}.
\tag{1}
\]

Consequently the explicitly normalized polynomials

\[
F_+(Z)=(-1)^{r(r-1)/2}B_nZ^{-r}L_+(Z),
\quad
F_-(Z)=(-1)^{r(r+1)/2}B_nZ^{-r}L_-(Z)
\tag{2}
\]

are monic, even polynomials in Z, both of degree Delta_n, and

\[
\boxed{F_+,F_-\in\mathbb Z[Z],\qquad
d_nP_n(Z)=(-1)^rZ^nF_+(Z)F_-(Z).}
\tag{3}
\]

The signs are normalization choices only; they do not affect magnitudes. The factor Z^n is removed before comparing degrees.

## 2. Orders at zero and parity

At Z=0, a plus-factor entry vanishes for every odd row index s, independently of the column a: its parity multiplier is 1+(-1)^s. A minus-factor entry vanishes for every even s. There are exactly r vanishing rows in either determinant, and every entry of each such row is divisible by Z. Hence ord_0 L_+>=r and ord_0 L_->=r. Since their product P_n has exact order 2r, both inequalities are equalities.

Changing Z to -Z and then replacing the coefficient variable z by -z multiplies row s of the plus matrix by (-1)^s. In the minus matrix it multiplies row s by (-1)^(s+1). For even n, both determinants therefore satisfy

\[
L_\pm(-Z)=(-1)^rL_\pm(Z).
\]

Division by Z^r leaves an even polynomial.

## 3. Exact leading coefficients

Group the columns first by even a=0,2,...,2r-2 and then by odd a=1,3,...,2r-1. The permutation has sign (-1)^(r(r-1)/2). Expand each column by powers of z:

\[
\sum_{j=0}^a\binom aj(Z/2)^{a-j}v_{j,p},
\qquad p=a\bmod2,
\]

where v_(j,p) is the coefficient vector of z^j(e^(z/2)+(-1)^p e^(-z/2)) in the plus case. Within each parity group, a nonzero determinant term must use distinct j. The least possible sum of these indices is r(r-1), attained exactly when both groups use j=0,...,r-1. This gives degree r^2, provided that coefficient is nonzero.

For either parity group, the binomial determinant is

\[
\det\left[\binom{2b+p}{j}\right]_{j,b=0}^{r-1}
=2^{r(r-1)/2}.
\]

This follows by writing binom(x,j) as a degree-j polynomial with leading coefficient 1/j!: its evaluation determinant is a Vandermonde divided by product j!, and the points 2b+p have twice the usual spacing. The two groups together contribute 2^(r(r-1)); the powers of Z/2 contribute 2^(-r^2), leaving 2^(-r).

Let W be the matrix whose first r columns are coefficients of z^j e^(z/2), and whose last r columns are coefficients of z^j e^(-z/2), j=0,...,r-1. The confluent Vandermonde formula, after dividing row s by s!, gives

\[
\det W=(-1)^{r^2}A_r^2/A_n.
\]

For completeness, the usual derivative form of that formula has columns d^j(x^s)/dx^j evaluated at the two nodes +1/2 and -1/2. Its determinant is (product_(j<r) j!)^2 times the node difference to power r^2. This can also be obtained by coalescing r distinct evaluation nodes at each endpoint in the ordinary Vandermonde determinant. The Taylor coefficient form divides by product_(s<n)s!=A_n.

Passing from the columns of W to their sums and differences multiplies its determinant by (-2)^r. Since r^2+r is even, the resulting vector determinant is 2^r A_r^2/A_n. Including the previous factors and the parity-grouping sign gives the stated leading coefficient of L_+.

For L_-, the sum and difference groups are interchanged. The interchange has sign (-1)^(r^2)=(-1)^r, giving the leading coefficient in (1). In particular neither leading coefficient is zero.

## 4. Why the normalized factors are integer polynomials

Equations (1)--(2) make F_+ and F_- monic polynomials over Q. The integer polynomial (-1)^r d_nP_n is monic, and (3) factors it over Q as the product of the monic polynomials Z^n, F_+, F_-. Gauss's lemma implies that every monic rational factor of a monic integer polynomial has integer coefficients. Thus both F_+ and F_- lie in Z[Z]. Their leading coefficient one also makes them primitive.

This argument uses the proved exact even denominator. It does not assume that arbitrary rational reflection factors separately admit equally small clearings.

## 5. One factor retains a small actual-period value

Put t=2*pi, omega=it, and

\[
U_n=t^{n^2}A_n^2/A_{2n},\qquad
V_n=\sqrt{d_nU_n}/t^{n/2}.
\]

The previous actual-period estimate |P_n(omega)|<=U_n and (3) imply

\[
\boxed{\min\{|F_+(\omega)|,|F_-(\omega)|\}\le V_n.}
\tag{4}
\]

This does not identify the smaller sign uniformly in n. It does suffice for the comparison below, because both factors have the same degree, the same derivative bound, and the same exact rational-value denominator property.

Using the established denominator asymptotic, with C=3/4+ln(pi/sqrt(2)), gives

\[
\ln V_n=-\tfrac14n^2\ln n+\tfrac12C n^2+O(n\ln n),
\quad
\frac{\ln V_n}{\Delta_n}=-\ln n+2C+o(1).
\tag{5}
\]

Both the gain and the rational-evaluation degree have been approximately halved. Their quotient still has the same leading threshold.

## 6. The direct rational comparison after splitting

Let p/q be reduced, q>=1, and write 2p/q=a/b in lowest terms, so b=q/gcd(2,q). Assume b>=2. Because F_+ and F_- are monic even integer polynomials of even degree Delta_n, their values at ia/b are rational with the exact reduced denominator

\[
\boxed{\operatorname{den}_{\mathrm{reduced}}F_\pm(ia/b)=b^{\Delta_n},
\quad |F_\pm(ia/b)|\ge b^{-\Delta_n}>0.}
\tag{6}
\]

Indeed, after multiplication by b^Delta_n, the numerator reduces modulo every prime dividing b to i^Delta_n a^Delta_n, which is nonzero. Rewriting the factors as polynomials in Z^2 does not halve this arithmetic denominator again: the new argument has denominator b^2.

For an explicit derivative bound, put S_n=n(n-1)/2 and

\[
M_n=B_n n!2^n n^2 5^{S_n}.
\]

On any segment joining omega to ia/b of length at most one, one has 1<|Z|<8. Each reflection entry has magnitude at most 2*5^a and derivative at most a*5^(a-1). Thus

\[
|L_\pm(Z)|\le n!2^n5^{S_n},\quad
|L_\pm'(Z)|\le n!2^{n-1}S_n5^{S_n-1}.
\]

Differentiating (2), and using |Z|>=1 and r+S_n/10<=n^2, proves |F_\pm'(Z)|<=M_n on the segment.

Choose any even n>=4 with V_n<=1/(2b^Delta_n), and choose the smaller factor in (4). Equations (4), (6), and the derivative bound then yield

\[
\boxed{\left|\pi-\frac pq\right|
\ge\frac{1}{4b^{\Delta_n}B_n n!2^n n^2 5^{S_n}}.}
\tag{7}
\]

If the segment has length greater than one, the left side exceeds 1/2 and (7) is automatic. No proximity assumption remains. The sufficient condition is exactly certifiable with rational arithmetic by squaring it:

\[
4b^{2\Delta_n}t^{n^2-n}\frac{d_nA_n^2}{A_{2n}}\le1.
\]

Replacing t by 44/7 gives a rational sufficient test.

Let kappa=e^(3/2)pi^2/2=e^(2C). For any fixed c>kappa, taking n to be the least even integer at least cb makes the sufficient true-t condition hold for all sufficiently large b, by (5). The logarithm of the right side of (7) is

\[
\boxed{-\tfrac12c^2b^2\ln b+O_c(b^2).}
\tag{8}
\]

Here ln B_n=(n^2/4)ln n+O(n^2), and the other derivative factors cost only O(n^2). This halves the leading logarithmic loss in the previous full-polynomial comparison, but its scale remains exp(-O(q^2 ln q)). No q^(-2) lower bound follows, and no universal obstruction to more sophisticated uses of the factors has been proved.

## Verification and remaining question

The companion standard-library verifier reconstructs L_+ and L_- exactly for even n=4,...,12, verifies their degrees, leading coefficients, orders, monic integer normalizations, and parity, and compares their product to the original raw determinant at additional rational arguments. It verifies exact rational-value denominators, derivative coefficient majorants, and finite rational sufficient-height tests. The general proofs are the arguments above.

Exploratory factorization of the previously known n=4,6,8,10 examples over Q found only the factor Z^n and these two nonzero factors; their degrees are respectively (2,2), (6,6), (12,12), and (20,20). This is finite evidence only. Neither general irreducibility nor the impossibility of an additional low-degree factor has been proved. Establishing such a factor, or a different arithmetic comparison with lower effective degree, remains a separate problem.


## Prior checkpoints and stronger companion bound

The [endpoint-content theorem](../../pi-primitive-period-2026-11-02/denominator-reflection-route/ENDPOINT_CONTENT_AND_DENOMINATORS.md) supplies the exact even denominator. The reflection identity and order at zero are also recorded in the [general near-period root proof](../../pi-primitive-period-2026-11-02/near-root-general/GENERAL_NEAR_PERIOD_ROOTS.md). The companion [normalized Gram disk argument](../next-local-period-sensitivity/LOCAL_PERIOD_DISC.md) gives a stronger exp(-O(q log q)) error bound without factoring. This note records the structural integer factors and their direct comparison, not the strongest bound of the present checkpoint.
