# A complete two-variable packet with an imbalanced raw minor and a new root

2026-11-02. This note gives one finite raw-column certificate beyond the balanced normal-block products. Both spatial derivative directions occur in the complete source. The normal row and column profiles differ, so a contributing determinant term must cross normal blocks. Its residual factor has roots different from the earlier consecutive one-dimensional family. The certified new root is not a better approximation to the first exponential period. No asymptotic coupled family or critical-weight theorem is claimed.

## 1. Model, complete packet, and raw columns

Use K=2, m=2, H=6, time weights w0=v0=1, theta=1/10, spatial column weights w=(2/5,1/2), and spatial source weights v=w/theta=(4,5). Thus K theta^m=1/50<1 and K(w0/v0)theta^m=1/50<1. These rational weights satisfy the basic numerical inequalities and all displayed budgets; they are not asserted to satisfy the original separated-weight theorem's successive lower thresholds or its arithmetic logarithmic-weight regime.

The complete source consists of all (j,s,beta1,beta2) with j=0,1 and

\[
s+4\beta_1+5\beta_2<6.
\]

For each center, these are exactly beta=(0,0), s=0,...,5; beta=(1,0), s=0,1; and beta=(0,1), s=0. There are 18 rows. Order centers first, then the three beta types in the preceding order, then s.

Choose 18 raw columns Y^h X_1^k X_2^a, ordered first by h=0,1 and then by the following ordered list (k,a):

\[
\mathcal A=\{(0,0),(0,1),(0,2),
(1,0),(1,1),(1,2),(1,3),(1,4),(1,6)\}.
\]

Every column satisfies h+(2/5)k+a/2<=6. In fact its maximum weight is 22/5, so the minimum column slack is 8/5. No column mixture or partial source packet is substituted.

Define the ordinary Taylor matrix

\[
B(Z)_{(j,s,\beta);(h,k,a)}
=[z^su_1^{\beta_1}u_2^{\beta_2}]
e^{hz}(jZ+z+u_1)^k(jZ+z+u_2)^a,
\qquad P(Z)=\det B(Z).
\]

At Z=2*pi*i the formal row equals the actual exponential endpoint row, because h is an integer. At general Z it is only the formal-period polynomial.

Take X_2 as the base direction and X_1 as the normal direction. The source normal profile is (14,4) at beta1=0,1, whereas the column profile is (6,12) at k=0,1. The matrix is block upper triangular in the support sense, but its diagonal blocks have sizes 14-by-6 and 4-by-12, not square. Any surviving determinant permutation must assign eight beta1=0 rows to k=1 columns. The total normal-degree surplus is 12-4=8. Thus the balanced-profile product theorem does not apply to this minor.

This certifies a finite mismatch in the chosen normal grading. It is not an assertion that no other change of coordinates or basis can simplify the polynomial.

## 2. Exact determinant and arithmetic content

The exact identity is

\[
\boxed{P(Z)=\frac{Z^{21}}{864}
S(Z),\qquad
S(Z)=13Z^6+720Z^4+6480Z^2-63936.}
\tag{1}
\]

The least positive integer coefficient denominator is 864: multiplication by 864 makes all coefficients integers, and the leading coefficient 13/864 is reduced. The resulting integer polynomial Q=864P is primitive, of degree 27 and leading coefficient 13.

The identity has a finite exact certificate. After multiplying row (s,beta) by s! beta1! beta2!, all entries are integer polynomials. In every determinant term, its Z degree is at most the total column spatial degree minus the total source spatial derivative degree, here 50-6=44. The verifier evaluates the integer raw determinant and the proposed right side at all 45 integers Z=0,...,44. Their equality at more points than this proved degree bound establishes polynomial identity (1). Additional rational evaluations are checked independently. This is an exact polynomial certificate, not a numerical extrapolation from nearby values.

Write f(x)=13x^3+720x^2+6480x-63936, so S(Z)=f(Z^2). Modulo 11, this degree-three polynomial has no root; the verifier checks all eleven residues. Thus f is irreducible over Q. In fact S is also irreducible over Q: if alpha is a root of f, then [Q(alpha):Q]=3 and its field norm is 63936/13, which is not a rational square because its 13-adic valuation is -1. Hence alpha cannot be a square in Q(alpha); adjoining its square root has degree two and yields a root of S of degree six.

In particular P is not a constant times a product of the earlier whole consecutive period polynomials P_n. After removing zero powers, S is irreducible of degree six. The nonzero degrees of those whole polynomials are zero for n=2,3, four for n=4, six with a reducible residual for n=5, and at least twelve for n>=6. This does not preclude relationships with more general one-dimensional lacunary constructions.

## 3. A certified new imaginary root, with an honest quality comparison

Put g(t)=S(it)=f(-t^2), a real polynomial. Exact rational arithmetic gives opposite signs at

\[
\ell=6302063/10^6,\qquad u=6302064/10^6:
\qquad g(\ell)>0,\quad g(u)<0.
\]

On this interval, 39<t^2<40 and

\[
g'(t)=-6t(13t^4-480t^2+2160)<0,
\]

because the quadratic 13y^2-480y+2160 is positive for 39<=y<=40. Thus there is exactly one root t_* in (ell,u), producing a simple nonzero imaginary root i*t_* of P. It lies above 44/7>2*pi, so its distance from the actual period is positive and is bounded below by ell-44/7.

For a consecutive-exponent comparator, replace only (1,6) in the raw column list by (1,5), keeping all rows, frequencies, weights, and profile counts. Its exact determinant is

\[
P_{\rm consecutive}(Z)
=\frac{Z^{19}}{864}(Z^2+24)(Z^4+36Z^2-144)
=-Z^{14}P_5(Z).
\tag{2}
\]

It has an imaginary root i*c with c=sqrt(18+6*sqrt(13)). Exact sign checks put c in (62954/10^4,62956/10^4). Both endpoints exceed 44/7, and c<ell. Therefore

\[
0<c-2\pi<t_*-2\pi.
\]

The new certified root is strictly farther from the first period than this comparator root. A changed zero location is not, by itself, an approximation improvement. Replacing the column leaves the profiles unchanged, so the comparator is not called a balanced-profile example.

## 4. Exact rational-value denominator

Let a,b be coprime integers with b>=2. Since S is even, S(ia/b) is rational. Its reduced denominator is exactly

\[
\boxed{\operatorname{den}_{\rm reduced} S(ia/b)
=\frac{b^6}{\gcd(b,13)}.}
\tag{3}
\]

To prove this, multiply by b^6. For a prime p dividing b other than 13, only the unit-leading term -13a^6 survives modulo p, so there is no p cancellation. If 13 divides b, that term has valuation exactly one while every other term contains b^2 and has valuation at least two. Thus exactly one factor 13 cancels, regardless of its multiplicity in b.

Likewise Q(ia/b) is purely imaginary, and the reduced denominator of its real coefficient after division by i is

\[
\boxed{\operatorname{den}_{\rm reduced}\bigl(Q(ia/b)/i\bigr)
=\frac{b^{27}}{\gcd(b,13)}.}
\tag{4}
\]

The extra numerator a^21 is coprime to b and does not change the cancellation argument. Both values are nonzero, including when b is divisible by 13. These exact denominators show that the new residual factor does not eliminate rational-specialization cost.

## 5. Positive surplus alone is not enough

A smaller complete control packet has H=4, v0=1, spatial v=(2,3), theta=1/3, w0=1/6 and spatial w=(2/3,1). Per center, its rows are beta00 with s=0,...,3, beta10 with s=0,1, and beta01 with s=0. Take raw columns h=0,1 with spatial exponents (a,0), a=0,...,4, and (a,1), a=0,1, ordered by h and then by this list. All columns are legal, with minimum slack 7/6. The two basic density inequalities are strict.

In the X_2 normal grading, row counts are (12,2), column counts are (10,4), and the normal-degree surplus is two. Nonetheless the exact ordinary Taylor determinant is only

\[
\boxed{P_{\rm control}(Z)=-2Z^{13}.}
\tag{5}
\]

It is nonzero at the actual period, but has no nonzero root at all. This gives an exact control against treating a positive profile surplus as a sufficient condition for a useful new zero.

## Verification and scope

The companion standard-library verifier enumerates both complete source packets and all chosen raw columns, checks the numerical constraints and minimum budget slacks, establishes (1), (2), and (5) through the proved degree bound, certifies the rational root brackets and cubic irreducibility modulo 11, and checks (3)--(4) at many reduced rational inputs including multiples of 13. Full row and column lists are included in its JSON output.

The construction establishes one nonzero complete-source raw minor with an imbalanced normal profile and a new algebraic root. It does not establish an infinite coupled family, a better root than the preceding one-dimensional examples, a small primitive value at the period, a critical half-density construction, or any statement deciding pi's bad approximability.


## Related proofs

The [normal-profile theorem](../normal-filtration-general/NORMAL_PROFILE_AND_STRAIGHTENING.md) gives the exact support and factorization criteria used to distinguish the profiles here. The old consecutive-family degree formulas used in the nonproduct argument are proved in the [endpoint-content theorem](../../pi-primitive-period-2026-11-02/denominator-reflection-route/ENDPOINT_CONTENT_AND_DENOMINATORS.md).
