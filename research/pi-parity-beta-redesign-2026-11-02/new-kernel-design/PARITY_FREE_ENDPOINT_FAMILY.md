# A parity-free primitive family with the same odd period bill

2026-11-02. This is a constructive variation of the integer-antiderivative coordinates. It removes the imposed odd parity of R and hence removes an unrequested zero at x=-1. All candidate costs below refer to the actually reduced fraction.

## 1. Exact primitive coordinates and the one missing numerator coefficient

Fix c>=2. Let R be any rational polynomial with R(0)=0 and deg R<=2c-1, and let alpha be rational. Define

\[
N=(1+x^2)R'-2(c-1)xR+\alpha(1+x^2)^{c-1}.
\]

Then

\[
I=\int_0^1\frac{N}{(1+x^2)^c}dx
=R(1)/2^{c-1}+\alpha\pi/4.
\]

No logarithmic period is introduced by allowing even powers in R. In fact the image of this rational coefficient map is exactly the space of polynomials N of degree<=2c for which

\[
\boxed{[x^{2c-1}]N=0.}
\]

The odd-R plus alpha columns give all even numerators, by the parent's invertible primitive matrix. The even-R columns (with zero constant term) give all odd numerators of degree<=2c-3: in ascending odd powers this is a triangular map with nonzero diagonal 2,4,...,2c-2. Its apparent top coefficient at x^(2c-1) vanishes because differentiation of x^(2c-2)/(1+x^2)^(c-1) loses that leading power. The two parity parts are independent, proving the assertion.

Equivalently, for deg N<=2c-2 the logarithm coefficient is automatically zero. The odd part is xP(x^2), deg P<=c-2; substituting t=x^2 reduces its integral to a sum of powers (1+t)^(-s) with s>=2. Thus one may generate N directly in that degree range and recover the unique primitive afterward.

## 2. Actual endpoint zeros without an extra reflected zero

Require

\[
N=x^A(1-x)^B H(x),\qquad A,B>=1,
\quad 2c>=A+B+2.
\]

Write D=2c-A-B. If deg H<=D, the single log-free condition above is one explicit linear condition. For the proposed fair parameters

\[
c=4m+1,\quad A=B=2m,\quad D=4m+2,
\]

it is

\[
H_{D-1}=2m H_D.
\]

An unconstrained direct subfamily is deg H<=4m. The full rational kernel has dimension 4m+2. The old odd-R design has dimension m+2 at the same c and actual endpoint orders 2m,2m; it additionally forces a zero of order 2m at -1. This dimension comparison is not by itself an approximation-quality theorem.

## 3. Odd-prime period image remains explicit

Let Lambda be the image of the saturated integer endpoint kernel under (R,alpha)->(alpha,R(1)). Let O_A be the largest odd integer<=A (O_A=1 for A=1), and M=max(O_A,B). For every odd prime p,

\[
\boxed{\Lambda\otimes\mathbb Z_p
=p^{\lfloor\log_p M\rfloor}\mathbb Z_p\times\mathbb Z_p.}
\]

At zero the condition N=O(x^A) means R/(1+x^2)^(c-1)=-alpha arctan(x)+O(x^(A+1)). Taylor integrality requires precisely alpha/(2j+1) to be p-integral for 2j+1<=A. At one, prescribe z=R(1). The expansion of arctan(1+u)-pi/4 has coefficient h_(k-1)/k, where 1/(2+2u+u^2)=sum h_j u^j. Every nonzero h_j is a dyadic unit, and h_(p^e-1) is nonzero for odd p. This gives the cutoff B. The constant z/2^(c-1) imposes no odd-prime condition.

Conversely these local data specify R modulo x^(A+1) and (x-1)^(B+1). Ordinary polynomial CRT over Z_p produces a polynomial of degree<=A+B+1<=2c-1; no parity restriction is imposed. This proves sufficiency and the full image. Localization of the integer kernel and image preserves this conclusion.

For A=B=2m the mandatory odd bill is therefore odd(lcm(1,...,2m)), exactly the earlier cutoff. A two-adic formula is not inferred from the odd-prime proof. Exact small period images computed by Hermite normal form are

| m | c | integer kernel rank | image in (alpha,R(1)) |
|---:|---:|---:|---|
| 1 | 5 | 6 | Z x Z |
| 2 | 9 | 10 | 3Z x Z |
| 3 | 13 | 14 | 15Z x 2Z |

## 4. A stronger candidate-wise analytic budget

For A=B=2m,c=4m+1 set h=||H||_1. The elementary inequality

\[
\frac{x(1-x)}{(1+x^2)^2}<\frac2{11}
\qquad(0<=x<=1)
\]

follows from 2x^4+15x^2-11x+2>0. The quadratic part is nonnegative outside [1/3,2/5]. Inside that interval it is at least -1/60, while 2x^4>=2/81>1/60. Consequently

\[
|I|\le\frac\pi4 h(4/121)^m,
\qquad
\boxed{Q<\frac{22}{7}\frac{q^2h}{|\alpha|}(4/121)^m.}
\]

This gives the sufficient construction budget q^2h/|alpha|=o((121/4)^m). It is a bound for the new candidate coordinates. It is not a free improvement for an unchanged old numerator: the embedding of an old multiplier is H_new=(1+x)^(2m)H_old(x^2), whose coefficient norm can grow by 4^m. The useful freedom is the ability to choose H without that extra factor.

## 5. Bounded generation controls

The companion full_parity_search.py builds the full saturated integer kernel by a Smith right transform, then uses a declared rational sample grid and the five fixed cost weights 2^(-e), e in {0,8m,12m,16m,20m}. No pi value or target fraction enters generation. Exact Machin intervals are used only afterward to evaluate candidate quality.

Recovering 355/113 at c=5 is a finite design benchmark. The previous fixed odd-R generator first found that benchmark at its tested c=9; this does not prove that the old mathematical family cannot represent it at c=5 with other coefficients. All primitive data and exact quality bounds are retained in the result JSON when a run completes. Sampling features are not integral upper certificates.

The initial nine-point grid was insufficient at the m=3 kernel rank14 and produced no positive nonzero-period LLL basis rows there; no completed result was claimed for it. The final grid uses max(9,c+2) sample points, providing more features than the integer-kernel dimension in the three declared scales. No broader scan is part of this experiment.

The completed final run gives:

| m | c | sample points | best retained basis-row fraction | Q (display) |
|---:|---:|---:|---|---:|
| 1 | 5 | 9 | 355/113 | 0.003406311930138... |
| 2 | 9 | 11 | 104703/33328 | 1.37175951497... |
| 3 | 13 | 15 | 229592869/73081680 | 8.82238318148... |

The m=1 retained data include alpha=339 and

R coefficients at x^1,...,x^9:
(-339,0,-1243,-2,-1630,-78,-745,-164,-59).

They have R(1)=-4260 and therefore -4R(1)/(2^4 alpha)=355/113. Their numerator is x^2(1-x)^2H(x), with

H=(0,-8,83,-286,376,-118,-59)

listed in ascending powers. The endpoint and no-log conditions hold exactly. This is a known benchmark produced at lower tested pole order, not a new best pi fraction or an optimality theorem.

`verify_parity_free_output.py` independently rechecks all 24 distinct retained candidates with standard-library integer/rational arithmetic, including the complete derivative numerator, endpoint quotient, gcd-reduced fraction, rational Machin quality interval and the proposed norm certificate. It also checks the three period images by containment and the gcd of all two-by-two minors of the saved integer-kernel period matrix.

The generation script requires SymPy and uses the same exact nearest-integer-floor correction as the parent's saturated search, avoiding PythonMPQ-to-float rounding in SymPy1.14. It does not alter the integer LLL size or Lovasz assertions. The output verifier has no SymPy dependency.
