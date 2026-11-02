# Exact norm design on the integer period quotient

2026-11-02. This continues the main constructive family with full parity freedom. It uses the rank-two period system as an actual design quotient, retaining the many integer representatives of each period. The purpose is to avoid spending the short-vector search entirely on zero-period total derivatives and to replace sampled integrand values by a proved rational upper objective.

## 1. Fixed family and exact integer period split

Fix m>=1,c=4m+1. The integer primitive coordinates are R(0)=0,deg R<=2c-1,alpha in Z, with

N=(1+x^2)R'-2(c-1)xR+alpha(1+x^2)^(c-1)
=x^(2m)(1-x)^(2m)H(x).

Write z=R(1). The integral is I=z/2^(c-1)+alpha*pi/4, and the reduced candidate is p/q=-4z/(2^(c-1)alpha), alpha!=0.

Let K be any integer basis of the saturated endpoint kernel, written as rows in primitive coordinates. Its period matrix P has columns (alpha,z). It has rank two, as already proved by the endpoint CRT construction. A Smith decomposition SPT=D, with T unimodular, changes the kernel basis to T^tK. The first two rows give integer lifts C_0,C_1 of a basis of the period image; the remaining rows Z_1,...,Z_r have both periods zero. They form the full integer kernel of the period map, not only a rational nullspace with individually cleared denominators.

Thus every integer candidate is represented uniquely as

H=Cv+Zu,

where v is an integer two-vector specifying the period and u is an integer vector choosing a representative with that same period. Here C and Z are coefficient matrices of the associated multiplier polynomials, with columns the two period lifts and the zero-period lifts respectively. The primitive R coordinates are carried alongside them, so integer reconstruction and actual gcd reduction remain exact.

## 2. A rational Hilbert objective gives a rigorous quality upper bound

For H(x)=sum h_j x^j define

||H||_2^2=integral_0^1 H(x)^2 dx=h^t G h,
G_(ij)=1/(i+j+1).

The Hilbert matrix G is rational positive definite. The proved pointwise kernel inequality x(1-x)/(1+x^2)^2<2/11 gives

I^2 <= (4/121)^(2m)||H||_2^2 integral_0^1(1+x^2)^(-2)dx
< (9/14)(4/121)^(2m)||H||_2^2,

since the last integral is 1/4+pi/8<9/14. Consequently

\[
\boxed{Q^2<\frac{72}{7}\frac{q^4}{\alpha^2}
\left(\frac4{121}\right)^{2m}h^tGh.}
\]

The q is the actually gcd-reduced denominator. Every quantity in the upper expression is rational. Unlike an l1 norm, the exact Hilbert quadratic form retains cancellation among the polynomial coefficients in the integrand norm. It still bounds the absolute integral by Cauchy--Schwarz, so it does not capture every possible signed-integral cancellation.

## 3. The exact quotient metric and a uniform integer lift bound

Put A=Z^tGZ, B=Z^tGC, D=C^tGC. For a fixed period vector v the real minimum is achieved at u_*=-A^(-1)Bv, and completing the square gives

\[
\boxed{||Cv+Zu||_2^2=v^t\mathcal Sv+(u-u_*)^tA(u-u_*),
\quad\mathcal S=D-B^tA^{-1}B.}
\]

The Schur complement S is rational positive definite. Coordinatewise rounding already proves the existence of an integer lift satisfying

||H||_2^2 <= v^t S v +(1/4)sum_(ij)|A_(ij)|.

A more useful constructive bound uses exact nearest-plane rounding in a zero-period basis. Let Z_i^* be its Gram--Schmidt orthogonal vectors in the Hilbert inner product, with d_i=||Z_i^*||_2^2. Reverse nearest-plane rounding has rounding errors e_i in [-1/2,1/2]. Orthogonality gives the exact identity and uniform bound

\[
\boxed{||H||_2^2=v^t S v+\sum_i e_i^2d_i
\le v^t S v+\frac14\sum_i d_i.}
\]

The radius R_m^2=(1/4)sum d_i is independent of v. It is computable from rational numbers, without square roots or pi values. It records the cost of reaching an integer representative; the real Schur minimum alone is not passed off as an achieved integer cost.

Combining the preceding two displays gives the all-scale charged construction certificate

\[
\boxed{Q^2<\frac{72}{7}\frac{q^4}{\alpha^2}
\left(\frac4{121}\right)^{2m}
\left(v^t S v+R_m^2\right).}
\]

No bound on the growth of R_m or the resulting reduced q is assumed. This is a computable cost rule, not an asymptotic success assertion.

## 4. Declared generation algorithm and bounded controls

The implemented rule does the following:

1. Split the saturated integer endpoint kernel by its Smith period quotient.
2. Integrally precondition the zero-period basis using Euclidean LLL on H coefficients, keeping the original primitive coordinates as exact guards. This only changes the basis; it is not claimed to solve a Hilbert closest-vector problem.
3. Compute the exact rational Hilbert Gram, its quotient Schur metric, and its zero-period orthogonal basis.
4. Add two declared rational cost coordinates for alpha and R(1), with the same five dyadic weights as the preceding bounded construction. Reduce the two-dimensional quotient by exact Gauss reduction, then take the two reduced basis directions and their sum/difference.
5. Lift each quotient direction by exact Hilbert nearest-plane rounding. Reject alpha=0; reconstruct R,N,H and fully reduce p/q. Compute the charged Hilbert certificate from the actual lift.
6. Only afterwards use a rational Machin interval to check achieved quality.

No numerical pi, continued fraction, or target fraction enters steps1--5. The quotient direction score is a rational metric; it is not described as a global minimizer of the nonlinear gcd-normalized quality.

The bounded m=1,2,3 run gives 12,14,11 distinct positive candidates respectively. Best post-check qualities were approximately0.00340631,9.34167,10.78224, at 355/113,3217/1024,1146053/364800. The first is a known benchmark; no new record or global optimum is claimed. All stored candidates satisfy the exact lift identity and the rational joint-quality upper bound.

The exact Gram objective also gives a concrete improvement in the certificate for the same multipliers. For the three retained best rows, the previous l1 quality upper bounds are approximately4766.55,670764.51,2.43452e10, while the square roots of the exact Hilbert quality-squared upper bounds are approximately9.06510,210.31091,14090.17621. These comparisons are recorded as exact rational squared ratios in hilbert_bound_comparison.json; the square roots and decimal figures are display-only. The achieved fractions and their actual qualities do not change when the bound is improved.

This control does not justify increasing the same search range. Its reusable result is the exact rank-two period design with a charged integer lift, which can be combined with a better kernel or with a period recurrence from the elliptic bridge. It makes explicit which part of the cost comes from the desired period and which comes from finding an integer multiplier realizing it.

## 5. Arrange an odd common period factor before lifting

There is a concrete congruence design inside this quotient. Set

\[
D_m=\operatorname{odd}\operatorname{lcm}(1,\ldots,2m).
\]

The endpoint-period theorem already proves that alpha is divisible by D_m for every integer endpoint vector. Impose the additional condition R(1)=0 modulo D_m on the two-dimensional period quotient. Every zero-period correction preserves both congruences, so the common factor is arranged before integrand minimization rather than searched for afterwards.

Let the second row of the two-period basis matrix be (z_0,z_1), and put g=gcd(z_0,z_1). The proved odd period image gives gcd(g,D_m)=1. Choose integers u,v with u z_0+v z_1=g and define

\[
T_{D_m}=\begin{pmatrix}uD_m&-z_1/g\\vD_m&z_0/g\end{pmatrix}.
\]

Then det(T_(D_m))=D_m, the second period row becomes (gD_m,0), and the first period row remains divisible by D_m. Its image is exactly the congruence kernel z_0v_0+z_1v_1=0 mod D_m, since that kernel has index D_m in Z^2.

Replace the two period lifts C by C T_(D_m), keep the zero-period basis Z, and perform the same Hilbert quotient/lift construction. Now

\[
D_m\mid\alpha,\quad D_m\mid R(1),\qquad
\gcd(4R(1),2^{c-1}\alpha)\ge D_m.
\]

The price is explicit: a quotient sublattice of index D_m, a changed Schur metric T_(D_m)^t S T_(D_m), and the resulting integer representative cost. The congruence does not provide a free quality gain. Overall coefficient content is recorded separately so that a global rescaling is not confused with a genuine common period factor.

The same bounded m=1,2,3 controls, with the same five cost weights, were repeated on this forced-period sublattice. There were 12,14,14 distinct positive candidates. At m=2 the best selected fraction changes from unrestricted 3217/1024 (quality about9.34167) to355/113 (quality about0.00340631). At m=3 the best selected fraction remains1146053/364800. The retained m=2 and m=3 best primitive coordinate vectors have content one, while gcd(alpha,R(1)) is exactly3 and15 respectively. Thus those common factors are not merely global coefficient scaling.

This is a verified design mechanism and a bounded selection comparison, not a proof that imposing this congruence improves every candidate or every scale. A known benchmark is not a new approximation record. The exact Hilbert and actual-gcd certificate is retained for every constructed vector.

## Reproduction

`hilbert_quotient_lift.py` requires SymPy for integer Smith/LLL basis operations and uses exact rational arithmetic afterwards. It reads the endpoint-kernel generator from the earlier parity-beta checkpoint. Run it with `--force-period` for the congruence variant. `verify_hilbert_output.py` independently checks both saved datasets: 77 candidates and60 basis vectors in total, including the period-sublattice determinant, forced gcd, coefficient content, Schur projection, nearest-plane error ledger and rational certificate. It uses the Python standard library only. The sampling-only prototype `quotient_lift.py` is a support module and earlier control; its results are not used as evidence for the Hilbert certificate.
