# Independent review of the multivariate full-packet escape certificate

2026-11-01. This review checks the parent researcher's fiber construction and finite complete-packet witness. It distinguishes a single nonzero transport term from an aggregated polynomial coefficient.

## 1. Fiber factorization

Let B be a finite lower set in nonnegative integer m-space. For each fixed y=(beta_2,...,beta_m), its first-coordinate fiber is {0,...,n_y-1}. Keep the same y labels in the target set and choose any n_y distinct nonnegative first-coordinate targets J_y. Order the y labels by a linear extension of their componentwise order.

The multivariate Pascal matrix is block upper triangular: its entry contains binom(y,beta'), which is zero unless beta'<=y. Each diagonal block is the one-variable Pascal matrix for J_y. Therefore its determinant is the product of the positive one-variable Vandermonde quotients. In particular, replacing n_y-1 by n_y in k eligible fibers gives total error exponent k e_1 and coefficient product_y n_y. Eligibility includes the homogeneous total-degree cap and column-type support; it is not supplied just by the Pascal determinant formula.

This refutes a dimension-free square of the total number of replacements as a transport-monomial lower bound. It does not invalidate the one-variable k^2 law within one fiber. At fixed dimension and weights, the number of available fibers is lower order than the total row count, so the finite first-order examples do not by themselves decide the asymptotic normalized cost.

## 2. The complete finite packet

Use m=2, K=2, H=20, v_0=w_0=1, theta=1/2, and transverse weights (2,5). The complete source packet has 156 rows and the weighted polynomial space has 242 raw monomial columns. In center j=0, total-order d=10, the source transverse types are

\[
(0,0),(1,0),(2,0),(3,0),(0,1).
\]

Replace the last point in each of its two first-coordinate fibers:

\[
(0,7,(3,0))\mapsto(0,6,(4,0)),\qquad
(0,9,(0,1))\mapsto(0,8,(1,1)).
\]

Both new rows have weighted cost 22, outside the strict source cutoff 20. Both preserve total order 10, have nonnegative time index, and fit the polynomial transverse support. The five-by-five Pascal determinant is 4. All other groups use identity transport. Thus the normalized full transport coefficient is exactly 4 zeta_1^2, although two rows are replaced.

## 3. Independent modular nonzero checks

`verify_full_packet_review.py` constructs matrix entries by multiplying univariate polynomial factors and convolving with the exponential Taylor series. This is independent of the originating probe's double Taylor-allocation sum. It uses reduced-row elimination for pivot selection and a separate determinant routine.

Take the prime 65537, with i represented by 256, since 256^2=-1 in this field. The source rational centers are 44i/7 and 710i/113. Reduction is a homomorphism from the localized ring Z[i,1/7,1/113] to this finite field; it is not a homomorphism from the entire characteristic-zero field Q(i).

Each ordinary time row is multiplied by s! to avoid factorial denominators. All factorials used are invertible modulo the prime. On the same 156 selected columns, the independent calculation gives:

| Matrix | Determinant modulo 65537, with s!-scaled rows |
|---|---:|
| Target period matrix, formal common period W=1 | 21680 |
| Original rational source matrix | 35592 |
| Original period source matrix, W=1 | 21281 |

The source row-factorial product is 26110 and the target product is 29542. Their ratio is 63. Thus the normalized target and rational-source determinants are 57507 and 36181, respectively. If using the scaled matrices in the transport formula itself, the displayed transport coefficient becomes 252 zeta_1^2; the factor 4 belongs to the original normalized coefficient convention.

A nonzero modular specialization proves that the target determinant is a nonzero polynomial in W with rational coefficients. Known transcendence of 2 pi i then proves its nonzero value at the actual period. The nonzero rational-source modular value independently certifies its original complex determinant. No unknown interpolation threshold is needed for this finite certificate.

## 4. Aggregation was checked separately

A single nonzero term would not rule out cancellation with other target row sets. We therefore also checked the fully aggregated coefficient of eta^2 zeta_1^2 in

\[
\det(B_L+\eta R(\zeta)B_O),
\]

before restoring the diagonal raw-column prefactor. Set the other formal error coordinate to zero. The first-order boundary matrix E_1 has 44 potentially active rows. In the s!-scaled convention its row is

\[
(E_1)_\rho=-s(\beta_1+1)B_{(j,s-1,\beta+e_1)}^{\rm scaled}
\]

when the shifted row lies outside the packet, and is zero otherwise. The extra s is necessary because source and target time orders differ by one.

At W=1 modulo 65537, let X=B_L^{-1}E_1. Independent computation gives tr X=-6 and tr(X^2)=6. Hence

\[
[\eta^2\zeta_1^2]\det(B_L+\eta R B_O)
=\frac{\det B_L}{2}\big((\operatorname{tr}X)^2-\operatorname{tr}(X^2)\big)
=57067\pmod{65537},
\]

which is nonzero. The aggregate coefficient is a polynomial in W with rational coefficients. It is therefore also nonzero at the actual period by transcendence. This disproves the universal multivariate eta-degree-two implies error-degree-at-least-four claim for these actual polynomial columns, including the sum over target row sets.

`verify_aggregate_coefficient.py` additionally verifies this coefficient without the trace formula: Sylvester reduction produces 44-by-44 determinants; evaluating at 45 field points and exact Newton interpolation gives degree-two coefficient 15 relative to det B_L, hence again 57067. The reduced calculation is cross-checked against the full 156-by-156 determinant at eta=2 and eta=7.

The identities tr X=-6 and tr(X^2)=6 do not make X a negative rank-six projector: its rank is 44 and X^2 is not -X in the finite-field specialization. No such structural claim should be inferred.

## 5. Scope

The nonzero aggregate statement concerns the polynomial coefficient in independent formal errors, at the actual period. It does not evaluate the entire mixed-boundary sum at the genuine error vector, prove domination over other error monomials, or estimate its high-dimensional size. The packet is finite and meets the displayed elementary weights and denominator choices; it does not certify the manuscript's unknown successive interpolation thresholds or an asymptotic counterexample in its critical family.

The separate full lower-grid packing proof in `boundary-blast-audit/BUCKET_BASE_PACKING.md` was also reviewed: the volume shadow estimate, footprint support estimate in Q^A, matroid rank condition, and stated bucket-cost bounds are consistent. The shell-only shortcut has correctly been removed. Retaining the low-total-order rows leaves a genuine global period-rank gap, as that note explicitly acknowledges.
