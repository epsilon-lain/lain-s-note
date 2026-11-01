# A complete finite bulk transport with compatible nonzero period and rational minors

2026-11-01. This is a finite exact certificate, with all low-total-order rows retained. It establishes neither an asymptotic family nor a quantitative period-minor bound. It complements the general positive-fraction transport law by checking that bulk transport and actual period independence can coexist in at least one complete packet.

## 1. The complete source and column spaces

Take m=2, K=2, H=20, theta=1/2, v_0=w_0=1 and w=(2,5). Thus v=(4,10) and c=(3,9). The source rows and raw columns are

\[
L=\{(j,s,\beta):j=0,1,\ s+4\beta_1+10\beta_2<20\},
\qquad
\mathcal C=\{(h,\alpha):h+2\alpha_1+5\alpha_2\le20\}.
\]

Their sizes are 156 and 242. All indices except j are nonnegative integers. Row ordering and selected columns are recorded explicitly in `bulk_transport_verification.json`.

For a formal common period Z define ordinary Taylor-coefficient rows

\[
B^Z_{j,s,\beta;h,\alpha}
=[z^su^\beta]e^{hz}
(jZ+z+u_1)^{\alpha_1}(jZ+z+u_2)^{\alpha_2}.
\]

The factorial-scaled matrix has entries s! times this expression, and belongs to Z[Z]. At the actual period Z=omega=2*pi*i these are the genuine local rows of exponential-polynomial functions; the omitted center exponential is one because h and j are integers. A rational-center source matrix A uses r_1=44i/7 and r_2=710i/113 in the two polynomial factors in place of Z. These are twice i times 22/7 and 355/113; the certificate claims no new approximation to pi.

## 2. Transport a full collection of buckets

Write N=s+|beta|. Within each center and each bucket with N=8,...,19, send

\[
\beta\longmapsto a=(2\beta_1+j,\beta_2),\qquad
(j,s,\beta)\longmapsto(j,N-|a|,a).
\]

For N=0,...,7 retain the original row. The algorithm selects these buckets by requiring every image at both centers to obey |a|<=N and 2a_1+5a_2<=20. Thus no row is discarded and no time order is negative. The output J contains 156 distinct row labels. Precisely 66 rows belong to the transported buckets, 51 labels change, and 33 rows of J lie outside L. A changed label need not lie outside the source packet; these counts should not be conflated.

In every transported bucket the source type set is lower closed. With source and image lists in corresponding order, its Pascal minor is

\[
\det\left[\binom a\beta\right]_{
  \beta\in B_N,\ a=(2\beta'_1+j,\beta'_2)}
=2^{\sum_{\beta\in B_N}\beta_1}.
\]

For example, this follows by triangular change from binomial polynomials to coordinate Newton polynomials: doubling the first coordinate multiplies each leading monomial by 2^(beta_1); translation by j contributes only lower terms. The untouched blocks have determinant one. All 24 transported block determinants are also checked exactly over Q.

Their product and total exponent increment are

\[
\prod_{j,N}\det P_{j,N}=2^{60},\qquad
\Gamma=\sum_{J}a-\sum_L\beta=(93,0).
\]

Consequently the complete grouped transport coefficient, for ordinary Taylor rows, is exactly

\[
\boxed{\det T_{L,J}=2^{60}(-\zeta_1)^{93}.}
\]

The source ordering and its paired image ordering yield the displayed sign convention. If a different canonical target order is chosen, the transport minor and the period minor acquire matching permutation signs, so their product is unchanged.

## 3. One common raw-column family works on both sides

Work modulo p=65537, with i mapped to 256; 256^2=-1 mod p. At formal period Z=1, lexicographic pivot selection chooses 156 of the 242 raw columns. On exactly these same columns the factorial-scaled determinants are

| Matrix | Determinant modulo 65537 |
|---|---:|
| Period target B^1_J | 9333 |
| Rational source A^r_L | 44706 |

The normalized ordinary-Taylor residues are 50847 and 41347 respectively. All row factorials are units modulo p. The two programs obtain the same selected columns and residues using different entry constructions (direct Taylor allocation versus polynomial convolution) and different pivot implementations.

Here is the exact characteristic-zero implication. The scaled period determinant is a polynomial P(Z) in Z[Z]. Its nonzero residue at Z=1 proves P is not the zero polynomial. Since 2*pi*i is transcendental, P(2*pi*i) is nonzero. For the rational source, reduction is a homomorphism from the localization Z[i,1/7,1/113] to F_p, with i mapped to 256. Nonzero reduction therefore proves the characteristic-zero determinant is nonzero. This is not a map from the whole field Q(i) to a finite field.

The transport coefficient in the factorial-scaled convention is multiplied by

\[
\frac{\prod_{(j,s,\beta)\in L}s!}
     {\prod_{(j,t,a)\in J}t!}.
\]

The exact integer ratio is stored in both JSON certificates. It must not be silently dropped when comparing scaled and unscaled matrices.

## 4. What this certifies, and what it does not

The coefficient and its accompanying actual-period minor are both nonzero, on a raw-column family where the rational source is also nonzero. This closes the finite gluing issue for this particular bulk example. It uses known transcendence only to transfer a polynomial nonzero certificate; it does not invoke or independently verify the uploaded manuscript's high-dimensional interpolation theorem.

If a formal marker eta counts outside rows, this target set contributes at eta^33*zeta_1^93 in the exterior expansion after the inside triangular factor is removed. This is one target-set contribution. Other target sets can contribute at the same degree and cancel it. No claim that the complete aggregate coefficient is nonzero is made here, and no lower bound for the magnitude at the actual period follows from a modular residue.

The open step is a family of such compatible selections in the required growing-parameter regime, together with quantitative bounds for their complete contribution. Neither a single finite example nor a positive-fraction transport coefficient alone supplies that step.

## Reproduce

Run `python verify_bulk_transport.py` and `python verify_bulk_review.py` from this directory. Both use only the Python standard library and write their JSON certificates next to the scripts. The review file is independently implemented rather than importing the originating checker.
