# Independent review: positive-fraction law and finite bulk compatibility

2026-11-01. Reviewed against the proofs and exact scripts available on this date. No conclusion about bad approximability of pi is inferred. The two results below have different scopes: one is asymptotic at the transport level, while the other is a complete finite compatible-rank certificate.

## 1. Positive-fraction transport law

Reviewed `POSITIVE_FRACTION_LAW.md`, including its nine displayed claims. No mathematical error was found under its stated fixed-parameter hypotheses.

The exact minimization is valid: for an unrestricted cost bucket with b points, replacing k of them minimizes displacement by removing its k most expensive points and inserting the k cheapest outside points. Thus its expression in terms of the first-n cost sums L(n) is exact, including ties.

The lattice estimates used uniformly in k are sufficient. Unit cubes sandwich the cost quantiles within an additive constant depending only on the fixed weights. Consequently L(n) has error O(n+1), and substituting b+k and b-k gives error O((h+1)^m). The endpoint k=b causes no singularity in this argument. The inequality F_m(rho)>=rho^2/m follows from the stated second derivative and remains valid at rho=1 by continuity.

The total-order cap is handled in the correct direction. A target outside the capped source bucket, but satisfying the same total-order cap, necessarily lies outside the unrestricted cost bucket. The k removed source points form a subset of the unrestricted inside points. Enlarging both choices gives the unrestricted minimum as a lower bound, so the capped comparison does not silently replace the source count by an unrestricted count.

The subsequent Cauchy--Schwarz constants check independently. With V=1/(m! product c_i),

\[
\sum_{j,N} h_N^{m-1}\sim \frac{KH^m}{mv_0},\qquad
M_H\sim\frac{KH^{m+1}}{(m+1)!v_0\prod_i v_i}.
\]

They yield exactly

\[
\liminf\frac{c\cdot\Gamma}{M_HH}
\ge \frac{\eta^2}{(m+1)\prod_i(v_i/c_i)}.
\]

The summed error is O(H^(m+1)), hence lower order than M_H H. This relies on fixed dimension and fixed weights before the H limit. The cost lattice bounds positive h_N away from zero.

The determinant bound follows from the full Pascal-column identity
sum_beta binom(a,beta)^2=product_i binom(2a_i,a_i)<=4^|a|. Restricting the rows can only decrease that column norm. The type cap w a<=H therefore gives the claimed log(2)/w_* term. The zeta_i are explicitly fixed before H tends to infinity; without such a hypothesis, multiplying an asymptotic error by an unbounded H-dependent quality parameter would require a new estimate.

For the construction, transverse fibers admit a common coordinatewise linear-extension ordering. The multivariate Pascal matrix is block triangular in that ordering, with the positive one-dimensional Vandermonde quotients on its diagonal. Thus the claimed nonzero product is justified in every fixed dimension. Each fiber's displacement is exactly c_1 k_y^2. Integrating n_y^2 gives the factor 2/((m+1)c_1), and summing over the selected complete-packet buckets gives the two constants in equation (9).

The selected-bucket cutoff guarantees the unrestricted source bucket already fits the total-order cap. The enlarged targets satisfy both |a|<=N and, for sufficiently large H, w a<=H. In particular, the displayed type-cap leading ratio is strictly below one even at theta(1+rho)=1, since rho>0 and v_0>0.

The companion script was rerun: 175 exhaustive finite target-set tests, 12 rational Pascal determinant tests, and 12 complete-packet support/count cases passed. These checks support the finite formulas; the asymptotic claims rest on the written argument above.

**Scope retained:** the nonzero Pascal coefficient does not supply a nonzero actual-period minor, a lower bound on that minor, or a common low-cost rank certificate. The subsequent large-weight limit and then dimension comparison must not be interchanged silently with the fixed-parameter H limit.

## 2. Finite bulk-compatible transport

Reviewed `verify_bulk_transport.py` using `verify_bulk_review.py`. The independent, self-contained verifier uses polynomial convolution for entries and reduced-row elimination for pivots, adapting the independent algorithm in `../pi-exterior-blast-2026-11-01/exterior/verify_full_packet_review.py`. It does not import any previous script or the probe's double Taylor-allocation entry formula.

The complete finite packet has K=2, m=2, H=20, theta=1/2, v_0=1, v=(4,10), and w=(2,5). Its 156 source rows satisfy

\[
s+4\beta_1+10\beta_2<20,
\]

and all 242 available raw columns satisfy h+2 alpha_1+5 alpha_2<=20. For N=8,...,19, each center's entire source bucket is transported by

\[
f_j(\beta)=(2\beta_1+j,\beta_2),\qquad
s'=N-|f_j(\beta)|.
\]

The remaining buckets use the identity. Every chosen target obeys s'>=0 and 2a_1+5a_2<=20. The parity of a_1 separates the two centers' transverse target sets inside each packed bucket. The low-order identity buckets remain in the 156-row test.

The exact integer affine-grid determinant in each packed bucket is 2^(sum_beta beta_1). Twenty-four independent rational block determinant checks give the complete transport coefficient

\[
\boxed{2^{60}(-\zeta_1)^{93}.}
\]

There are 66 rows in packed buckets, 51 changed rows, and 33 targets outside the original source packet. These three counts measure different things and must not be identified.

The independently computed common-column certificates are:

| Quantity | Exact result |
|---|---:|
| Prime | 65537 |
| Image of i | 256 |
| Target rank at formal period Z=1 | 156 |
| Factorial-scaled target determinant modulo the prime | 9333 |
| Factorial-scaled rational-source determinant modulo the prime | 44706 |

Both residues use the same 156 raw columns. Their complete index list and raw-column descriptions are in `bulk_review_verification.json`. The rational source centers are 44i/7 and 710i/113. Modular reduction takes place in the localized coefficient ring Z[i,1/7,1/113], where the denominators are invertible modulo 65537. It is not a homomorphism on all of Q(i).

The period-side entries, after replacing the true exponential period factor by its exact value one, are polynomials in the formal period Z. Evaluation at Z=1 is a polynomial specialization used solely to prove that the determinant polynomial is nonzero. It does not assert e^1=1 or reinterpret Z=1 as an exponential period. Known transcendence of 2 pi i then certifies that this particular polynomial is nonzero at the actual period. Row factors s! are all nonzero and invertible at this prime, so normalization does not affect rank or nonvanishing.

No general interpolation threshold is used in this finite computation: full rank is verified directly on the stated 242-column space. The test therefore does not establish that the manuscript's separated-weight interpolation hypotheses hold in a growing parameter family.

**Scope retained:** this proves one complete finite packet has a nonzero grouped transport coefficient and compatible nonzero source and actual-period minors on common raw columns. It provides neither an asymptotic family of such packets nor a determinant magnitude estimate. In the auxiliary expansion where outside rows carry eta, this individual target-set contribution has degree eta^33 zeta_1^93. Other complete target sets may contribute to the same coefficient and cancel it. This result is not a certificate that the aggregate coefficient of eta^33 zeta_1^93 is nonzero.

## 3. Relation to the boundary-rank filtration

The positive-fraction law and this folder's rank-filtration result are compatible. Fixed total error degree sees only a thin boundary strip, but a fixed positive fraction of replaced rows incurs an H-scaled total error degree. The thin-strip reduction cannot be extrapolated to say that the whole error operator has negligible rank. Conversely, a large coefficient-transport rank need not survive multiplication by the compatible period matrix.

The missing global step remains a quantitative estimate on actual compatible period minors, or an exact cancellation statement for the full sum over target sets. The new finite bulk example narrows what a universal rank-obstruction claim could say, but does not settle that asymptotic step.

## 4. Review of the proposed exterior-sum envelope

Reviewed the proposed `EXTERIOR_SUM_ENVELOPE` statement communicated by the lead author. It is valid under the same hypotheses as the positive-fraction estimate: m>=1, fixed K and positive rational weights, c_i=v_i-v_0>0, admissible H on a common cost lattice, and fixed errors 0<|zeta_i|<1. The determinant and its minors must use the same normalized shear-transport convention.

Every nonzero transport target preserves N=s+|beta|. Since c_i>0, a source row has N<H/v_0. After removing target types that are identically zero on the weighted column space, the finite ambient row set can therefore be taken to be

\[
\Omega_H=\{(j,t,a):0\le j<K,\quad t+|a|<H/v_0,\quad wa\le H\}.
\]

Its cardinal D_H satisfies

\[
D_H\le K\lceil H/v_0\rceil\prod_i(\lfloor H/w_i\rfloor+1)
=O(H^{m+1}).
\]

Thus the number of M-element target sets is at most 2^D_H, whose logarithm is o(MH). This count is for the exact grouped exterior expansion. It need not count individual row permutations, which have already been included in each Pascal determinant.

Here is an explicit justification of uniformity in the target set. Put

\[
V=\frac1{m!\prod_i c_i},\qquad
S_H=\sum_{0\le j<K,\ 0\le N<H/v_0}(H-v_0N)^{m-1}.
\]

The finite uniform bucket bound and Cauchy--Schwarz give, for every admissible target set J,

\[
c\cdot\Gamma_J\ge
\frac{k_J^2}{mVS_H}-C H^{m+1},\qquad
k_J=|J\setminus L|,
\]

with one constant C independent of J. Consequently

\[
\frac{\log|\det T_{L,J}|}{MH}
\le\frac{\log2}{w_*}
-\tau\frac{M}{mVS_HH}\left(\frac{k_J}{M}\right)^2
+O(H^{-1}).
\]

The coefficient M/(mVS_HH) tends to 1/((m+1) A), where A=product_i(v_i/c_i). Since 0<=k_J/M<=1, replacing that coefficient by its limit creates one o(1) error uniform in J. The fixed-zeta hypothesis is used here to keep tau fixed.

For any compatible D_H-by-M matrix B, Cauchy--Binet and the triangle inequality therefore imply

\[
\frac{\log|\det(TB)|}{MH}
\le\frac{\log2}{w_*}
+\max_J\left\{
\frac{\log|\det B_J|}{MH}
-\frac{\tau}{(m+1)A}\left(\frac{|J\setminus L|}{M}\right)^2
\right\}+o(1).
\]

The maximum may be restricted to sets with det(T_{L,J}) det(B_J) nonzero; if there are none, both the determinant and the sum are zero, and the logarithmic assertion is interpreted using log(0)=-infinity. No size estimate or asymptotic regularity of B is needed for this upper-bound implication, because the transport remainder and set-count error are uniform independently of B.

The exact identity contribution J=L can be split from the sum when L is contained in Omega_H; the standing theta<=1/2 hypothesis ensures this (theta<=1 is already enough for the column-type inclusion). A bound on the sum excluding J=L is not automatically a bound on the full determinant unless the source contribution is also controlled.

**Scope retained:** this is a sufficient upper envelope obtained by discarding possible cancellation. It is not a lower bound, not a proof that a maximizing term survives the sum, and not an estimate of the unknown period-minor maximum. No issue was found in the proposed statement once these hypotheses and conventions are explicit.
