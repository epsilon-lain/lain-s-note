# Independent review of active spans and normal profiles

2026-11-02. Reviewer: near-total-cost-coupling agent. Reviewed `ACTIVE_SPAN_AND_PROFILE_CAPACITY.md` and `EXACT_BIVARIATE_BLOCKS.md` independently.

**Proof review: passed; no substantive mathematical correction required.**

- Equal transverse weights preserve every fixed-time, fixed-total-degree derivative space under the adapted linear coordinates. Normal coordinates vanish at all centers and along the time direction, giving the exact Kronecker matching of normal profiles. The resulting rank sum is exact; replacing block ranks by their row/column minima is correctly stated only as an upper bound.
- The aggregate capacity integral uses the correct row-to-column factor A_d. The visible-column normal degree distribution is Beta(m-d,d+2), including the unused budget coordinate. Both critical rank-fraction formulas and the m=10 values check exactly.
- The full source matrix at displaced centers is also covered when its active span is a plane. Therefore the equal-weight rank loss cannot be interpreted as a period-only collapse with an automatically available nonzero arithmetic source minor.
- For unequal weights, the increasing greedy basis generates exactly the filtered column restriction space. In the decreasing greedy basis, a dual direction annihilates every coordinate of weight above its pivot weight. Products of these directional derivatives belong to the original source packet, and their label combinations are independent. This proves the stated necessary quotient inequality; it does not prove sufficient interpolation.
- The raw/adapted Gram distinction is correctly retained. The metric Q=C^(-1)C^(-*) generally couples profiles even when the basis change has determinant one.
- The appendix's condition v_0>=v is precisely what makes its time/transverse shear preserve the anisotropic packet. Under that hypothesis the direct Hermite decomposition proves the exact block-rank formula. The 15-by-15 counterexample outside that hypothesis has the displayed polynomial kernel and a nonzero 14-by-14 minor, so the warning against automatic maximal block rank is justified.

The bivariate verifier was rerun successfully: 36 exact models, 14,592 coordinate-identity entries, rank 14 with minor 1/640 in the counterexample, and eight invisible-polynomial cases.

The main verifier was also inspected and rerun successfully: 5,206 coordinate entries across five exact profile models, four adapted Gram products, 138 rational beta identities, 15 finite profile counts, five greedy-basis models with 70 filtration checks, and 15 exact quotient-image ranks. Its rational center specializations are correctly treated as formal matrix tests rather than exponential periods.

The results are rank obstructions and necessary weight-separation conditions. No conclusion about pi or sufficient interpolation follows.
