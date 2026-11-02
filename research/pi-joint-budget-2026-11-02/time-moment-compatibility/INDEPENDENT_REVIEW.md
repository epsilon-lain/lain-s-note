# Independent review of the phase-cycled dilation

2026-11-02. Reviewer: near-total-cost-coupling agent. Independently reviewed `PHASED_DILATION_COMPATIBILITY.md`, including the near-total-escape paragraph, and inspected and reran `verify_phased_dilation.py`.

**Result: passed; no mathematical correction required.**

1. The low-order threshold guarantees nonnegative target time. Above it, the coordinate dilation is increasing and injective within each conserved bucket; below it the identity map applies. Distinct bucket and center labels then guarantee distinct target rows. The global transverse budget follows from A<1 and a fixed additive phase error.

2. The Newton determinant proof is valid on every source lower set. Each Newton polynomial has unit leading coefficient in the normalized binomial basis, its lower terms stay inside the lower set, and evaluation vanishes unless the node index dominates the row index coordinatewise. The diagonal product is positive and phases cancel. Thus the nonzero transport determinant is proved, not inferred from matching alone. The harmonic-sum estimate for its extra coefficient is also correct.

3. For arbitrary raw subsets, the intrinsic filtration formula uses the rank of consecutive Taylor tests on distinct h values, so downward closure is unnecessary. The target low-time count is bounded by r unchanged source times plus an interval of length r+S for changed source times. Condition (7) therefore proves all r<=T_H cuts simultaneously; larger cuts are trivial. This remains a necessary-support statement, not actual-period nonvanishing.

4. The exact phase count is Q per coordinate: the condition that a has a preimage with a phase in 0,...,L-1 is equivalent to an interval of length Q for the integer preimage. Away from the fixed coordinate boundary it therefore contains exactly Q integers. Each full phase vector occupies one residue class modulo L^m. This proves the stated uniform cone density and its square cost with all parameters fixed before H.

5. The low-bucket correction affects at most K N_0 times the number of source types. Every type occupancy is bounded by K ceil(H/v_0), for both the valid and auxiliary maps. Hence the stated l1 and square-cost errors, including the factor 4K/v_0, are justified. Taking transverse weights large at each fixed dimension controls these errors without interchanging the phase-period and height limits.

6. The escape calculation uses the correct coefficients lambda_i=lambda-(lambda-1)v_0/v_i>1. Its inner simplex is contained in the source simplex and has relative volume product(lambda_i^(-1)). The low-order correction and auxiliary negative times are both covered by the same low-bucket mass error. After separation this gives the claimed outside fraction 1-lambda^(-m).

7. In the fixed-exponential critical regime, the rational dilation numerator can be chosen O(m), so the sufficient time-cut factor grows only polynomially while w_0/v_0 grows exponentially when C theta<1. The construction therefore realizes the cone mean and square cost, coordinatewise matching and a nonzero transport minor on the stated order of limits. Actual-period rank on the selected columns, determinant size, and complete-sum cancellation remain unproved, as the note states.

The main verifier passed with 77,360 distinct source/target rows, 101 nontrivial time cuts, 217 exact bucket Pascal matrices, and positive residual coordinate degrees on the selected raw columns. An additional independent enumeration checked 1,680 one-coordinate phase counts for all integers 1<=Q<L<=9; every count was exactly Q.

No conclusion about pi or an actual nonzero period minor follows from this review.
