# Independent review of the mean-degree collision relaxation

2026-11-02. Reviewed `MEAN_DEGREE_AND_COLLISIONS.md`, its prerequisites on the same-column residual degree, and `verify_moment_collision.py`. The formulas and the stated fixed-parameter relaxation obstruction check. The concentration assertion uses the bounded normalized arithmetic/analytic budget explicitly assumed in Section 5; it is not unconditional for arbitrary B,C.

## Variational bound and sharpness

The low-mean optimizer is the normalized positive part of an affine function, with cutoff A=(m+2)t/m. Direct integration gives mass one, mean t, and the displayed squared norm. The square-completion dual inequality applies to every nonnegative profile, without compactness assumptions. When the affine optimizer has full support, orthogonal projection onto span{1,x} gives the middle branch; positivity holds exactly at t>=m/(m+2). Above the unweighted mean m/(m+1), the constraint is a mean upper bound and the constant profile is optimal. The note correctly distinguishes this from prescribing an exact mean above that threshold.

Sampling, normalizing, and rounding changes each multiplicity by at most one. At fixed parameters there are Q_H types and each real multiplicity is O(H), so the square-cost error is O(Q_H H); division by M_H H makes it O(1/H). Mean errors also tend to zero. A strict mean margin followed by a diagonal limit establishes one-sided feasibility. This proves sharpness for the stated mass/mean relaxation, without constructing actual row labels or nonzero minors.

The stronger witness in Section 6 also checks: under the change of variables y_i=w_i a_i/H, both the target profile and cheapest-transverse-column selection are radial in sum_i y_i. Weighted lattice sums therefore converge to coordinate-symmetric integrals even when the original weights differ. Each coordinate receives t/m in the target and u_m/m in the columns; the source receives s_m/m. For fixed sufficiently large m, the strict inequalities s_m<t<u_m dominate all lattice, tie-selection, and rounding errors. Thus Gamma_i and E_i are eventually nonnegative for every coordinate, using the same selected columns. This does not supply compatible time grades or determinant rank, which remain explicitly omitted.

As an additional independent finite check, `verify_independent_relaxed_witness.py` constructs an integer target profile of total mass 11,904 in the two-dimensional setup with H=60 and selects the same number of cheapest transverse columns from 20,336 raw columns. It finds

\[
\Gamma=(54891,27322),\qquad E=(24092,11992),\qquad D=\Gamma+E,
\]

with source mean 7/62, target mean 38035/142848, and selected-column mean 79417/238080. All three means have the required strict ordering. This verifies the coordinatewise mechanism in a finite instance; it is not an actual exterior term or a finite verification of the large-m critical-family limit.

## Parameter obstruction and its scope

For fixed theta<B, the collision factor has its exponential transition at B. Uniform divergence below B-epsilon and a feasible low-mean witness above B+epsilon give the claimed supremum limit 1+theta-2B. When B<=theta, the corresponding limit is 1-theta. The witness shows that the retained relaxation really permits positive values; an unattained estimate alone is not being mistaken for an impossibility theorem.

The selected-column time-degree estimates have the correct directions. Their upper bound makes the particular Cauchy growth term vanish when CB=1, while their lower bound forces that term to diverge when CB>1. This yields CB<=1 only for suppressing that term in this analytic estimate. The LCM lower bound is valid: each prime-adic carry contributes at most one to the central binomial valuation, so the central binomial coefficient divides the LCM. The consequent lower bound on kappa forces v_0 to grow in the stated large-weight normalization, hence C theta>B along the fixed exponential family.

Combining these two restrictions gives B^2<theta. Optimizing the relaxation then leaves the positive infimum (1-sqrt(theta))^2 for every fixed 0<theta<1. This argument does not prove a lower bound for a determinant, forbid cancellation, or rule out a different denominator normalization, analytic estimate, or a joint theta_m-to-one regime. These restrictions are correctly preserved in the note.

## Verification

The supplied verifier passed on rerun: 72 exact continuum integral checks, 27 finite optimization problems with 565 feasible integer profiles, 36 packet-normalization records, four rational concentration thresholds, 199 LCM divisibility checks, and 23,118 admissible parameter triples. The independent rounded-profile verifier also passed. The universal lattice and parameter limits rely on the written proofs, not those finite tests.
