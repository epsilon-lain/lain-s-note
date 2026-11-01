# Independent review of the exterior-transport regrouping

2026-11-01. Reviewed `exterior-boundary-cancellation/EXTERIOR_TRANSPORT.md` directly. The conclusions below concern its algebraic identities, not the unproved high-dimensional period-minor estimate.

1. **Conserved groups are necessary.** The transport preserves center j and total derivative order d=s+|beta|. A full-size minor must select exactly the original number of targets from every such group. Otherwise a block has too few or too many columns and every determinant permutation contains a zero entry.

2. **Common error monomial is correct.** In a group, factoring (-zeta)^a from each column and (-zeta)^(-beta) from each row leaves the integer Pascal evaluation matrix. Every nonzero term therefore has the same multiexponent Gamma=sum a-sum beta. A nonzero permutation proves Gamma>=0 coordinatewise. Zero coordinate shifts should be interpreted by the resulting polynomial identity; zero Pascal minors are simply zero and must not be written as zero times an undefined negative power.

3. **Cauchy--Binet really groups the prior boundary terms.** Internal transport U has determinant one. Thus multiplying by U^{-1} preserves every full-source-size minor of the transport matrix, including those selecting a mixture of inside and outside columns. Grouping by target sets before taking absolute values retains cancellations between assignments which separate row coefficient estimates lose.

4. **The two-dimensional zero example is correct.** With v0=1, v=(3,4), H=8 and d=4, the source lower set is precisely {(0,0),(1,0),(0,1)}. The target nodes {(0,0),(1,1),(2,2)} satisfy the homogeneous cutoff. Their coordinate evaluation rows coincide, so the Pascal determinant is zero. The displayed 2 by 2 boundary remainder has two nonzero products which cancel exactly. This is a genuine coefficient cancellation, not a consequence of a vanishing individual error.

5. **Rounded lower-grid determinant is correct.** The multivariate Newton polynomials have the same leading monomial x^beta/beta! as the binomial basis and involve only gamma<=beta. The lower-set condition makes the basis change unit triangular. Evaluation on a coordinatewise increasing grid is triangular in a common linear extension, with the stated positive diagonal. The integer-dilation specialization and the exact error-exponent formula follow. Target cutoff, weighted column support, and actual period-minor nonvanishing remain separate requirements, as the note states.

6. **The shell warning is valid.** A lower-grid footprint bound does not transfer merely by removing a lower subset. The polynomial x-2 on shell {2,3} has support one, whereas its shell upward shadow has size two. This refutes only that proposed shortcut, not a different packing theorem.

The one-dimensional full-group formula agrees with the independently derived Cauchy boundary formula in `ALTERNATING_BOUNDARY_SUPPORT.md`, Section 6. For consecutive outside targets it gives exactly a unit coefficient multiplying an error power. Consequently, determinant-level grouping can annihilate some multivariate assignments, but cannot be assumed to annihilate every relevant boundary target set.

The remaining unproved step is quantitative control of the sum of compatible period minors with these exact integer coefficients. The accepted algebraic identities do not establish a contradiction against the arithmetic lower bound or bad approximability of pi.

## Added review: the partial-boundary quadratic power

The new one-dimensional partial-replacement law is correct. If the source group is {0,...,b-1} and k target indices lie outside it, the smallest possible target sum retains {0,...,b-k-1} and adds {b,...,b+k-1}. The difference from the source sum is exactly k squared, and the Pascal Vandermonde determinant of that target set is nonzero when the homogeneous cutoff admits it. Thus the exponent is sharp at the transport-coefficient level.

Multiplying the outside columns of [I,R] by a formal variable eta makes the eta^k coefficient divisible by zeta^(k squared), because every relevant Cauchy--Binet target set has exactly k outside columns. This statement concerns the determinant of the corrected period-row expression; the separate raw-column diagonal prefactor b^D is not part of that polynomial-divisibility statement. For prescribed counts k_{j,d}, the exponent adds as sum k_{j,d} squared. If only a total replacement count is specified, the exponent must be minimized over its distributions among groups; the square of the total count cannot be substituted.

A minor cross-reference typo in the rounded-grid section was reported to the author: its integer-dilation specialization should refer to the newly renumbered Newton product equation (8).
