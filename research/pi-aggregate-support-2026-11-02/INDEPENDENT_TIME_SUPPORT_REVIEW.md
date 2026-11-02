# Independent review of the filtered time support

2026-11-02. The cut inequalities, intrinsic filtration, feedback rank bound and marker-degree bounds are valid with the stated qualifications. This audit checks the proof independently of the finite symbolic examples.

For a nonzero transport matching, every time moves downward or stays fixed. Hence E_r counts the number of matched intervals crossing cut r and is nonnegative. A nonzero period matching requires target times to dominate column valuations in sorted order, which is exactly E_r<=C_r. Summation recovers the total-surplus condition, but cannot recover the individual cuts from that sum.

The intrinsic statement requires no chosen basis: the subspace W intersect (Y-1)^r C[Y,X] has dimension d_r and is annihilated by all low-time rows. Injectivity of a square target evaluation therefore requires at least d_r high-time rows. In a fixed transverse raw-column block, the normalized jets at Y=1 are polynomials binom(h,i) in the distinct raw exponents. Their Vandermonde determinant is nonzero, proving filtration orders 0,...,k-1 even when the exponents have gaps. This argument changes the basis inside the existing span, and correctly makes no untracked arithmetic-normalization claim.

For the feedback bound, the preimage under the invertible source evaluation B of its S_r-dimensional high-time coordinate space contains the d_r-dimensional vanishing subspace. Low-time outside evaluations factor through the quotient, so their restriction has rank at most S_r-d_r. The orientation in the claimed submatrix X_(outside time<r,source time>=r) is correct. Multiplication by R from source times <=r uses only outside times <r and proves the stated bound on N. In particular the time-r diagonal is covered; positive C_r does not imply nilpotence.

For the marker coefficient, a nonzero exterior contribution selects omitted source rows and outside rows. A nonzero minor of R supplies a strictly downward matching between them. Each omitted source row at time s crosses cut s, and each selected outside row at time t crosses cut t+1. Bounding their counts separately by the available rows and the cut capacities gives both degree estimates. The target-period minor imposes the same intrinsic filtration restriction. Cancellation can lower the degree further, so the bounds are not asserted to be equalities in general.

The polynomial-only family in section 6 attains the source-side bound: an outside derivative (t,a) is binom(t+a,a) times the source derivative at time t+a. Thus RX is diagonal with entry (1-zeta)^s-1. The product marker formula and its two different degrees follow exactly. The source-only-beta=0 hypothesis and the finite column budget are essential; this family is not silently extended to fixed-parameter H infinity after transverse rows enter.

The six-row example correctly separates total feasibility from cut feasibility. The chosen forbidden replacement crosses the saturated cut 2, while the surviving replacement crosses cut 1. This is a genuine additional support obstruction, not a numerical near-zero test.

The accompanying verifier uses SymPy for exact symbolic determinants. The final publication runs it and records its output separately. These support and rank results do not provide the unresolved uniform size estimate for critical near-total escape.
