# Independent review of the anisotropic joint column quantile

2026-11-02. Reviewed `JOINT_COLUMN_QUANTILE.md` and `verify_joint_quantile.py`. The cost coefficients, quantile formula, constrained dual, and explicit degree frontier are correct in their stated moment-selection scope. The fast-theta frontier is correctly excluded from the original K>=2 interpolation parameter regime.

## Cost and quantile

Expanding the earlier residual expression gives column coefficients d_i+ln(2R), target coefficients -(ell_i+ln R-ln 2), and source coefficients ell_i-d_i. Thus (1)--(2) retain the same columns and have the correct signs. This is a decomposition of an upper estimate; a lower bound on its positive part is not a determinant lower bound.

The signs of the divided-difference denominators in (6) are correct. Differentiating the displayed F and integrating its exponential transform gives the stated uniform-simplex transform; polynomial interpolation cancels the apparent tails. This proves the formula also when some knots are negative. The stop-loss integral J starts at minus infinity in that signed case. The identity Q=q-J/lambda follows from the truncated first moment and is unchanged by signed costs. Repeated knots must be treated by continuity or confluent divided differences, as stated.

The variational formula (8) has derivative 1-F/lambda, so its maximizer is the lower quantile. The resulting value is concave, **positively** homogeneous, and 1-Lipschitz for the coefficient sup norm. Positive homogeneity is the correct convention here: negative scaling exchanges lower and upper cost tails.

## Source-degree constraints and the closed frontier

Subtracting eta_i times each lower moment constraint yields the dual lower bound in (10). For theta<1, the constant density f=lambda has all spatial means strictly above the bounds while preserving mass and lying strictly between zero and one. Compactness of the finite moment set and separation therefore give the stated strong duality. The lattice approximation argument correctly mixes in this strictly feasible density before rounding, so the necessary finite coordinate inequalities can be enforced with a positive margin. It does not claim full determinant rank.

For equal spatial costs, averaging at fixed time and total spatial degree preserves the objective. Every averaged spatial mean becomes the average of the previous coordinate means, which remains at least the common source bound. Thus isotropic reduction is valid for this constrained moment problem.

In the explicit frontier, the selection a z_0+sum_i z_i<=theta has mass theta^(m+1)/a=lambda and coordinate means theta/(m+2). Its time mean is lambda theta^(-m)/(m+2). The multiplier eta=kappa-A/a is nonnegative under the stated condition; the same selection minimizes the modified cost and makes every source constraint an equality. Primal and dual values therefore coincide exactly at (12).

At half density and theta<1, theta^(m+1)>=1/2 implies 2 theta^m>1. Hence every K>=2 violates K theta^m<1 on this explicit branch. The note correctly treats the formula as a moment diagnostic, rather than an eligible nonzero-minor construction.

The separate exact D_i=0 warning is also correct. Every potentially nonzero determinant permutation must then match alpha=beta coordinatewise. Entries off these equal-type blocks may therefore be set to zero without changing the determinant. Each remaining block has identical rows across the repeated centers at the same time order, so its determinant vanishes when K>=2. An asymptotically zero normalized degree gap does not imply this exact finite condition.

## Verification

The supplied verifier passed on rerun: 52 independent polygon probability/moment checks, 50 higher-dimensional moment identities, four rational quantile enclosures, 36 finite weak-dual checks, 63 exact frontier identities, and 21 K=2 parameter-ineligibility checks. I additionally compared the general distinct-knot CDF and stop-loss expression with the independently integrated special-time formulas in 35 exact rational two-dimensional cases; all agreed.

These checks support the formulas, not a rank assertion for the selected columns. The general sharpness and duality claims rest on the written continuum and lattice arguments.
