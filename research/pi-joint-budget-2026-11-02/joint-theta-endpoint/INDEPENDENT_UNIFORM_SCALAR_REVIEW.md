# Independent review of the uniform scalar obstruction

2026-11-02. Reviewed `UNIFORM_SCALAR_OBSTRUCTION.md` and its verifier independently of their author. The argument is correct for the explicitly stated scalar mass/mean relaxation and Cauchy/LCM normalization. It does not exclude improvements using actual determinant compatibility or other normalizations.

The new column link is exact before any limit. If N_alpha selected columns share a transverse type, their distinct nonnegative time exponents have sum at least N_alpha(N_alpha-1)/2. Summing and applying the established collision theorem to these same multiplicities gives

\[
y\ge I_m(u)/[4(m+1)].
\]

The factor follows from the half-density mass M/(QH)->1/[2(m+1)w_0]; the discarded linear term is w_0/(2H). Thus the selected spatial mean and time cost cannot be adjusted independently. No assertion that the moment optimizer respects every time-exponent capacity is needed for this lower bound.

For s<=u<=a, I independently obtained

\[
P(\sqrt{su})=P(u)(u/s)^{m/2},\qquad
(\kappa_0-1)(u-s)=P(s)(m+2)(u-s)/(K\theta).
\]

These identities give exactly (10). If x>=1/4, completing the square bounds xr^2/K-r below by -K; if x<1/4, r<2. In both ranges the remaining 7K/c term is strictly larger. The base term is the nonnegative square (sqrt(u)-sqrt(s))^2.

For u>=a, the chosen t=a sqrt(theta) lies in the required scalar interval. Splitting the unit arithmetic term gives (u-a)+a(1-sqrt(theta))^2. The remaining analytic, LCM and collision terms have precisely the bracket in (11). For x>=1, its first two positive terms dominate c by the arithmetic-geometric mean inequality. For x<1, p>=exp(-x/2)>1/2 makes its first term alone exceed c. The additional arithmetic contribution omitted by replacing u with a is nonnegative. Boundary cases theta=1 and u=a are included.

The strict positivity is pointwise in each fixed parameter set. Its margin need not be bounded below along a joint sequence. Accordingly the finite-weight error statement only excludes a negative margin bounded away from zero after those errors vanish; it does not silently replace the stated order of limits by a uniform estimate.

The scope distinction concerning targets is essential and is correctly present. A radial scalar optimizer need not fit the coordinatewise degree totals of an arbitrary imbalanced original column set. Conditional averaging on the weighted simplex preserves a scalar selection's time and spatial total means and supplies a symmetric representative; it does not preserve an individual minor or its time grades. The theorem therefore obstructs the scalar envelope, not every coordinate-constrained determinant argument.

I reran the exact verifier: 1,728 low-mean and 1,728 high-mean rational endpoint checks, 317 finite distinct-column subsets, 884 finite moment dual checks, 27 selected-column moment records and 24 joint-theta families all passed. The enclosure 1/3<c<7/20 and affine dependence of each tested score on c correctly justify checking both rational endpoints. The general conclusion follows from the written inequalities, not extrapolation of these finite tests.
