# Independent review of diagonal profiles and the weighted-column kernel

2026-11-01. The completed `DIAGONAL_PROFILES_AND_COLUMN_KERNEL.md` and its standard-library verifier were reviewed independently. No error was found in the coordinate identities, fixed-type capacity bound, or kernel comparison. No endpoint determinant estimate follows merely from these identities.

## Column and derivative filtrations

For the least-weight anchor, both substitutions X_i=T+U_i and U_i=X_i-T preserve the column budget, since w_1<=w_i. The new monomials have the same indexing set as the old weighted monomials. Expansion is triangular in the total exponent of coordinates 2,...,m, with diagonal coefficient one, which proves the asserted integer determinant-one basis change.

Extracting the old transverse coefficient from

\[
(jZ+z+u_1)^a\prod_{i=2}^m(u_i-u_1)^{\nu_i}
\]

gives exactly the sign, binomial factors and tangent order ell=|beta|-|nu| in equation (3). It generally mixes multiple normal profiles; treating these profiles as independent choices for each original row would be invalid.

For the greatest-weight anchor, partial_Xm=partial_T-sum_i<m partial_Ui. Expanding its beta_m-th power and dividing by beta! gives the binomial product in equation (4), including all factorial factors. Each replacement moves derivative weight from the greatest-weight coordinate toward coordinates of no greater weight. The inverse expansion of partial_T=sum_i partial_Xi has the same budget property. Both matrices are triangular with diagonal one when ordered by the greatest-coordinate derivative. Thus the complete source derivative space is preserved. The corresponding column change need not preserve the assigned column weights, exactly as stated in the note. Neither claim establishes impossibility for every alternative adapted coordinate system.

## Fixed-type differential capacity

For a fixed original type beta, every permitted column has remaining degree d=|alpha-beta| with w_* d<=w dot (alpha-beta). Conversely, placing all d units in one least-weight coordinate realizes each permitted pair (h,d). This proves equality with the stated exponential-polynomial space, not only inclusion.

Its distinct h-blocks are independent. On h>0, differentiation acts as h+D on the finite polynomial factor and is invertible; on h=0, D^r kills exactly min(r,d_0+1) polynomial dimensions. Hence the displayed upper bound for all time orders at least r is correct. The image is itself an exponential-polynomial solution space of the reduced constant-coefficient differential equation, so consecutive initial derivatives give the stated abstract sharpness. The final note correctly qualifies that these consecutive time orders need not fit the separate conserved-order budget of a transport problem.

For a formal period Z, the entries in equation (1) omit the center exponential e^(hjZ). This does not invalidate the rank or recurrence statements: blockwise phase correction and translation commute with differentiation and preserve its polynomial kernel. At the actual period the phase correction equals one. The final note explicitly retains this distinction.

The example with image span{1,e^z} has capacity one after one derivative, explaining the earlier dependent type-2 rows exactly. Away from the weighted-column boundary the profile dimension grows quadratically in H, while the relevant count of rows per type grows only linearly. Thus this constraint alone does not settle the asymptotic high-escape sector.

## Genuine weighted-column kernel

After clearing rational weight denominators, each raw monomial of cost c contributes t^c/(1-t) to the cumulative generating series. Summing the independent geometric series in h and alpha proves equation (8).

For fixed 0<t<1, strict convergence at all centers gives a neighborhood where all finite jet extractions are legitimate. The infinite weighted Gram matrix assigns a column of cost c the positive weight t^c. Multiplication by t^(-H) assigns every permitted column c<=H a weight at least one; columns above H add positive semidefinite terms. This proves the Loewner comparison in equation (9), independent of numerical conditioning. Cauchy--Binet bounds any selected raw-column minor by the full finite Gram determinant, and positive-semidefinite determinant monotonicity supplies the final inequality.

The normalized determinant charge -ln(t)/2 is necessary and is included. This coupled product kernel keeps the common monomial budget and coefficient relations, but its determinant still needs a new uniform estimate. It cannot be replaced by a direct sum of unconstrained type-function spaces without losing precisely the compatibility that this construction retains.

## Verification rerun

`python verify_diagonal_profiles.py` passed independently: 10,980 transformed-entry comparisons, three integer column basis changes, six derivative-space changes, twelve sharp capacity cases, twenty-six generating-function coefficients, and fifteen exact positive principal minors of a Loewner difference. The finite kernel instance uses a rational formal center and does not assert that it is an exponential period. The infinite determinant comparison rests on the positive-semidefinite proof above, not on extrapolating that finite instance.
