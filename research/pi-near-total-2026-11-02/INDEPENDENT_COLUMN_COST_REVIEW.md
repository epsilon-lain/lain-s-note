# Independent review of the raw-column and residual-degree budget

2026-11-02. Reviewed `COLUMN_COST_AND_RESIDUAL_DEGREE.md` and `verify_column_cost.py` independently of their author. The stated bounds are correct within their explicit raw-column, fixed-parameter scope. No critical near-total determinant margin or approximation property of pi follows from them.

## Quantile and arithmetic cost

Integrating out the time exponent in the weighted column simplex gives the transverse-cost density

\[
f_m(x)=m(m+1)x^{m-1}(1-x).
\]

Its mass and first moment up to A are respectively

\[
(m+1)A^m-mA^{m+1},\qquad
mA^{m+1}-\frac{m(m+1)}{m+2}A^{m+2}.
\]

Their quotient is exactly the note's q_m(lambda). The conditional mean is increasing in A, since its derivative is f_m(A)(A-q_m)/F_m(A)>0 for 0<A<1. Thus the rational enclosing intervals used by the verifier also enclose the conditional mean in the stated direction. The limiting density is continuous; a cutoff shell and lattice-count errors are lower order. Selecting the cheapest M distinct columns gives the uniform lower bound for every raw-column selection, independently of any rank hypothesis.

The source simplex has mean weighted transverse cost m theta/(m+2). A nonzero source determinant has at least one nonzero permutation, and every matched entry requires alpha_i>=beta_i. This separately supplies D_i>=0 and justifies the maximum with zero in (8). It also justifies (10): the arithmetic-normalization error is bounded by delta times w dot D/(MH), which is at most delta.

For Y=m(1-x), the limiting density is y exp(-y). Selecting the smallest x means retaining the tail Y>=y, not its lower tail. The tail mass is exp(-y)(1+y), and its conditional mean is

\[
\frac{y^2+2y+2}{y+1}=y+1+\frac1{y+1}.
\]

Together with the bounded second moments, this proves the direction and coefficient in (9a). The stated theta=1/2 correction follows by subtracting m/[2(m+2)]=1/2-1/m+O(m^{-2}). The order of limits H first, then m, is essential and is correctly stated.

## Residual degree and the collision estimate

For a nonzero transport minor, a nonzero determinant permutation gives Gamma_i>=0. For a nonzero period minor, the corresponding matching gives E_i>=0. The coordinate identity D=Gamma+E is exact on the same row and column sets; neither nonnegativity condition is sufficient for rank.

Expand each type function in powers of z/R. Cauchy's estimate bounds each coefficient in type a and column (h,alpha) by

\[
e^{Rh}2^{|alpha|}R^{|alpha|-|a|}.
\]

In each coefficient determinant, factor R^{-|a|} from every row and e^{Rh}(2R)^{|alpha|} from every column. Every remaining entry has modulus at most one, so its determinant is bounded by M!. The residual power is exactly R^{sum_i E_i}; it is independent of the Taylor-degree labels being summed. The ordinary-jet evaluation bound contributes r^d for each Taylor degree d. Within a repeated type, distinct degrees are necessary for a nonzero coefficient determinant, and their sum is

\[
\frac{n!r^{n(n-1)/2}}{\prod_{j=1}^n(1-r^j)}.
\]

Multiplying over types, using product n_a!<=M!, and bounding the denominator below by (1-r)^M produces precisely the remainder 2 ln(M!)+c_R M-M ln(1-r) in (12). It is o(MH) at fixed parameters. This derivation applies to sparse ordinary time rows and does not import an extra logarithmic-jet factorial or time-extraction cost.

After factoring the common error monomial from a grouped transport determinant, each Pascal column has Euclidean norm at most its nonnegative full column sum 2^{|a|}. This proves (13). Combining the two estimates and using D=Gamma+E gives (14), including the signs of both the residual and transport terms. The text correctly describes (14) as a sufficient upper estimate, not a lower bound or evidence that any one term dominates.

## Verification and limits

I reran the supplied verifier successfully: 317 exhaustive column subsets at 20 cardinalities, 9 finite packet counts, 1120 target/column degree checks and exact residual factorizations, 66 jointly nonzero terms, and 3 rational half-density bounds. I additionally integrated the mass and first-moment polynomials independently using rational arithmetic for 32 pairs (m,A), m=1,...,8, and checked the formula for E[Y^2] in those eight dimensions. All checks passed.

The finite packet records illustrate the limit but do not themselves certify an asymptotic rate. The verifier's formal period specializations certify polynomial coefficient identities, not exponential periodicity at rational surrogate centers. Both limitations are respected in the note. The arithmetic conclusion also remains conditional on the displayed normalization ratios d_i/w_i and says nothing stronger about the actual value of a determinant than the stated arithmetic bound.
