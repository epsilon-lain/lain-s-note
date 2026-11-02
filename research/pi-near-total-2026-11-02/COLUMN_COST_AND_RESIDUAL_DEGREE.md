# Unavoidable raw-column cost and an exact residual-degree budget

2026-11-02. This note couples the scalar transport cost to the degrees of the same selected raw columns. It proves that a cheaper raw-column choice cannot remove the leading arithmetic cost in the critical half-density family. It also retains the residual degrees in the analytic estimate instead of assigning every column the full budget H. No conclusion about pi's bad approximability follows.

## 1. Setup

Fix m>=1, K>=1, positive rational weights, and 0<theta<=1. Source rows are

\[
L_H=\{(j,s,\beta):0\le j<K,\quad
v_0s+\sum_i v_i\beta_i<H\},\qquad v_i=w_i/\theta.
\]

Let M_H be their number. Raw columns are distinct monomials Y^h X^alpha with

\[
w_0h+w\cdot\alpha\le H.
\]

Let N_H count all such columns and assume the limiting row-to-column ratio

\[
\lambda=\lim\frac{M_H}{N_H}
=\frac{K\theta^m w_0}{v_0}
\tag{1}
\]

lies in (0,1]. For lambda=1, the discussion applies along heights for which M_H<=N_H. All parameters are fixed before H tends to infinity through the admissible weight lattice. We select M_H raw columns, denoted C_H.

For an exterior target set J of M_H rows, define coordinate totals

\[
D_i=\sum_{C_H}\alpha_i-\sum_{L_H}\beta_i,
\qquad \Gamma_i=\sum_J a_i-\sum_{L_H}\beta_i,
\qquad E_i=\sum_{C_H}\alpha_i-\sum_Ja_i.
\tag{2}
\]

Thus D=Gamma+E exactly. A nonzero Pascal transport determinant implies Gamma_i>=0. A nonzero period determinant on these columns implies E_i>=0, because every nonzero determinant permutation pairs each target a with a column alpha>=a coordinatewise. In particular

\[
\boxed{0\le\Gamma_i\le D_i,\qquad E_i=D_i-\Gamma_i\ge0.}
\tag{3}
\]

These are necessary conditions, not sufficient rank certificates.

## 2. Minimum transverse degree of any M_H distinct raw columns

For the uniform counting measure on all raw columns, let x=(w dot alpha)/H. Lattice integration over h and alpha gives the limiting density

\[
f_m(x)=m(m+1)x^{m-1}(1-x),\qquad 0\le x\le1.
\tag{4}
\]

Indeed, type cost x has simplex-shell volume proportional to x^(m-1), and the available h interval has length proportional to 1-x. The weights cancel after normalization. Its cumulative distribution is

\[
F_m(A)=(m+1)A^m-mA^{m+1}.
\tag{5}
\]

Let A=A_m(lambda) be the unique point in (0,1] with F_m(A)=lambda. Selecting the M_H cheapest transverse costs minimizes their total over every M_H-column choice. The resulting limiting minimum mean is

\[
\boxed{q_m(\lambda)=
\frac{mA}{m+2}
\frac{(m+2)-(m+1)A}{(m+1)-mA}.}
\tag{6}
\]

To obtain (6), integrate x f_m(x) from zero to A and divide by lambda. Weighted-lattice cube estimates imply convergence of the empirical distribution and truncated first moments. The limiting distribution is continuous and strictly increasing on (0,1), so any finite tie at the selection cutoff is lower order. Therefore, uniformly over all selected raw-column families,

\[
\boxed{\frac1{M_HH}\sum_{C_H}w\cdot\alpha
\ge q_m(\lambda)-o(1).}
\tag{7}
\]

This statement uses distinct columns only. Imposing nonzero source or period rank can only increase the minimum.

The complete source simplex has mean

\[
\frac1{M_HH}\sum_{L_H}w\cdot\beta
\longrightarrow\frac{m\theta}{m+2}.
\]

Consequently every nonzero source raw minor obeys

\[
\boxed{\liminf\frac{w\cdot D}{M_HH}
\ge\max\left\{0,q_m(\lambda)-\frac{m\theta}{m+2}\right\}.}
\tag{8}
\]

The zero in (8) also uses D_i>=0 for a nonzero source minor. It matters in some small-dimensional parameter choices.

## 3. The half-column-density family keeps its arithmetic cost

In the previous critical construction, v_0=2K theta^m w_0, so lambda=1/2 exactly in the H limit. For every fixed lambda>0,

\[
q_m(\lambda)\longrightarrow1\quad(m\to\infty).
\tag{9}
\]

For a short proof, for each fixed epsilon>0 the entire-column fraction with x<=1-epsilon is F_m(1-epsilon), which tends to zero. Among any fixed positive fraction lambda of columns, all but a vanishing fraction therefore have x>=1-epsilon. Let epsilon decrease to zero. The upper bound x<=1 gives the converse inequality.

There is also a first-order correction. For fixed 0<lambda<1, let y>0 solve

\[
e^{-y}(1+y)=\lambda.
\]

The rescaled boundary distance Y=m(1-x) has density
((m+1)/m)y(1-y/m)^(m-1) on [0,m], which converges to y exp(-y). Its second moment is 6m^2/((m+2)(m+3)), so the first moments are uniformly integrable. Consequently the cutoff m(1-A_m) tends to y, and the truncated first moments give

\[
q_m(\lambda)=1-\frac{y+1+1/(y+1)}m+o(m^{-1}).
\tag{9a}
\]

In particular, at theta=1/2 the lower bound in (8) is
1/2-[y+1/(y+1)]/m+o(m^(-1)). The finite-dimension reduction in this degree cost is only of order 1/m. This does not compare it to any favorable analytic estimate of the same order.

The available upper bound w dot D/(M_HH)<=1-m theta/(m+2)+o(1), together with (8), shows that in this order of limits every raw-column choice has the same leading transverse gap 1-theta.

To transfer this to the actual arithmetic cost, write

\[
d_i=\ln(q_iL_i|b_i|),\qquad \kappa_i=d_i/w_i,
\qquad \mathfrak A_H=\frac1{M_HH}\sum_i d_iD_i.
\]

If max_i|kappa_i-1|<=delta and D_i>=0, then

\[
\left|\mathfrak A_H-\frac{w\cdot D}{M_HH}\right|\le\delta,
\tag{10}
\]

because w dot D<=M_HH. More precisely, take a family indexed by m with fixed lambda>0 and theta, and suppose delta_m=max_i|kappa_i-1| tends to zero. First take the H limit separately for each fixed member of this family; then take m to infinity. The arithmetic exponent converges to 1-theta uniformly over nonzero selected raw minors. This statement does not require delta to tend to zero at fixed m, which a persistent lcm cost could prevent. At theta=1/2 the exponent tends to 1/2. This is a lower bound on the size of the *denominator exponent in the proved arithmetic bound*, not a lower bound on the actual logarithmic size of every determinant.

For scale, exact rational bisection of (5) at lambda=1/2 proves the following weaker rational bounds, before the delta error in (10):

| Dimension m | Lower bound for q_m(1/2)-m/(2(m+2)) |
|---:|---:|
| 40 | > 0.451 |
| 80 | > 0.475 |
| 120 | > 0.483 |

The quantile cutoff is an algebraic number; the verifier uses rational enclosing intervals. It does not round a floating-point root to obtain these claims.

## 4. Analytic growth with the actual residual degrees

At omega=2*pi*i use ordinary Taylor rows and the same raw columns. Fix R>7K, let r=7K/R and c_R=ln(1/r)/2. Write

\[
S_h=\sum_{C_H}h,\quad
S_\alpha=\sum_{C_H}|\alpha|,\quad
S_a=\sum_J|a|.
\]

For a column P=Y^hX^alpha, its type-a function is

\[
f_{a,P}(z)=\binom\alpha a e^{hz}z^{|\alpha|-|a|}.
\]

On |z|=R, its norm is at most

\[
e^{Rh}\,2^{|\alpha|}\,R^{|\alpha|-|a|}.
\tag{11}
\]

In the coefficient determinant used by the distinct-Taylor-degree collision argument, factor R^(-|a|) from every row and e^(Rh)(2R)^|alpha| from every selected column. The remaining entries have modulus at most one. Therefore the same summation proof yields

\[
\boxed{
\ln|\det B_J|
\le R S_h+(\ln2)S_\alpha+(\ln R)\sum_iE_i
-c_R\sum_a n_a(J)^2+\mathcal R_H,}
\tag{12}
\]

where

\[
\mathcal R_H=2\ln(M_H!)+c_RM_H-M_H\ln(1-r)=o(M_HH).
\]

The estimate is uniform over all target sets and raw-column choices at fixed parameters. It retains the selected column degrees and the target degree increase. In particular it does not separately charge R times H/w_0 to every column and ln(R) times its full transverse degree while forgetting the target-row degree.

For the grouped Pascal transport, Hadamard's column estimate gives

\[
\ln|\det T_{L,J}|
\le(\ln2)S_a-\sum_i\ell_i\Gamma_i,
\qquad \ell_i=-\ln|\zeta_i|.
\tag{13}
\]

This is the unrounded form of the previous uniform ln(2)/w_* transport bound. Combining (12), (13), and the arithmetic exponent on the same raw columns gives

\[
\boxed{
\ln|\det T_{L,J}\det B_J|+\sum_i d_iD_i
\le R S_h+(\ln2)(S_\alpha+S_a)
+\sum_i\bigl[(d_i+\ln R)E_i-(\ell_i-d_i)\Gamma_i\bigr]
-c_R\sum_a n_a(J)^2+\mathcal R_H.}
\tag{14}
\]

The relation D=Gamma+E is used exactly once. In an endpoint-like regime with d_i/w_i->1 and ell_i/w_i->2, the main weighted degree competition is w dot E-w dot Gamma, plus the separately displayed exponential-column and collision terms. It is therefore necessary to control the residual degree as well as the amount moved; Gamma<=D alone goes in the wrong direction to force a favorable contradiction.

Formula (14) is a sharpened same-column sufficient upper estimate. It is not a lower bound on a term, a claim that one term dominates, or an assertion that the surviving spectral/collision contribution is favorable near total escape. The number of exterior target sets still has logarithm o(M_HH), so a strict uniform negative margin in (14) would suffice for the absolute-value approach. Such a margin has not been established.

## 5. Scope of diagonal rescaling

The two nonnegative degree vectors Gamma and E measure different steps: source-to-target transport and target-to-column residual degree. Multiplying coordinate i by any nonzero scalar rescales these two factors by the corresponding powers and rescales the full determinant by the power D_i=Gamma_i+E_i. A common diagonal change therefore cannot erase a residual-degree cost while retaining the original arithmetic normalization. More flexible norms may improve (14), but their column and row factors must be carried through both sides.

## Verification

`verify_column_cost.py` checks finite cheapest-column minimization against exhaustive enumeration; the source/column density counts and predicted quantile moments; rational enclosures for the half-density constants; nonzero determinant matching degree conditions; and the residual-degree factorization in the actual-period coefficient formulas at exact formal specializations. Finite tests do not prove the asymptotic claims, which follow from the written lattice and integral arguments.
