# Joint raw-column cost, exact quantiles, and source-degree constraints

2026-11-02. This note optimizes the positive analytic and arithmetic cost on the same family of distinct raw columns. It gives an exact weighted-simplex quantile, a constrained variant retaining the necessary coordinate degree inequalities, and a closed form in a fast theta-to-one regime. The optimization does not certify a nonzero source or period minor.

## 1. The combined cost that actually occurs

For the selected raw columns, use normalized means

\[
y_H=\frac{w_0\sum_C h}{M_HH},\qquad
u_{i,H}=\frac{w_i\sum_C\alpha_i}{M_HH},\qquad
t_{i,H}=\frac{w_i\sum_J a_i}{M_HH},\qquad
s_{i,H}=\frac{w_i\sum_L\beta_i}{M_HH}.
\]

At fixed parameters, the complete source packet has
\(s_{i,H}\longrightarrow s_i=\theta/(m+2)\). The exact same-column residual-degree estimate can be written

\[
\begin{split}
\frac{\ln|\det T_{L,J}\det B_{J,C}|+\sum_i d_iD_i}{M_HH}
\le{}&A y_H+\sum_i c_i u_{i,H}
-\sum_i b_i t_{i,H}+\sum_i e_i s_{i,H}\\
&-c_R\frac{\sum_a n_a^2}{M_HH}+o(1),
\end{split}
\tag{1}
\]

where

\[
A=R/w_0,\quad c_i=\frac{d_i+\ln(2R)}{w_i},\quad
b_i=\frac{\ell_i+\ln R-\ln2}{w_i},\quad
e_i=\frac{\ell_i-d_i}{w_i}.
\tag{2}
\]

Here d_i is the arithmetic denominator exponent and ell_i=-ln|zeta_i|. Formula (1) follows by expanding E_i=sum_C alpha_i-sum_J a_i and Gamma_i=sum_J a_i-sum_L beta_i in the earlier residual estimate. In particular, the positive column cost is **one** linear functional

\[
A y_H+\sum_i c_i u_{i,H}.
\tag{3}
\]

Minimizing its time and spatial parts on independently chosen column sets need not minimize (3). The target, source and collision terms in (1) remain necessary; a lower bound on (3) is not a lower bound on an actual determinant.

## 2. The sharp anisotropic column quantile

Let d=m+1. Under uniform counting of all raw columns with w_0h+w dot alpha<=H, the normalized exponent vector

\[
z=(z_0,z_1,\ldots,z_m)
=(w_0h,w_1\alpha_1,\ldots,w_m\alpha_m)/H
\]

converges to uniform probability measure mu on the standard d-dimensional simplex
\(\Delta_d=\{z\ge0:\sum z_i\le1\}\). The remaining barycentric coordinate is the unused budget. Suppose M_H/N_H tends to lambda in (0,1).

For a coefficient vector c=(A,c_1,...,c_m), write L_c(z)=A z_0+sum_i c_i z_i. The spatial c_i may be signed in the dual problem below. Define

\[
Q_\lambda(c)=\frac1\lambda\int_{L_c\le q_\lambda}L_c\,d\mu,
\qquad \mu(L_c\le q_\lambda)=\lambda.
\tag{4}
\]

For a nonconstant linear functional, its distribution has no atoms and the quantile is unique in the interior of its range. Ordering the actual raw columns by this cost proves

\[
\boxed{\min_{|C_H|=M_H}
\frac1{M_HH}\sum_{(h,\alpha)\in C_H}
\left(Rh+\sum_i(d_i+\ln(2R))\alpha_i\right)
\longrightarrow Q_\lambda(c).}
\tag{5}
\]

The statement is uniform over choices of C_H in the lower-bound direction. It follows from lattice Riemann sums and convergence of the truncated first moment at the quantile; ties at a finite cutoff have a lower-order effect. All weights and coefficients are fixed for this H limit. A subsequent dimension-dependent parameter choice does not silently provide an error estimate uniform in H and m.

There is an explicit one-variable formula. Let the d+1 vertex costs be

\[
k_0=0,\quad k_1=A,\quad k_{i+1}=c_i\quad(1\le i\le m).
\]

When these knots are pairwise distinct, put

\[
F(q)=\sum_{i=0}^{d}
\frac{(q-k_i)_+^d}{\prod_{j\ne i}(k_j-k_i)},\qquad
J(q)=\frac1{d+1}\sum_{i=0}^{d}
\frac{(q-k_i)_+^{d+1}}{\prod_{j\ne i}(k_j-k_i)}.
\tag{6}
\]

Then F is the distribution function of L_c and J=integral_{-infinity}^q F. In particular,

\[
\boxed{F(q_\lambda)=\lambda,\qquad
Q_\lambda(c)=q_\lambda-\frac{J(q_\lambda)}\lambda.}
\tag{7}
\]

One proof of (6) integrates exp(t L_c) on the simplex by induction on dimension, obtaining

\[
\int e^{tL_c}\,d\mu
=\frac{d!}{t^d}\sum_{i=0}^{d}
\frac{e^{tk_i}}{\prod_{j\ne i}(k_i-k_j)}.
\]

The derivative of F in (6) has the same transform. The polynomial interpolation identities cancel its tails beyond the largest knot and give total mass one. Integrating F gives J, and integration by parts gives (7). Coincident knots are handled by continuity, equivalently confluent divided differences; one must not divide by a zero knot difference. The special-time-cost companion note gives formulas without this issue when all spatial coefficients coincide.

There is also a useful variational representation, valid with repeated knots:

\[
\boxed{Q_\lambda(c)=\sup_q\left[q-
\lambda^{-1}\int(q-L_c)_+\,d\mu\right].}
\tag{8}
\]

Its derivative is 1-F(q)/lambda, so the maximum occurs at the quantile. The function Q_lambda is concave and positively homogeneous in the cost vector, since it is the infimum of linear costs over selections of fixed mass. If two cost vectors differ by at most epsilon in every coordinate, their Q values differ by at most epsilon, because sum z_i<=1. To resolve a joint-theta gap tending to zero, coefficient errors must be smaller than that gap; convergence merely to zero is not enough.

## 3. Retaining the necessary source-degree inequalities

A nonzero source minor must have D_i=sum_C alpha_i-sum_L beta_i>=0 for every i. At the continuum level this requires the selected means u_i>=s_i. Define the relaxed sharp cost

\[
V_\lambda(c;s)=
\inf\left\{\lambda^{-1}\int L_c f\,d\mu:
0\le f\le1,\ \int f\,d\mu=\lambda,
\ \lambda^{-1}\int z_i f\,d\mu\ge s_i\ (1\le i\le m)\right\}.
\tag{9}
\]

For the source s_i=theta/(m+2) and theta<1, the constant selection density f=lambda is strictly feasible in all coordinate inequalities. The compact convex set of feasible moment vectors and the separating-hyperplane argument for a linear program therefore give

\[
\boxed{V_\lambda(c;s)=
\sup_{\eta_i\ge0}\left[
\sum_i\eta_i s_i+
Q_\lambda(A,c_1-\eta_1,\ldots,c_m-\eta_m)\right].}
\tag{10}
\]

The lower-bound direction also follows immediately by subtracting the nonnegative multipliers times the coordinate constraints. Strict feasibility rules out a duality gap; a maximizing multiplier need not be presumed before taking the supremum.

This relaxation is sharp for asymptotic *distinct-column moment selection*. Approximate a density f by a simple density on small simplex cells; choose the required fraction of lattice points in each cell and adjust the total to M_H. Lattice counts give the prescribed limiting moments. To enforce the source inequalities despite rounding, first mix f with a small amount of the strictly feasible constant density. Then take the H limit, followed by the mixing and cell-size limits. This argument supplies neither independent interpolation rows nor a nonzero determinant.

If all normalized spatial costs equal kappa and all source bounds equal theta/(m+2), conditional averaging at fixed (z_0,sum_i z_i) preserves mass, time mean, spatial total and objective. It gives equal spatial coordinate means. Thus this isotropic moment problem loses nothing by using a density radial in the weighted spatial coordinates. At fixed parameters, distinct columns approximate that averaged density as above. This is not a claim that a radial target construction is compatible with every original unaveraged or time-labeled column set.

## 4. An exact relaxed degree frontier when theta approaches one quickly

Let n=m+1, c_i=kappa>0, and assume

\[
\theta^n\ge\lambda,\qquad
0<A\le\frac{\kappa\theta^n}{\lambda}.
\tag{11}
\]

Then the constrained optimum has the closed form

\[
\boxed{V_\lambda(A,\kappa;s)
=\frac{\kappa m\theta+A\lambda\theta^{-m}}{m+2},
\qquad s_i=\frac\theta{m+2}.}
\tag{12}
\]

Indeed, put a=theta^n/lambda>=1 and select the simplex

\[
S=\{z\ge0:a z_0+\sum_i z_i\le\theta\}.
\]

Because a>=1 and theta<=1, it lies inside the original column simplex. Its relative volume is theta^n/a=lambda. Its conditional means are

\[
u_i=\frac\theta{m+2},\qquad
y=\frac{\theta}{a(m+2)}
=\frac{\lambda\theta^{-m}}{m+2}.
\tag{13}
\]

For optimality, use the common multiplier eta=kappa-A/a>=0 in (10). The modified linear cost is (A/a)(a z_0+sum_i z_i), whose cheapest lambda fraction is exactly S. The source constraints are equalities there, so the dual bound is attained and equals (12).

The positive time growth plus the arithmetic degree gap is therefore at least

\[
\boxed{V_\lambda(A,\kappa;s)-\kappa\sum_i s_i
=\frac{A\lambda\theta^{-m}}{m+2}.}
\tag{14}
\]

For lambda=1/2, the condition theta^(m+1)>=1/2 describes a fast approach to one, roughly (m+1)(1-theta)<=ln2. At the limiting value theta=1, the formula gives a time mean lambda/(m+2), rather than the smaller unrestricted minimum time mean asymptotic to (1-ln2)/m at half density. The companion note quantifies the unrestricted joint time/spatial cost at all time coefficients.

This particular fast regime is not a parameter-eligible construction for the original multi-center argument. At lambda=1/2 and theta<1, theta^(m+1)>=1/2 implies theta^m>1/2. For every K>=2 it therefore violates the earlier necessary count condition K theta^m<1 used by that interpolation route. Formula (12) diagnoses the degree-only relaxation and must not be presented as a candidate nonzero source minor in the original multi-center regime. The general quantile and constrained dual remain available outside this explicit branch.

The frontier (12) is only a moment frontier. For K>=2 and an exact finite raw-column selection with D_i=0 for every i, every nonzero determinant permutation would have to match alpha=beta coordinatewise. One may therefore set all entries with alpha!=beta to zero without changing the determinant. The remaining type blocks have entries independent of the center j, giving duplicate source rows and determinant zero. The asymptotic condition D_i/(M_HH)->0 does not imply D_i=0 at finite H, so this elementary observation does not by itself invalidate approximation of the continuum frontier or prove that a nonzero minor attains it. Additional rank-compatible degree conditions, such as the separate center-coalescence estimate, must be imposed separately.

## 5. What the optimization settles

The positive cost of selecting actual distinct columns can be computed jointly rather than by unrelated time and spatial extrema. Formula (10) also prevents an optimizer from borrowing a negative arithmetic degree gap forbidden to a nonzero source minor. In the fast joint-theta regime, (12) makes that correction explicit.

The formulas are building blocks for the same-column upper envelope (1). They do not optimize the target collision and transport terms, remove their time-grade or rank constraints, or establish a negative endpoint margin. In particular, an optimized positive upper-bound cost must not be reported as a lower bound on a determinant.
