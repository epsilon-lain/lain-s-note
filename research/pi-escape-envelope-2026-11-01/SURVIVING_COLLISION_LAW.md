# Surviving-source collisions and concentration near total escape

2026-11-01. The source packet, target row set, and period matrix below are the same ones used in the preceding positive-fraction transport note. The main result is a uniform upper bound for compatible period minors as a function of the escape fraction. It does not close the endpoint inequality. A conditional consequence locates the remaining contribution in a thin region near total escape.

## 1. Objects and normalization

Fix an integer dimension m>=1, K>=1 centers, positive rational weights v_0,v_1,...,v_m, w_0,w_1,...,w_m, and 0<theta<=1, with v_i=w_i/theta>v_0. All these data are fixed before H tends to infinity on the admissible weight lattice. Let

\[
L_H=\{(j,s,\beta):0\le j<K,\quad v_0s+v\cdot\beta<H\},
\qquad M_H=|L_H|.
\]

All derivative indices are nonnegative integers. A target set J contains M_H distinct rows (j,t,a), with t>=0 and w dot a<=H. Write k=|J\setminus L_H| and rho_H=k/M_H. For each transverse type a, let n_a(J) count its target rows, including all centers and all time orders.

The source type multiplicities are

\[
n_\beta=K\left\lceil\frac{H-v\cdot\beta}{v_0}\right\rceil
\quad(v\cdot\beta<H).
\]

At the actual period omega=2*pi*i use ordinary normalized time derivatives:

\[
B_{j,t,a;h,\alpha}
=[z^tu^a]e^{hz}\prod_i(j\omega+z+u_i)^{\alpha_i},
\qquad w_0h+w\cdot\alpha\le H.
\]

Choose any M_H raw columns satisfying this budget. The bounds below hold uniformly over those column choices. The ordinary time convention matters in Section 3.

## 2. The exact surviving-count optimization

The M_H-k source rows retained by J have multiplicities r_beta with 0<=r_beta<=n_beta and sum r_beta=M_H-k. Since n_beta(J)>=r_beta,

\[
\sum_a n_a(J)^2\ge
\min\left\{\sum_\beta r_\beta^2:
0\le r_\beta\le n_\beta,\quad\sum_\beta r_\beta=M_H-k\right\}.
\tag{1}
\]

Allowing real r_beta only lowers this minimum. Its minimizer is

\[
r_\beta=\min(n_\beta,\lambda),\qquad
\sum_\beta\min(n_\beta,\lambda)=M_H-k.
\tag{2}
\]

Indeed, two unsaturated positive coordinates with different values can be moved toward each other, preserving their sum and reducing the sum of squares. The boundary cases k=0 and k=M_H follow by continuity. The integer minimum is obtained by allocating the cheapest marginal costs 1,3,5,... in each capped coordinate; it differs from (2) by at most O(sum n_beta + number of types), which is O(H^(m+1)).

Under the uniform type count on the simplex, y=1-v dot beta/H has limiting density m(1-y)^(m-1), 0<=y<=1. In units KH/v_0, the cap n_beta tends to y. Writing lambda=(KH/v_0)ell, the retained-count constraint gives

\[
\frac{\int_0^1\min(y,\ell)m(1-y)^{m-1}\,dy}
     {\int_0^1y m(1-y)^{m-1}\,dy}
=1-(1-\ell)^{m+1}=1-\rho.
\]

Thus ell=1-rho^(1/(m+1)). Integrating the squared caps proves

\[
\boxed{\frac{\sum_a n_a(J)^2}{M_HH}
\ge\frac{2K}{v_0(m+2)}F_m(\rho_H)-\varepsilon_H,}
\tag{3}
\]

where one epsilon_H tends to zero uniformly over all target sets, and

\[
F_m(\rho)=1-(m+2)\rho+(m+1)\rho^{(m+2)/(m+1)}.
\tag{4}
\]

For uniformity, the functions min(y,ell) and min(y,ell)^2 are uniformly bounded and uniformly Lipschitz in y for 0<=ell<=1. Weighted-lattice cube estimates therefore give uniform Riemann-sum errors. Their integrated retained-count function is continuous and strictly increasing in ell. Its inverse is uniformly continuous on [0,1], including the endpoints. These facts transfer the count constraint to the cost uniformly in k. The ceiling errors are lower order after division by M_H H.

Useful consequences, with continuous endpoint interpretations, are

\[
F_m(0)=1,\quad F_m(1)=0,\quad F_m(\rho)>0\ (\rho<1),
\]

\[
\boxed{F_m(\rho)\ge\frac{m+2}{2(m+1)}(1-\rho)^2.}
\tag{5}
\]

For (5), F_m(1)=F_m'(1)=0 and F_m''(rho)=(m+2)/(m+1) rho^(-m/(m+1)) is at least (m+2)/(m+1). Integrate backward twice. For fixed rho, the separate limit m->infinity is 1-rho+rho ln rho, with 0 ln 0=0. The minimization is sharp as a deletion-only count problem; it does not assert that every optimizer can be realized by a nonzero transport and period minor simultaneously.

## 3. Collision estimate for ordinary time jets

For each type a and raw column P let

\[
f_{a,P}(z)=[u^a]P(e^z,z+u_1,\ldots,z+u_m)
=\sum_{d\ge0}c_{a,d,P}(z/R)^d,
\]

where R>7K is fixed. Put

\[
E_R=\max\{R/w_0,\ln(2R)/w_*\},\quad
w_* =\min_i w_i,\quad r=7K/R,\quad
c_R=\tfrac12\ln(1/r).
\]

The raw column budget gives |c_{a,d,P}|<=exp(H E_R). For the ordinary normalized time test, Cauchy's formula on |z|=1 gives

\[
\left|[z^t]\left(\frac{j\omega+z}{R}\right)^d\right|
\le r^d,
\tag{6}
\]

because |j omega+z|<=2*pi*(K-1)+1<7K. There is no factor 2^t here. Such a factor was needed in the earlier logarithmic-time tests on |t|=1/2; it must not be carried into this ordinary-time estimate. Conversely, a sparse logarithmic-jet packet cannot simply be substituted into (6).

Expand the determinant in these coefficient rows. Repeated (a,d) labels give zero. For a type appearing n times the absolutely convergent distinct-degree sum is

\[
\sum_{d_1,\ldots,d_n\text{ distinct}}r^{\sum_i d_i}
=\frac{n!\,r^{n(n-1)/2}}{\prod_{j=1}^n(1-r^j)}.
\]

The coefficient determinant is bounded by M_H! exp(M_H H E_R). Using product n_a!<=M_H! and product(1-r^j)>= (1-r)^M_H gives the uniform finite inequality

\[
\frac{\ln|\det B_J|}{M_HH}
\le E_R-c_R\frac{\sum_a n_a(J)^2}{M_HH}
+\frac{2\ln M_H+c_R-\ln(1-r)}H.
\tag{7}
\]

Zero determinants satisfy this with ln 0=-infinity. Combining (3) and (7) gives the period-minor envelope

\[
\boxed{\frac{\ln|\det B_J|}{M_HH}
\le E_R-S_m F_m(\rho_H)+o(1),\qquad
S_m=\frac{2c_RK}{v_0(m+2)}.}
\tag{8}
\]

This is uniform over all allowed targets and raw columns. In particular, rho_H->0 retains the original collision saving S_m. It does not compare a chosen minor to its own possibly much smaller source-period value.

## 4. Complete transport sum and the region it leaves open

For the same normalized shear transport T and fixed errors 0<|zeta_i|<1, the previous positive-fraction theorem gives uniformly

\[
\frac{\ln|\det T_{L,J}|}{M_HH}
\le \frac{\ln2}{w_*}
-a_0\rho_H^2+o(1),\qquad
a_0=\frac{\tau}{(m+1)\mathcal A},
\quad \tau=\min_i\frac{-\ln|\zeta_i|}{v_i-v_0},
\quad \mathcal A=\prod_i\frac{v_i}{v_i-v_0}.
\]

There are only exp(O(M_H)) possible target sets at fixed parameters, so the complete exterior sum satisfies

\[
\boxed{\limsup_{H\to\infty}\frac{\ln|\det(TB)|}{M_HH}
\le E_R+\frac{\ln2}{w_*}
-\min_{0\le\rho\le1}\{S_mF_m(\rho)+a_0\rho^2\}.}
\tag{9}
\]

This is a proved sufficient upper estimate, not a favorable endpoint inequality. The displayed minimum is at most a_0, by evaluating at rho=1. In the slow-quality endpoint regime tau->2 theta and subsequent large-weight limit A->1, a_0 is of order 1/m. Thus (9) alone does not beat an arithmetic cost tending to 1-theta. Allowing rho=1 in this relaxation is not a proof that such a compatible nonzero target exists.

Nevertheless the low-escape region can be removed quantitatively. Let S_low be the sum over rho_H<=1-epsilon, and S_high the remaining sum. If the rational arithmetic bound on the same raw minor gives

\[
|\det(TB)|\ge\exp(-M_HH\mathfrak A_H),
\]

and S_m F_m(1-epsilon)>E_R+ln(2)/w_*+limsup A_H+delta for some delta>0, then, for all sufficiently large H,

\[
|S_{\rm low}|\le e^{-M_HH(\mathfrak A_H+\delta/2)},
\qquad
|S_{\rm high}|\ge e^{-M_HH\mathfrak A_H}
(1-e^{-M_HH\delta/2}).
\tag{10}
\]

The factor in parentheses on the second line multiplies the preceding exponential. This follows from the triangle inequality applied to det(TB)=S_low+S_high. It is a statement about the grouped high-escape sum, not the sign or dominance of a particular term.

## 5. A shrinking surviving fraction in the earlier high-dimensional family

Take theta=1/2, B_0=3/5, C=13/10, K=floor(C^m), w_0=B_0^(-m), v_0=2K theta^m w_0 and R=14K. Each m is fixed before the transverse weights and then H are chosen. Here c_R=ln(2)/2 and

\[
S_m=\frac{c_R}{m+2}(6/5)^m.
\]

For sufficiently large m define

\[
\epsilon_m=\sqrt{12(m+1)(5/6)^m}<1.
\]

Equation (5) gives the exact uniform lower bound

\[
\boxed{S_mF_m(1-\epsilon_m)\ge6c_R=3\ln2>2.}
\tag{11}
\]

Under the previously stated raw-rank and arithmetic assumptions, one can choose the transverse weights so that E_R+ln(2)/w_*+limsup A_H<1. For precision, if w_i=ceil(ln q_i), 1/2<=|b_i|<=2, and L_i=lcm(1,...,ceil(w_i/(theta v_0))-1), set kappa_i=ln(q_i L_i |b_i|)/w_i. The elementary lcm bound yields

\[
1-2/w_*\le\kappa_i
\le1+4\ln2/(\theta v_0)+1/w_*.
\]

The column budget and the source simplex mean then give

\[
\limsup\mathfrak A_H
\le1-\frac{m\theta}{m+2}
+\frac{4\ln2}{\theta v_0}
+\frac{1+2m\theta/(m+2)}{w_*}.
\tag{12}
\]

The sum of the right side with E_R+ln(2)/w_* tends to 1/2 along this order of choices. The source-rank assumption is still required to obtain a nonzero rational minor; this note does not prove the manuscript's interpolation theorem or supply actual approximation tuples meeting arbitrary thresholds.

Consequently (10) applies to epsilon_m. In this conditional regime, the part of the sum able to supply the arithmetic lower bound must have escape fraction greater than 1-epsilon_m. Exact rational bounds on epsilon_m^2 give, for example, epsilon_80<0.022 and epsilon_120<0.0007. These percentages describe this asymptotic determinant construction, not a density of rational approximations to pi.

The remaining problem is now concentrated near total escape: estimate or cancel those terms using the true column compatibility. Neither the deletion-count optimization nor qualitative period nonvanishing controls that region.

## Verification scope

`verify_surviving_collision.py` checks the integer optimization against exhaustive finite enumeration, the continuous capped minimizer against the integer optimum, complete-packet counts and limiting formulas, ordinary time coefficient bounds, and rational high-dimensional constants. Asymptotic statements are proved by the preceding lattice and integral arguments, not inferred from those finite checks.
