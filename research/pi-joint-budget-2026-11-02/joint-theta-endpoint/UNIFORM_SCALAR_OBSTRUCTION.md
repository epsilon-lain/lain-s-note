# A uniform obstruction to the scalar mean/collision endpoint bound

2026-11-02. The earlier fixed-theta obstruction left open a joint limit theta=theta_m->1. This note closes that loophole for the **scalar mass/mean collision relaxation with the displayed raw-column Cauchy and LCM costs**. For every allowed scalar column mean, an explicit scalar target mean makes its endpoint upper-bound expression strictly positive. The result holds at every fixed dimension after the stated H and large-weight limits; it needs no fixed exponential rates B or C.

This is a statement about a relaxed upper estimate, not a positive lower bound for an actual determinant. An admissible scalar target mean need not satisfy all coordinatewise degree inequalities for an arbitrary imbalanced column set. Exact time grades, actual-period rank, complete-sum cancellation, other norms and reduced denominator arguments are not covered.

## 1. Parameters and the endpoint expression

Fix m>=1, 0<theta<=1, K>=1 and w_0>0. The critical row/column count relation is

\[
\lambda=\tfrac12,\qquad
v_0=2K\theta^m w_0,\qquad
R=14K,\qquad c=\tfrac12\ln2.
\tag{1}
\]

Keep every weight fixed while H tends to infinity. The leading scalar source and selected-column means are

\[
a=\frac m{m+2},\qquad s=\theta a,\qquad
u=\lim\frac{\sum_{C_H}w\cdot\alpha}{M_HH},\qquad
y=\lim\frac{w_0\sum_{C_H}h}{M_HH}.
\tag{2}
\]

Pass to a convergent subsequence when needed. A nonzero joint term requires s<=u; the raw budget gives u+y<=1. All columns are distinct.

Let I_m(t) be the sharp scalar mass/mean collision value proved in [the preceding mean/collision note](../../pi-aggregate-support-2026-11-02/MEAN_DEGREE_AND_COLLISIONS.md). Define

\[
P(t)=\frac{c}{2(m+1)w_0}I_m(t).
\tag{3}
\]

For 0<t<=a,

\[
P(t)=\frac{c}{w_0(m+2)}\left(\frac at\right)^m.
\tag{4}
\]

For all t>0, I_m(t)>=1. For t<=a the collision optimization has an attaining radial density, with integer lattice approximations after all additional row constraints are discarded.

The actual LCM normalization, with w_i=ceil(ln q_i) and |b_i| bounded above and away from zero, gives in the later large-weight limit

\[
\kappa_i=\frac{\ln(q_iL_i|b_i|)}{w_i}
\ge\kappa_0:=1+\frac{\ln2}{\theta v_0}
=1+\frac{c}{K\theta^{m+1}w_0}.
\tag{5}
\]

More precisely, the finite-weight lower bound is kappa_i>=kappa_0-epsilon_w with epsilon_w->0 for each fixed parameter set. This follows from the central-binomial divisibility bound for the LCM in the preceding note. It does not assert an intrinsic lower bound for a different, reduced denominator normalization.

At endpoint shear quality ell_i/w_i=2, the optimistic scalar upper-bound expression is

\[
\boxed{\mathcal F(u,y,t)=
\frac R{w_0}y+\kappa_0(u-s)-2(t-s)-P(t),
\qquad s\le t\le u.}
\tag{6}
\]

The omitted binomial and residual-degree Cauchy terms are nonnegative. The word "optimistic" means that making (6) nonnegative already prevents this particular upper-bound calculation from proving a uniform negative margin. It does not give any sign information about the actual determinant.

## 2. Distinct time exponents couple column cost to collisions

For each transverse type alpha, let N_alpha count the selected columns Y^h X^alpha. Within a type, the h values are distinct nonnegative integers. Therefore the exact finite inequality is

\[
\boxed{\sum_{C_H}h\ge\frac12\sum_\alpha N_\alpha(N_\alpha-1).}
\tag{7}
\]

Apply the sharp mean/collision theorem to the column multiplicities N_alpha themselves. Their total mass is M_H and their transverse mean is u. After H tends to infinity,

\[
\boxed{y\ge\frac1{4(m+1)}I_m(u).}
\tag{8}
\]

The discarded term in (7), after multiplication by w_0/(M_HH), is w_0/(2H), which tends to zero at fixed weights. Continuity of I_m handles a sequence of means converging to u. Thus (8) holds for every selected column set, not only a specially chosen cutoff or parity family.

Combining (3) and (8),

\[
\boxed{\frac R{w_0}y\ge\frac R{2c}P(u)
=\frac{7K}{c}P(u).}
\tag{9}
\]

Equation (9) is the missing same-column link: concentrating transverse column degrees also forces distinct time exponents and their analytic cost. The selected column mean u and the analytic time cost y cannot be optimized independently.

## 3. Universal positive scalar witness when s<=u<=a

Choose

\[
t=\sqrt{su},\qquad
z=P(u)>0,\qquad
x=\frac{(m+2)(u-s)}{\theta},\qquad
r=\left(\frac us\right)^{m/2}
=\left(1+\frac xm\right)^{m/2}.
\]

This satisfies s<=t<=u<=a. The power law (4) gives

\[
P(t)=zr,\qquad P(s)=zr^2,
\]

while the LCM excess has the exact relation

\[
(\kappa_0-1)(u-s)=\frac{x}{K}P(s)=\frac{zr^2x}{K}.
\]

Using (9) in (6),

\[
\boxed{\mathcal F(u,y,t)
\ge(\sqrt u-\sqrt s)^2
+z\left(\frac{7K}{c}+\frac{r^2x}{K}-r\right)>0.}
\tag{10}
\]

The last inequality is elementary and uniform in all parameters. If x>=1/4, completing the square gives xr^2/K-r>=-K/(4x)>=-K, so the bracket is at least K(7/c-1)>0. If 0<=x<1/4, then r<exp(1/8)<2, while 7K/c>14. No limiting assumption on m or theta was used.

## 4. Universal positive scalar witness when u>=a

Choose t=a sqrt(theta)<=a<=u, and put

\[
p=K\theta^{m/2},\qquad
\rho=\frac{\theta^{-m/2}}{w_0},\qquad
x=\frac{m(1-\theta)}{\theta}.
\]

Now (8) and I_m(u)>=1 imply y>=1/[4(m+1)]. Also u-s>=a(1-theta). Substitution in (6) gives

\[
\boxed{
\mathcal F(u,y,t)
\ge (u-a)+a(1-\sqrt\theta)^2
+\frac{\rho}{m+2}
\left[
\frac72\frac{m+2}{m+1}p+\frac{cx}{p}-c
\right]>0.}
\tag{11}
\]

For x>=1, the first two terms in the bracket have sum at least 2sqrt(7cx/2)>c. For x<1,

\[
p\ge\theta^{m/2}
=\left(1+\frac xm\right)^{-m/2}
\ge e^{-x/2}>\tfrac12,
\]

so the first term alone is greater than 7/4>c. This includes theta=1 and the boundary u=a. Again the strict positivity is dimension-independent as an assertion, although its numerical margin may tend to zero along a parameter sequence.

## 5. What is and is not ruled out

For every scalar pair (u,y) arising from any selected raw-column family with u>=s,

\[
\boxed{\sup_{s\le t\le u}\mathcal F(u,y,t)>0.}
\tag{12}
\]

In particular, optimizing the raw-column choice cannot make this **scalar** envelope strictly negative. The result covers arbitrary joint sequences theta_m->1, B_m, C_m or K_m once the critical count relation and the specified Cauchy/LCM normalization are retained. There is no fixed-theta assumption and no condition such as C B<1 or C theta>B. The finite-weight error in (5) changes the conclusion by at most epsilon_w(u-s); in a joint large-weight regime where these errors vanish, it cannot create a uniform negative margin bounded away from zero.

The scalar target means in (10)--(11) are feasible for the mass/mean collision optimization, whose radial optimizer is known. They do not automatically satisfy every inequality Gamma_i>=0 and E_i>=0 for a particular imbalanced column set. For scalar column-cost optimization, conditionally average a continuum selection density over spatial directions at each fixed time coordinate and total transverse cost. This preserves its mass, u, y, and every isotropic scalar cost, and makes the weighted coordinate means all equal to u/m. The radial target profiles then obey the coordinatewise degree constraints. Such selection densities admit distinct-column lattice approximations; use a small strict-feasibility perturbation before rounding at an equality boundary.

This symmetrization identifies feasible representatives of the scalar relaxation. It does **not** prove that an arbitrary imbalanced nonzero minor has a compatible radial target, or that imposing all coordinatewise and time-grade constraints leaves the same minimax value. Those constraints, stronger actual-period estimates, complete-sum cancellation, different norms and reduced denominator arguments remain possible avenues. No claim about pi's bad approximability follows.

## Verification

`verify_uniform_scalar_obstruction.py` checks the exact finite column inequality (7), its finite dual moment consequences, both scalar witness formulas with rational square choices of theta and u, joint theta_m->1 families, and finite selected-column moments. It encloses c=ln2/2 between rational endpoints and checks the endpoint expressions at both, so its positivity tests do not rely on floating-point roots or logarithms. Infinite-lattice and parameter statements follow from the proofs above.
