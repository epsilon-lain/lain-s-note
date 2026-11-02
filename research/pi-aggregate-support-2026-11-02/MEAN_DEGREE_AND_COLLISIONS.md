# Target mean degree, collisions, and a parameter obstruction

2026-11-02. This note gives a sharp optimization of target-type collisions under a bound on their mean weighted degree. It adds a second concentration statement to the earlier near-total-escape result: in the critical family, assuming the stated nonzero minor and a bounded normalized analytic/arithmetic budget, a surviving exterior contribution must also have mean target degree approaching at least B. The resulting mean/collision relaxation still cannot close the endpoint for any fixed 0<theta<1 under the current parameter and arithmetic normalization. This is an obstruction to that relaxation, not a counterexample to the actual determinant problem or a conclusion about pi.

## 1. Fixed-parameter setup

Fix dimension m>=1 and positive rational transverse weights w=(w_1,...,w_m). First let H tend to infinity, keeping every parameter fixed. Write

\[
\mathcal A_H=\{a\in\mathbb Z_{\ge0}^m:w\cdot a\le H\},\qquad
Q_H=|\mathcal A_H|,\qquad x_a=\frac{w\cdot a}{H}.
\]

Let M_H be the number of target rows, and assume

\[
\frac{M_H}{Q_HH}\longrightarrow c>0.
\tag{1}
\]

A target set determines nonnegative integer multiplicities n_a, satisfying sum_a n_a=M_H. Its mean degree is

\[
t_H=\frac1{M_H}\sum_a x_a n_a.
\tag{2}
\]

The source/raw-column density lambda from the preceding notes gives

\[
c=\frac{\lambda}{(m+1)w_0}.
\tag{3}
\]

The empirical distribution of x_a on all types tends to the probability measure

\[
d\mu_m(x)=m x^{m-1}\,dx,\quad 0\le x\le1.
\]

In this note, "sharp" refers to the relaxation retaining only total mass and the mean-degree bound. It does not assert that its minimizing profiles satisfy the exact time grades, transport minors, or actual-period rank conditions.

## 2. The sharp continuum problem

Put

\[
\mu=\frac m{m+1},\quad
\sigma^2=\frac{m}{(m+1)^2(m+2)},\quad
a_m=\frac m{m+2}.
\]

For 0<t_0<=1, minimize integral g^2 dmu_m over g>=0 with integral g dmu_m=1 and integral xg dmu_m<=t_0. The exact value is

\[
\boxed{I_m(t_0)=
\begin{cases}
\displaystyle\frac{2(m+1)}{m+2}
\left(\frac{m}{(m+2)t_0}\right)^m,
&0<t_0\le a_m,\\[6pt]
\displaystyle1+\frac{(\mu-t_0)^2}{\sigma^2},
&a_m\le t_0\le\mu,\\[6pt]
1,&\mu\le t_0\le1.
\end{cases}}
\tag{4}
\]

For the first branch, set A=(m+2)t_0/m<=1. The unique optimizer up to null sets is

\[
g(x)=\frac{m+1}{A^{m+1}}(A-x)_+.
\tag{5}
\]

Its mass is one, its mean is mA/(m+2)=t_0, and its squared norm is the first expression in (4). To prove optimality directly, for any u>=0 and feasible g,

\[
\int g^2d\mu_m
\ge2u(A-t_0)-u^2\int(A-x)_+^2d\mu_m.
\]

The pointwise inequality follows by completing the square when x<A, and because the right-hand linear term is nonpositive when x>=A. Using

\[
\int(A-x)_+^2d\mu_m
=\frac{2A^{m+2}}{(m+1)(m+2)}
\]

and maximizing over u gives equality at (5).

For the middle branch the optimizer is the nonnegative affine density

\[
g(x)=1+\frac{\mu-t_0}{\sigma^2}(\mu-x).
\tag{6}
\]

Its nonnegativity is exactly t_0>=a_m in the relevant decreasing case. Orthogonal projection onto the span of 1 and x proves the stated norm and optimality. For t_0>=mu the constant g=1 is feasible and Cauchy--Schwarz gives the lower bound one. These last formulas use a mean **upper bound**; they are not the solution for an exact prescribed mean above mu. At t_0=0 the continuum infimum is infinite.

## 3. Transfer to target multiplicities, including sharpness

Uniformly over all target profiles with t_H<=t_0,

\[
\boxed{\liminf_{H\to\infty}
\frac{\sum_a n_a^2}{M_HH}\ge c I_m(t_0).}
\tag{7}
\]

One can avoid any compactness assumption on the profiles. In the low-mean branch, the exact finite inequality is

\[
\boxed{\sum_a n_a^2\ge
\frac{M_H^2(A-t_0)^2}
{\sum_{a\in\mathcal A_H}(A-x_a)_+^2},
\qquad A>t_0.}
\tag{8}
\]

It follows from the same pointwise square completion, followed by optimization of u. Choose A=(m+2)t_0/m. The denominator divided by Q_H converges to the integral above. In the middle branch the corresponding affine dual test proves (7), and the last branch follows from sum n_a^2>=M_H^2/Q_H. The estimates are uniform over the target profiles because the only asymptotic errors are in sums of fixed test functions over the weighted lattice.

The lower value is asymptotically attained even by integer profiles when all other row constraints are discarded. Sample a bounded optimizing density, normalize its lattice values to total mass M_H, round each multiplicity down, and distribute the remaining fewer than Q_H units among distinct types. The rounding changes each multiplicity by at most one; since the original multiplicities are O(H), its normalized square-cost error is O(1/H). To ensure a one-sided mean constraint, first use a density with mean strictly below t_0 and then let that margin decrease after the H limit. A diagonal choice gives a single sequence. These profiles need not be supported on any nonzero exterior minor.

In particular, for t_0<=a_m and the source density (3),

\[
\boxed{\frac{\sum_a n_a^2}{M_HH}
\ge \frac{2\lambda}{w_0(m+2)}
\left(\frac{m}{(m+2)t_0}\right)^m-o(1).}
\tag{9}
\]

This counts collisions of **all** target rows. It does not rely on a positive fraction of source rows surviving.

## 4. Same-column arithmetic and transport coupling

Retain the ordinary actual-period collision estimate and the same raw columns as in [the preceding same-column estimate](../pi-near-total-2026-11-02/COLUMN_COST_AND_RESIDUAL_DEGREE.md). Write normalized transverse means

\[
s_H=\frac{\sum_{L_H}w\cdot\beta}{M_HH},\quad
t_H=\frac{\sum_Jw\cdot a}{M_HH},\quad
u_H=\frac{\sum_{C_H}w\cdot\alpha}{M_HH}.
\]

The exact degree identity D=Gamma+E implies, for a nonzero joint term,

\[
s_H\le t_H\le u_H\le1,\qquad
\frac{w\cdot\Gamma}{M_HH}=t_H-s_H,
\quad \frac{w\cdot D}{M_HH}=u_H-s_H.
\tag{10}
\]

Let d_i/w_i<=1+delta and ell_i/w_i>=2-epsilon, where d_i is the arithmetic denominator exponent and ell_i=-ln|zeta_i|. With R>7K and c_R=ln(R/(7K))/2, the established residual-degree bound yields

\[
\frac{\ln|\det T_{L,J}\det B_{J,C}|+\sum_i d_iD_i}{M_HH}
\le G_H+(1+\delta)(u_H-s_H)
-(2-\epsilon)(t_H-s_H)
-c_R\frac{\sum_an_a^2}{M_HH}+o(1),
\tag{11}
\]

where one may use the uniform bound

\[
0\le G_H\le\frac R{w_0}
+\frac{2\ln2+\ln R}{w_*},\qquad w_* =\min_iw_i.
\tag{12}
\]

All factors in (11) refer to the same columns. The arithmetic exponent has not been replaced by a source-independent norm. Moreover, the distinct-column theorem proves u_H->1 uniformly after the subsequent large-m limit at any fixed positive column density. Thus replacing u_H by one loses no leading constant in that limit.

The exponential-column part of G_H can also retain its actual mean. Summing h over all raw columns and dividing by the selected fraction gives

\[
\limsup_{H\to\infty}\frac{R\sum_{C_H}h}{M_HH}
\le\frac{R}{\lambda w_0(m+2)}.
\tag{12a}
\]

For a lower comparison, the fraction of all raw columns with w_0h/H<=x tends to 1-(1-x)^(m+1)<=(m+1)x. Set x=lambda/[2(m+1)]. At least half of any selected lambda fraction has cost at least x, so

\[
\liminf_{H\to\infty}\frac{w_0\sum_{C_H}h}{M_HH}
\ge\frac{\lambda}{4(m+1)}.
\tag{12b}
\]

These bounds apply uniformly to every selected raw-column set.

## 5. Critical family and a new mean-degree concentration

Fix 0<theta<1, 0<B<1 and C>1. For each m set

\[
\lambda=\tfrac12,\quad w_0=B^{-m},\quad
K=\lfloor C^m\rfloor,\quad
v_0=2K\theta^m w_0,\quad R=14K.
\tag{13}
\]

For each fixed member choose transverse weights sufficiently large, and then take H to infinity. Subsequent m limits below mean the limits of those fixed-parameter H statements. Here c_R=ln2/2 and s_H->s_m=m theta/(m+2). The collision saving forced by a mean bound t_H<=t_0<=a_m is

\[
\boxed{P_m(t_0)=
\frac{c_R}{m+2}
\left(\frac{B m}{(m+2)t_0}\right)^m.}
\tag{14}
\]

For fixed t_0<B this tends to infinity; for fixed t_0>B it tends to zero. Thus the moment condition produces an exponential threshold at B.

For a quantitative shrinking threshold, choose any fixed T>0 and put

\[
\tau_m(T)=B\frac m{m+2}
\left(\frac{c_R}{T(m+2)}\right)^{1/m}.
\tag{15}
\]

For sufficiently large m this is in (0,a_m), P_m(tau_m)=T, and

\[
\tau_m(T)=B-O(\log m/m).
\]

Assume the same nonzero rational minor and arithmetic lower bound as in the preceding notes. If the analytic growth, binomial factors and that arithmetic cost have total normalized upper bound strictly less than T, the sum over all targets with t_H<=tau_m(T) is exponentially smaller than the arithmetic lower bound. The number of possible target sets has logarithm o(M_HH), so the assertion follows by summing the absolute-value estimate. Consequently the complementary grouped sum with t_H>tau_m(T) must supply the arithmetic lower bound, up to an exponentially small relative error. This statement concerns the grouped exterior sum, not dominance of a single term.

The same argument can be combined with the prior surviving-source concentration: subtract both negligible sums and the remaining contribution lies in the intersection of their complementary regions.

For B=3/5, theta=1/2, C=13/10, the following exact rational thresholds give collision saving greater than two, using c_R>1/3:

| m | Threshold t_0 | Certified rational lower bound for saving |
|---:|---:|---|
| 40 | 497/1000 | (1/[3(m+2)]) [3m/(5(m+2)t_0)]^m > 2 |
| 80 | 541/1000 | same expression > 2 |
| 120 | 558/1000 | same expression > 2 |
| 200 | 573/1000 | same expression > 2 |

The verifier also checks that, at these four dimensions, the preceding arithmetic/analytic upper budget can be made less than one by taking transverse weights large. Thus under the nonzero-minor hypotheses the table gives actual sufficient concentration thresholds for this determinant construction. It does not describe a distribution of rational approximations to pi.

## 6. The sharp relaxation still has a positive endpoint gap

Ignore errors tending to zero in (11), and use the exact source mean s_m. The sharp mass/mean collision relaxation has the upper-envelope expression

\[
F_m(t)=1+s_m-2t
-\frac{c_R}{2(m+1)w_0}I_m(t),
\qquad s_m\le t\le1.
\tag{16}
\]

For fixed theta<B<1, its maximum has the limit

\[
\boxed{\lim_{m\to\infty}\sup_{s_m\le t\le1}F_m(t)
=1+\theta-2B.}
\tag{17}
\]

For a direct proof, on every fixed interval t<=B-epsilon the collision term tends uniformly to infinity. Above B-epsilon, discard its nonnegative contribution to bound F_m by 1+s_m-2(B-epsilon). Conversely, choose any fixed t=B+epsilon<1; the collision term tends to zero, giving the matching lower limit after epsilon decreases to zero. The same argument for B<=theta gives limit 1-theta, with the maximum at the source boundary in the limiting problem.

In the low-mean branch the exact interior critical point is

\[
t_m^*=\left[
\frac{m c_R}{2(m+2)}
\left(\frac{Bm}{m+2}\right)^m
\right]^{1/(m+1)}\longrightarrow B.
\]

This describes the optimizer of the relaxed upper estimate, not an actual nonzero period minor. Its positivity proves only that these inequalities do not yield a uniform negative bound.

There is a feasible witness for the retained relaxation, so this is not merely a failure to derive a stronger lower bound. Fix theta<B and choose t=B+epsilon<1. For each sufficiently large m, use the integer profiles from Section 3 that approach the density (5) with mean t. Select M_H raw columns by the cheapest-transverse-degree cutoff from the preceding note, whose mean u_m=q_m(1/2) tends to one. For all sufficiently large m, s_m<t<u_m. Radial symmetry in the rescaled coordinates w_i a_i gives the limiting coordinate identities

\[
\frac{w_i\Gamma_i}{M_HH}\longrightarrow\frac{t-s_m}{m}>0,
\qquad
\frac{w_iE_i}{M_HH}\longrightarrow\frac{u_m-t}{m}>0.
\]

Thus the exact integer profiles satisfy Gamma_i>=0 and E_i>=0 at sufficiently large H, while retaining the same selected columns and D=Gamma+E. Their collision cost tends to zero in the later m limit because t>B. Their arithmetic degree cost tends to 1-theta. The right-hand expression in (11), evaluated with the endpoint qualities d_i/w_i->1 and ell_i/w_i->2 and its stated Cauchy factors, therefore tends to 1+theta-2(B+epsilon). Letting epsilon decrease proves the lower direction in (17) with feasible mass, degree, and selected-column data. The witness still omits exact time-grade matching and determinant nonvanishing, and supplies no actual approximation tuple; those omitted constraints are possible sources of improvement. If B<=theta, the same construction with t=theta+epsilon yields the limit 1-theta.

For the current normalization and this particular Cauchy estimate, suppressing the analytic column-growth term is possible at CB=1 by (12a), while CB>1 makes it diverge by (12b). Thus the sharp fixed-exponential parameter condition is

\[
CB\le1.
\tag{18}
\]

Controlling the arithmetic qualities kappa_i=d_i/w_i near one by taking large weights requires v_0->infinity, hence

\[
C\theta>B.
\tag{19}
\]

Section 7 supplies an elementary lower bound showing that this requirement is necessary for the specific LCM denominator exponent used here, not merely sufficient for its previous upper estimate. Conditions (18)--(19) imply

\[
\boxed{B^2<\theta.}
\tag{20}
\]

Yet crossing the arithmetic/transport endpoint in (17) would require B>(1+theta)/2. The arithmetic-geometric mean inequality gives (1+theta)/2>=sqrt(theta), so the requirements are incompatible. For any fixed 0<theta<1, optimizing B and C within (18)--(19) leaves the strictly positive infimum

\[
\boxed{(1-\sqrt\theta)^2.}
\tag{21}
\]

At theta=1/2 this is 3/2-sqrt(2), between 0.0857 and 0.0858. The existing B=3/5 family leaves the larger gap 3/10. Varying fixed theta toward one makes (21) tend to zero, but no fixed theta<1 closes the relaxation. A joint choice theta=theta_m->1 is not ruled out by this fixed-parameter result: lower-order arithmetic, degree, and analytic terms would have to be controlled anew.

These statements preserve the available column-density and arithmetic normalization. They do not exclude improved period-rank information, cancellation, tighter reduced denominators, or a different analytic estimate.

## 7. Why bounded v_0 cannot make this LCM cost disappear

Use the existing normalization w_i=ceil(ln q_i), with |b_i| bounded above and bounded away from zero, and

\[
L_i=\operatorname{lcm}(1,\ldots,N_i),\qquad
N_i=\left\lceil\frac{w_i}{\theta v_0}\right\rceil-1,
\quad d_i=\ln(q_iL_i|b_i|).
\]

For every N>=2, put r=floor(N/2). The central binomial coefficient binom(2r,r) divides lcm(1,...,N). Indeed, for each prime p, its valuation is a sum of terms floor(2r/p^j)-2floor(r/p^j), each equal to zero or one, and there are at most floor(log_p N) nonzero terms. Therefore

\[
\ln\operatorname{lcm}(1,\ldots,N)
\ge\ln\binom{2r}{r}
\ge(N-1)\ln2-\ln(N+1).
\tag{22}
\]

For each fixed m, theta and v_0, as the weights tend to infinity this gives

\[
\boxed{\liminf_{w_i\to\infty}\kappa_i
\ge1+\frac{\ln2}{\theta v_0}.}
\tag{23}
\]

Thus a later fixed-theta sequence with qualities approaching one needs v_0->infinity. This is a statement about the displayed proved denominator exponent; cancellation or a new denominator argument could change it.

## Verification and scope

`verify_moment_collision.py` uses exact rational arithmetic to check the continuum optimizer integrals and branch matching, finite lattice quadratic-program certificates against exhaustive integer profiles, finite packet normalizations, the rational concentration thresholds and budget conditions, and central-binomial divisibility in the LCM argument. It also checks the parameter incompatibility on a rational grid. Infinite limits and the general parameter obstruction follow from the proofs above, not the finite examples.
