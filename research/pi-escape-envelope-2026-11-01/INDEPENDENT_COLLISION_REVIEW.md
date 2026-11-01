# Independent review of the surviving-source occupancy law

2026-11-01. This review independently checks the proposed capped-occupancy minimization, its uniform asymptotic constants, the ordinary-jet collision estimate, and the conditional concentration of the exterior expansion near total escape. No error was found under the stated order of limits. These are transport and period-matrix estimates; they do not prove a Diophantine property of pi.

## 1. The finite minimization has the claimed direction

Let the full source packet be

\[
L_H=\{(j,s,\beta):0\le j<K,\quad v_0s+v\cdot\beta<H\},
\]

with all weights positive and fixed. Write M_H=|L_H|. Its multiplicity at type beta is

\[
n_\beta=K\left\lceil\frac{H-v\cdot\beta}{v_0}\right\rceil
\quad(v\cdot\beta<H).
\]

For any target set J with M_H distinct row labels, put k=|J\setminus L_H|. Then |J\cap L_H|=M_H-k. If r_beta counts the rows in this intersection of type beta, then

\[
0\le r_\beta\le n_\beta,\qquad
\sum_\beta r_\beta=M_H-k.
\]

The target's total multiplicity N_a(J) is at least r_a at every source type. Hence

\[
\sum_a N_a(J)^2\ge\sum_\beta r_\beta^2
\ge\min\left\{\sum_\beta u_\beta^2:
0\le u_\beta\le n_\beta,\quad\sum_\beta u_\beta=M_H-k\right\}.
\tag{1}
\]

The minimum in (1) allows real u_beta, so it is a lower bound for the integer problem. Its minimizer is

\[
u_\beta=\min(n_\beta,\lambda),
\]

where lambda is chosen to obtain the required sum. To check this without an optimization theorem, any two unsaturated entries with unequal values can be moved closer together while preserving their sum, decreasing their squared sum. Saturated entries must have caps no larger than the common unsaturated value. This yields precisely the displayed clipped allocation.

The resulting profile is sharp for this surviving-source minimization. It is not asserted sharp for the total target occupancy: the added outside rows, conserved-order restrictions, Pascal nonvanishing and actual period rank can all impose additional cost.

## 2. Continuum calculation and uniformity in the escape fraction

Let rho=k/M_H and scale multiplicities by A_H=KH/v_0. The number of source types is asymptotic to H^m/(m! product v_i), and

\[
\frac{n_\beta}{A_H}=1-\frac{v\cdot\beta}{H}+O(H^{-1}).
\]

Under normalized counting of source types, y=1-v dot beta/H has limiting density m(1-y)^(m-1) on [0,1]. In particular its mean is 1/(m+1). If a=lambda/A_H, the retained fraction tends to

\[
\frac{\mathbb E\min(y,a)}{\mathbb E y}
=1-(1-a)^{m+1}.
\]

Thus the clipping level corresponding to escape fraction rho is

\[
a=1-\rho^{1/(m+1)}.
\]

The second clipped moment is

\[
\mathbb E\min(y,a)^2
=2\int_0^a t(1-t)^m\,dt.
\]

Since M_H is asymptotic to K H^(m+1)/((m+1)! v_0 product v_i), substitution gives

\[
\boxed{
\frac1{M_HH}\sum_a N_a(J)^2
\ge\frac{2K}{v_0(m+2)}F_m(\rho)-o(1),
}
\tag{2}
\]

where

\[
F_m(\rho)=1-(m+2)\rho+(m+1)\rho^{(m+2)/(m+1)}.
\]

The o(1) may be chosen uniformly for all target sets and all rho in [0,1], with dimension, weights and K fixed. Here is a justification that also covers the endpoint fractions. The functions min(y,a) and min(y,a)^2 form uniformly bounded, equicontinuous families on a common compact interval containing the normalized finite multiplicities. Lattice Riemann sums therefore converge uniformly in a. The limiting first-moment function 1-(1-a)^(m+1) is continuous and strictly increasing on [0,1], so its inverse is uniformly continuous, even though its derivative is unbounded at one endpoint. The clipped second moments are uniformly continuous as well. These facts give uniform convergence of the optimized value as a function of the retained fraction. A Lipschitz inverse is neither true at every endpoint nor needed.

The ceiling error and the lattice boundary error are lower order. No uniformity as the dimension or weights vary with H follows from this argument; the prescribed fixed-parameter H limit must be taken first.

## 3. Useful exact inequalities for the profile

One has F_m(1)=F_m'(1)=0, and, for 0<rho<=1,

\[
F_m''(\rho)=\frac{m+2}{m+1}\rho^{-m/(m+1)}
\ge\frac{m+2}{m+1}.
\]

Integrating from rho to 1 twice proves

\[
\boxed{F_m(\rho)\ge
\frac{m+2}{2(m+1)}(1-\rho)^2.}
\tag{3}
\]

Also F_m decreases on [0,1]. For each fixed rho>0,

\[
F_m(\rho)\longrightarrow1-\rho+\rho\ln\rho
\quad(m\to\infty).
\]

The endpoint value at rho=0 is one, agreeing with the continuous convention rho ln rho=0. These are properties of the scalar profile, separate from any simultaneous limit in the matrix problem.

## 4. Ordinary normalized jets remove the earlier extraction charge

For a raw column P=Y^h X^alpha, define the entire type function

\[
f_{a,P}(z)=\binom\alpha a e^{hz}z^{|\alpha|-|a|},
\]

with zero value if a is not coordinatewise at most alpha. At omega=2*pi*i, its ordinary local coefficient is

\[
[z^s]f_{a,P}(j\omega+z),
\]

which equals the actual period row, because e^(hj omega)=1. No analogous assertion is made for an arbitrary substituted value of omega.

Fix R>7K and put r=7K/R. Expand f=sum_d c_d(z/R)^d. Cauchy's estimate on the local circle |z|=1 gives

\[
\left|[z^s]\left(\frac{j\omega+z}{R}\right)^d\right|
\le r^d,
\]

since |j omega+z|<=2*pi*(K-1)+1<7K. The radius contributes 1^(-s), so there is no factor 2^s. This applies to arbitrary nonnegative ordinary jet orders, including sparse selections. The earlier ln(2)/v_0 extraction charge belonged to logarithmic jets estimated on |t|=1/2; it cannot be carried into this ordinary-jet estimate as an unavoidable cost.

For columns with w_0 h+w dot alpha<=H, the uniform function bound is

\[
|f_{a,P}(z)|\le e^{H E_R}\quad(|z|\le R),
\qquad E_R=\max\{R/w_0,\ln(2R)/w_*\}.
\]

The established sum over distinct Taylor degrees then gives

\[
\frac{\ln|\det B_J|}{M_HH}
\le E_R-c_R\frac{\sum_aN_a(J)^2}{M_HH}+o(1),
\qquad c_R=\tfrac12\ln(R/(7K)).
\tag{4}
\]

Its remainder is uniform over target sets and raw-column selections: two factorial factors contribute O(ln M_H/H), and the geometric-product denominator contributes O(1/H). There is no need to assume the rows are independent; zero determinants only strengthen the upper bound. Combining (2) and (4) is therefore valid uniformly in rho.

## 5. Conditional concentration near total escape

For theta=1/2, B_0=3/5, C=13/10, let

\[
K=\lfloor C^m\rfloor,\qquad w_0=B_0^{-m},\qquad
v_0=2K\theta^m w_0,\qquad R=14K.
\]

The collision coefficient in (2)--(4) becomes exactly

\[
S_m=\frac{2c_RK}{v_0(m+2)}
=\frac{c_R}{m+2}(6/5)^m,
\qquad c_R=\frac{\ln2}{2}.
\]

Put epsilon_m^2=12(m+1)(5/6)^m, and restrict to dimensions large enough that epsilon_m<1. For every rho<=1-epsilon_m, (3) yields

\[
S_m F_m(\rho)\ge6c_R=3\ln2>2.
\tag{5}
\]

The strict inequality follows, for example, from ln2>2/3. The floor in K causes no error in this calculation because K cancels exactly against its factor in v_0.

Now use the same normalization on both sides of the exterior identity. Suppose its full sum D_H has an arithmetic lower bound |D_H|>=exp(-M_HH A_H), and the total normalized analytic growth, transport-column bound and limsup A_H are strictly below one. The transport bound alone supplies the uniform growth ln2/w_*; its additional error suppression may be discarded for this step. The number of exterior target sets has logarithm o(M_HH). Equations (4)--(5) then show that the sum restricted to rho<=1-epsilon_m is exponentially smaller than the arithmetic lower bound, with a fixed positive margin.

If D_low denotes that restricted sum and D_high the remaining sum, then for some delta>0 and all sufficiently large H,

\[
|D_{\rm low}|\le e^{-M_HH(A_H+\delta)},\qquad
|D_{\rm high}|\ge e^{-M_HH A_H}
\left(1-e^{-M_HH\delta}\right).
\]

This is a lower bound for the complete high-escape sum, obtained by the reverse triangle inequality. It does not say that its terms have common phase, that one prescribed term dominates, or that internal cancellation is absent. The arithmetic lower bound and existence of a nonzero rational source minor remain separate inputs; in the intended asymptotic application their original interpolation hypotheses have not been independently proved here.

The concentration narrows the unresolved sector to almost total escape in this parameter family. It does not eliminate that sector. In particular, the separate positive-fraction displacement law has only O(1/m) cost near rho=1 under its endpoint-quality assumptions, and cannot be treated as a finished contradiction.

## 6. Review of the completed note and its raw arithmetic constants

The completed `SURVIVING_COLLISION_LAW.md` was subsequently read, including its complete-sum inequality (9), concentration statement (10), and raw arithmetic bounds (12). They agree with the argument reviewed above.

For w_i=ceil(ln q_i) one has w_i-1<ln q_i<=w_i. The assumed interval 1/2<=|b_i|<=2 gives -ln2<=ln|b_i|<=ln2. If T_i=ceil(w_i/(theta v_0)), then T_i-1<w_i/(theta v_0), so the stated elementary least-common-multiple bound gives

\[
0\le\ln L_i\le\frac{4\ln2}{\theta v_0}w_i.
\]

Consequently the note's bounds

\[
1-2/w_*\le\kappa_i
\le1+4\ln2/(\theta v_0)+1/w_*
\]

hold. For sufficiently large transverse weights, kappa_max is positive. On a nonzero raw source minor, D_i=sum_columns alpha_i-sum_source beta_i, and the column budget gives

\[
\mathfrak A_H\le\kappa_{\max}
-\frac1{M_HH}\sum_{\rm source,i}w_i\beta_i\kappa_i.
\]

Each normalized source mean w_i beta_i/H tends to theta/(m+2). Substituting the uniform bounds on the kappas yields exactly equation (12). The interval assumption on b_i here is only an arithmetic estimate; it does not replace the strict 0<|zeta_i|<1 hypothesis of the transport estimate.

The rational finite-dimensional cost certificates are also correctly directed. They use ln2<25/36, so 4ln2/(theta v_0)<50/(9v_0) when theta=1/2. Replacing the holomorphic maximum by the sum and using ln(28K)<28K gives the script's larger cost bound. The extra coefficient (28K+2+m/(m+2))/w_* accounts for the logarithmic holomorphic term, ln2/w_*, and the row/column arithmetic errors. No cost is silently omitted.

The completed verifier was rerun successfully: 231 exhaustive survivor allocations, 39 optimization cases, 27 complete-packet cases, 3,575 ordinary-coefficient checks and four conditional high-dimensional rational parameter bounds. Its asymptotic comparisons are labelled as proved limits rather than inferred finite trends. These checks do not construct an approximation tuple or verify the manuscript's separated-weight interpolation threshold.
