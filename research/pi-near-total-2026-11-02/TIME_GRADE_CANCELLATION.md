# A same-column time-grade filter and exact boundary cancellation

2026-11-02. This note adds a necessary support condition to the exterior transport expansion. It produces exact vanishing of every nonidentity transport-supported term on a genuine Hermite staircase, and shows how diagonal Gram weights recover that structure. It also proves that the same time-saturation certificate is unavailable in the previously selected critical high-dimensional parameter family. No conclusion about pi's bad approximability follows.

## 1. Fix the columns before comparing transport and period minors

Use ordinary coefficient jets, meaning derivatives divided by their factorials. Let the complete source packet be

\[
L=\{(j,s,\beta):0\le j<K,\ v_0s+v\cdot\beta<H\},
\qquad v_i=w_i/\theta,
\]

and write M=|L|. Select M distinct shifted monomial columns

\[
P_{h,\alpha}(Y,X)=(Y-1)^hX^\alpha,
\qquad w_0h+w\cdot\alpha\le H.
\tag{1}
\]

At the true period omega=2*pi*i their target-row entries are

\[
B_{j,t,a;h,\alpha}
=[z^tu^a](e^z-1)^h
\prod_i(j\omega+z+u_i)^{\alpha_i}.
\tag{2}
\]

Two entrywise zeros will be essential:

\[
B_{j,t,a;h,\alpha}=0
\quad\hbox{unless }h\le t\hbox{ and }a_i\le\alpha_i\ \forall i.
\tag{3}
\]

The full shifted basis is an integral unit-triangular change of the full raw basis Y^h X^alpha, preserving the weighted column budget. A general partial raw-column subset need not be closed under this change. We are choosing the columns in (1), rather than silently replacing an arbitrary existing partial minor. The single transverse degree alpha of each shifted column is unchanged, so the coordinatewise arithmetic scaling still applies.

The complete shear transport matrix T has a nonzero entry from (j,s,beta) to (j,t,a) only when

\[
a\ge\beta\text{ coordinatewise},\qquad
t+|a|=s+|\beta|.
\tag{4}
\]

On that support its entry is binom(a,beta)(-zeta)^(a-beta). The inside block U=T_{L,L} is unit triangular. These are the same exact transport conventions as the previous boundary checkpoint.

For a target row set J with M rows and the fixed selected column set C, define

\[
\Gamma_i=\sum_Ja_i-\sum_L\beta_i,
\quad D_i=\sum_C\alpha_i-\sum_L\beta_i,
\quad \Delta_{\rm time}=\sum_Ls-\sum_Ch.
\tag{5}
\]

**Time-grade filter.** If both det T_{L,J} and det B_{J,C} are nonzero, then

\[
\boxed{0\le\Gamma_i\le D_i\quad\forall i,
\qquad |\Gamma|\le\Delta_{\rm time}.}
\tag{6}
\]

Proof: a nonzero determinant contains at least one nonzero permutation term. A nonzero transport permutation gives Gamma_i>=0 and, by (4),

\[
\sum_Jt=\sum_Ls-|\Gamma|.
\]

A nonzero period permutation gives sum_C h<=sum_J t and sum_J a_i<=sum_C alpha_i by (3). Combining the two proves (6). The two permutations need not coincide; only their sums are used.

Moreover det T_{L,J}!=0 and Gamma=0 force J=L. All coordinate increments in a nonzero transport permutation are nonnegative; if their total is zero, each increment is zero. The center and conserved total in (4) then force each matched row to be identical.

Consequently, if Delta_time=0, every nonidentity transport-supported target has zero period minor on these same columns. In Cauchy--Binet only J=L can remain. This is termwise rank vanishing, not a claim about cancellation of many nonzero scalar terms. Unrelated target sets with zero transport minor are not covered by the assertion about their period minors.

The time part of (6) extends to polynomial columns divisible by (Y-1)^h. The coordinatewise D statement requires a genuine single transverse degree or a separately justified multidegree bound; a general mixed column must not be assigned a fictitious alpha.

## 2. A complete Hermite staircase saturates the time budget

Consider m=1, K>=2, v_1=w_1/theta, with

\[
K\theta<1,\qquad w_0<v_0.
\tag{7}
\]

For each integer 0<=s<H/v_0 put

\[
d_s=\left\lceil\frac{H-v_0s}{v_1}\right\rceil.
\]

There are Kd_s source rows at time s, namely beta=0,...,d_s-1 at each center. Select precisely the Kd_s columns

\[
(Y-1)^sX^\alpha,\qquad 0\le\alpha<Kd_s.
\tag{8}
\]

Their budget is valid for all sufficiently large H. Indeed

\[
w_0s+w_1\alpha
<H\max(K\theta,w_0/v_0)+w_1(K-1),
\]

so it suffices that H>w_1(K-1)/(1-max(K theta,w_0/v_0)). Row and column multiplicities agree at every time, hence Delta_time=0.

With increasing-time block order the source matrix B is block lower triangular. Its time-s diagonal block is the ordinary-coefficient confluent Vandermonde

\[
\left[\binom\alpha\beta(j\omega)^{\alpha-\beta}\right]_{(j,\beta),\alpha}.
\]

Therefore, with center-major order inside a time block,

\[
\boxed{\det B=
\prod_s\prod_{0\le i<j<K}((j-i)\omega)^{d_s^2}\ne0.}
\tag{9}
\]

This staircase determinant formula was already established in `pi-nonzero-determinants-2026-11-01/support_rank_constraints.md`; the new point here is its time-grade consequence for every transported target and for the boundary operator. The degree difference is

\[
D=\frac{K(K-1)}2\sum_sd_s^2.
\]

Because d_s decreases with s, lowering h in a column of (8) stays in the raw staircase. Thus in this particular selection the integral shifted-to-raw column change is internal and has determinant one. The same formula with omega replaced by a nonzero rational center difference r proves the arithmetic source minor nonzero as well.

In the original rational/period normalization the prefactor b^D is retained:

\[
\det A(r)=b^D\det(B+RB_O),\qquad b=r/\omega.
\]

The boundary cancellation below does not discard b^D or the arithmetic denominator cost. Indeed (9) makes the identity det A(r)=b^D det B explicit.

## 3. The relative boundary operator is nilpotent

Let B_O be the ordinary period matrix on the outside target rows, and let R=U^{-1}T_{L,O} be the exact boundary transport. The finite outside set may be restricted to rows reached by transport. Set

\[
X=B_OB^{-1},\qquad N=RX.
\]

Write s for a source time, t for an outside time, and s' for another source time. Transport support and the triangular inverse of U give

\[
R_{s,t}=0\quad\hbox{unless }t<s.
\tag{10}
\]

The inverse of the block lower triangular B has (column-time h, source-time s') entry zero unless s'<=h. By (3), B_O has (outside-time t, column-time h) entry zero unless h<=t. Hence

\[
X_{t,s'}=0\quad\hbox{unless }s'\le t.
\tag{11}
\]

Combining (10)--(11), N has an entry only if s'<s. It is strictly block lower triangular and nilpotent, with nilpotence index at most the number of source time levels. In particular

\[
\boxed{\det(B+\eta RB_O)=\det B\quad\text{for every complex }\eta.}
\tag{12}
\]

The outside matrix and RB_O can be nonzero and have substantial rank. Nevertheless all positive-degree coefficients of this determinant's boundary marker vanish. Equation (6) gives the stronger explanation at the level of each transport-supported exterior target.

## 4. Support-aware diagonal Gram weights recover the exact bound

The previous escape-Gram checkpoint established, for arbitrary positive diagonal source and outside matrices E,F,

\[
|\det(I+RX)|^2
\le\det(I+\widetilde R\widetilde R^*)
       \det(I+\widetilde X^*\widetilde X),
\tag{13}
\]

where R_tilde=E^(1/2) R F^(-1/2) and X_tilde=F^(1/2) X E^(-1/2). Use the midpoint time grading

\[
E_{s,s}=\varepsilon^{4s},\qquad
F_{t,t}=\varepsilon^{4t+2},\qquad 0<\varepsilon<1.
\tag{14}
\]

For each potentially nonzero entry,

\[
\widetilde R_{s,t}=\varepsilon^{2(s-t)-1}R_{s,t},\qquad
\widetilde X_{t,s'}=\varepsilon^{2(t-s')+1}X_{t,s'}.
\]

Both exponents are positive integers, by (10)--(11). At fixed H both scaled matrices tend to zero as epsilon decreases to zero. The right side of (13) tends to 1, while the left side is exactly 1. Therefore the infimum of this valid diagonal product-Gram bound is exactly 1 in the staircase family.

The midpoint shift matters. Weights E_ss=epsilon^(2s), F_tt=epsilon^(2t) would leave some entries of X_tilde constant when t=s'; they do not suffice for this conclusion. This is a finite-matrix optimization statement, not a uniform bound for a prescribed epsilon as H or the dimension varies. The exact cancellation (12) does not require a limit.

## 5. Why this certificate fails in the chosen critical family

For the general complete source simplex, at fixed dimension and weights as H tends to infinity,

\[
\frac{\sum_Ls}{MH}\longrightarrow\frac1{(m+2)v_0}.
\]

Every admissible shifted column has h<=H/w_0. Therefore every M-column choice obeys

\[
\boxed{\frac{\Delta_{\rm time}}{MH}
\ge\frac1{(m+2)v_0}-\frac1{w_0}+o(1).}
\tag{15}
\]

If w_0>(m+2)v_0, time saturation is impossible for large H. The previously studied critical family has

\[
v_0/w_0=2K\theta^m,\quad \theta=1/2,\quad K=\lfloor(1.3)^m\rfloor,
\]

so (m+2)v_0/w_0 tends to zero. Equation (15) rules out this Delta_time=0 certificate there. The limit is taken separately for each fixed parameter member before comparing large dimensions; no uniform H cutoff is claimed.

Any polynomial in the allowed raw space has Y-degree at most floor(H/w_0), so changing to a different time-adapted polynomial basis cannot exceed that individual vanishing-order cap. This does not exclude a different cancellation mechanism with positive Delta_time, nor an estimate using more than the time grade. It identifies the concrete budget mismatch in this particular mechanism.

## 6. Reproducible exact checks

Run `python verify_time_grade.py`. The standard-library verifier uses K=2, w_0=1/4, w_1=1, v_0=1, v_1=4, theta=1/4 and exact rational arithmetic. Its formal center values Z=1,2 are polynomial specializations of (2); it never asserts that exp(Z)=1 at those values.

At H=4 it checks all 576 transport-supported target choices. Every transport determinant is nonzero, the source period minor is nonzero, and all 575 nonidentity period minors vanish. At H=8 it checks the complete 24-row staircase. Across the two examples it verifies 1248 boundary entries against the independent alternating-cutoff formula, the inverse supports, nonzero boundary operators of ranks 3 and 11, nilpotence indices 2 and 6, and the constant marker determinant. It also evaluates the midpoint Gram bound at three exact rational epsilons.

The finite checks detect normalization and support mistakes; the proofs above establish the infinite families. The independent audit is in `INDEPENDENT_TIME_GRADE_REVIEW.md` in the published folder.
