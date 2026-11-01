# Exact diagonal profiles, time capacity, and the genuine weighted-column kernel

2026-11-01. This note isolates the compatibility that is lost when every transverse coefficient is treated as an independent analytic function. It supplies exact coordinate transformations, an exact time-capacity formula, and a weighted-column Gram bound. No asymptotic improvement of the endpoint comparison is claimed.

Let w_0,w_1,...,w_m be positive weights, and let V_H be the span of raw monomials Y^h X^alpha with w_0 h+w·alpha<=H. At a common formal period Z the ordinary coefficient rows are

\[
B_{j,s,\beta;h,\alpha}(Z)
=\binom\alpha\beta[z^s]e^{hz}(jZ+z)^{|\alpha|-|\beta|}.
\tag{1}
\]

An entry is zero unless alpha>=beta coordinatewise. At Z=2*pi*i these are the genuine period rows, because the center exponential is one. At other Z, (1) is a polynomial specialization, not an assertion that e^Z=1. Time rows multiplied by s! have the same rank.

## 1. A budget-exact column decomposition

Order weights w_1<=...<=w_m. Put T=X_1 and U_i=X_i-X_1 for i=2,...,m, assigning T weight w_1 and U_i weight w_i. Both directions of this coordinate change preserve weighted degree: X_i=T+U_i and U_i=X_i-T involve only terms of weight at most w_i. Consequently

\[
\boxed{V_H=\bigoplus_{\nu\in\mathbf Z_{\ge0}^{m-1}}
U^\nu\,\operatorname{span}\{Y^hT^a:
w_0h+w_1a\le H-\sum_{i=2}^m w_i\nu_i\}.}
\tag{2}
\]

Only nonempty summands are included. This is equality of the full weighted spaces, not a rectangular enlargement of the column budget. The basis change from Y^h T^a U^nu to raw monomials has integer entries and determinant one: its leading term is Y^h X_1^a product X_i^(nu_i), and every other term has strictly lower total U degree when indexed by raw exponents in coordinates 2,...,m.

For a new column Y^hT^aU^nu, the old row beta is

\[
\boxed{
B_{j,s,\beta}(Y^hT^aU^\nu)
=(-1)^{|\nu|-|\beta'|}\binom\nu{\beta'}\binom a\ell
[z^s]e^{hz}(jZ+z)^{a-\ell},
\quad \ell=|\beta|-|\nu|,}
\tag{3}
\]

provided nu>=beta'=(beta_2,...,beta_m), ell>=0 and a>=ell; otherwise it is zero. To prove (3), substitute T=jZ+z+u_1 and U_i=u_i-u_1. Extracting u_i^(beta_i) from U_i^(nu_i) leaves (-u_1)^(nu_i-beta_i), so the remaining tangent derivative is exactly ell.

Thus the actual row matrix can be organized in normal profiles nu, with explicitly coupled old rows. It is not valid to allocate each old row to one profile without accounting for the sum in (3) and possible determinant cancellation.

## 2. The two natural anchors preserve different filtrations

The least-weight anchor preserves the column space in (2), but generally does not preserve a source derivative packet. Indeed the new tangent derivative is partial_T=sum_i partial_Xi. A permitted cheap derivative can therefore involve an old derivative whose weight is too large. For example, with transverse derivative weights (1,3), a cutoff between 1 and 3 contains partial_X1 but excludes partial_X2, whereas partial_T contains both. This statement concerns the coordinate change at fixed center and fixed time order; it does not silently add other time tests.

There is a complementary exact statement. Choose the greatest-weight anchor T=X_m and U_i=X_i-X_m, i<m. For any complete source packet

\[
v_0s+\sum_i v_i\beta_i<H,
\qquad v_1\le\cdots\le v_m,
\]

the normalized transverse-derivative space is invariant under this coordinate change, at each center and time order. Explicitly,

\[
\frac{\partial_X^\beta}{\beta!}
=\sum_{t+|k|=\beta_m}
(-1)^{|k|}\prod_{i<m}\binom{\beta_i+k_i}{\beta_i}
\frac{\partial_T^t\partial_U^{\beta'+k}}
{t!(\beta'+k)!}.
\tag{4}
\]

Every new index has weight

\[
v_mt+\sum_{i<m}v_i(\beta_i+k_i)
=v\beta-\sum_{i<m}(v_m-v_i)k_i\le v\beta.
\]

The inverse expansion comes from partial_T=sum_i partial_Xi and partial_Ui=partial_Xi. It obeys the same weight inequality. Ordering by the highest-coordinate derivative gives an integer unipotent transformation, hence determinant one. The factorial normalization in (4) is essential.

For this greatest-weight anchor the nominal weighted column space need not be preserved: a cheap X_i becomes U_i+T, but T has the larger weight w_m. These two anchor-coordinate choices therefore cannot in general be used to preserve both filtrations at once for unequal weights. This is not an impossibility theorem for every conceivable adapted coordinate transformation.

## 3. Exact capacity of one transverse type

Fix an original type beta with w·beta<=H, and put L_beta=H-w·beta and w_*=min_i w_i. Its diagonal image is exactly

\[
\mathcal E_\beta
=\operatorname{span}\{e^{hz}z^d:
h,d\ge0,\ w_0h+w_*d\le L_\beta\}.
\tag{5}
\]

Inclusion follows from |alpha-beta|<= (H-w·beta-w_0h)/w_*. Equality follows by taking alpha=beta+d e_i at a least-weight coordinate. The coefficient binomial(alpha,beta) is nonzero and does not affect the image space.

Write d_h=floor((L_beta-w_0h)/w_*), for 0<=h<=floor(L_beta/w_0), and N_beta=sum_h(d_h+1). Distinct exponential-polynomial blocks are independent, so dim E_beta=N_beta. They are annihilated by the exact constant-coefficient differential operator

\[
Q_\beta(D)=\prod_h(D-h)^{d_h+1}.
\tag{6}
\]

Here D is differentiation in z. If all selected rows of this fixed type have time orders s>=r, their rank on the entire V_H is at most

\[
\boxed{N_\beta-\min(r,d_0+1).}
\tag{7}
\]

Indeed D^r is invertible on every h>0 exponential-polynomial block, while on the h=0 polynomial block its kernel consists of degrees below r. Applying any further time derivatives and evaluations cannot restore this kernel. For a formal Z that is not a true period, the row formula omits the factor e^(hjZ); this is a separate nonzero rescaling within each h-block and commutes with the translation and differentiation used here. The kernel-count argument therefore still gives (7) for the formal rows.

The bound is sharp as a capacity in the unrestricted family of time tests: at the zero center, the consecutive orders r,...,r+C-1, where C is (7), have rank C. The image D^r E_beta is again an exponential-polynomial solution space of an order-C constant-coefficient ODE, and its first C Taylor data determine it uniquely. These orders may exceed the source or transport time cap; this sharpness assertion does not guarantee that an admissible target selection can attain the capacity.

Equation (6) also gives exact linear recurrences among the factorial-scaled time rows at every fixed center. If Q_beta(x)=sum_q a_q x^q, then sum_q a_q R_{j,s+q,beta}=0. In ordinary coefficient normalization the factors (s+q)! must be retained.

For the previous gluing obstruction H=5,w_0=1,w_1=2,beta=2, the image is span{1,e^z}. Thus N_beta=2 and r=1 gives capacity one: the first and second derivative rows cannot be independent. This explains the exact dependence without a numerical rank calculation.

This capacity constraint alone is not an asymptotic solution. If H-w·beta>=epsilon H for fixed epsilon>0, then N_beta is of order H^2, while the usual fixed-parameter row counts per type are of order H. The constraint becomes restrictive near the column-budget edge; further relations between different types can still matter elsewhere.

## 4. A kernel that keeps the entire weighted column budget

Clear rational weight denominators so that all w_i and H are integers. Define the actual full-column kernel

\[
K_H((Y,X),(V,W))
=\sum_{w_0h+w\alpha\le H}
(Y\overline V)^h\prod_i(X_i\overline W_i)^{\alpha_i}.
\]

Its cumulative generating function is exactly

\[
\boxed{\sum_{H\ge0}K_Ht^H
=\frac1{1-t}\frac1{1-t^{w_0}Y\overline V}
\prod_i\frac1{1-t^{w_i}X_i\overline W_i}.}
\tag{8}
\]

Pull back Y=e^z, X_i=j omega+z+u_i, and the conjugate variables at another center. Extract the requested time and transverse coefficients. This gives the full raw-column Gram matrix G_H[J]=B_J B_J^* on any row selection J. All correlations coming from shared monomial columns and their common weighted budget remain present.

There is a useful exact upper bound. Choose 0<t<1 so that t^(w_i)|(K-1)omega|^2<1 for every transverse coordinate. Let G_infinity(t)[J] be the jet Gram matrix of the product kernel in (8) with the cumulative factor 1/(1-t) removed. Then

\[
\boxed{G_H[J]\preceq t^{-H}G_\infty(t)[J],
\qquad
|\det B_{J,C}|^2\le\det G_H[J]
\le t^{-M H}\det G_\infty(t)[J],}
\tag{9}
\]

for every family C of M=|J| permitted raw columns. The first inequality is termwise positive semidefinite comparison: each column of cost c<=H receives coefficient t^(c-H)>=1, and all additional columns add positive semidefinite terms. Strict convergence at the centers also ensures convergence of every finite jet. Cauchy--Binet gives the middle bound, and determinant monotonicity on positive semidefinite matrices gives the last bound.

Equation (9) is a finite inequality for any t satisfying its convergence conditions, including an H-dependent choice. It supplies a genuine weighted-column alternative to replacing all coefficient functions by a generic norm. Its determinant is still a nontrivial object; no estimate uniform in H or exit fraction, and no improvement from optimizing t, has been proved here. The cost -log(t)/2 after MH normalization must be included, and one cannot replace this coupled kernel by the earlier direct-sum zero-error kernel without a new argument. General column mixing or arithmetic column scaling must likewise be charged separately.

## Verification

`verify_diagonal_profiles.py` uses exact rational arithmetic to check the budget-preserving column transformation and row formula, greatest-weight derivative transformations and their inverses, the sharp fixed-type capacities, the cumulative kernel generating function, and a finite Loewner/determinant instance of (9). The finite checks do not assert a new high-dimensional interpolation theorem or a Diophantine consequence.
