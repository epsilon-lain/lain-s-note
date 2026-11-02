# A weighted center-coalescence cutoff for the complete aggregate

2026-11-02. This is a finite support theorem for the **complete signed aggregate**, followed by its fixed-parameter height asymptotics. It does not bound the remaining coefficients or imply that any individually selected target minor vanishes. It strengthens the available information beyond a scalar mean/collision estimate without deciding pi's bad approximability.

## 1. Fixed objects and the cheapest spatial weights

Fix integers m>=1 and K>=1, height H>0, positive rational weights w_0,w_1,...,w_m,v_0, and 0<theta<=1. The complete source is

\[
L_H=\{(j,s,\beta):j\in\{0,\ldots,K-1\},\ s\in\mathbb Z_{\ge0},\ \beta\in\mathbb Z_{\ge0}^m,
\ v_0s+w\cdot\beta/\theta<H\},\qquad M_H=|L_H|.
\]

Columns are fixed polynomials, independent of the formal center vector r, in the space spanned by raw monomials Y^hX^alpha with w_0h+w dot alpha<=H. Their square source determinant is

\[
P_C(r)=\det\bigl([z^su^\beta]P_c(e^z,jr+z\mathbf1+u)\bigr)_{L_H,C}.
\tag{1}
\]

The formal removal of the center exponential is part of this definition. Only r=omega times the vector of ones, omega=2pi i, is being identified with the actual common-period matrix. An arbitrary rational r used for a polynomial check is not asserted to be an exponential period.

For each nonempty time layer define

\[
B_s=\{\beta:w\cdot\beta<\theta(H-v_0s)\},\quad
b_s=|B_s|,\quad S_s=\sum_{\beta\in B_s}w\cdot\beta.
\]

Let C_w(n) be the sum of the n smallest numbers in the multiset

\[
\{w\cdot\gamma:\gamma\in\mathbb Z_{\ge0}^m\},\qquad C_w(0)=0.
\]

Distinct exponent vectors count separately even when their weights coincide. Any ordering within a tie gives the same sum. Positive weights ensure only finitely many vectors below a bounded weight. The complete strict-cutoff packet B_s is an initial weight segment, so C_w(b_s)=S_s exactly.

Put

\[
\boxed{L_H^{(w)}=\sum_{s:b_s>0}\bigl[C_w(Kb_s)-K S_s\bigr]\ge0.}
\tag{2}
\]

Nonnegativity follows because the mean of the first n weights is nondecreasing; hence C_w(Kb)>=K C_w(b). No generic-position assumption enters this definition.

## 2. Every center monomial pays this weight

For each column define the row vector of zero-center jets

\[
\ell_{s,\gamma;c}=[z^su^\gamma]P_c(e^z,z\mathbf1+u).
\]

Translation in the spatial variables gives the exact finite identity

\[
\boxed{\mathcal R_{j,s,\beta}(r)
=\sum_{\gamma\ge\beta}\binom\gamma\beta
j^{|\gamma-\beta|}r^{\gamma-\beta}\ell_{s,\gamma}.}
\tag{3}
\]

The sum is finite because the columns are polynomials in X. The convention 0^0=1 handles j=0. Time differentiation does not change the expansion: the center translation is independent of z, so the same vector ell_(s,gamma) occurs for every center at fixed s.

Take the exterior product of the Kb_s rows at that time. A surviving summand must use Kb_s **distinct gamma vectors**: a repeated gamma repeats the same row vector. Its center monomial has weighted degree

\[
\sum_{\text{chosen }\gamma}w\cdot\gamma-K S_s
\ge C_w(Kb_s)-K S_s.
\]

Wedge multiplication over the time layers adds these degrees. Dependence between different layers can annihilate additional terms; it cannot create a term of smaller weight. Thus

\[
\boxed{[r^\eta]P_C(r)\ne0\quad\Longrightarrow\quad
w\cdot\eta\ge L_H^{(w)}.}
\tag{4}
\]

This is a lower bound on each surviving center monomial's weight, not merely the order after a radial specialization. It remains valid for fixed mixed polynomial columns. It is a lower bound on support and does not assert positivity of coefficients or a lower bound for P_C at any nonzero center.

## 3. Consequence for the complete shear aggregate

For the rest of this section each column has one spatial exponent:

\[
P_c(Y,X)=g_c(Y)X^{\alpha_c},\qquad
D_i=\sum_c\alpha_{c,i}-\sum_{L_H}\beta_i\ge0.
\]

If some D_i is negative, the source determinant is identically zero. In the nonzero case the earlier reciprocal identity defines the ordinary polynomial

\[
Q_C(x)=x^D P_C(x_1^{-1},\ldots,x_m^{-1}),
\]

and the **entire signed** normalized shear aggregate is

\[
\mathcal F_C(\zeta)
=\omega^{|D|}Q_C((1-\zeta_1)/\omega,\ldots,(1-\zeta_m)/\omega).
\tag{5}
\]

Coefficient extraction bounds each center exponent eta by D coordinatewise. A monomial r^eta becomes x^(D-eta) in Q; (4) therefore gives weight at most w dot D minus L_H^(w). Substituting x_i=(1-zeta_i)/omega can only decrease each exponent. Consequently

\[
\boxed{[\zeta^\Gamma]\mathcal F_C=0
\quad\text{whenever}\quad
w\cdot\Gamma>w\cdot D-L_H^{(w)}.}
\tag{6}
\]

This is an exact identity before an asymptotic limit. Contributions from individual targets having the same Gamma can be nonzero and must then cancel. It does not say those targets lack a matching or have zero period determinant. It also cannot be obtained by separately taking absolute values of all terms.

The time-filtration and ordinary center-coalescence cutoffs remain available simultaneously. Equation (6) supplies a weighted constraint that need not become negligible when w_max/w_min is very large.

For arbitrary mixed columns (4) still holds. Formula (6), with a single prefactor D, is asserted only in the single-spatial-exponent setting where (5) has been proved.

## 4. An exact two-variable equality case

Take K=2,H=3,w_0=1,w=(1,2),theta=1/2,v_0=2. At each center the complete source consists of

\[
(s,\beta)=(0,(0,0)),\ (0,(1,0)),\ (1,(0,0)).
\]

Use the six allowed raw columns

\[
1,\ Y,\ X_1,\ X_2,\ X_1^2,\ X_1X_2.
\]

With center-major rows in the displayed order, exact calculation gives

\[
P_C(r_1,r_2)=-2r_1^2r_2.
\]

Here b_0=2,b_1=1,S_0=1,S_1=0. The four cheapest weights are 0,1,2,2; hence C_w(4)=5 and C_w(2)=1. Thus

\[
L_H^{(w)}=(5-2)+1=4=w\cdot(2,1).
\]

Also D=(2,2), w dot D=6, Q_C(x)=-2x_2, and

\[
\mathcal F_C(\zeta)=-2\omega^3(1-\zeta_2).
\]

The weighted aggregate degree is exactly 2=6-4. This verifies that the finite bound can be attained with nontrivial spatial and time derivatives on the same raw columns. It is not an example in the critical high-dimensional regime.

For m=1, C_w(n)=w_1 n(n-1)/2. Equation (2) reduces exactly to w_1 binom(K,2) sum_s b_s^2, the previously derived one-dimensional confluent order. Thus the new weight-count formulation is consistent with that sharp case.

## 5. Height asymptotics and the remaining high-dimensional window

Keep m,K,theta and all weights fixed before H tends to infinity. Weighted lattice counting gives

\[
\#\{\gamma:w\cdot\gamma\le T\}
\sim\frac{T^m}{m!\prod_iw_i},
\quad
C_w(n)\sim\frac{m}{m+1}(m!\prod_iw_i)^{1/m}n^{1+1/m}.
\]

At a time slice with positive rescaled remaining height, B_s is a weighted simplex. Replacing b_s by Kb_s increases its leading cutoff by K^(1/m). Applying Riemann sums in s, with the boundary layers lower order, yields

\[
\boxed{\frac{L_H^{(w)}}{M_HH}\longrightarrow
\frac{m\theta}{m+2}\bigl(K^{1/m}-1\bigr).}
\tag{7}
\]

One can justify the boundary uniformly by separating H-v_0s>=epsilon H and its complement: the lattice estimates are uniform on the first region at fixed weights, while the latter contributes O(epsilon^(m+2) H^(m+2))+o(H^(m+2)). Then let epsilon tend to zero.

A necessary condition for a nonzero source determinant follows from its upper weighted center degree. Write

\[
u_H=\frac{\sum_C w\cdot\alpha}{M_HH},\quad
s_H=\frac{\sum_{L_H}w\cdot\beta}{M_HH}
\longrightarrow\frac{m\theta}{m+2}.
\]

Then u_H-s_H>=L_H^(w)/(M_HH). In particular any limiting column mean satisfies

\[
u\ge K^{1/m}\frac{m\theta}{m+2}.
\tag{8}
\]

More significantly for complete summation, a nonzero aggregate coefficient indexed by Gamma must obey (6). A target set contributing to that coefficient has weighted mean t_H=s_H+w dot Gamma/(M_HH), so

\[
\boxed{t_H\le u_H-\frac{L_H^{(w)}}{M_HH}.}
\tag{9}
\]

The word “nonzero” here concerns the **combined coefficient**. A target violating (9) can still contribute an individual term that cancels with other targets.

Along the previously studied subsequent dimension sequence theta=1/2, K=floor((13/10)^m), (7) tends to 3/20. Since u_H<=1 for allowed raw columns, the complete aggregate contains no surviving coefficients whose target mean stays above 17/20=0.85 after the height limits. The earlier collision estimate conditionally focused contributions near or above B=0.6; these results can be imposed together, leaving a smaller mean window rather than a contradiction.

In particular, the scalar witness near t=sqrt(1/2) is below 0.85. The new cutoff does not remove it or resolve the previous scalar obstruction. Coefficient sizes, actual target ranks on suitable columns, and cancellation inside the remaining weighted range are still uncontrolled. No simultaneous-limit error bound in H,m and the weights is asserted.

## Verification

The companion exact verifier checks finite cheapest-weight sums, the one-dimensional identity, the displayed complete six-row two-variable determinant, and its reciprocal aggregate. Fixed-weight lattice records illustrate (7); the limit is proved above, not inferred from finite experiments. No numerical value of pi is needed for the support theorem.
