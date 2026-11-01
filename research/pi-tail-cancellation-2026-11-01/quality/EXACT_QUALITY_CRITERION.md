# Exact coordinate quality after eliminating logarithmic tails

2026-11-01. This is a conditional extension of the uploaded manuscript's determinant framework. It uses its separated-weight interpolation theorem to obtain a nonzero determinant. It does not establish the irrationality measure or bad approximability of π. The numerical example below is synthetic parameter arithmetic, not a list of approximations to π.

## 1. The inputs and why tails need not consume the quality gain

Fix the original full row packet

\[
\mathcal R_H=\{(j,s,\beta):0\le j<K,\quad v_0s+w\beta/\theta<H\},
\]

with \(K\ge2\), \(0<\theta<1\), and positive weights. As in the manuscript, \(w_i=\lceil\ln q_i\rceil\) for integers \(q_i\ge2\), and \(r_i=2ip_i/q_i\). All parameters are fixed before admissible \(H\to\infty\). Assume the original interpolation hypotheses.

Choose the integer truncations

\[
T_i=\left\lceil\frac{w_i}{\theta v_0}\right\rceil,
\qquad G_i(t)=\sum_{1\le k<T_i}\frac{(-1)^{k+1}t^k}{k},
\qquad L_i=\operatorname{lcm}(1,\ldots,T_i-1),
\]

where an empty least common multiple equals one. In this note \(\ln(1+t)\) means its formal Taylor series at zero.

There is an exact full-row determinant identity between the matrix using \(G_i\) and the matrix using \(\ln(1+t)\). Indeed, writing \(\tau_i=G_i-\ln(1+t)\), the former row \((j,s,\beta)\) is

\[
\sum_{d\ge0}\binom{\beta+d}{\beta}
\sum_{k\le s}[t^k]\tau(t)^d\,
\operatorname{row}_{\ln}(j,s-k,\beta+d).
\]

A nonzero coefficient has \(k\ge\sum_i d_iT_i\), and hence

\[
v_0(s-k)+w(\beta+d)/\theta
\le v_0s+w\beta/\theta<H.
\]

Thus all rows used remain in the same packet. Every nonidentity contribution has \(k>0\), and therefore a strictly smaller time order. Ordering by time order gives a unitriangular row transformation, with determinant one. Equality at \(v_0T_i=w_i/\theta\) is harmless: the original row cutoff is strict and time order still decreases.

This identity holds for every chosen common column family. It allows us to retain polynomial truncations for arithmetic integrality and to use only constant translations for the analytic determinant expansion. There is no tail index \(d\) and no tail-quality loss in that expansion.

## 2. Exact arithmetic and approximation coordinates

Define

\[
\lambda_i=-\ln\bigl(q_i^2|\pi-p_i/q_i|\bigr),
\quad \alpha_i=2\ln q_i+\lambda_i
=-\ln|\pi-p_i/q_i|,
\]

\[
\nu_i=\frac{\alpha_i}{w_i}>0,
\qquad u_i=\frac{\ln(q_iL_i)}{w_i}>0.
\]

Positivity of every \(\nu_i\) is an explicit hypothesis for the anisotropic simplex count. It holds, in particular, when every selected approximation satisfies \(|\pi-p_i/q_i|\le q_i^{-2}\). The symbol \(\nu_i\) is an exact quality-to-weight ratio; it is not a common exponent hypothesis.

Scale a monomial column \(Y^hX^\gamma\) by \(\prod_i(q_iL_i)^{\gamma_i}\) and a row \((j,s,\beta)\) by \(\prod_i(q_iL_i)^{-\beta_i}\). Every entry is then a Gaussian integer because

\[
q_iL_i\bigl(jr_i+G_i(t)\bigr)\in\mathbb Z[i][t].
\]

A nonzero minor consequently satisfies

\[
\frac{\ln|\Delta_H|}{M_HH}
\ge-u_{\max}+b_{u,H},
\quad
b_{u,H}=\frac1{M_HH}\sum_{\rho,i}u_iw_i\beta_{\rho,i}.
\tag{1}
\]

Here \(u_{\max}=\max_i u_i\), and the column bound used only \(w\gamma\le H\). The row-simplex moment is

\[
b_{u,H}\longrightarrow\frac{\theta}{m+2}\sum_i u_i.
\]

Put

\[
R_{\rm ar}=u_{\max}-\frac{\theta}{m+2}\sum_i u_i>0,
\qquad D_\nu=\frac{\theta}{m+2}\sum_i\nu_i.
\tag{2}
\]

The least-common-multiple cost is already contained in \(u_i\); it must not be charged again as a separate arithmetic error. For reference, the elementary bound \(\ln\operatorname{lcm}(1,\ldots,n)\le\Lambda n\), \(\Lambda=4\ln2\), gives

\[
u_i\le1+\frac{\Lambda}{\theta v_0},
\]

because \(T_i-1<w_i/(\theta v_0)\). This is a maximum-coordinate charge, not the sum of \(m\) separate column charges.

The common bounded-column-mixing extension has the same threshold: after pre-scaling every monomial column by \(\prod_i(q_iL_i)^{\gamma_i}\), the arithmetic lower bound is \(b_{u,H}\); the analytic bound pays at most \(u_{\max}\) for that scaling. Bounded mixing contributes only \(o(1)\) after division by \(M_HH\). Thus (2) remains the comparison threshold, provided the common-mixing lemma is explicitly included.

## 3. A scalar estimate retaining every quality coordinate

After the determinant identity in §1, a row is expanded only by constant shifts from \(jr_i\) to \(j\omega\), where \(\omega=2\pi i\). Its scalar for \(a\ge\beta\) is

\[
\xi_{j,\beta,a}
=\binom a\beta\prod_i[j(r_i-\omega)]^{a_i-\beta_i}.
\tag{3}
\]

At \(j=0\), only \(a=\beta\) survives and its scalar is one. For \(j\ge1\),

\[
|\xi_{j,\beta,a}|=
\binom a\beta(2j)^{|a-\beta|}
\exp\!\left(-\sum_i\nu_iw_i(a_i-\beta_i)\right).
\tag{4}
\]

Since \(\binom a\beta\le2^{|a|}\), \(wa\le H\), and \(w_*=\min_iw_i\),

\[
\ln|\xi_{j,\beta,a}|
\le-\sum_i\nu_iw_i(a_i-\beta_i)
+H E_{\rm tr},
\quad
E_{\rm tr}=\frac{\ln(4(K-1))}{w_*}.
\tag{5}
\]

No ceiling loss proportional to \(\nu_i\), maximum quality penalty, or omitted-tail term appears. In the same nonzero class the reverse inequality without the error is also valid, because the prefactor in (4) is at least one.

## 4. An anisotropic determinant comparison

Let a valid test-determinant estimate be

\[
\frac{\ln|\det B_H|}{M_HH}
\le-\frac{c\sum_a n_a^2}{M_HH}+E_{\rm hol}+o(1),
\tag{6}
\]

uniform in all test choices, with fixed \(c>0\). The manuscript supplies \(c=\ln2/4\) and its displayed \(E_{\rm hol}\); a separately proved improvement can be substituted here without changing the argument. Put \(E=E_{\rm tr}+E_{\rm hol}\).

Choose a real threshold \(A>0\) and \(0<\eta<1\). Partition test types by

\[
\sum_i\nu_iw_i a_i\le AH.
\]

Ignoring the additional domain restriction \(wa\le H\) only enlarges this low-type set. Its count is

\[
D_A(H)\sim\frac{(AH)^m}{m!\prod_i(\nu_iw_i)}.
\]

Together with the original row count, this yields

\[
L_\nu=
\frac{\eta^2K\theta^m}{(m+1)v_0 A^m}\prod_i\nu_i.
\tag{7}
\]

If at least \(\eta M_H\) rows have low types, the collision contribution is at least \(cL_\nu+o(1)\). Otherwise the scalar contribution in (5) is at least

\[
A(1-\eta)-D_\nu+o(1).
\]

The number of choices per row is polynomial in \(H\), so its normalized logarithm tends to zero. Consequently

\[
\limsup_{H\to\infty}\frac{\ln|\Delta_H|}{M_HH}
\le E+\max\{-cL_\nu,-A(1-\eta)+D_\nu\}.
\tag{8}
\]

Combining (1) and (8) gives the following explicit conditional criterion:

\[
\boxed{\quad cL_\nu>R_{\rm ar}+E,
\qquad A(1-\eta)>R_{\rm ar}+D_\nu+E.\quad}
\tag{9}
\]

If the actual approximation tuple and interpolation hypotheses satisfied (9), the resulting nonzero determinant would have incompatible lower and upper bounds. The inequalities are strict, which absorbs the fixed-parameter \(o(1)\) terms.

This is a sufficient contradiction criterion. Failure of (9) proves no approximation statement.

## 5. What nonuniform quality buys, and its limitation

With the arithmetic/background data fixed, \(D_\nu\) depends on the arithmetic mean

\[
\bar\nu=\frac1m\sum_i\nu_i,
\]

whereas \(L_\nu\) depends on the product. The arithmetic-geometric mean inequality gives

\[
\prod_i\nu_i\le\bar\nu^m.
\tag{10}
\]

Thus the displayed relaxed-cutoff criterion cannot benefit from heterogeneity at fixed mean: its high-type threshold is unchanged, and its collision quantity can only decrease compared with formal equal ratios \(\nu_i=\bar\nu\). This is a comparison of the criterion, not a claim that actual rational approximations can be changed to have prescribed equal quality. Intersecting the low-type simplex with \(wa\le H\), using row-specific constraints, or exploiting determinant cancellations may improve the criterion and is not excluded by (10).

A synthetic exact example makes the distinction concrete. Let

\[
m=500,\quad\theta=1/2,\quad B=17/25,\quad
K=\lfloor(139/100)^m\rfloor,
\quad w_0=B^{-m},\quad v_0=2K\theta^mw_0,
\]

and formally set \(u_i=1\), \(E=1/100\), \(\eta=1/100\), \(A=13/5\). For all \(\nu_i=4\),

\[
L_\nu=\frac{1}{10^4\cdot1002}\left(\frac{68}{65}\right)^{500}
\approx626.48447.
\]

For 250 ratios equal to 3 and 250 equal to 5, the mean remains 4 but

\[
L_\nu^{\rm heterogeneous}=L_\nu(15/16)^{250}
\approx0.0000616211.
\]

Both have the identical high-type margin

\[
A(1-\eta)-R_{\rm ar}-D_\nu-E
=\frac{8782}{125500}>0.
\]

Using only the rational enclosure \(14/81<c<25/144\) for \(c=\ln2/4\), the first passes the collision inequality and the second fails it. The verifier checks all inequalities by rational arithmetic. These \(u_i,\nu_i,E\) are formal illustrative inputs; no actual \(p_i/q_i\), source weight thresholds, or analytic error bound is certified by this example.

## 6. Actual quality in the nonzero shell witness

Assume the preceding sparse-shell construction and its common-column step. Its parameters obey

\[
\theta<T<\theta'<1,
\quad t=T/\theta,
\quad t^m>K,
\quad v_0'<v_0,
\quad w_0(\theta')^m<v_0'.
\]

Every injection produced there has independent test rows. Define its exact-quality score

\[
\mathcal C_H^{\rm quality}
=\frac{\sum_{\rho,i}\nu_iw_i(a_{\rho,i}-\beta_{\rho,i})
+c\sum_a n_a^2}{M_HH}.
\tag{11}
\]

The fractional shell matching can be decomposed into integral injections with the same edge marginals. In weighted coordinates \(y_i=w_ix_i\), the limiting simplex and radial shell transport are permutation symmetric. Therefore the expected weighted displacement in each coordinate has the same leading term. Rounding and cube-boundary differences contribute only lower-order terms for fixed parameters. Summing the different coefficients \(\nu_i\) replaces the uniform coefficient by \(\bar\nu\), not by its minimum or maximum.

As every injection has \(n_a\le\lceil H/v_0\rceil\), the collision bound is deterministic. Selecting an injection with no greater than expected linear cost gives

\[
\boxed{\limsup_{H\to\infty}\mathcal C_H^{\rm quality}
\le\bar\nu\,J+\frac c{v_0},\quad
J=\frac{K-1}{K}\frac m{m+2}
\frac{t^m}{t^m-1}(T-\theta).}
\tag{12}
\]

This is an existence statement for genuinely independent choices, conditional on the same interpolation theorem as the shell result. No intersection of unrelated good-probability events is assumed, because every injection is already independent.

For the explicit family \(\theta=1/2,B=3/5,C=13/10\), with the previous rational choices of \(t,T,\theta',v_0'\), one has \(J\to3/20\) and \(v_0\to\infty\). If \(\bar\nu_m\to\bar\nu_\infty\), then

\[
\limsup_{m\to\infty}\limsup_{H\to\infty}
\mathcal C_H^{\rm quality}\le\frac3{20}\bar\nu_\infty.
\tag{13}
\]

If \(w_*\to\infty\) and \(m^{-1}\sum_i\lambda_i/w_i\to0\), then \(\bar\nu_m\to2\), even if the individual unnormalized gains \(\lambda_i\) are unbounded. The bound in (13) remains 0.3. The bounds \(1-1/w_*\le u_i\le1+\Lambda/(\theta v_0)\) are uniform in the coordinates, so \(w_*\to\infty\) and \(v_0\to\infty\) imply \(R_{\rm ar}\to1-\theta=0.5\). More generally a limiting mean ratio below \(10/3\) leaves this fixed parameter family below that threshold. This is conditional on having approximation tuples with the stated quality/weight behavior; it is not an assertion about the actual sequence of approximations to π.

For this constant-translation witness, the scalar product itself is controlled on both sides. If \(S_H\) is the first term of (11), then

\[
-S_H\le\frac1{M_HH}\ln\left|\prod_\rho\xi_\rho\right|
\le-S_H+E_{\rm tr}.
\tag{14}
\]

Thus the exact quality has been fully accounted for in these scalars. Equation (14) is still not a lower bound for the accompanying test determinant. A nonzero test determinant may be much smaller than its collision upper bound, and different terms may cancel. Neither stronger function estimates nor cancellation arguments are ruled out.

## 7. Verification and scope

`verify_quality_criterion.py` checks the synthetic criterion with exact fractions; the equal-mean product penalty; finite anisotropic lattice capacities; and the exact scalar prefactor bounds over a complete small index set. It does not numerically test π's bad approximability and does not certify Theorem 2.1.

Sources: uploaded manuscript §3, together with the preceding conditional sparse-shell construction. The tail-removal and column-row denominator calculations are reproduced in §§1–2, rather than treated as unverified black boxes. All infinite-dimensional language is avoided: dimensions and weights are fixed for each determinant limit.
