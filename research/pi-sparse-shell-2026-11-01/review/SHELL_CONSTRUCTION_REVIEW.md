# Independent review: a sparse shell construction with a genuine period minor

Date: 2026-11-01. Status: the construction below is valid **conditional on the uploaded manuscript's Theorem 2.1**. It does not independently establish that interpolation theorem, an irrationality measure, or bounded partial quotients of π.

This note independently reviews the sparse shell construction. It includes a direct collision bound, requiring no density integration, and an endpoint-compatible refinement of the manuscript's parameter choice.

## 1. Data and statement

Fix positive rational weights \(w_0,v_0,w_1,\ldots,w_m\), an integer \(K\ge2\), and \(0<\theta<1\). Write \(wa=\sum_iw_ia_i\). Original rows are

\[
\mathcal R_H=\{(j,s,\beta):0\le j<K,\quad v_0s+w\beta/\theta<H\}.
\]

Assume that parameters \(T,\theta',v_0'\) can be chosen with

\[
\theta<T<\theta'<1,\qquad v_0'<v_0,
\qquad K<(T/\theta)^m,
\qquad w_0(\theta')^m<v_0'.                 \tag{1}
\]

For the original matrix's rank certificate also assume its numerical conditions \(K\theta^m<1\) and \(Kw_0\theta^m/v_0<1\). Choose the separated weights \(w_i\) to satisfy both the original \(K\)-center interpolation theorem and its auxiliary one-center version with parameters \((w_0,v_0',\theta')\). Successive thresholds can be combined by taking their maximum at each choice of \(w_i\).

For all sufficiently large admissible \(H\), there exists a selection of test rows

\[
(j,s,\beta)\longmapsto(j,s,a),\qquad a\ge\beta,
\]

which fixes \(a=\beta\) at \(j=0\), has all pairs \((s,a)\) distinct, lies in the one-center packet \(J_{v_0',w/\theta'}(H)\), and has full row rank at the exact period \(\omega=2\pi i\).

Here “full row rank” concerns the actual test-row matrix on the polynomial domain \(P_W(H)\). The proof is not merely a support matching argument or generic nonvanishing assertion.

## 2. Shell transport and a joint injection at each time order

Fix \(s\) for which \(L=\theta(H-v_0s)>0\), and set

\[
\sigma=\sum_iw_i,\qquad S=L+\sigma,\qquad t=T/\theta,
\qquad\lambda=t^m-1>K-1.
\]

The source cubes \(\beta+[0,1)^m\), with \(w\beta<L\), are disjoint and entirely contained in \(\{x\ge0:wx<S\}\). On the positive orthant away from the origin put

\[
F(x)=\frac{(S^m+\lambda(wx)^m)^{1/m}}{wx}\,x.                 \tag{2}
\]

The value at zero is irrelevant to all measures. If \(r=wx\), its image radius is

\[
r'=(S^m+\lambda r^m)^{1/m}.
\]

For \(0<r<S\), one has \(S<r'<tS\) and \(r'>r\). Thus \(F(x)\ge x\) coordinatewise. The radial map is injective, and its Jacobian is

\[
\det DF=(r'/r)^{m-1}\frac{dr'}{dr}=\lambda.                 \tag{3}
\]

The same calculation follows after the diagonal change of coordinates \(y_i=w_ix_i\); its determinant factors cancel.

For each source cube choose \(x\) uniformly and use \(a=\lfloor F(x)\rfloor\) as a fractional target assignment. It has the following properties:

\[
a\ge\beta,\qquad L<wa<tS=T(H-v_0s)+t\sigma.                \tag{4}
\]

The strict lower bound follows because rounding down loses strictly less than \(\sigma\), while \(wF(x)>S\). This separates all new types from the anchor types, whose weights satisfy \(w\beta<L\).

For one copy of the source cubes, the pushforward density is at most \(1/\lambda\), since their union has density at most one and \(F\) is injective with constant Jacobian. Across the \(K-1\) nonzero centers, the fractional load in any target unit cube is therefore at most

\[
(K-1)/\lambda<1.                                         \tag{5}
\]

There are finitely many source and target vertices. The resulting fractional bipartite matching has source degree exactly one and target degree at most one, so integrality of the bipartite matching polytope supplies an injection of **all** pairs \((j,\beta)\), \(1\le j<K\), into types \(a\), at this fixed \(s\). Each chosen edge has positive fractional weight and therefore obeys (4). This is a joint injection across centers; separate injections at each center would not suffice.

Keep the anchor assignment \(a=\beta\) for \(j=0\). Its targets are disjoint from the shell by (4). Consequently all selected \((s,a)\) pairs are distinct across every center.

## 3. Uniform containment in one auxiliary packet

For the shell assignments,

\[
v_0's+wa/\theta'
\le (T/\theta')H+
       \bigl(v_0'-(T/\theta')v_0\bigr)s+t\sigma/\theta'.
\]

Since \(0\le s<H/v_0\), the right side is at most

\[
H\max\{T/\theta',v_0'/v_0\}+t\sigma/\theta'.              \tag{6}
\]

Both ratios are strictly below one. Thus (6) is strictly less than \(H\) for all sufficiently large \(H\), uniformly over the selected rows. The anchors satisfy the same estimate with \(T\) replaced by \(\theta\). This handles all boundary layers without an additional injection argument.

## 4. From a collapsed center to the actual period

For a raw monomial column \(Y^hX^\alpha\), define the period-polynomial matrix by

\[
B_{(j,s,a),(h,\alpha)}(W)
=\binom\alpha a[t^s](1+t)^h
       (jW+\log(1+t))^{|\alpha|-|a|},                     \tag{7}
\]

with the entry zero when \(a\not\le\alpha\). Each entry lies in \(\mathbb Q[W]\). Equation (7) is the polynomial matrix obtained after using the period identity \(e^{hj\omega}=1\); it is not the full original analytic expression with \(e^{hjW}\) at arbitrary \(W\).

At \(W=0\), the label \(j\) disappears, and (7) is precisely the coefficient functional

\[
[t^su^a]P(1+t,u_1+\log(1+t),\ldots,u_m+\log(1+t)).
\]

The selected \((s,a)\) pairs are distinct and contained in the auxiliary one-center packet. Theorem 2.1 with \(K=1\), whose numerical conditions follow from (1), makes the whole packet independent for sufficiently large admissible \(H\). Hence the selected rows are independent at \(W=0\).

Some maximal minor \(D(W)\in\mathbb Q[W]\) therefore satisfies \(D(0)\ne0\). It cannot vanish at \(W=2\pi i\): otherwise the transcendental number \(2\pi i\) would be a root of a nonzero rational polynomial. Thus the actual period test matrix has full row rank.

This proves actual nonvanishing at π, not merely nonvanishing for a generic center.

### Common columns require a separate step

This argument gives a nonzero raw monomial minor for the test matrix. It does **not**, on its own, prove that the same raw monomial minor is nonzero for the rational original-row matrix. Two full-rank matrices can have disjoint sets of nonzero maximal coordinate minors.

To obtain an actual nonzero term in one determinant expansion using a common column family, invoke the previously proved bounded mixed-columns/common-basis lemma for the original and test matrices. Its bounded rational or integral column combinations must also be carried through the arithmetic and analytic estimates. This is an explicit dependency, not a consequence of (7) alone.

For the selected rows, the \(d=0,k=0\) summand of the manuscript's row expansion has scalar

\[
\binom a\beta\prod_i\bigl(j(r_i-\omega)\bigr)^{a_i-\beta_i}.
\]

It equals one for an anchor and is nonzero for a nonzero center, because the rational approximation parameters \(r_i\) differ from \(\omega\). Thus the row choices are realized by nonzero scalar summands, rather than only by formal support labels.

## 5. A direct collision bound

Because all pairs \((s,a)\) are distinct, a fixed type \(a\) occurs at most once at each time order. The original time range is

\[
s=0,\ldots,\lceil H/v_0\rceil-1.
\]

Writing \(n_a\) for the number of rows of type \(a\), one obtains the exact bound

\[
\sum_a n_a^2\le\lceil H/v_0\rceil\sum_a n_a
=\lceil H/v_0\rceil M_H.
\]

Therefore

\[
\limsup_{H\to\infty}\frac{c\sum_a n_a^2}{M_HH}
\le\frac c{v_0}.                                         \tag{8}
\]

No randomized matching choice or limiting density computation is needed for this estimate.

## 6. Linear cost and the high-dimensional limit

At each nonzero center, (4) implies

\[
w(a-\beta)\le T(H-v_0s)-w\beta+O(1),
\]

where the fixed \(O(1)\) is independent of \(H\). At the anchor the cost is zero. The elementary simplex averages of the original rows are

\[
\frac{\sum v_0s}{M_HH}\longrightarrow\frac1{m+2},
\qquad
\frac{\sum w\beta}{M_HH}\longrightarrow\frac{m\theta}{m+2}.
\]

The distributions are identical at all \(K\) centers, so

\[
\limsup_{H\to\infty} C_H
\le \nu\frac{K-1}{K}
       \frac{T(m+1)-m\theta}{m+2}+\frac c{v_0}.             \tag{9}
\]

Here \(C_H=[\nu\sum w(a-\beta)+c\sum_a n_a^2]/(M_HH)\) is exactly the two-part cost under study. Equation (9) is a bound on that cost, not on the full determinant including every auxiliary error term.

On the manuscript's core path

\[
K=\lfloor C^m\rfloor,\quad w_0=B^{-m},\quad
v_0=2K\theta^mw_0,\quad B<C\theta<1,
\]

choose any fixed \(1<\kappa<2\) and take \(T\) with

\[
T^m=\kappa K\theta^m.
\]

If rational parameters are required, use an arbitrarily small rational perturbation inside the strict interval \(K\theta^m<T^m<2K\theta^m\). One can then choose rational \(\theta'>T\) and \(v_0'<v_0\) such that \(w_0(\theta')^m<v_0'\). For large \(m\), \(T<1\), because \(T\to C\theta<1\).

Moreover \(v_0\to\infty\), since \(C\theta/B>1\), and (9) yields

\[
\limsup_{m\to\infty}\limsup_{H\to\infty}C_H
\le\nu\theta(C-1).                                      \tag{10}
\]

The order of limits matters: choose \(m\), the parameters and separated weights first, then let admissible \(H\) tend to infinity.

## 7. Endpoint-compatible refinement of the original parameter choice

In the manuscript's Lemma 4.1, \(A<\sqrt\theta\), and the required interval for \(C\) is \(A/\theta<C<1/A\). It can be narrowed to

\[
A/\theta<C<1/\sqrt\theta,
\]

which is nonempty because \(A<\sqrt\theta\). Choose

\[
A<B<C\theta.
\]

Then \(B<C\theta<\sqrt\theta<1/C\), so the original constraint \(B<1/C\) is also retained. The endpoint value \(\nu=2\) in (10) obeys

\[
2\theta(C-1)<2(\sqrt\theta-\theta)<1-\theta,
\]

where the final inequality is strict because \(0<\theta<1\). This retains the original parameters selected for a fixed exponent \(\nu_0>2\) while leaving a strict low-cost margin at the endpoint cost coefficient two. These are different roles for the exponent. The source gap \(\nu_0(A-\theta)>1-\theta\) is **not** asserted with \(\nu_0=2\); indeed, it cannot hold at two when \(A<\sqrt\theta\). Thus the endpoint cost conclusion is not a counterexample to the manuscript's full \(\nu_0>2\) analytic bound. The numerical example below certifies geometric and endpoint-cost feasibility only unless its separate \(A,\nu_0\) source-gap choices are supplied.

For the simple numerical choice \(\theta=0.5,B=0.6,C=1.3\), (10) gives endpoint cost at most \(0.3\), below \(1-\theta=0.5\). The parameter argument above shows the mechanism is not restricted to that example.

## 8. Scope of the conclusion

The shell construction bypasses the earlier whole-\(K\)-packet containment obstruction by collapsing distinct center labels to a single center at the algebraic specialization \(W=0\). It still uses an auxiliary **one-center** full packet and the source interpolation theorem. Anchor rigidity is respected exactly.

Conditional on that theorem, it establishes a genuine independent sparse test selection and a low two-part cost in the critical separated-weight regime \(v_0/w_0\to0\). It does not by itself establish a noncanceling determinant sum, refute every possible strengthened analytic estimate, or prove a statement about π's irrationality measure.
