# Exact cancellation at the actual period is diagonal-ideal vanishing

2026-11-01. This continues the full-row tail cancellation and constant-translation expansion. The algebraic statements below do not depend on the manuscript's interpolation theorem. Claims that the required nonzero high-dimensional minors exist still require a separate rank input. No claim about π's bad approximability or a new rational approximation is made.

## 1. What can cancel at the actual period

Let $m\ge2$, let $D(r_1,\ldots,r_m)$ be a nonzero polynomial over $\mathbb Q$ (or any algebraic number field), and put $\omega=2\pi i$. In the application, $D$ is a fixed full-row determinant after replacing the truncated logarithms by $\ln(1+t)$ exactly. Its columns have algebraic coefficients and its center coordinates are $jr_i$. Finite coefficient extraction makes $D$ a polynomial in the $r_i$.

Introduce the ordinary polynomial coordinates

\[
 z=r_1,\qquad \delta_i=r_i-r_1\quad(2\le i\le m).
\]

There is a unique finite expansion

\[
 D(z,z+\delta_2,\ldots,z+\delta_m)
 =\sum_{\gamma\in\mathbb Z_{\ge0}^{m-1}}
       F_\gamma(z)\delta^\gamma.                       \tag{1}
\]

Define $k$ to be the smallest $|\gamma|$ for which $F_\gamma$ is not the zero polynomial. Equivalently, if

\[
 I=(r_2-r_1,\ldots,r_m-r_1),
\]

then $D\in I^k\setminus I^{k+1}$. This is the order of vanishing along the full diagonal, calculated entirely over the coefficient field.

**Actual-period order theorem.** The least total degree of a nonzero term in the Taylor polynomial

\[
 D(\omega+\varepsilon_1,\ldots,\omega+\varepsilon_m)
\]

is exactly $k$. Its leading homogeneous part is

\[
 \boxed{\sum_{|\gamma|=k}F_\gamma(\omega)
       \prod_{i=2}^m(\varepsilon_i-\varepsilon_1)^{\gamma_i}.} \tag{2}
\]

**Proof.** Substitute $z=\omega+\varepsilon_1$ and $\delta_i=\varepsilon_i-\varepsilon_1$ in (1). No degree below $k$ can occur. Known transcendence of π implies that every nonzero polynomial $F_\gamma$ with algebraic coefficients satisfies $F_\gamma(\omega)\ne0$. The linear coordinates $\varepsilon_1,\varepsilon_2-\varepsilon_1,\ldots,\varepsilon_m-\varepsilon_1$ are independent, so the expression in (2) is a nonzero homogeneous polynomial of degree $k$. ∎

In particular,

\[
 D(\omega,\ldots,\omega)=0
 \quad\Longleftrightarrow\quad
 D(Z,\ldots,Z)\equiv0
 \quad\Longleftrightarrow\quad D\in I.                 \tag{3}
\]

Thus exact vanishing at this particular transcendental period is not an accidental algebraic coincidence peculiar to π. It is the algebraic property of vanishing along the entire diagonal. This does not preclude very small nonzero values there.

For $m=1$, the ideal has no normal variables and every nonzero polynomial has order zero at $\omega$; this recovers the earlier one-variable distinction without confusing center coalescence at zero with approximation error at the period.

### Exact regrouping of the remaining row terms

For a full-row minor, write the constant-translation expansion at a formal common center $W$ and group terms by the total coordinate increment

\[
 \kappa_i=\sum_{\rho}(a_{\rho,i}-\beta_{\rho,i}).
\]

Then the grouped identity is

\[
 D(W+\varepsilon_1,\ldots,W+\varepsilon_m)
 =\sum_\kappa G_\kappa(W)\varepsilon^\kappa,
 \qquad G_\kappa(W)\in\overline{\mathbb Q}[W].           \tag{4}
\]

More explicitly, $G_\kappa(W)$ is the sum over those row selections of their test determinants times the rational factors

\[
 \prod_\rho\binom{a_\rho}{\beta_\rho}
                 j_\rho^{|a_\rho-\beta_\rho|}.
\]

This follows directly from determinant multilinearity applied to the exact constant translations. At $j=0$ any positive increment has zero coefficient. It uses the polynomial row formula $P(1+t,jW+\ln(1+t)+u)$; no assertion $e^{jW}=1$ at arbitrary $W$ is needed.

For each fixed $\kappa$, $G_\kappa(\omega)=0$ if and only if $G_\kappa$ is identically zero. Therefore all exact cancellations of a specified total increment at the actual period can be certified by rational polynomial arithmetic. A merely small value requires additional quantitative estimates; this equivalence does not supply them.

## 2. A weighted column budget limits normal vanishing order

Let the original source rows be $(j,s,\beta)$ and let the fixed polynomial columns be $P_1,\ldots,P_M$. Give the center variables $r_i$ the positive weights $w_i$. Let

\[
 b_\ell=\max\{w\alpha:Y^hX^\alpha\text{ occurs in }P_\ell\},
 \qquad E=\sum_{\ell=1}^M b_\ell-\sum_{\rho}w\beta_\rho.
\]

For a single entry, differentiation by $u^\beta$ reduces the available weighted center degree by $w\beta$, and coefficient extraction in $t$ cannot increase center degree. Determinant expansion consequently gives

\[
 \deg_w D\le E.                                         \tag{5}
\]

This remains true for fixed mixed polynomial columns; their $b_\ell$ is the maximum transverse degree of their actual monomial support. In the original domain $w_0h+w\alpha\le H$, every $b_\ell\le H$, so

\[
 E\le MH(1-\bar b_H),\qquad
 \bar b_H=\frac1{MH}\sum_\rho w\beta_\rho.               \tag{6}
\]

Assume the coordinates are ordered so that $w_1\le w_i$ for every $i$. Give $z$ weight $w_1$ and $\delta_i$ weight $w_i$. The substitution $r_1=z$, $r_i=z+\delta_i$ does not increase weighted degree. Hence every nonzero monomial $z^a\delta^\gamma$ in (1) satisfies

\[
 w_1a+\sum_{i=2}^m w_i\gamma_i\le E.
\]

For a nonzero $D\in I^k$ this implies

\[
 \boxed{k\min_{i\ge2}w_i\le E\le MH(1-\bar b_H).}       \tag{7}
\]

This argument transforms only the center polynomial. It does **not** assume that replacing transverse jet coordinates by $u_i-u_1$ preserves the original weighted row packet; for unequal weights that assertion would require separate proof and can fail.

### A quantitative upper bound, with its missing input visible

Let $L$ be the sum of the absolute coefficients in the normal-coordinate polynomial (1), and let $d$ be its ordinary total degree. If $|\varepsilon_i|\le\eta$ with $0<\eta\le1/2$, then $|\delta_i|\le2\eta\le1$ and $|z|\le|\omega|+1$. Thus

\[
 |D(\omega+\varepsilon)|
 \le L\max(1,|\omega|+1)^d(2\eta)^k.                    \tag{8}
\]

The order $k$ alone gives no lower bound on $D$, and the leading polynomial (2) may vanish along special directions of the actual error vector, causing further cancellation.

For monomial columns the coefficient norm in (8) also has an explicit bound. On the formal $t$ circle of radius $1/2$, replace the coefficients of $\ln(1+t)$ by their absolute values, whose sum is $\ln2$. The normal center $jr_i$ has coefficient norm at most $2j$ in $(z,\delta)$, so a single entry has coefficient norm at most

\[
 2^s(3/2)^h\binom\alpha\beta(2j+\ln2)^{|\alpha|-|\beta|}
 \le2^s(3/2)^h(4K)^{|\alpha|}.
\]

The center-zero case satisfies the same upper bound. The coefficient norm is submultiplicative, so for a fixed raw monomial minor,

\[
 \boxed{L\le M!\,2^{\sum_\rho s_\rho}
            (3/2)^{\sum_{\rm cols}h}
            (4K)^{\sum_{\rm cols}|\alpha|}.}            \tag{8a}
\]

Consequently, with $w_* = \min_iw_i$ and the common weighted column budget,

\[
 \frac{\ln L}{MH}
 \le\frac{\ln M}{H}+\frac{\ln2}{v_0}
 +\max\left\{\frac{\ln(3/2)}{w_0},\frac{\ln(4K)}{w_*}\right\}.
\]

Since $d\le MH/w_*$ and $|\omega|+1<8$, (8) yields the concrete direct normal-order estimate

\[
\boxed{
 \frac{\ln|D(\omega+\varepsilon)|}{MH}
 \le -\frac{k}{MH}\ln\frac1{2\eta}
 +\frac{\ln2}{v_0}
 +\max\left\{\frac{\ln(3/2)}{w_0},\frac{\ln(4K)}{w_*}\right\}
 +\frac{\ln8}{w_*}+\frac{\ln M}{H}.
}                                                        \tag{8b}
\]

This estimate is useful only for a nonzero minor with a suitable positive normal order. It has no collision saving when $k=0$. For general fixed polynomial columns, replace the corresponding factor in (8a) by

\[
 \sum_{h,\alpha}|[Y^hX^\alpha]P_\ell|(3/2)^h(4K)^{|\alpha|}.
\]

Thus column scaling and mixing must still be charged through their actual coefficients. In particular, pre-scaling by $\prod_i(q_iL_i)^{\alpha_i}$ can add a main exponent and cannot be ignored when comparing with the arithmetic lower bound.

A direct limitation of this particular bound is already visible. Let $w_{\rm normal}=\min_{i\ge2}w_i$ and $L_\eta=\ln(1/(2\eta))\ge0$. If $L_\eta\le w_{\rm normal}$, then (7) gives

\[
 \frac{kL_\eta}{MH}
 \le (1-\bar b_H)\frac{L_\eta}{w_{\rm normal}}
 \le1-\bar b_H.                                        \tag{8c}
\]

The existing uniform arithmetic comparison would require this saving to exceed $1-\bar b_H$ plus its positive arithmetic and analytic errors. Hence that comparison cannot close using (8b) under the stated condition. This is a limitation of these displayed bounds, not a lower bound on the actual determinant; sharper column-specific arithmetic, directional rather than common-$\eta$ information, or further cancellation are not excluded.

## 3. Complete finite packets with positive order at the actual diagonal

These examples concern exact finite matrices. They are not asserted to meet the manuscript's unknown successive weight thresholds.

### A raw monomial minor

Take

\[
 m=2,\quad K=2,\quad H=2,\quad v_0=w_0=w_1=w_2=1,
 \quad\theta=1/3.
\]

The complete strict source row packet consists of $\beta=(0,0)$, $s=0,1$ at each center. Order the rows by $(j,s)$ and choose the columns

\[
 1,\quad X_1,\quad YX_1,\quad X_2.
\]

All columns satisfy the source budget. Exact coefficient extraction gives

\[
 \boxed{D(r_1,r_2)=r_1(r_1-r_2)=-z\delta_2.}            \tag{9}
\]

Thus its actual-period order is one:

\[
 D(\omega+\varepsilon_1,\omega+\varepsilon_2)
 =-\omega(\varepsilon_2-\varepsilon_1)
  -\varepsilon_1(\varepsilon_2-\varepsilon_1).
\]

The original determinant is nonzero whenever $r_1\ne0$ and $r_1\ne r_2$, whereas the same minor at the diagonal period is zero. This illustrates why full rank of the entire diagonal matrix does not imply that every fixed original minor has diagonal order zero.

The same four-row minor also fits a separated finite budget $H=4,v_0=2,w_0=1,(w_1,w_2)=(1,3),\theta=1/5$. Its weighted column budget in (5) is $E=5$, and the center polynomial has weighted degree four. No row-coordinate change is used.

### A pure normal factor of order four

Take

\[
 H=3,\quad v_0=2,\quad w_0=1/4,\quad w_1=w_2=1,
 \quad\theta=1/4,\quad K=2.
\]

Again the complete rows are $\beta=(0,0),s=0,1$ at each center, ordered by $(j,s)$. Put $Z=X_2-X_1$ and choose the columns

\[
 1,\quad Y-1,\quad Z^2,\quad (Y-1)Z^2.
\]

The common logarithms cancel in $Z$, leaving $j(r_2-r_1)+u_2-u_1$. Therefore, and also by direct verification,

\[
 \boxed{D(r_1,r_2)=(r_2-r_1)^4.}                        \tag{10}
\]

Here $k=4$ and $E=4$, so (7) is sharp. The same columns and row profile also fit $H=7,v_0=4,w_0=1/4,(w_1,w_2)=(1,3),\theta=1/8$; now $E=12$ and $kw_2=12$, so the bound is sharp with unequal weights as well. At the true diagonal period the full expression is exactly $(\varepsilon_2-\varepsilon_1)^4$, not a power of the center separation from zero.

The construction extends to any polynomial power $Z^a$ and a finite number of time blocks $(Y-1)^s$ whenever all chosen rows and columns obey the stated budget. Its determinant is a signed power $(r_2-r_1)^{a\cdot\text{number of time blocks}}$. This observation alone does not build a fixed-parameter, all-$H$ complete row family: transverse derivative rows eventually enter the source packet and must then be handled.

## 4. Difference factors have a calculable gain and an arithmetic cost

Write the actual real approximation errors as

\[
 e_i=p_i/q_i-\pi,\qquad r_i=2ip_i/q_i,
 \qquad \alpha_i=-\ln|e_i|.
\]

If $|e_j|\le|e_i|/2$, the triangle inequalities give the two-sided estimate

\[
 \boxed{|e_i|\le|r_j-r_i|=2|e_j-e_i|\le3|e_i|.}         \tag{11}
\]

Thus a difference anchored at a much weaker approximation has a logarithmic gain determined by that weaker error, even when its other coordinate is extremely accurate.

There is also a weighted-degree restriction for a whole explicit product. Suppose

\[
 D(r)=F(r)\prod_{i<j}(r_i-r_j)^{N_{ij}},\qquad F\ne0,
\]

with nonnegative integers $N_{ij}$. Weighted degree is additive under multiplication of nonzero polynomials, so (5) implies

\[
 \sum_{i<j}N_{ij}\max(w_i,w_j)\le E.                    \tag{12}
\]

If the error separation in (11) holds on each used edge, then the logarithmic gain supplied by this explicit factor product satisfies

\[
 \frac{-\ln\left|\prod_{i<j}(r_i-r_j)^{N_{ij}}\right|}{MH}
 \le\frac{E}{MH}
     \max\!\left\{0,\max_{N_{ij}>0}\frac{\alpha_i}{\max(w_i,w_j)}\right\}. \tag{13}
\]

The outer maximum with zero makes the bound valid even when an error exceeds one. If no factors occur, the gain is zero. For example, with increasing weights and with the later error on every used edge at most half the earlier one, the maximum is at most $\max(0,\max_{i<m}\alpha_i/w_{i+1})$. If these ratios tend to zero, this explicit factor mechanism supplies no fixed positive normalized gain. This is a condition on the chosen approximation tuple; it is not claimed for every tuple for π.

Neither (12) nor (13) bounds the residual $F$. They also do not control cancellations between sums of different normal monomials. They constrain the stated product-factor construction only.

For the pure fourth-power example (10), let $Q=\operatorname{lcm}(q_1,q_2)$. If the two rational approximations are distinct, then

\[
 Q\left(\frac{p_2}{q_2}-\frac{p_1}{q_1}\right)\in\mathbb Z\setminus\{0\},
 \qquad |D|\ge\left(\frac2Q\right)^4.
\]

On the other hand $|D|\le[2(|e_1|+|e_2|)]^4$. Comparing these two estimates gives exactly

\[
 |e_1|+|e_2|\ge\frac1Q,                                 \tag{14}
\]

the ordinary separation bound for distinct rationals. Raising the normal difference to a high power amplifies both the apparent smallness and its denominator cost; in this example the powers cancel out of the final comparison. Any stronger consequence requires structure beyond that factor alone.

## 5. Verification

Run `python3 verify_diagonal_order.py`. The standard-library verifier reconstructs both four-row determinants by exact multivariate polynomial arithmetic; verifies their complete row packets and column budgets, including a sharp unequal-weight example; independently groups 469 exact constant-translation row choices by total increment; checks normal and error-coordinate orders on 728 nonzero rational polynomials; verifies the weak-anchor bound in 336 rational cases; checks 32 entry coefficient norms and the coefficient upper bound in 54 cases. It also checks exact denominator clearing for the known fractions 3, 22/7 and 355/113 without making numerical claims about their π errors.

The general order theorem rests on the proof and known transcendence of π. The computations certify the finite examples and algebraic transformations, not the manuscript's interpolation theorem or π's bad approximability.
