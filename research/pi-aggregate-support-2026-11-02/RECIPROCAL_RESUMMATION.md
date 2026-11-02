# Exact resummation by the reciprocal center polynomial

2026-11-02. This note rewrites the complete normalized shear sum as a single finite polynomial attached to the same selected columns. It also proves a restriction on exact radial cancellation at the actual transcendental period: a nonconstant radial aggregate cannot lose its first-order term, or any intermediate error degree. This is qualitative; it supplies no useful lower bound on a nonzero coefficient and no conclusion about pi's bad approximability.

## 1. Columns, centers, and their finite polynomial

Fix the complete source packet

\[
L=\{(j,s,\beta):0\le j<K,\quad v_0s+v\cdot\beta<H\},
\qquad M=|L|.
\]

Select M columns of the form

\[
g_c(Y)X^{\alpha_c},\qquad g_c\in\overline{\mathbb Q}[Y],
\quad\alpha_c\in\mathbb Z_{\ge0}^m,
\]

with their actual polynomial support inside the permitted weighted column space. Each column has a single transverse multidegree. Both raw and shifted monomials are included. Put

\[
\nu_c=\operatorname{ord}_{Y=1}g_c(Y),\qquad
D_i=\sum_c\alpha_{c,i}-\sum_L\beta_i,
\qquad \Delta=\sum_Ls-\sum_c\nu_c.
\tag{1}
\]

Zero columns can be discarded as giving a zero determinant. Define the period-normalized center polynomial

\[
P_C(r_1,\ldots,r_m)
=\det\left([z^su^\beta]g_c(e^z)
 \prod_i(jr_i+z+u_i)^{\alpha_{c,i}}\right)_{L,c}.
\tag{2}
\]

This is a polynomial with algebraic coefficients. At a common center r_i=omega=2*pi*i its entries are actual-period jets. At arbitrary r, (2) is the specified formal-period polynomial, not a claim that exp(jr_i)=1.

If P_C is nonzero, every D_i is nonnegative, and it has the support form

\[
\boxed{P_C(r)=\sum_{\delta}c_\delta r^{D-\delta},
\qquad 0\le\delta_i\le D_i,\quad |\delta|\le\Delta.}
\tag{3}
\]

Indeed, one entry uses some nonnegative number d_i of its remaining spatial powers as z powers, and leaves center powers alpha_i-beta_i-d_i. The factor g_c(e^z) has vanishing order nu_c, so sum_i d_i<=s-nu_c. Adding these quantities over a determinant permutation gives (3). Determinant cancellation may delete support but cannot add it.

The reversed polynomial

\[
Q_C(x)=x^D P_C(x_1^{-1},\ldots,x_m^{-1})
=\sum_\delta c_\delta x^\delta
\tag{4}
\]

is therefore an ordinary polynomial of total degree at most Delta, and coordinate degree at most D_i. For a raw monomial Y^h the valuation nu is zero, not h. The intrinsic column-space filtration developed in `FILTERED_TIME_SUPPORT.md` can give better valuations after an explicitly tracked basis change within a fixed transverse type. No arbitrary raw subset is silently replaced by shifted monomials.

## 2. The complete aggregate is an affine substitution in Q_C

Let omega=2*pi*i and define the normalized slope matrix

\[
A_C(\zeta)_{j,s,\beta;c}
=[z^su^\beta]g_c(e^z)
 \prod_i\bigl(j\omega+(1-\zeta_i)z+u_i\bigr)^{\alpha_{c,i}}.
\]

Set F_C(zeta)=det A_C(zeta). Then the exact identity is

\[
\boxed{F_C(\zeta)=\omega^{|D|}
Q_C\left(\frac{1-\zeta_1}{\omega},\ldots,
         \frac{1-\zeta_m}{\omega}\right)
=\sum_\delta c_\delta\omega^{|D|-|\delta|}(1-\zeta)^\delta.}
\tag{5}
\]

For a direct proof set b_i=(1-zeta_i)^(-1) and r_i=b_i omega. In the entry of (2), factor b_i^(alpha_i) from the column and b_i^(-beta_i) from coefficient extraction after rescaling u_i. The resulting matrix is A_C(zeta), giving

\[
F_C(\zeta)=b^{-D}P_C(b_1\omega,\ldots,b_m\omega).
\tag{6}
\]

Substitute (3) to get (5). Initially the argument uses zeta_i!=1; both sides of (5) are polynomials, so the identity extends everywhere.

The exact shear formula also gives

\[
A_C(\zeta)=T_{L,L}B_{L,C}+T_{L,O}B_{O,C}.
\]

For the complete source packet U=T_{L,L} is unit triangular. Consequently

\[
\boxed{F_C(\zeta)=\det(B_{L,C}+R(\zeta)B_{O,C}),
\qquad R=U^{-1}T_{L,O}.}
\tag{7}
\]

Thus (5) represents the full signed exterior sum, including its phase cancellation, rather than a bound by the largest single target. It does not need B_{L,C} to be invertible. Returning to the original arithmetic center determinant still requires the prefactor b^D in (6).

In particular every complete coefficient has the finite expression

\[
[\zeta^\gamma]F_C
=(-1)^{|\gamma|}
\sum_{\delta\ge\gamma}
 \binom\delta\gamma c_\delta\omega^{|D|-|\delta|}.
\tag{8}
\]

This groups all target selections with the same error multidegree. A nonzero individual target contribution does not establish that (8) is nonzero.

There is an exact zero test by homogeneous slices. Let Q_d(x)=sum_(|delta|=d) c_delta x^delta. Then (8) vanishes at omega if and only if, for every d,

\[
\boxed{\sum_{|\delta|=d,\,\delta\ge\gamma}
 \binom\delta\gamma c_\delta
=\frac{\partial^\gamma Q_d(1,\ldots,1)}{\gamma!}=0.}
\tag{8a}
\]

Indeed different d multiply different powers of the transcendental omega. With algebraic column coefficients these powers cannot cancel to zero unless every algebraic coefficient is zero. This criterion does allow exact cancellations within the same homogeneous slice; it retains the distinction between general error directions and the radial direction below.

## 3. Exact radial cancellation has no isolated missing degree

Consider the common radial direction zeta_1=...=zeta_m=t. This keeps all center coordinates equal in (6). Put

\[
c_d=\sum_{|\delta|=d}c_\delta,
\qquad F_{\rm rad}(t)=\sum_d c_d\omega^{|D|-d}(1-t)^d.
\tag{9}
\]

If every c_d is zero then P_C(Z,...,Z) is identically zero and F_rad is identically zero. Otherwise let d_* be the largest d with c_d!=0.

**Radial coefficient theorem.** Every coefficient of t^k for 0<=k<=d_* is nonzero at omega, and no higher coefficient occurs.

For each such k,

\[
[t^k]F_{\rm rad}
=(-1)^k\sum_{d=k}^{d_*}\binom dk c_d\omega^{|D|-d}.
\tag{10}
\]

The polynomial in omega on the right is nonzero: its lowest power is |D|-d_*, with nonzero coefficient (-1)^k binom(d_*,k)c_(d_*). All other d give different powers. Known transcendence of omega ensures that a nonzero polynomial with algebraic coefficients cannot vanish there.

In particular

\[
\boxed{F_{\rm rad}'(0)=0
\quad\Longleftrightarrow\quad F_{\rm rad}\text{ is constant}.}
\tag{11}
\]

The constant may be zero. In one transverse dimension this is the entire error variable: a nonconstant normalized aggregate cannot begin with a quadratic or higher error term while its linear coefficient vanishes. Equivalently, first-order vanishing forces P_C(r)=a r^D, including a=0.

For several transverse dimensions this is only a radial statement. For example the genuine four-row raw minor with columns 1,X_1,YX_1,X_2 has

\[
P_C(r_1,r_2)=r_1(r_1-r_2),\qquad D=(2,1),
\]

so

\[
F_C(\zeta_1,\zeta_2)=\omega^2(\zeta_1-\zeta_2).
\]

Its radial restriction vanishes identically although its normal variation is nonzero. Independent approximation errors generally do not follow a common radial line. Equation (11) cannot rule out positive normal-diagonal order or higher-order cancellation in those other directions.

## 4. Exact finite aggregate tests

The companion `AGGREGATE_PERIOD_POLYNOMIAL.md` computes full source polynomials and compares them with direct complete exterior sums. In its six-row example,

\[
P_C(Z)=Z^3(Z+3),\qquad D=4,
\]

and (5) gives exactly

\[
F_C(\zeta)=\omega^4+3\omega^3-3\omega^3\zeta.
\]

There are 36 target choices and 11 nonzero target minors. Despite these individual contributions, every aggregate coefficient from degree two through six is zero. Degrees five and six are already excluded by D=4; the lower cancellations must be checked in the signed sum. The radial theorem predicts the nonzero constant and linear coefficients once P_C is known.

For the previously certified 32-row packet, the full source polynomials in both recorded bases have D=62 and center valuation 27. Their complete error degree is therefore 35, and every coefficient of degrees 0 through 35 is nonzero at omega. In particular the full coefficient of degree 30 is nonzero. This does not identify any individual target term as surviving unchanged, nor bound the size of that coefficient.

`verify_reciprocal_resummation.py` checks small multivariate source determinants and the full coefficient identity by exact polynomial arithmetic, as well as the diagonal grouping and radial support statement at the formal-polynomial level. The aggregate companion verifier independently checks the complete six-row expansion and the larger exact interpolation certificates.

## 5. What this does and does not remove

The new identity replaces uncertainty about whether a complete coefficient is zero with a finite algebraic question about one source polynomial and its reversed coefficients. It also prevents an incorrect inference that a positive time surplus necessarily allows a leading error of arbitrarily high order in the one-dimensional or common-radial direction.

It leaves the central quantitative issue intact: a nonzero algebraic-coefficient polynomial can have a very small value at omega, and its coefficient height and degree both grow with the packet. Transcendence alone gives no uniform lower bound on (8), no upper bound sufficient for the arithmetic contradiction, and no noncancellation of the aggregate evaluated at the actual error tuple. The critical high-dimensional case and pi's Diophantine conclusion remain unresolved.
