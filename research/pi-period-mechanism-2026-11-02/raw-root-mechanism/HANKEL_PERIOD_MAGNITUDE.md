# A Fourier–Hankel formula and a two-sided actual-period estimate

2026-11-02. This note treats a particular genuine raw-column family in one transverse dimension. It proves a uniform two-sided estimate at the actual period, not only a numerical near-root phenomenon. The estimate also quantifies exactly how the standard row-factorial normalization changes the scale. It does not prove a root-convergence theorem, a rational approximation result, or an obstruction for other determinant normalizations.

## 1. The complete packet and the two different determinants

For an integer n>=2 use K=2, H=2n, w0=1, w1=2, v0=2, v1=2n and theta=1/n. The complete source packet is

\[
L_n=\{(j,s,0):j=0,1,\ 0\le s<n\}.
\]

Use the 2n raw columns Y^h X^a, h=0,1 and 0<=a<n, ordered first by h and then by a. Every column obeys h+2a<=2n. Order rows first by j and then by s. Define the formal-period polynomial

\[
P_n(Z)=\det\left[[z^s]e^{hz}(jZ+z)^a\right]_{(j,s),(h,a)}.
\tag{1}
\]

For comparison define the actual exponential endpoint determinant

\[
F_n(Z)=\det\left[[z^s]e^{h(jZ+z)}(jZ+z)^a\right]_{(j,s),(h,a)}.
\tag{2}
\]

These are different functions at general Z. At omega=2*pi*i they agree because exp(omega)=1. No arbitrary formal Z is called an exponential period.

Write

\[
A_n=\prod_{a=0}^{n-1}a!,\qquad B_n=\prod_{k=0}^{2n-1}k!.
\]

The main conclusion is

\[
\boxed{\frac16 U_n<|P_n(\omega)|\le U_n,
\qquad U_n=(2\pi)^{n^2}\frac{A_n^2}{B_n},\qquad n\ge2.}
\tag{3}
\]

Thus these actual-period determinants are nonzero for an explicit analytic reason; transcendence is not required for their nonvanishing.

## 2. Taylor remainders turn the determinant into moments

Let E be the n-by-n matrix of multiplication by exp(z), truncated to degrees below n, and T(Z) the exact translation matrix on polynomials of degree below n:

\[
E_{s,a}=\begin{cases}1/(s-a)!,&s\ge a,\\0,&s<a,\end{cases}
\qquad T(Z)_{s,a}=\binom as Z^{a-s}.
\]

The block matrix of (2) is [I,E;T,e^Z ET]. Eliminating its first column block gives

\[
F_n(Z)=\det(e^Z ET-TE).
\]

For g_a(x)=e^x x^a, the (s,a) entry is its order-n Taylor remainder differentiated s times and divided by s!:

\[
\frac{1}{s!(n-s-1)!}
\int_0^Z (Z-t)^{n-s-1}g_a^{(n)}(t)\,dt.
\tag{4}
\]

The contour is the straight segment from 0 to Z; the integrand is entire. Since

\[
g_a^{(n)}(t)=e^t(1+\partial_t)^n t^a,
\]

and (1+partial_t)^n is a unit triangular operator on polynomials of degree below n, a determinant-one column change replaces g_a^(n) by e^t t^a. The row polynomials (Z-t)^(n-s-1) have reversed degrees. The sign of this reversal cancels the product of their leading signs. Their normalization contributes A_n^(-2). Substituting t=Zu therefore proves the entire identity

\[
\boxed{F_n(Z)=\frac{Z^{n^2}}{A_n^2}
\det\left[\int_0^1 u^{i+j}e^{Zu}\,du\right]_{i,j=0}^{n-1}.}
\tag{5}
\]

This is a direct Taylor-remainder proof, with no matrix-integral identity assumed.

There is also an exact algebraic version for (1). Set

\[
\widehat\mu_0(Z)=0,\qquad
\widehat\mu_k(Z)=Z^{-1}-kZ^{-1}\widehat\mu_{k-1}(Z),\quad k\ge1.
\]

Then

\[
\boxed{P_n(Z)=\frac{Z^{n^2}}{A_n^2}
\det[\widehat\mu_{i+j}(Z)]_{i,j=0}^{n-1}.}
\tag{6}
\]

To justify the formal substitution, retain a symbol E for exp(Z). Integration by parts gives mu0=(E-1)/Z and muk=E/Z-k*mu_(k-1)/Z. Both sides of (5) are then rational in Z and polynomial in E. An identity after E=exp(Z) is an identity in E: a finite sum of polynomials in Z times distinct exp(kZ) cannot vanish identically, as seen by its largest exponential at positive real infinity. Setting E=1 yields (6). The apparent negative powers cancel because its left side is a polynomial.

## 3. A positive density controls the absolute upper bound

Expanding the two determinants and integrating termwise gives

\[
\det\left[\int_0^1 u^{i+j}e^{Zu}\,du\right]
=\frac1{n!}\int_{[0,1]^n}
e^{Z\sum u_j}\prod_{i<j}(u_j-u_i)^2\,du.
\tag{7}
\]

At Z=0 the Cauchy determinant evaluates to

\[
H_n(0):=\det[1/(i+j+1)]_{i,j=0}^{n-1}
=\frac{A_n^4}{B_n}.
\tag{8}
\]

Normalize the nonnegative integrand in (7) at Z=0 to a probability density. Let S=sum u_j and X=S-n/2 under that density. Reflection u_j->1-u_j shows that X has a symmetric distribution. Hence

\[
\frac{H_n(it)}{H_n(0)}
=e^{int/2}\,\mathbb E\cos(tX),
\qquad t\in\mathbb R.
\tag{9}
\]

The absolute upper bound in (3) follows at once from |E cos(tX)|<=1 and (5)--(8).

## 4. Three exact moments give a uniform positive lower bound

The centered moments needed below are

\[
m_2=\mathbb E X^2=\frac{n^2}{4(4n^2-1)},
\]

\[
m_4=\mathbb E X^4
=\frac{3n^2(n^2-2)}{16(4n^2-9)(4n^2-1)},
\]

\[
m_6=\mathbb E X^6
=\frac{15n^2(4n^6-33n^4+54n^2-16)}
{64(4n^2-25)(4n^2-9)(4n^2-1)^2}.
\tag{10}
\]

Here is a finite algebraic derivation of these moment formulas. For a partition lambda of k, let its cells have zero-based coordinates (i,j), and let h_ij be their hook lengths. The alternant definition of Schur polynomials, followed by (7) and the Cauchy determinant, gives

\[
\mathbb E s_\lambda(u)
=\prod_{(i,j)\in\lambda}
\frac{(n+j-i)^2}{h_{ij}(2n+j-i)}.
\]

Multiplication of an alternant by sum u_j adds one cell; iteration gives the coefficient k!/product h_ij. Consequently

\[
\boxed{\mathbb E S^k
=k!\sum_{\lambda\vdash k,\ \ell(\lambda)\le n}
\prod_{(i,j)\in\lambda}
\frac{(n+j-i)^2}{h_{ij}^2(2n+j-i)}.}
\tag{11}
\]

Only partitions of k<=6 are needed. Centering these six raw moments gives (10). For n>=6 all partitions occur; the accompanying verifier checks the resulting identities as polynomial identities in n, not merely at sampled n. For n=1,...,5 it uses the length-restricted formula (11) directly. It also recursively counts cell-addition paths for every partition of size at most six and checks their agreement with k!/product h_ij, so that finite multiplicity calculation does not rely on an untested general hook-length formula.

For every real x,

\[
\cos x\ge1-x^2/2+x^4/24-x^6/720.
\]

For example, the sixth derivative of the difference is 1-cos x>=0, and its first six initial values vanish; integrate six times for x>=0 and use evenness. Thus

\[
\mathbb E\cos(tX)\ge L_n(t):=
1-t^2m_2/2+t^4m_4/24-t^6m_6/720.
\tag{12}
\]

For t>0, L_n is strictly decreasing. Indeed m4^2<=m2*m6 by Cauchy--Schwarz, so the quadratic

\[
m_2-(t^2/6)m_4+(t^4/120)m_6
\]

has negative discriminant and is positive. Its negative times t is L_n'(t).

Use 2*pi<44/7. For n=2,

\[
L_2(44/7)=113553277/555891525>1/6.
\]

For n>=3 put x=n-3>=0 and

\[
D(n)=(2n-5)(2n-3)(2n-1)^2(2n+1)^2(2n+3)(2n+5)>0.
\]

Direct substitution of (10) gives

\[
L_n(44/7)-\frac16=\frac{R(n-3)}{235298 D(n)},
\]

where, in descending powers,

\[
\begin{aligned}
R(x)={}&2841056x^8+68185344x^7+680567864x^6\\
&+3658868208x^5+11464693580x^4\\
&+21140470608x^3+21947378508x^2\\
&+11262098088x+1908713103.
\end{aligned}
\]

Every coefficient is positive. Therefore E cos(2*pi*X)>1/6 for every n>=2, proving (3). In this column order the phase is also fixed:

\[
P_n(2\pi i)=i^{n^2}(-1)^n U_n\,c_n,
\qquad 1/6<c_n\le1.
\tag{13}
\]

## 5. The scale changes after the exact row-factorial normalization

Elementary summation of log(a!) gives

\[
\log A_n=\tfrac12n^2\log n-\tfrac34n^2+O(n\log n).
\]

Since the factor c_n in (13) stays between 1/6 and 1, (3) implies the actual asymptotic statement

\[
\boxed{\log|P_n(2\pi i)|
=-n^2\log n+
\bigl[\tfrac32+\log(\pi/2)\bigr]n^2+O(n\log n).}
\tag{14}
\]

The normalized Taylor rows in (1) contain 1/s!. Multiplying every row by s! clears these denominators, with exact determinant factor A_n^2. Accordingly

\[
\boxed{\log|A_n^2P_n(2\pi i)|
=\log(\pi/2)n^2+O(n\log n).}
\tag{15}
\]

In particular, the row-factorial-cleared determinant grows exponentially on the n^2 scale, even though the coefficient-normalized determinant tends to zero faster than exp(-c*n^2) for each fixed c.

This concerns this specified clearing operation. A_n^2 need not be the least common denominator of the already expanded polynomial P_n; an asymptotic theorem for the least denominator has not been established here. At the actual period the cleared value is not an integer. At a rational imaginary formal center ip/q, the usual additional factor q^(n(n-1)) does give a Gaussian integer. No comparison with a varying rational approximant is supplied by (15).

The family is also outside the critical high-dimensional half-density regime: m=1, theta=1/n, M=2n and the ambient raw-column count is (n+1)^2, so the density tends to zero. The weights have not been tied to arbitrarily large approximation denominators. Its near-period algebraic roots and actual determinant scale cannot be presented as a rational approximation record or a proof about pi's bad approximability.

## 6. Finite least-denominator check

The companion near-root certificates provide the full rational coefficients of P_n for n=4,...,10. Let d_n be their least common denominator. In every listed case d_n P_n has integer coefficients with gcd one. Re-enclosing pi by a separate high-precision Machin interval gives the following exact outward integer enclosures:

| n | d_n | A_n^2/d_n | Interval containing the cleared absolute value |
|---:|---:|---:|---:|
| 4 | 144 | 1 | [279078, 279079] |
| 5 | 864 | 96 | [1012149, 1012150] |
| 6 | 74649600 | 16 | [5256609833, 5256609834] |
| 7 | 22394880000 | 27648 | [6533585621, 6533585622] |
| 8 | 758487711744000000 | 20736 | [46368157078953, 46368157078954] |
| 9 | 26759446470328320000000 | 955514880 | [13261296588019, 13261296588020] |
| 10 | 489407583931698134384640000000000 | 6879707136 | [60061286897514303, 60061286897514304] |

Thus row factorials do overclear substantially in some cases. Even the primitive integer polynomial has actual-period absolute value greater than one in all seven certified cases. These finite values are not monotone and are not used to assert any asymptotic law for d_n.

## Verification

`verify_hankel_magnitude.py` uses only the Python standard library. It checks the full endpoint and moment determinants at exact rational formal values of Z and E, the centered moment identities as polynomial identities, the finite small-n exceptions, the positive-coefficient lower-bound certificate and the factorial normalizations. A finite Machin interval verifies the numerical upper bound on pi used in the proof. The finite supplement reads `../root-separation-next/raw_near_period_verification.json`, records its SHA-256, computes the least denominators, and independently re-encloses the actual-period values. Run the companion `verify_raw_near_period.py` first to regenerate and certify those polynomial coefficients. Infinite-size conclusions rely on the displayed identities and positive-polynomial argument, not numerical extrapolation.
