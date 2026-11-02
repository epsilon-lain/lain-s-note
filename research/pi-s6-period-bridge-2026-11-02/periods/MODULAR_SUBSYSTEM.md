# A computable modular subsystem of the S6 period construction

2026-11-02. This note carries out a concrete tool-transfer test. It uses the modular coordinate in the recent S6 construction and identifies an explicit classical route to a rational hypergeometric series for 1/pi. These are known identities used as a bridge and a baseline, not a new pi-approximation record.

## 1. The precise geometric input and its scope

The original manuscript at https://alpo.ge/s6.pdf has, in Section 3.1, a modular coordinate tau satisfying j(tau)=1728t. Remark 3.12 displays the elliptic curve

\[
y^2=4x^3-3t^3(t-1)x-t^4(t-1)^2.
\]

Its invariant is directly j=1728t: with g2=3t^3(t-1), g3=t^4(t-1)^2, one has g2^3-27g3^2=27t^8(t-1)^3. We use this explicitly singled-out modular coordinate. We do not assert that the generic complex two-torus in the S6 construction splits as an elliptic product or has an elliptic quotient.

Set z=1/t=1728/j(tau). The following normalized/gauged period is classical:

\[
f(z)={}_2F_1(1/12,5/12;1;z)=E_4(\tau)^{1/4},\qquad f(0)=1.
\]

The fourth root is chosen by its cusp expansion. This is not an unscaled integral of dx/y; factors involving pi and algebraic gauge factors cannot be omitted when comparing raw geometric periods. The modular identity is recorded, for example, by Don Zagier in *Arithmetic and Topology of Differential Equations*, Example 2 on page 24 of the PDF:
https://people.mpim-bonn.mpg.de/zagier/files/tex/ECMHirzebruchLecture/ECMHirzebruchLectureCorrected.pdf

## 2. Exact differential equation, recurrence and integer coefficients

For theta=z d/dz,

\[
\boxed{[\theta^2-z(\theta+1/12)(\theta+5/12)]f=0.}
\]

Writing f=sum a_n z^n gives

\[
a_0=1,\qquad (n+1)^2a_{n+1}=(n+1/12)(n+5/12)a_n.
\]

A convenient all-index integer normalization is

\[
u_n=1728^n a_n
=\frac{12^n\prod_{k=0}^{n-1}(12k+1)(12k+5)}{(n!)^2}\in\mathbb Z.
\]

For primes p>=5, each complete block of p^e consecutive k gives at least one factor divisible by p^e in each arithmetic progression, cancelling the two copies of v_p(n!). For p=2 and 3, the factor 12^n pays the entire remaining denominator. More precisely

\[
v_2(u_n)=2s_2(n),\qquad v_3(u_n)=s_3(n),
\]

where s_b is the base-b digit sum. The exact integer recurrence is

\[
(n+1)^2u_{n+1}=12(12n+1)(12n+5)u_n.
\]

The first values are 1,60,39780,38454000,43751038500,54538294552560. Thus one obtains an explicit rational-coefficient period system with a single exponential coefficient denominator, rather than merely a statement that its monodromy is integral.

## 3. The logarithmic companion and the missing normalization constant

Put

\[
g=f\log z+h,\qquad h(z)=\sum_{n\ge1}a_n S_n z^n,
\]
\[
S_n=\sum_{k=0}^{n-1}\left(\frac{12}{12k+1}+\frac{12}{12k+5}-\frac2{k+1}\right).
\]

This is a second solution of the same equation. Its chosen cusp normalization is

\[
\boxed{2\pi i\tau=\frac gf-\log1728.}
\]

In particular h_1=31/72. The constant log1728 is required by j=q^{-1}+744+... and q=e^{2pi i tau}; it cannot be dropped to pretend that the logarithmic period already supplies a pi linear form. With L_m=lcm(1,...,m),

\[
L_{12n}\,1728^n[z^n]h\in\mathbb Z.
\]

Near the cusp, a positive counterclockwise loop around z=0 gives g -> g+2pi i f and tau -> tau+1. The original manuscript names an oppositely oriented cusp generator g0 with tau -> tau-1; these are inverse loops, not contradictory monodromy normalizations. For the Wronskian,

\[
\boxed{fg'-f'g=\frac1{z\sqrt{1-z}}.}
\]

Hence d tau/dz=1/(2pi i z sqrt(1-z)f^2). For an all-parameter justification of the normalization, the modular derivative identities give dz/dtau=2pi i z E6/E4=2pi i z sqrt(1-z) f^2. Comparing with the derivative of g/f and then using q~z/1728 fixes the displayed constant. This does not make tau at an arbitrary Gaussian rational point a known algebraic number.

The larger S6 period matrix also includes functions mu and beta. Theorem 3.4 and Section 3.3 explicitly allow beta=beta_part+c0 with a free complex constant. We avoid that additional period: integral monodromy does not by itself control these constants or their coefficient denominators.

## 4. The symmetric square produces integer hypergeometric data

Clausen's identity gives

\[
G(z)=f(z)^2={}_3F_2(1/6,1/2,5/6;1,1;z)
=\sum_{n\ge0}M_n(z/1728)^n,
\]
\[
\boxed{M_n=\frac{(6n)!}{(3n)!(n!)^3}\in\mathbb Z.}
\]

These are multinomial integers and obey

\[
(n+1)^3M_{n+1}=24(2n+1)(6n+1)(6n+5)M_n.
\]

They satisfy the third-order equation

\[
[\theta^3-z(\theta+1/6)(\theta+1/2)(\theta+5/6)]G=0.
\]

This is an exact scalar system on which coefficient costs and rational transformations can be tested.

## 5. An explicit pullback to a rational 1/pi formula

Let lambda be the modular lambda parameter and u=4lambda(1-lambda). Then

\[
j=\frac{64(4-u)^3}{u^2},\qquad
z=\frac{27u^2}{(4-u)^3}.
\]

Set

\[
H(u)={}_3F_2(1/2,1/2,1/2;1,1;u)
=\sum_{n\ge0}\binom{2n}{n}^3(u/64)^n.
\]

The theta identities and Clausen transformation give the actual bridge

\[
\boxed{G\left(\frac{27u^2}{(4-u)^3}\right)
=\sqrt{1-u/4}\,H(u).}
\]

Both sides have constant term 1, and the equality holds near u=0. It is also checked by exact rational coefficients independently of any numerical value of pi. The accompanying independent operator audit checks the full pulled-back differential operator.

At u0=1/64, j=255^3 is the CM point tau=i sqrt7, discriminant -28. Ramanujan's classical formula is

\[
\boxed{5H(u_0)+42u_0 H'(u_0)
=\sum_{n\ge0}(42n+5)\binom{2n}{n}^3/4096^n=\frac{16}{\pi}.}
\]

Differentiating the bridge shows that this is the same explicit modular subsystem as

\[
(432+7182\theta)G((4/85)^3)=85\sqrt{255}/\pi,
\]

or

\[
\sum_{n\ge0}(133n+8)\frac{(1/6)_n(1/2)_n(5/6)_n}{(n!)^3}
(4/85)^{3n}=\frac{85\sqrt{255}}{54\pi}.
\]

These are Ramanujan's formulas (29) and (34), respectively, in *Modular equations and approximations to pi*, Quarterly Journal of Mathematics 45 (1914), 350-372:
https://ramanujan.sirinudi.org/Volumes/published/ram06.pdf
https://ramanujan.sirinudi.org/Volumes/published/ram06.html

The algebraic sqrt255 factor in the level-one coordinate is absorbed by the rational pullback into the H coordinate, where the 1/pi prefactor is rational. This is a concrete reason to examine arithmetic coordinates, not just generic period analogies.

## 6. What has actually been achieved by this bridge

The test does reach an explicit rational 1/pi series with integer coefficients from the modular coordinate already present in the S6 construction. It does not use the full six-sphere topology, nor establish that these classical identities are new. Truncating a fast pi-computation formula is not automatically a good q-squared rational approximation; the reduced denominator and the series tail still need joint accounting.

The companion program checks the differential/Frobenius relations, integrality and Wronskian through 40 coefficients; the pullback through 30 coefficients; and the modular normalization using an independently formed E4 and Delta q-expansion through 12 coefficients. All checks are exact rational arithmetic. The general displayed formulas are justified by the derivations and classical modular identities, not by finite verification alone.
