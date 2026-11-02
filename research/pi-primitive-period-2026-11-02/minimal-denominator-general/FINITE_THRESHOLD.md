# An explicit threshold after balanced coefficient clearing

2026-11-02. This is a supplementary exact estimate for the complete two-frequency raw-column family. The denominator theorem is proved in `../denominator-reflection-route/ENDPOINT_CONTENT_AND_DENOMINATORS.md`; the actual-period estimate is proved in `../../pi-period-mechanism-2026-11-02/raw-root-mechanism/HANKEL_PERIOD_MAGNITUDE.md`. No empirical formula for odd polynomial content is needed.

## Statement

Write

\[
A_n=\prod_{j=0}^{n-1}j!,\qquad
r=\lfloor n/2\rfloor,\quad s=\lceil n/2\rceil,
\qquad D_n=\frac{A_n^2}{A_r^2A_s^2}.
\]

The balanced coefficient-clearing theorem gives

\[
Q_n(Z):=D_nP_n(Z)\in\mathbb Z[Z].
\]

For even n, D_n is exactly the minimal positive integer coefficient denominator; for odd n it is a sufficient denominator. In either parity:

\[
\boxed{0<|Q_n(2\pi i)|<1\qquad\text{for every }n\ge23.}
\]

The value at an irrational argument is not an integer. This statement therefore gives no contradiction and does not by itself prove a rational approximation result for pi. It establishes that the large scale seen after row-factorial clearing is not intrinsic to integer coefficient clearing.

## A rational majorant

The actual-period estimate gives, for n>=2,

\[
0<|P_n(2\pi i)|\le (2\pi)^{n^2}\frac{A_n^2}{A_{2n}}.
\]

Put

\[
T=44/7>2\pi,\qquad a=T^2/16=121/49<3,
\]

and define the positive rational number

\[
W_n=T^{n^2}\frac{A_n^4}{A_{2n}A_r^2A_s^2}.
\]

Then \(0<|Q_n(2\pi i)|<W_n\). The ratio of two terms of the same parity is exactly

\[
\frac{W_{n+2}}{W_n}
=T^{4n+4}
\frac{[n!(n+1)!]^4}
{(2n)!(2n+1)!(2n+2)!(2n+3)!\,(r!s!)^2}.
\tag{1}
\]

Since the central binomial coefficient is the largest of the 2n+1 coefficients in \((1+1)^{2n}\),

\[
\binom{2n}{n}\ge\frac{4^n}{2n+1}.
\]

Substituting this into (1), and collecting the four consecutive factorials in its denominator, gives

\[
\frac{W_{n+2}}{W_n}
\le \frac{T^4}{4}\frac{2n+1}{2n+3}(n+1)^2
\left(\frac{a^n}{r!s!}\right)^2
<b_n,
\]

where

\[
b_n=\frac{T^4}{4}(n+1)^2
\left(\frac{a^n}{r!s!}\right)^2.
\]

## Uniform decrease, established without floating point

The auxiliary ratio satisfies

\[
\frac{b_{n+2}}{b_n}
=\left(\frac{n+3}{n+1}\right)^2
\frac{a^4}{[(r+1)(s+1)]^2}.
\]

For n>=16, the three elementary estimates

\[
\frac{n+3}{n+1}<2,\qquad a<3,\qquad (r+1)(s+1)\ge81
\]

imply

\[
\frac{b_{n+2}}{b_n}<\frac{4\cdot3^4}{81^2}=\frac4{81}.
\]

The exact rational base checks are

\[
b_{16}<\frac12,\qquad b_{17}<\frac1{10}.
\]

Induction in the two parities yields \(b_n<1/2\) for every n>=16, hence

\[
\frac{W_{n+2}}{W_n}<\frac12\qquad(n\ge16).
\]

The two remaining rational comparisons are

\[
W_{23}<\frac1{60000},\qquad W_{24}<\frac1{10^{10}}.
\]

They prove the displayed theorem by parity induction. They also give the explicit bounds

\[
W_{23+2j}<\frac{2^{-j}}{60000},\qquad
W_{24+2j}<10^{-10}2^{-j},\qquad j\ge0.
\]

The code `verify_finite_threshold.py` uses only standard-library integer arithmetic and `Fraction` for all proof checks. It checks the four finite premises by cross multiplication, audits (1) and the auxiliary ratio independently for n=2,...,40, and records `finite_threshold_verification.json`. The uniform conclusion comes from the inequalities and induction above, not from the finite audit range. Floating-point logarithms in the JSON are diagnostic only.

## Scope

The sufficient threshold 23 is not asserted to be the first n for which the actual cleared value drops below 1. For example, the rational upper bound is still above 1 at n=22, but an upper bound above 1 does not establish that the actual value is above 1.

Even minimal coefficient clearing now has superexponentially small actual-period values. To use this for rational approximation, one must also compare the value at a rational argument, its denominator after substitution, and the derivative cost. Integer coefficients alone do not turn the small value at 2*pi*i into an integer or a linear form in pi.
