# Closed base periods and a linear truncation schedule

2026-11-02. This supports the complete-truncation generator in ../constructive-review/COMPLETE_TRUNCATION_GENERATOR.md. It gives explicit base periods, a uniform nonzero-period rule K=13m, and ordinary convergence to pi with reduced denominators tending to infinity. It does not prove quality tending to zero.

For m>=1 put c=4m+1 and

\[
J_j=\int_0^1\frac{x^{2m+2j}(1-x^2)^{2m}}{(1+x^2)^{4m+1}}\,dx
=\beta_j+\alpha_j\pi/4,\quad j=0,1.
\]

Write

\[
A=\binom{2m}{m},\quad C=\binom{4m}{2m},\quad
S=2m-1,\quad R=6m+1,\quad B=1/(SC).
\]

Then the exact periods are

\[
\boxed{\alpha_0=A/64^m,\quad\beta_0=0,
\quad\alpha_1=\frac RS\alpha_0,\quad\beta_1=-B.}
\tag{1}
\]

## A beta-integral proof

Put x=tan(theta), then v=sin^2(2theta). The angular kernel identity yields

\[
J_0=\frac{4^{-m}}4\operatorname B(m+1/2,m+1/2).
\]

Also tan^2(theta)=(2-v-2sqrt(1-v))/v, giving

\[
J_1=\frac{4^{-m}}4\left[
2\operatorname B(m-1/2,m+1/2)
-\operatorname B(m+1/2,m+1/2)
-2\operatorname B(m-1/2,m+1)\right].
\]

The first beta ratio is 4m/(2m-1), while the final half-integer beta value is rational. The gamma factorial identities simplify these expressions to (1). Equivalently J_1=(R/S)J_0-B.

Consequently

\[
\boxed{\Delta=\alpha_0\beta_1-\alpha_1\beta_0
=-\frac{A}{64^mSC}\ne0.}
\tag{2}
\]

A common denominator for all four base periods therefore divides

\[
\boxed{D_0\mid64^mSC,}
\tag{3}
\]

which is sharper here than the general Fourier/lcm bound.

## Uniform schedule without a stopping search

Let mu_j,epsilon_K be the complete rational truncation quantities in the companion note. Since R/S<=7, M_1<=M_0, and

\[
M_0=\int_0^1x^{2m}(1-x^2)^{2m}\,dx\le(4/27)^m,
\]

the nonzero test is implied by

\[
\varepsilon_K\,8(4/27)^m\le B/2.
\tag{4}
\]

Set K=13m. The rational remainder estimate, using binom(c+K,K+1)<=2^(c+K) and the tail-ratio factor<=3, gives

\[
\varepsilon_{13m}\le\frac43\left(\frac{2^{21}}{3^{17}}\right)^m.
\tag{5}
\]

Since C<=16^m, condition (4) follows if

\[
\frac{64}{3}(2m-1)\left(\frac{2^{27}}{3^{20}}\right)^m<1.
\]

The left side at m=1 is 2^33/3^21<1. Each next-to-previous ratio is at most 3*(2^27/3^20)<1. This proves (4) for every m>=1. In particular alpha_K is nonzero at every scale.

The same bound proves epsilon_K<2^(-c), because the ratio of the right side of (5) to 2^(-4m-1) is (8/3)*(2^25/3^17)^m<1 already at m=1 and decreases. Thus S_K(t)>0 on [0,1] and mu_0,mu_1>0.

## Exact reduced fractions and a simplified denominator bill

Let L=lcm(den(mu_0),den(mu_1)) and u_j=L mu_j be integers. Then

\[
\boxed{\frac{p_m}{q_m}
=\frac{2^{6m+2}u_0}{AC(Ru_0-Su_1)},}
\tag{6}
\]

followed by complete gcd reduction. The denominator in (6) is positive by the nonzero-period test: alpha_K=alpha_0(mu_1-(R/S)mu_0)<0. In fact

\[
\frac1{2C}\le\frac{Ru_0-Su_1}{L}\le\frac3{2C},
\]

so a useful bound on the fully reduced denominator is

\[
\boxed{q_m\le\frac32 A L.}
\tag{7}
\]

This retains cancellations among the closed base-period factors. The universal truncation moment bound gives L dividing 3^(17m+1)*odd(lcm(1,...,32m+3)). Neither the factor L nor the final gcd can be dropped from a quality argument.

## Ordinary convergence is proved, quality-zero is not

The complete-truncation error obeys

\[
|I_K|\le2\varepsilon_K(1+\varepsilon_K)M_0M_1,
\quad |\alpha_K|\ge|\Delta|/2.
\]

Using A>=4^m/(2m+1), C<=16^m gives

\[
|\Delta|\ge\frac{256^{-m}}{(2m-1)(2m+1)}.
\]

Combining these with (5), M_0M_1<=(4/27)^(2m), and epsilon_K<1 proves

\[
\boxed{\left|\pi-\frac{p_m}{q_m}\right|
\le\frac{128}{3}(4m^2-1)
\left(\frac{2^{33}}{3^{23}}\right)^m\longrightarrow0.}
\tag{8}
\]

Here 2^33/3^23<1. Since pi is irrational, convergence in (8) forces the reduced denominators q_m to tend to infinity. Multiplying (8) by q_m^2 is still costly; no quality-zero conclusion follows. A charged joint bound follows from the complete-truncation certificate and (3):

\[
\boxed{Q_m\le12\,64^m A S C L^2\varepsilon_K(1+\varepsilon_K)M_0M_1.}
\tag{9}
\]

The actual-gcd bound in the companion note remains preferable when computing a particular scale. Formula (9) is an unconditional a priori bill, not a proof of net improvement.

The exact period formulas were independently checked against the triangular rational-antiderivative extractor for m=1,...,6. All statements (1)--(9) above have all-scale proofs; finite checks are only algebra audits.
