# A complete-truncation cancellation rule with a nonzero-period certificate

2026-11-02. This is a constructive alternative to cancelling every low moment separately. It gives an explicit generator inside the user's new family, together with coefficient, reduced-denominator and signed-integral bounds. It does **not** prove a quality-zero limit.

## 1. Fixed data and exact rational base periods

Fix positive integers a,b and an integer c>=a+b+1. Write

\[
w(x)=x^{2a}(1-x^2)^b,\quad f(t)=(1+t)^{-c},
\quad M_j=\int_0^1w(x)x^{2j}\,dx
=\frac{2^b b!}{\prod_{r=0}^{b}(2a+2j+2r+1)}.
\]

All M_j are positive rational numbers. Define two rational-period pairs by

\[
J_j=\int_0^1w(x)x^{2j}f(x^2)\,dx
=\beta_j+\frac{\alpha_j\pi}{4},\qquad j=0,1.
\]

The coefficients alpha_j,beta_j are rational and computable from a,b,c without a numerical pi value. One convenient exact procedure is as follows. For N(t)=sum n_j t^j, expand N(y-1)=sum d_l y^l. Let

\[
K_s=\int_0^1(1+x^2)^{-s}\,dx=u_s+v_s\pi/4,
\]
\[
u_1=0,\quad v_1=1,\qquad
u_s=\frac{2s-3}{2s-2}u_{s-1}+\frac1{(s-1)2^s},
\quad v_s=\frac{2s-3}{2s-2}v_{s-1}.
\]

Then beta=d_c+sum_(s=1)^c d_(c-s)u_s and alpha=sum_(s=1)^c d_(c-s)v_s. Apply this to N=t^a(1-t)^b and tN for j=0,1.

## 2. The two base period vectors are independent

Put

\[
\boxed{\Delta=\alpha_0\beta_1-\alpha_1\beta_0\ne0.}
\tag{1}
\]

Here is a proof covering the permitted parameters, rather than a numerical nonzero test. If the determinant vanished, some nonzero linear multiplier H(t)=h_0+h_1t would have both rational periods zero. The invertible rational antiderivative map from the parent note would then give alpha=0 and R(1)=0 for

\[
\frac{w(x)H(x^2)}{(1+x^2)^c}
=\frac{d}{dx}\frac{R(x)}{(1+x^2)^{c-1}}.
\]

The zero of order 2a at zero forces x^(2a+1)|R. The zero of order b at one, together with R(1)=0, forces (x-1)^(b+1)|R. Odd parity also forces (x+1)^(b+1)|R. Therefore

\[
R=x^{2a+1}(1-x^2)^{b+1}S(x^2),\qquad
\deg S\le d-2,\quad d=c-a-b.
\]

If d=1 this already forces R=0. Otherwise, if S has degree e, the leading coefficient of the differential numerator is multiplied by 2(e-d)+5, an odd nonzero integer. Thus the quotient H has degree exactly e+2>=2. This contradicts deg H<=1. Hence (1).

## 3. Cancel one entire rational truncation

For K>=c define the rational polynomial

\[
S_K(t)=\left(\frac23\right)^c
\sum_{k=0}^{K}(-1)^k\binom{c+k-1}{k}
\left(\frac{2t-1}{3}\right)^k.
\]

The parent's uniform remainder estimate is

\[
\|f-S_K\|_{[0,1]}\le\varepsilon_K,
\quad
\varepsilon_K=
\left(\frac23\right)^c
\frac{\binom{c+K}{K+1}}{3^{K+1}}
\left(1-\frac{c+K+1}{3(K+2)}\right)^{-1}.
\tag{2}
\]

This is positive, rational, and tends to zero for each fixed a,b,c. Define exact rational numbers

\[
\mu_j=\int_0^1w(x)x^{2j}S_K(x^2)\,dx,\quad j=0,1,
\qquad
\boxed{H_K(t)=\mu_1-\mu_0t.}
\tag{3}
\]

Then

\[
\boxed{\int_0^1w(x)H_K(x^2)S_K(x^2)\,dx=0.}
\]

This cancels the whole truncated approximation in one rational equation. It does not demand that each monomial moment vanish, so it is not the maximal Rodrigues rank-one specialization.

The resulting exact periods are

\[
\alpha_K=\mu_1\alpha_0-\mu_0\alpha_1,
\qquad\beta_K=\mu_1\beta_0-\mu_0\beta_1,
\qquad I_K=\beta_K+\alpha_K\pi/4.
\tag{4}
\]

No target fraction, numerical pi, or continued fraction appears in this generator.

## 4. A computable rule guarantees alpha_K is nonzero

The bounds |mu_j-J_j|<=epsilon_K M_j and (1) give

\[
|\alpha_K-\Delta|
\le\varepsilon_K(|\alpha_0|M_1+|\alpha_1|M_0).
\]

Choose K to be the least integer K>=c satisfying the exact rational test

\[
\boxed{\varepsilon_K(|\alpha_0|M_1+|\alpha_1|M_0)
\le |\Delta|/2.}
\tag{5}
\]

This search terminates for every fixed a,b,c, and every chosen K, or every larger K meeting the same test, has

\[
\frac{|\Delta|}{2}\le|\alpha_K|\le\frac{3|\Delta|}{2}.
\tag{6}
\]

The rule can thus be repeated at a=m,b=2m,c=4m+1. It has an explicit nonzero-period certificate at every scale. The fixed schedule below additionally proves convergence to pi and divergence of the reduced denominators, but not monotonic denominator growth or improving quality.

### Closed base periods and a fixed linear schedule for the user's family

For a=m,b=2m,c=4m+1, the two base period pairs simplify exactly to

\[
\boxed{\beta_0=0,\qquad
\alpha_0=\frac{\binom{2m}{m}}{64^m},\qquad
\alpha_1=\alpha_0\frac{6m+1}{2m-1},\qquad
\beta_1=-\frac1{(2m-1)\binom{4m}{2m}}.}
\tag{5a}
\]

A direct verification uses x=tan(theta). The first integrand becomes 2^(-4m) sin^(2m)(4theta), giving beta_0=0 and alpha_0 as stated. On the full interval 0<theta<pi/2, replace tan^2(theta) by the average of tan^2(theta) and cot^2(theta); the beta integral then gives the alpha ratio (6m+1)/(2m-1). The antisymmetric quarter-interval difference gives

\[
\beta_1=-4^{-m}\int_0^1u^{2m-2}(1-u^2)^m\,du
=-\frac1{(2m-1)\binom{4m}{2m}}.
\]

Thus the base period extraction requires only binomial coefficients. Write B=|beta_1| and r=alpha_1/alpha_0=(6m+1)/(2m-1)<=7. Since M_j<=(4/27)^m and B>=16^(-m)/(2m-1), equation (5) follows whenever

\[
8\varepsilon_K(4/27)^m\le B/2.
\]

The completely explicit choice

\[
\boxed{K_m=13m}
\tag{5b}
\]

always satisfies this. Indeed (2), using K>=c and binom(c+K,K+1)<=2^(c+K), gives

\[
\varepsilon_{13m}\le(4/3)\left(\frac{2^{21}}{3^{17}}\right)^m.
\]

The sufficient inequality's left side divided by its right side is consequently at most

\[
\frac{64}{3}(2m-1)\left(\frac{2^{27}}{3^{20}}\right)^m<1.
\]

For m=1 this follows from 2^33<3^21, and the ratio of each consecutive expression to its predecessor is at most 3*(2^27/3^20)<1. Hence alpha_K never vanishes for this fixed, linear truncation schedule.

This rule also proves absolute-error convergence. The central-binomial estimates binom(2m,m)>=4^m/(2m+1) and binom(4m,2m)<=16^m imply

\[
|\Delta|=\alpha_0 B\ge
\frac{256^{-m}}{(2m-1)(2m+1)}.
\]

Since M_0M_1<=(4/27)^(2m), epsilon_K<1, and the remainder bound below gives |I_K|<=4epsilon_K M_0M_1, we obtain

\[
\boxed{\left|\pi-\frac{p_m}{q_m}\right|
\le\frac{128}{3}(4m^2-1)
\left(\frac{2^{33}}{3^{23}}\right)^m\longrightarrow0.}
\tag{5c}
\]

The exponential base is less than one (indeed less than 1/9). Therefore the actually reduced denominators q_m tend to infinity: a bounded-denominator subsequence in a bounded interval would contain only finitely many rational numbers, none equal to the irrational pi. This is an absolute-error convergence theorem; it does not imply that q_m^2 times the error tends to zero.

There is even an explicit reduced-fraction formula. If u_j=L_K mu_j are integers, reduce the ratio

\[
\boxed{\frac{p_m}{q_m}=
\frac{2^{6m+2}u_0}
{\binom{4m}{2m}\binom{2m}{m}
\bigl((6m+1)u_0-(2m-1)u_1\bigr)}.}
\tag{5d}
\]

The displayed numerator and denominator must still be divided by their gcd; the nonzero-period theorem certifies the denominator is nonzero.

## 5. Integer realization and its cost

Let L_K=lcm(den(mu_0),den(mu_1)). Then L_K H_K is an integer polynomial. The matrix M_c in the parent construction has determinant magnitude

\[
T_c=2^{c-1}(c-1)!.
\]

Multiplying the numerator by T_c L_K gives integer R and integer alpha, by the adjugate formula for M_c inverse. This places the construction literally in the parent's integer-antiderivative family. Its multiplier satisfies

\[
\|T_cL_KH_K\|_1
\le T_cL_K(1+\varepsilon_K)(M_0+M_1),
\]

because |mu_j|<=(1+epsilon_K)M_j. The large common scaling is not called a denominator saving; it cancels from the exact ratio and quality identity.

A universal moment-denominator bound is

\[
\boxed{L_K\mid3^{c+K}\operatorname{odd}\operatorname{lcm}
(1,\ldots,2(a+b+K+1)+1).}
\tag{7}
\]

Indeed coefficients of S_K have denominator dividing 3^(c+K); each integrated monomial contributes an odd denominator 2r+1 with r<=a+b+K+1.

For a sharper reduced-denominator bill, let D_0 be the least common multiple of the denominators of alpha_0,alpha_1,beta_0,beta_1. Put

\[
A_K=D_0L_K\alpha_K\in\mathbb Z,\quad
B_K=D_0L_K\beta_K\in\mathbb Z.
\]

The fully reduced candidate and denominator are exactly

\[
\boxed{\frac{p_K}{q_K}=-\frac{4B_K}{A_K},
\qquad q_K=\frac{|A_K|}{\gcd(|A_K|,4|B_K|)}
\le D_0L_K|\alpha_K|.}
\tag{8}
\]

The harmless scaling T_c in the integer realization does not need to be charged again in (8); this is an exact reduction of the period ratio.

For the user's family, the closed periods yield the particularly simple common denominator D_0|64^m(2m-1)binom(4m,2m). Together with (6), this gives

\[
\boxed{q_m\le\frac32\binom{2m}{m}L_{13m}.}
\tag{8a}
\]

Thus the large fixed base-period denominator cancels from this upper bill. The growing moment denominator L_(13m) must still be charged; no smallness theorem for that integer is assumed.

One explicit bound on D_0 is

\[
\boxed{D_0\mid4^{c-1}\operatorname{odd}\operatorname{lcm}(1,\ldots,c-1).}
\tag{9}
\]

To see this, put x=tan(theta). The Fourier expansion of cos^(2s-2)(theta), integrated from 0 to pi/4, gives

\[
v_s=4^{-(s-1)}\binom{2s-2}{s-1},\quad
u_s=4^{-(s-1)}\sum_{r=1}^{s-1}
\binom{2s-2}{s-1-r}\frac{\sin(r\pi/2)}r.
\]

Only odd r contribute to u_s. The shift coefficients d_l are integral, proving (9).

## 6. Joint upper bounds, with no omitted denominator cost

The cancellation in (3) proves the exact rational certificate

\[
|I_K|\le\varepsilon_K(|\mu_1|M_0+|\mu_0|M_1).
\]

Hence the achieved quality obeys

\[
\boxed{q_K^2\left|\pi-\frac{p_K}{q_K}\right|
\le\frac{4q_K^2}{|\alpha_K|}\varepsilon_K
(|\mu_1|M_0+|\mu_0|M_1).}
\tag{10}
\]

This uses the actually reduced denominator. Combining (6),(8) and the elementary mu bounds also gives a wholly a priori bill

\[
\boxed{Q_K\le12(D_0L_K)^2|\Delta|\varepsilon_K
(1+\varepsilon_K)M_0M_1.}
\tag{11}
\]

These quality bounds are not asserted to tend to zero. For fixed a,b,c, replacing L_K by the explicit universal right-hand side of (7) gives a denominator envelope with logarithm (ln3+2)K+o(K), while epsilon_K has rate 3^(-K) up to a polynomial. The resulting substituted upper certificate grows rather than decays. This does not prove that the actual L_K or the actual certificate must grow, and it does not exclude improvements after exact denominator cancellation. It shows that this available uniform bill does not prove quality tending to zero. At growing m the costs of D_0,Delta,M_j and K must also be retained.

This generator is therefore a precise next test with an explicit nonzero-period rule, not a completed pi breakthrough. It separates the genuine signed cancellation from the arithmetic condition still needed for quality to improve.
