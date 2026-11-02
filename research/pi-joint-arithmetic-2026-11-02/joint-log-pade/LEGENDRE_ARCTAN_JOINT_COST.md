# Exact reduced-denominator control for the Legendre–arctangent Padé family

2026-11-02. This note supplies a concrete joint error/denominator result for a classical fixed-parameter Padé family. It is proved for the **actual reduced rational fractions**, without substituting an unreduced common denominator. The conclusion obstructs convergence of the entire natural sequence to quality zero. It does not exclude a good subsequence and does not decide whether pi is badly approximable.

## 1. Rational sequences and their exact target

Fix rational numbers a>0 and c>0. Let L_n(t) be the Legendre polynomial normalized by L_n(1)=1; equivalently

\[
L_0(t)=1,\quad L_1(t)=t,\quad
(n+1)L_{n+1}(t)=(2n+1)tL_n(t)-nL_{n-1}(t).
\]

The Cauchy-transform numerator is the polynomial

\[
W_n(z)=\frac12\int_{-1}^{1}\frac{L_n(z)-L_n(t)}{z-t}\,dt.
\]

Thus W_0=0 and W_1=1. For n>=1 it obeys the same three-term recurrence as L_n. To see this, subtract the Legendre recurrence at z and t, divide by z-t, and integrate; the extra term is a multiple of the zero integral of L_n(t).

Define real rational numbers

\[
\ell_n(a)=i^{-n}L_n(ia),\qquad
w_n(a)=i^{1-n}W_n(ia).
\]

Their initial conditions and recurrence have no complex arithmetic:

\[
\ell_0=1,\ \ell_1=a,\quad w_0=0,\ w_1=1,
\]
\[
\boxed{(n+1)x_{n+1}=(2n+1)a x_n+n x_{n-1}\quad(n\ge1),}
\tag{1}
\]

for either sequence x. In particular every ell_n is positive. The approximants and their reduced fractions are

\[
R_n(a)=\frac{w_n(a)}{\ell_n(a)},\qquad
r_n=cR_n(a)=\frac{p_n}{q_n},\quad
\gcd(p_n,q_n)=1,\quad q_n>0.
\tag{2}
\]

Their target is alpha=c arctan(1/a). For a=1,c=4 this is pi.

## 2. A positive exact error identity

Set

\[
C(z)=\frac12\int_{-1}^{1}\frac{dt}{z-t}.
\]

For z outside [-1,1], subtracting W_n from L_n C gives

\[
L_n(z)C(z)-W_n(z)=\frac12\int_{-1}^{1}\frac{L_n(t)}{z-t}\,dt.
\]

Orthogonality of L_n to polynomials of degree less than n further gives

\[
L_n(z)\int_{-1}^{1}\frac{L_n(t)}{z-t}\,dt
=\int_{-1}^{1}\frac{L_n(t)^2}{z-t}\,dt.
\]

Indeed their difference is the integral of L_n(t) times the polynomial
(L_n(z)-L_n(t))/(z-t), whose degree in t is at most n-1. Therefore

\[
C(z)-\frac{W_n(z)}{L_n(z)}
=\frac1{2L_n(z)^2}\int_{-1}^{1}\frac{L_n(t)^2}{z-t}\,dt.
\tag{3}
\]

At z=ia, symmetry makes the real part of the integral zero, while C(ia)=-i arctan(1/a). Since L_n(ia)^2=(-1)^n ell_n^2, equation (3) becomes

\[
\boxed{(-1)^n\left(\arctan(1/a)-R_n(a)\right)
=\frac{a}{2\ell_n(a)^2}
\int_{-1}^{1}\frac{L_n(t)^2}{a^2+t^2}\,dt>0.}
\tag{4}
\]

This holds for n=0 as well. It proves alternation and nonzero error; there is no unproved cancellation assumption. The exact squared norm is

\[
\int_{-1}^{1}L_n(t)^2\,dt=\frac2{2n+1}.
\]

Both orthogonality and this norm can be established directly from Rodrigues' formula
L_n(t)=(2^n n!)^{-1}(d/dt)^n(t^2-1)^n. Integrating by parts n times proves orthogonality. For the norm use the leading coefficient (2n)!/[2^n(n!)^2] and the beta integral
integral_(−1)^1 (1−t^2)^n dt=2^(2n+1)(n!)^2/(2n+1)!.

Since a^2<a^2+t^2<a^2+1 except on sets of measure zero, (4) gives the strict bounds

\[
\boxed{\frac{a}{(a^2+1)(2n+1)\ell_n^2}
<\left|\arctan(1/a)-R_n(a)\right|
<\frac1{a(2n+1)\ell_n^2}.}
\tag{5}
\]

Thus the rational approximants converge exponentially in n for fixed a. One elementary growth justification is to choose 0<epsilon<min(a,1). Eventually the two coefficients in (1), after division by n+1, exceed 2a−epsilon and 1−epsilon. The positive root lambda of lambda^2=(2a−epsilon)lambda+(1−epsilon) exceeds one, since substituting lambda=1 makes the left side minus the right side equal 2epsilon−2a<0. Positive initial conditions then give ell_n>=C lambda^n after discarding finitely many initial terms. No denominator estimate is needed for convergence to the target itself.

## 3. The Wronskian forces a denominator product after full reduction

Define D_n=ell_n w_(n−1)−ell_(n−1)w_n. Equation (1) gives

\[
D_1=-1,\qquad D_{n+1}=-\frac n{n+1}D_n,
\qquad\boxed{D_n=\frac{(-1)^n}{n}.}
\tag{6}
\]

Hence adjacent approximants are distinct and satisfy the exact identity

\[
\boxed{|r_n-r_{n-1}|=\frac{c}{n\ell_n\ell_{n-1}}.}
\tag{7}
\]

For reduced fractions p_n/q_n and p_(n−1)/q_(n−1), their nonzero cross-difference is a nonzero integer. It follows that

\[
\frac1{q_nq_{n-1}}\le |r_n-r_{n-1}|,
\qquad
\boxed{q_nq_{n-1}\ge\frac{n\ell_n\ell_{n-1}}c.}
\tag{8}
\]

This is a statement about the **reduced denominators themselves**. Arbitrarily large gcd cancellation in either approximant cannot evade the neighboring-product constraint. Neither individual denominator is given a lower bound by (8).

## 4. A uniform joint-quality obstruction

Let

\[
Q_n=q_n^2\left|\alpha-\frac{p_n}{q_n}\right|,
\qquad \alpha=c\arctan(1/a).
\]

Multiplying the two lower error bounds (5) by the squared denominator product in (8) gives, for every n>=1,

\[
\boxed{Q_n Q_{n-1}>
\frac{a^2 n^2}{(a^2+1)^2(4n^2-1)}
>\frac{a^2}{4(a^2+1)^2}.}
\tag{9}
\]

The factors c cancel exactly. Therefore the complete sequence Q_n cannot tend to zero, even when every rational approximant is fully reduced. Every adjacent pair contains an index j with

\[
Q_j>\frac{a}{2(a^2+1)}.
\]

For pi, a=1,c=4, this specializes to

\[
\boxed{Q_nQ_{n-1}>\frac{n^2}{4(4n^2-1)}>\frac1{16}.}
\tag{10}
\]

Thus at least one member of every adjacent pair has normalized error greater than 1/4.

There is also an exact same-parity constraint. Add the two signed adjacent differences and use the ell recurrence to obtain, for n>=2,

\[
|r_n-r_{n-2}|=
\frac{c(2n-1)a}{n(n-1)\ell_n\ell_{n-2}}.
\]

The same integer cross-difference and strict error bounds now give

\[
\boxed{Q_nQ_{n-2}>
\frac{n^2(n-1)^2}{(2n-1)^2(a^2+1)^2(2n+1)(2n-3)}
\longrightarrow\frac1{16(a^2+1)^2}.}
\tag{11}
\]

For pi the limiting lower constant is 1/64. Thus neither the complete even-index nor the complete odd-index sequence has Q_n tending to zero. A **sparse** subsequence of either parity remains possible. For that stronger question one must still control actual gcds or the distribution of the individual q_n. Equations (9) and (11) also do not imply a lower bound for all rational approximations to alpha.

The [companion direct-integral investigation](../joint-arctan-integral/DIRECT_INTEGRAL_GCD_THRESHOLD.md) identifies its positive arctangent integral with this same a=1 family. Consequently its natural sequence is covered by (10); these are two representations of one candidate, not independent approximation mechanisms.

That companion also proves a sharper asymptotic arithmetic criterion. For a=1 write A_n=n! ell_n and B_n=n! w_n, which are integers, and g_n=gcd(A_n,4B_n). Then q_n=A_n/g_n exactly, and

\[
Q_n\sim\frac{\sqrt2}{n}\left(\frac{n!}{g_n}\right)^2.
\]

Thus the possible sparse-subsequence improvement is precisely a large-gcd question: along chosen indices, Q_n tends to zero if and only if n!/(g_n sqrt(n)) tends to zero. Its fixed-index-gap corollary strengthens the asymptotic neighboring constants to 1/8 for gap one and 1/32 for gap two, and requires the consecutive index gaps of any successful subsequence to tend to infinity. These sharper asymptotic bounds are compatible with the uniform finite lower bounds (10)--(11) above.

## 5. What exact experiments measure

The companion verifier constructs all reduced fractions through n=100 for a=1,5,239. It encloses arctan(1/5) and arctan(1/239) by alternating series, and pi by the exact Machin identity

\[
\pi=16\arctan(1/5)-4\arctan(1/239).
\]

For completeness, tan(2 arctan(1/5))=5/12 and tan(4 arctan(1/5))=120/119. The tangent subtraction with 1/239 equals 1. The angle 4 arctan(1/5)−arctan(1/239) lies in (0,pi/2), which selects pi/4, proving the identity without a numerical pi value.

For each fraction the program records its reduced denominator and a strict rational interval for Q_n. It also computes the same-index Machin combination

\[
s_n=16R_n(5)-4R_n(239),
\]

reduces that rational number, and encloses its normalized error relative to pi. The general adjacent-product theorem (9) applies to each individual arctangent sequence; it is **not** silently applied to this combination. Any lower/upper bounds for the combination use the independently certified pi interval and its actual reduced fraction.

The finite output is a diagnostic of arithmetic cost, not a new pi approximation record or an infinite-subsequence theorem. The all-n conclusion (9) follows from (4), (6), and the integer cross-difference, not from the table.

For the pi family c=4,a=1, the unique least Q_n in the checked range 0<=n<=100 occurs at n=2 and gives p_n/q_n=3. Its quality lies strictly between 0.14159265358 and 0.14159265359. Selected later values, with decimal endpoints interpreted as exact rationals, are:

| n | Number of decimal digits in reduced q_n | Strict interval enclosing Q_n |
|---:|---:|---:|
| 10 | 7 | (2.1370771243, 2.1370771244) times 10^5 |
| 20 | 16 | (4.2529684335, 4.2529684336) times 10^15 |
| 50 | 45 | (1.1551236157, 1.1551236158) times 10^51 |
| 100 | 90 | (1.3230775974, 1.3230775975) times 10^103 |

For the same-index Machin combination, the unique minimum for 1<=n<=100 is n=1 with reduced fraction 3804/1195 and Q_n in (59507.150857,59507.150858). At n=100 its reduced denominator has 426 decimal digits and Q_n lies in (5.9883695662,5.9883695663) times 10^650. These are measurements of these specified candidate families only.
