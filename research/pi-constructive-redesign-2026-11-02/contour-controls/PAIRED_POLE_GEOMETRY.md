# A constructive paired-pole family with arithmetic endpoint design

2026-11-02. This is a new project experiment using a classical contour-integral mechanism. It does not claim a new pi record or an infinite quality-zero sequence. No old prime-pool exclusion argument is advanced here.

## Primary-source mechanism and scope

Zeilberger and Zudilin, *The irrationality measure of pi is at most 7.103205334137...*, Moscow Journal of Combinatorics and Number Theory 9 (2020), 407–419, provide a Salikhov-type symmetric rational contour construction, arithmetic coefficient lemmas, and exact-recurrence software. Official author page and paper:

- https://sites.math.rutgers.edu/~zeilberg/mamarim/mamarimhtml/pimeas.html
- https://sites.math.rutgers.edu/~zeilberg/mamarim/mamarimPDF/pimeas.pdf

Their Part II uses a quartic numerator, a central even power, and paired real poles. This is a primary source for the mechanism, not a theorem that this family produces quality tending to zero. The source's exponent 7.1032 is an upper bound on irrationality measure, not the critical constructive estimate needed here. The general endpoint parametrization and elementary universal denominator bound below are derived directly in this project.

## 1. Exact family: choose the arithmetic geometry first

Let s,t be positive coprime integers, and put

\[
L=2s^2+2st+t^2,\quad u=s(2s+t),\quad v=st,\quad z=u+iv.
\]

The two real poles will be 0,L. Direct calculation gives

\[
\frac z{L-z}=\frac st(1+i),\qquad
|z|^2=2s^2L,\quad |L-z|^2=t^2L.
\tag{1}
\]

Thus the norm factors that occur in the endpoint evaluations are explicitly known before constructing any approximant.

For integers a≥1,b≥0,c≥0 define

\[
H(x)=((x-u)^2+v^2)((x-(L-u))^2+v^2),
\]
\[
R(x)=\frac{H(x)^a(2x-L)^{2b}}
{x^{c+1}(L-x)^{c+1}},\qquad
I=\frac1i\int_{u-iv}^{u+iv}R(x)\,dx.
\tag{2}
\]

Use the straight vertical path. Since 0<u<L, it does not meet either pole. All coefficients of R are rational, R(L-x)=R(x), and conjugating the oriented integral proves I is real.

The numerator degree is 4a+2b, the denominator degree is 2c+2. The rational function is proper precisely when c≥2a+b. Otherwise the polynomial part has degree d=4a+2b−2c−2≥0. No properness assumption is needed for the formulas below; polynomial terms are integrated explicitly. Fixed-ratio scaling a=αn,b=βn,c=γn keeps all degrees O(n).

## 2. The pi period and the branch check

Symmetry gives the partial-fraction expansion

\[
R(x)=P(x)+\sum_{j=0}^{c}r_j
\left(\frac1{x^{j+1}}+\frac1{(L-x)^{j+1}}\right),
\tag{3}
\]

where P is an integer polynomial (the denominator has leading coefficient ±1). Let N(x)=H(x)^a(2x-L)^{2b}. Then

\[
r_j=[x^{c-j}]\frac{N(x)}{(L-x)^{c+1}}.
\tag{4}
\]

For the simple-pole pair the primitive is Log x−Log(L-x). Both logarithms are continuous on the path because their arguments have positive real part. At the upper endpoint their imaginary-part difference is π/4 by (1); the lower endpoint is its conjugate. Therefore the change is exactly iπ/2, with no unspecified multiple of 2πi. Equivalently, the quotient has positive real part everywhere along this path, so its principal Log has no crossing.

Higher poles and P have rational endpoint contributions after division by i, because z and L−z are Gaussian integers and conjugation pairs the endpoints. Consequently

\[
\boxed{I=A\pi+B,\qquad A=r_0/2,\quad A,B\in\mathbb Q.}
\tag{5}
\]

An explicit formula used in the code is

\[
B=\sum_k P_k\frac{2\operatorname{Im}(z^{k+1})}{k+1}
+\sum_{j=1}^{c}\frac{2r_j}{j}
\operatorname{Im}\bigl((L-z)^{-j}-z^{-j}\bigr).
\tag{6}
\]

## 3. Precharged arithmetic cost and the exact objective

Put K=lcm(2s²,t²), M=max(c,d+1,1), and ell=lcm(1,...,M), treating the polynomial degree term as absent if P=0. The residue formula gives

\[
\operatorname{den}(r_j)\mid L^{2c+1-j}.
\]

By (1), an elementary common denominator is

\[
\boxed{D_0=2L^{2c+1}K^c\ell,\qquad D_0A,D_0B\in\mathbb Z.}
\tag{7}
\]

This is deliberately coarse. It makes the geometry's cost explicit; actual residue divisibility or the source's more specific arithmetic lemmas can lower the cost considerably. Merely making an endpoint closer to a useful saddle does not establish an arithmetic benefit.

For either D=D0 or the actual joint denominator, set α=DA,β=−DB,g=gcd(|α|,|β|). If A≠0, the reduced approximant is β/α, with denominator q=|α|/g. Exactly

\[
\boxed{q^2|\pi-\beta/\alpha|
=\frac{D^2|A|}{g^2}|I|.}
\tag{8}
\]

No comparison is based only on the size of I.

## 4. Remove false geometric improvements

A common dilation of s,t merely dilates the geometry and scales I; it does not change β/α. Coprimality removes this redundancy.

There is another exact equivalence:

\[
(s,t)\mapsto(t,2s),\qquad
L' =2L,\quad z'=2(L-\bar z),
\]
\[
\boxed{I' =2^{4a+2b-2c-1}I.}
\tag{9}
\]

Hence the approximant is unchanged. In particular, (s,t)=(1,1) and (1,2) give the same reduced fraction for every a,b,c. A smaller displayed raw denominator for one of those representations is not an arithmetic improvement.

A genuine geometry parameter is

\[
\kappa=t/s+2s/t\ge2\sqrt2,\qquad v/L=1/(\kappa+2).
\tag{10}
\]

The old geometry has κ=3. The coprime pair (2,3) has κ=17/6, while (5,7) has κ=99/35. These are inequivalent to κ=3. Pell pairs t²−2s²=±1 approach the limiting geometry, but their larger norm factors must be paid for in (7)–(8). This is a construction parameter worth testing, not a proof that it helps.

## 5. Finite executed controls

`endpoint_design.py` scanned s,t≤4 coprime, a≤3,b≤3,c≤6,a+b+c≤10: 842 evaluable positive rational candidates. The best normalized quality was again the known 355/113, not a new record. The classic fraction emerged without being input as a target.

`geometry_compare.py` compared 3 exponent ratios at n=1,2,4,8 and four endpoint choices. All coefficients and reduced fractions are exact. Machin alternating-series rational intervals certify errors; the displayed logarithms below are floating summaries, not rigorous log enclosures.

For (a,b,c)=(n,n,2n), n=8:

| (s,t) | κ | log(error)/n | log(q)/n | log(Q)/n |
|---|---:|---:|---:|---:|
| (1,2) | 3 | −7.90918 | 5.51056 | 3.11194 |
| (1,1) | 3 | −7.90918 | 5.51056 | 3.11194 |
| (2,3) | 17/6 | −9.04791 | 9.45875 | 9.86959 |
| (5,7) | 99/35 | −9.08316 | 13.32325 | 17.56333 |

This isolates an actual geometric improvement in error, and also measures the larger denominator cost that overwhelms it in these tested cases. It supplies no reason to promote these naive Pell substitutions as a successful new route.

The original fixed geometry also received 175 small-parameter controls and 13 fixed-ratio scale controls through n=12. None of those scale controls is asserted asymptotically successful. They are comparison data, not an additional exclusion theorem.

Eight independent 100-decimal numerical contour quadratures agree with the exact Aπ+B formulas. Exact checks verify (9) and the old/new normalization. Numerical quadrature supplements the algebraic proof; it is not its foundation.

## 6. Concrete remaining constructive step

The useful delivered change is that the geometry and its norm factorization are now construction variables, with an exact joint objective. The next test must add a divisor mechanism tied to these chosen Gaussian norms, or a multi-parameter coefficient identity that cancels them, before scaling the Pell geometry. Repeating only the exponent or endpoint scans has not yet supplied that mechanism.

No uniform prime-pool bound is pursued here. No new quality-zero subsequence has been found. The reproducible modules and finite controls make those limits explicit.
