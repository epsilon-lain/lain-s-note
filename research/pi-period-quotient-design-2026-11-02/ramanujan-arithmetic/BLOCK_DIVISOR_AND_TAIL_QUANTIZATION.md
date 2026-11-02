# A real Ramanujan block divisor and a charged tail construction

2026-11-02. This note examines an arithmetic structure in the Ramanujan bridge. It constructs a factor-preserving rational tail correction and records why the most immediate short-vector interpretation overstates its benefit. It does not replace the parent's main integer-kernel construction route.

## 1. Fixed series and an exponential integer divisor

Put b=4096, B_n=binom(2n,n), d_n=42n+5, and

\[
a_n=d_nB_n^3/b^n,\qquad
S_N=\sum_{n=0}^N a_n=C_N/b^N,
\quad C_N=bC_{N-1}+d_NB_N^3.
\]

Ramanujan's identity, in Guillera2018 formula(2), is S_infinity=16/pi. The primary publisher abstract for BingHe2017 states, for every integer n>1,

\[
\boxed{2nB_n\mid C_{n-1}.} \tag{1}
\]

Consequently B_n divides both C_(n-1) and C_n. This is a genuine exponential-size divisor, since B_n is asymptotic to 4^n/sqrt(pi n). It is not a cancellation in the original approximant 16b^N/C_N: that numerator is a power of two and has no odd B_n factor.

Sources checked: JesusGuillera, *A method for proving Ramanujan series for1/pi*, arXiv1807.07394, formula(2), https://arxiv.org/pdf/1807.07394 . BingHe, *On the divisibility properties concerning sums of binomial coefficients*, RamanujanJournal43(2017),313-326, DOI10.1007/s11139-016-9780-6, official statement https://link.springer.com/article/10.1007/s11139-016-9780-6 . The publisher's theorem statement was accessible; its full WZ proof was not accessed, so (1) is a cited theorem, not a newly supplied WZ certificate. The accompanying script separately checks (1) in its bounded examples.

## 2. A common factor in both coefficients of a new linear form

Let F(z)=sum d_nB_n^3 z^n. For integers u,v take Q(z)=u+vz and let P_N be QF truncated at degreeN. Evaluation at z=1/b gives

\[
\frac{16Q(1/b)}{P_N(1/b)}
=\frac{16b^{N-1}(bu+v)}{uC_N+vC_{N-1}}. \tag{2}
\]

Impose the integer congruence bu+v=B_N k. Equation(1) ensures that both integer coefficients in(2) are divisible by B_N. Using the recurrence for C_N, the exact quotient is

\[
\boxed{r_{N,u,k}
=\frac{16b^{N-1}k}
{kC_{N-1}+u d_NB_N^2}.} \tag{3}
\]

This common factor is arranged before integrating or using a numerical target. It is not merely dividing the numerator coefficient C_N by B_N while leaving the other coefficient unchanged.

In fact the common B_N in(3) follows already from the recurrence and bu+v=B_N k; it does not need the stronger theorem(1). The imported divisor additionally shows that, for k=1, the quotient denominator C_(N-1)+u d_NB_N^2 still contains B_N. Its odd part cannot cancel the final pure-two numerator. Thus this tail construction is not being presented as a way of converting all of He's divisor into denominator savings.

There is an important height charge. The congruence lattice has determinant B_N, but it contains the fixed zero-evaluation vector (u,v)=(1,-b). If k is nonzero, then

\[
\max(|u|,|v|)\ge B_N/(b+1).
\]

Indeed B_N<=|bu+v|<=(b+1)max(|u|,|v|). Thus a generic square-root-of-determinant short-vector estimate does not yield useful nonzero Q(1/b) at square-root cost. Short vectors in the fixed zero-evaluation direction supply no pi approximant. No asymptotic square-root denominator gain is claimed here.

## 3. An explicit noncircular tail rule

Set k=1 and define the exact hypergeometric tail factor

\[
R_N=\frac{\sum_{n=N}^\infty a_n}{a_N}.
\]

It is approximated rationally, not used as an unknown input. For n>=1, direct cross multiplication gives

\[
\frac{a_{n+1}}{a_n}<\frac1{64},
\]

because the difference is 168n^3+312n^2+132n-7>0. Take the least L>=1 satisfying

\[
4B_N\le b\cdot63\cdot64^{L-1}.
\]

Put

\[
G=\sum_{j=0}^{L-1}a_{N+j}/a_N,
\quad U=\frac1{63\cdot64^{L-1}}.
\]

Then G<R_N<G+U, and the interval [(B_N/b)G,(B_N/b)(G+U)] has width at most1/4. Choose u to be the nearest integer to its midpoint, using floor(midpoint+1/2). This guarantees

\[
\left|u-\frac{B_N}{b}R_N\right|\le\frac58. \tag{4}
\]

The generated candidate is explicitly

\[
\boxed{r_N=\frac{16b^{N-1}}{D_N},
\qquad D_N=C_{N-1}+u d_NB_N^2,} \tag{5}
\]

followed by complete gcd reduction. All input arithmetic is rational or integer and uses source coefficients only through index N+L-1. The upper width U is the proved geometric envelope, so it does not evaluate any additional source coefficient. The JSON records that exact maximum index, and the direct-truncation comparison consumes the same coefficients.

## 4. All-scale error and denominator costs

Let S_tilde=D_N/b^(N-1). Equations(4)-(5) imply

\[
|S_tilde-16/\pi|
\le\frac58\frac{d_NB_N^2}{b^{N-1}}.
\]

For N>=2, S_tilde>=5. Also S_infinity<5+47/504, while d_NB_N^2/b^(N-1) is decreasing from its N=2 value 89*36/4096<4/5. The successive ratio is less than32/b. Hence S_tilde<6. The resulting explicit ordinary-error certificate is

\[
\boxed{|\pi-r_N|\le
\frac25\frac{d_NB_N^2}{b^{N-1}}
\le\frac{8192}{5}(42N+5)256^{-N}.} \tag{6}
\]

Since the numerator in(5) is 2^(12N-8) and D_N<6b^(N-1)<2^(12N-8), its actual gcd is exactly 2^v_2(D_N). Thus

\[
\boxed{q_N=D_N/2^{v_2(D_N)}<6b^{N-1},}
\]
\[
\boxed{q_N^2|\pi-r_N|\le
\frac25 q_N^2\frac{d_NB_N^2}{b^{N-1}}.} \tag{7}
\]

The coarser substitution q_N<6b^(N-1) makes(7) grow. Actual gcds and exact tail intervals are retained in the script. This construction buys additional tail precision on a controlled coefficient grid; its quoted common factor does not by itself give a quadratic-quality theorem.

The same calculation with k=B_N or k=B_N^2 increases the correction grid precision and cancels further introduced B_N factors. At k=B_N^2 the denominator becomes C_(N-1)+u d_N, an integer quantization of the reciprocal-pi period. This explains why merely enlarging the forced block divisor can degenerate into ordinary reciprocal-digit rounding rather than a new approximation mechanism.

## 5. Bounded controls and recommendation

`block_tail_quantization.py` runs the fixed controls N=2,...,12 and exact rational error comparisons. The selected u is zero for N=2,...,6, so those outputs are exactly the earlier prefix and are not new candidates. Nonzero corrections occur for N=7,...,12. Their actual quadratic qualities are poor. Comparisons against direct truncation use the same maximum coefficient index used for G; U is a closed geometric bound requiring no additional term.

The useful arithmetic input is (1), together with the requirement that both coefficients of the final pi linear form share any proposed divisor. To turn it into a stronger construction, a paired-period or higher-degree numerator needs to preserve additional structure beyond the zero-evaluation congruence direction. The parent should not enlarge this elementary tail-quantization scan.
