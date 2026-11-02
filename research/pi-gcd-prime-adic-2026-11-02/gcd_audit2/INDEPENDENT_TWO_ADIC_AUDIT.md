# Independent audit of the exact 2-adic valuation

2026-11-02. Reviewed `../gcd_arithmetic/TWO_ADIC_VALUATION.md` for the exact recurrence
\(A_0=A_1=1\), \(A_{n+1}=(2n+1)A_n+n^2A_{n-1}\).

**Verdict:** no mathematical gap found. In particular the cancellation identity used for the residue class \(n\equiv3\pmod4\) is justified as a 2-adic identity, without importing a real integral value.

## Exact valuation statement

For every integer \(k\ge0\), the proof establishes
\[
v_2(A_{4k})=2k,\qquad v_2(A_{4k+1})=2k,
\]
\[
v_2(A_{4k+2})=2k+2,\qquad
v_2(A_{4k+3})=2k+3+v_2(k+1).
\]
The binomial reversal formulas agree with the coefficient formula for
\((1-2z-z^2)^{-1/2}\), including all powers of two and the odd-index prefactor \(2m+1\).

## The potentially delicate infinite identity

Let
\[
T(x)=\sum_{\ell\ge0}\frac{\ell!}{(2\ell+1)!!}\binom{x}{\ell}^2.
\]
The generalized binomial coefficients are 2-adic integers for \(x\in\mathbb Z_2\), so the terms tend uniformly to zero by \(v_2(\ell!)\to\infty\). At \(x=-1\), the binomial square is one.

To prove \(T(-1)=0\), algebraically integrate
\[
\sum_{\ell\ge0}2^\ell x^\ell(1-x)^\ell
=\frac1{1-2x+2x^2}
\]
between 0 and 1. Each summand is a degree-\(2\ell\) polynomial whose coefficients have valuation at least \(\ell\). Integration followed by evaluation at 1 loses at most \(\lfloor\log_2(2\ell+1)\rfloor\). Hence the resulting series converges, and termwise integration is valid; its term is
\(2^\ell(\ell!)^2/(2\ell+1)!=\ell!/(2\ell+1)!!\).

The factorization in \(\mathbb Q_2(i)\) is
\[
1-2x+2x^2=(1-(1+i)x)(1-(1-i)x).
\]
Both \(1+i\) and \(1-i\) have valuation \(1/2\), so the logarithm primitives converge at \(x=1\). Their endpoint logarithms are \(\log(i)\) and \(\log(-i)\). The convergent 2-adic logarithm is additive on \(1+\mathfrak m\), and these fourth roots of unity belong to that domain. Therefore \(4\log(\pm i)=\log1=0\), and both logarithms vanish. The proposed primitive has the correct derivative and value zero at zero. This proves the identity and resolves the possible interchange objection.

## Divided differences and finite remainder

For odd positive \(m\), set \(r=v_2(m+1)\). The binomial addition formula with upper arguments \(-1\) and \(m+1\), together with
\[
\binom{m+1}{j}=\frac{m+1}{j}\binom m{j-1},
\]
gives
\[
v_2\left(\binom m\ell-\binom{-1}\ell\right)
\ge r-\lfloor\log_2\ell\rfloor.
\]
The sum of the two binomials is also integral, so this bound applies to the squared difference. Division by \(m+1\) and multiplication by \(\ell!/(2\ell+1)!!\) gives the tail bound in the author note. For \(\ell\ge6\), its valuation is at least two and tends to infinity. Thus the whole omitted tail lies in \(4\mathbb Z_2\), including the terms \(\ell>m\) originating from \(T(-1)\).

I independently expanded the first five divided differences at \(m=2u-1\) and verified the displayed rational polynomial identity exactly with SymPy. Its denominator 155925 is odd and every coefficient after subtracting 2 has a factor 4. Hence
\(T_m/(m+1)\equiv2\pmod4\), proving the exact odd-case valuation.

The parity/modulo-four proof of the even-index sum and the translation through Legendre's factorial formula are also correct.

## Arithmetic consequence and its boundary

For \(g_n=\gcd(A_n,4B_n)\), define the positive odd-prime surplus
\[
E_n=\prod_{p\text{ odd prime}}
p^{\max(0,v_p(g_n)-v_p(n!))}.
\]
Since \(g_n\mid A_n\), the theorem yields
\[
\frac{g_n\sqrt n}{n!}
\le2(n+1)^2\sqrt n\,2^{-n/2}E_n.
\]
All negative odd-prime valuation differences were discarded only to obtain this upper bound. Thus any successful subsequence, on which the left side tends to infinity, necessarily satisfies
\[
\log E_n\ge\frac{\log2}{2}n-O(\log n)
\]
along that subsequence. An estimate \(E_n=\exp(o(n))\) would exclude success for this specified family. The proof does not establish that missing estimate, nor a lower bound for all rational approximations to \(\pi\).
