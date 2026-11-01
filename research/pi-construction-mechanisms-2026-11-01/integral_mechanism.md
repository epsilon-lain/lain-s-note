# How the contour integrals manufacture usable linear forms in π

Research note, 1 November 2026. The purpose is a working model of the construction, with its arithmetic and nonvanishing requirements visible. It is not a certification of every cited proof.

## Sources and actual reading coverage

| Source | Primary text and precise location | Coverage in this investigation |
|---|---|---|
| Hata (1993), *Rational approximations to π and some other numbers* | [Publisher record](https://www.impan.pl/en/publishing-house/journals-and-series/acta-arithmetica/all/63/4/107787/rational-approximations-to-and-some-other-numbers); [full PDF](https://matwbn.icm.edu.pl/ksiazki/aa/aa63/aa6344.pdf). Printed pp. 335–349; §2, Lemmas 2.1–2.4; §3, equation (3.1). | All 15 pages read, including proofs. No independent validation of the drawn steepest-descent contours. |
| Salikhov (2008), *On the irrationality measure of π* | [Journal record](https://www.mathnet.ru/eng/rm9175); [indexed primary PDF](https://www.mathnet.ru/links/67585709df7b89fa677d7c21b02e2e64/rm9175_eng.pdf), printed p. 570, equations (2)–(4). | Search extraction of p. 570 read. Direct PDF requests failed. Pages 571–572 were **not** read. |
| Salikhov (2010), *On the measure of irrationality of the number π* | [MathNet record](https://www.mathnet.ru/eng/mzm8688); [publisher](https://link.springer.com/article/10.1134/S0001434610090294), pp. 563–573. | Bibliographic record and abstract only. MathNet full-text request returned 403; Springer displayed subscription access. The integral below is independently visible in the 2008 announcement and reproduced in ZZ. **No full-reading claim.** |
| Zeilberger–Zudilin (2020), abbreviated ZZ | [Author PDF](https://sites.math.rutgers.edu/~zeilberg/mamarim/mamarimPDF/pimeas.pdf#page=4); [arXiv](https://arxiv.org/abs/1912.06345). Author-PDF pp. 4–12: Part II, equations (1)–(18), Lemmas 1–6 and Propositions 1–2; p. 13 final rates. | Entire 13-page paper read. The external Maple recurrence/certificate files were not audited. The paper itself reports correction of a false earlier termwise divisibility claim, pp. 2 and 9. |
| Bai (2026), v2 | [Full HTML](https://arxiv.org/html/2609.11276v2), §§1, 2, 4, 5; [abstract/version record](https://arxiv.org/abs/2609.11276v2). | Existence, author, version date and actual construction verified. §§1, 2, 4 and 5 read; selected inputs/outputs of §3 read. §6 and Appendix A **not** read in full or independently checked. The many-digit enclosures and local-minimum claim remain source claims here. |

Full-text caches, where obtained, are in `source_context/pi_literature/`; they are separate from these original notes. The reproducible script `integral_lab.py` uses no cached source text.

## 1. First design principle: poles decide which constants can appear

Here is the model we actually compute, a convenient sign normalization of the ZZ integral:

\[
R_n(t)=\frac{5t^{2n}(t^4+6t^2+25)^{2n}}{(25-t^2)^{3n+1}},\qquad
J_n=i\int_{-1-2i}^{-1+2i}R_n(t)\,dt.
\]

The source's `I_n` is \((-1)^{n+1}J_n\); absolute-value asymptotics agree. Put \(a=-1-2i\), \(b=\bar a\). All poles are at the rational points \(t=\pm5\), safely outside the integration region. The rational function is even. Polynomial division and uniqueness of principal parts therefore force

\[
R_n(t)=P_n(t)+\sum_{j=0}^{3n}A_{j,n}
 \left((5+t)^{-j-1}+(5-t)^{-j-1}\right),
\quad P_n\in\mathbb Z[t^2].
\]

The two equal pole coefficients are the mechanism, rather than a cosmetic symmetry. Terms with \(j\ge1\) have rational-function primitives. Their endpoint values are Gaussian rationals; conjugacy makes `i` times their difference rational. The same is true of the polynomial primitive.

Only the simple poles can produce logarithms. Their sum has integral

\[
\int_a^b\left(\frac1{5+t}+\frac1{5-t}\right)dt
=\log\frac{(4+2i)(6+2i)}{(4-2i)(6-2i)}
=\frac{\pi i}{2}.
\]

This identity uses the continuous logarithms along the segment. Indeed, the two numerator factors multiply to \(20+20i\), and the denominator is its conjugate. Their arguments differ by exactly \(\pi/2\). Consequently

\[
J_n=u_n+v_n\pi,\qquad u_n,v_n\in\mathbb Q,
\qquad v_n=-A_{0,n}/2.
\]

Thus the denominator locates the logarithms, the involution pairs their coefficients, and the endpoints make the surviving angle a rational multiple of π. None of these steps establishes smallness or integer coefficients yet.

An exact ablation exposes the role of evenness. With \(n=0\),

\[
i\int_a^b\frac5{25-t^2}dt=-\frac\pi4.
\]

Multiply the integrand by \(t\), breaking evenness. The same poles now give

\[
i\int_a^b\frac{5t}{25-t^2}dt
=-\frac{5i}{2}\log\frac{28+4i}{28-4i}
=5\arctan(1/7).
\]

That unwanted angle is not removed by merely keeping rational poles. The equal-residue pairing was doing essential work.

## 2. An exact small instance

For \(n=1\), direct rational-function algebra gives

\[
P_1(t)=5t^2+560,\qquad
(A_{0,1},A_{1,1},A_{2,1},A_{3,1})=(-2392,6880,-11200,8000).
\]

Integrate each term using

\[
\int\left((5+t)^{-j-1}+(5-t)^{-j-1}\right)dt
=\frac{(5-t)^{-j}-(5+t)^{-j}}j\quad(j\ge1).
\]

The result is

\[
J_1=-\frac{11272}{3}+1196\pi
=0.01148036005937986798815307\ldots.
\]

The paper's normalization at this index is \(M_1=12\), so

\[
12J_1=-45088+14352\pi\in\mathbb Z+\mathbb Z\pi.
\]

The script verifies the full decomposition and rationality exactly, then displays the decimal residual. It also checks indices 0–4 and the adjacent determinants for those indices only. Small numerical residuals alone are not an asymptotic proof.

## 3. Why these numerator zeros are useful twice

The quartic factors as

\[
t^4+6t^2+25=((t+1)^2+4)((t-1)^2+4).
\]

It puts zeros at both integration endpoints and at their reflected partners. Raising it to \(2n\) gives endpoint zeros of order \(2n\). The factor \(t^{2n}\) adds a central zero. There are two separate benefits.

First, these zeros make the exponential kernel small on a carefully deformed contour through complex stationary points. Pole orders rise simultaneously, so the claim needs an actual extremal-contour or recurrence calculation; endpoint zeros by themselves do not guarantee a useful rate.

Second, they control the arithmetic of the polynomial quotient. Expand the quotient at \(t_0=-1-2i\):

\[
P_n(t)=\sum_{k=0}^{4n-2}B_k(t-t_0)^k.
\]

For \(k<2n\), the corresponding derivative of \(R_n\) vanishes. Differentiating the partial-fraction identity expresses \(B_k\) entirely through the pole coefficients and the Gaussian endpoint distances. Their factors include

\[
4-2i=2(2-i),\qquad 6+2i=2(1+i)(2-i).
\]

Large powers of 2 in the pole coefficients therefore force powers of 2 in \(B_k\). On integrating, the endpoint difference is \(b-t_0=4i\), which contributes another \((4i)^{k+1}\). These powers offset the divisions by \(k+1\). This is why the endpoint-zero condition contributes to denominator removal as well as to analytic decay.

## 4. Division costs are part of the approximation, not bookkeeping

The exact arithmetic theorem in ZZ Part II is

\[
M_nJ_n\in\mathbb Z+\mathbb Z\pi,\qquad
M_n=2^{-\lfloor5n/2\rfloor+2}
\frac{\operatorname{lcm}(1,\ldots,4n)}{\Phi_n},
\]

where

\[
\Phi_n=\prod_{p\in\mathcal P_n}p,\qquad
\mathcal P_n=\left\{p>\max(5,\sqrt{3n}):
\tfrac12\le\{n/p\}<\tfrac23\right\}.
\]

The prime symbol in this set means primes. `M_n` can be rational: the requirement is that it normalizes both coefficients, not that the multiplier itself be an integer.

There are three sources of division: derivatives/Laurent coefficients; integration of pole terms by \(j\); integration of the polynomial by \(k+1\). A proof must cover all three. A common divisor of the π coefficient alone is insufficient.

The removable primes are justified at the coefficient indices that would otherwise require a factor \(p\) in an integration denominator. The corrected ZZ Lemma 2 proves the divisibility for the **summed coefficient**, using a transformed binomial sum. It does not assert that every original Leibniz summand has that divisor. Negative Laurent indices matter because they also govern the polynomial quotient. Two expressions for that quotient then combine its prime divisibility with its 2-adic divisibility.

Here is the resulting rate balance, recomputed numerically by the script from the source's algebraic saddle values. Write

\[
r=\log(21851.691396\ldots),\quad
s=\log(0.029458495928\ldots),\quad
\phi=\frac\pi{2\sqrt3}-\log\frac{3\sqrt3}{4}=0.64527561\ldots.
\]

The raw estimates are coefficient growth `r` and integral decay `s`. The remaining normalization cost is

\[
c=4-\frac52\log2-\phi.
\]

Then \(\sigma=r+c\), \(\tau=-s-c\). The cost raises the coefficient growth and consumes part of the decay simultaneously.

| Arithmetic retained | Normalization exponent `c` | Normalized decay `τ` | Consequence |
|---|---:|---:|---|
| Least common multiple alone | 4 | −0.4752270760… | These estimates give no shrinking integer forms. |
| Also divide the proved powers of 2 | \(4-\tfrac52\log2\) | 1.2576408754… | Finite bound 10.74774746567… |
| Also remove the proved prime product | \(4-\tfrac52\log2-\phi\) | 1.9029164856… | \(1+\sigma/\tau=7.103205334137\ldots\). |

This is a quantitative ablation: deleting the arithmetic saving can destroy the argument even when the raw integral is exponentially small.

## 5. Nonvanishing and coefficient growth: a further exact model

The statement `v_n≠0` can be proved directly, without inferring it from a few numerical approximants. Adapt the residue substitution used in Bai §4.1:

\[
z=\frac{t+5}{t-5},\qquad t=-5\frac{1+z}{1-z},\qquad
P(z)=2+6z+9z^2+6z^3+2z^4,
\]

\[
S(z)=\frac{(1+z)^2P(z)^2}{(1-z)^4}.
\]

Direct substitution in our ZZ model yields

\[
R_n(t)dt=\frac{(-1)^n4^n}{2}
z^{-3n-1}S(z)^n\,dz.
\]

Residues and \(v_n=-A_{0,n}/2\) now give the exact formula

\[
\boxed{v_n=(-1)^{n+1}4^{n-1}[z^{3n}]S(z)^n.}
\]

All coefficients of `S` are positive. Therefore \(v_n\ne0\) for every \(n\), and since π is irrational, \(J_n\ne0\). This uses the already known irrationality of π; it is not a new standalone irrationality proof.

The same substitution explains the ordinary growth limit. Choose the unique \(\zeta\in(0,1)\) for which \(\zeta S'(\zeta)/S(\zeta)=3\). The tilted probabilities \(\Pr(X=k)=[z^k]S(z)\zeta^k/S(\zeta)\) have mean 3, positive variance and span one. Coefficient extraction gives

\[
[z^{3n}]S(z)^n=S(\zeta)^n\zeta^{-3n}
\Pr(X_1+\cdots+X_n=3n).
\]

A local central limit estimate contributes only a factor of order \(n^{-1/2}\). Hence the exponential growth is exactly \(\log4+\log S(\zeta)-3\log\zeta\). This positivity route prevents a hidden cancellation in the dominant coefficient asymptotic. The source's integral rate is still a separate analytic input.

## 6. Exactly what converts forms into a finite irrationality bound

Here is a sufficient criterion with every required assumption explicit. Suppose rational forms \(j_n=a_n+b_n\theta\) have a normalization \(m_n\) such that

\[
\ell_n=m_nj_n=U_n+V_n\theta,\qquad U_n,V_n\in\mathbb Z.
\]

Assume, for every small \(\varepsilon>0\) and all sufficiently large `n`,

\[
|V_n|\le e^{(\sigma+\varepsilon)n},\qquad
|\ell_n|\le e^{-(\tau-\varepsilon)n},\qquad\tau>0,
\]

with \(\sigma>0\), and additionally

\[
D_n=U_nV_{n+1}-U_{n+1}V_n\ne0.
\]

For an arbitrary rational \(p/q\), put \(T_n=qU_n+pV_n\in\mathbb Z\). Consecutive `T` values cannot both vanish, since that would annihilate the nonzero determinant. Choose `n` just large enough that \(q|\ell_n|<1/2\), and use whichever of `n,n+1` has \(T\ne0\). Then

\[
1\le|T_j|\le q|\ell_j|+|V_j|\,|q\theta-p|,
\]

so

\[
|\theta-p/q|\ge\frac1{2q|V_j|}
\ge C_\varepsilon q^{-1-(\sigma+\varepsilon)/(\tau-\varepsilon)}.
\]

Therefore \(\mu(\theta)\le1+\sigma/\tau\). Determinants are a convenient sufficient hypothesis, not an automatic consequence of nonzero forms. The computation checks finitely many ZZ determinants and makes no all-index assertion.

There is another route, actually used by Hata's lemma and Bai §5: require an **ordinary** limit \(n^{-1}\log|V_n|\to\sigma>0\), eventual \(V_n\ne0\), and a decay limsup. If the integer `T` happens to vanish at the first eligible index, select its next nonzero index. The previous zero relates \(|q\theta-p|\) to the ratio \(|\ell|/|V|\); the two-sided coefficient-growth estimate bounds the possible gap. This yields the same exponent. Simply claiming that the first integer is nonzero leaves a real logical gap.

If \(n^{-1}\log|m_n|\to c\), \(n^{-1}\log|b_n|\to r\), and \(\limsup n^{-1}\log|j_n|\le s\), the usable rates are

\[
\sigma=r+c,\qquad \tau=-s-c>0,\qquad
\boxed{\mu(\theta)\le1+\frac{r+c}{-s-c}.}
\]

An infinite set of good approximations without index control does not by itself imply this upper bound.

## 7. What changes along the cited chain

Hata's starting kernel is \(F(z)=\prod_{j=1}^3(z-a_j)^2/z^3\), integrated as \(F(z)^n dz/z\). Laurent expansion leaves just one logarithm. For \((a_1,a_2,a_3)=(1,2,1+i)\), the endpoint logs are \(\log2\) and \(\tfrac12\log2+i\pi/4\); a combination of paths cancels \(\log2\). His binomial-divisibility lemma lowers the integration denominator cost, while distinct saddles control errors and coefficients. This yields his stated π bound through a linear-independence estimate.

Salikhov's announced kernel places five numerator zeros symmetrically about 5, with exponents `3n` against pole orders `5n+1` and even `n`. This is the source of the paired-pole architecture above. ZZ changes the exponent ratio to `2n` against `3n+1`, proves the required common-divisor statement and obtains its stated bound. Those are parameter changes coupled to new arithmetic proofs, rather than a consequence of numerator decay alone.

Bai v2 introduces independent exponents into the same rational kernel:

\[
R_n(t)=\frac{5t^{2an}(t^4+6t^2+25)^{bn}}{(25-t^2)^{cn+1}},\qquad
(a,b,c)=(1857,3714,5570).
\]

Its normalization is \(2^{4-4645n}\operatorname{lcm}(1,\ldots,7430n)/\Phi_n\). An especially useful mechanism is the finite-field proof of removable primes: reduce \(R_n=U^pH\), show the polynomial `H` has no monomials of degree \(-1\bmod p\), integrate it in characteristic `p`, and conclude \(R_n=(U^pV)'\). Its Laurent coefficient of degree \(-j-1\) is then divisible by `p` whenever `p|j`. This covers the negative indices too. The paper claims rates \(r\approx3.3319439862\), \(s\approx-1.1745439413\), \(c\approx0.5399938192\), giving \(\mu(\pi)<7.101862832357\). Its strict two-dimensional local-minimum claim concerns its explicit auxiliary function in a stated parameter chamber; it does not certify a global optimum.

## Verification status

**Exact checks supplied here:** paired-pole decompositions and integer normalizations for indices 0–4; the displayed `n=1` integral; the elementary symmetry ablation; the algebraic residue substitution and its positive coefficients. **Theorems used from sources:** the all-index arithmetic savings and analytic decay rates for ZZ; the prime-number-theorem conversion of prime products. **Numerical evidence only:** script decimal residuals, finite adjacent determinants and recomputed rate balances. **Not certified here:** the external ZZ recurrence certificate, every steepest-descent detail of Hata, Salikhov's full 2010 proof, Bai's Appendix A interval calculations or its §6 local-minimum proof.
