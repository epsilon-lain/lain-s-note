# Sources and dependency boundary

This checkpoint studies only the fixed recurrence stated in README.md. Its new divisibility and formal-window statements have proofs in this directory; no external source is being credited with the unproved moving-prime bound.

## Earlier repository inputs

- [Odd prime blocks and excess](../pi-gcd-blast-2026-11-02/gcd-prime-blocks/ODD_PRIME_BLOCKS_AND_EXCESS.md): Lucas congruences, the leading-pole formula and odd-prime bounds.
- [Local poles and the excess sieve](../pi-gcd-blast-2026-11-02/gcd-integrality-obstruction/LOCAL_POLES_AND_EXCESS_SIEVE.md): the sharper computable gain divisor and whole-prime deletion rule.
- [Binary numerator valuation](../pi-gcd-blast-2026-11-02/gcd-valuation-theory/two-adic-proof/TWO_ADIC_NUMERATOR_VALUATION.md) and [binary gcd saturation](../pi-error-exponent-2026-11-02/exponent-audit/TWO_ADIC_GCD_SATURATION.md): binary main term and the single exceptional index.
- [Reduced denominator exponent map](../pi-error-exponent-2026-11-02/exponent-error-map/REDUCED_DENOMINATOR_EXPONENT_MAP.md): the quality normalization and asymptotic error relations.
- [Global first-layer budget](../pi-global-anomaly-2026-11-02/global-s-budget/GLOBAL_FIRST_LAYER_BUDGET.md): the preceding upper-half S(n) formulation.
- [Previous checkpoint](../pi-global-anomaly-2026-11-02/README.md): reflection and the already completed finite anomaly scan. That scan is not newly repeated or extended here.

These files were present in the parent checkpoint 370339b97a9705f612740e0717d3936ac58e502b. The new finite verifiers also compare the claimed integer divisibilities against directly generated A_n, B_n and their gcds.

## Primary external references checked on 2026-11-02

1. NIST Digital Library of Mathematical Functions, §27.12, Asymptotic Formulas: Primes. https://dlmf.nist.gov/27.12 . The prime number theorem gives theta(x)=x+o(x), equivalently by partial summation from the prime-counting form. It supplies total candidate prime weight, not any bound on primes dividing C_n.
2. Andrew Sutherland, MIT 18.783, Spring 2021, Problem Set 4, Problem 1. https://math.mit.edu/classes/18.783/2021/ProblemSet4.pdf . Hasse-invariant coefficient and Frobenius trace congruence used only in the midpoint appendix.
3. Andrew Sutherland, MIT 18.783, Spring 2021, Lecture 7, Theorem 7.3. https://math.mit.edu/classes/18.783/2021/LectureNotes7.pdf . Hasse's bound used to turn the trace congruence into exact trace zero for p>=5. The necessary p=7 mod 8 restriction is proved in the appendix using rational torsion and the duplication formula.

No prime-independence hypothesis, short-interval prime number theorem, or uniform Legendre-zero-distribution theorem is used. In particular, neither the external references nor the finite verification prove the missing S(n) or W(n) upper bound.
