# Sources and exact proof dependencies

2026-11-02. Starting repository head: `db233082d7ee811cccc8d2cbd376a3daec6b89fb`.

The current request proposes the pointwise first-layer bound S(n)≤n/6 and supplies the n≤10000 finite maximum ln17/21. The budget implication was rederived; two independently implemented modular scans confirm that supplied finite range. The proposed eventual bound is not a premise of any unconditional result here.

## Fixed-sequence proof inputs

1. [Odd-prime Wronskian caps and local pole identity](../pi-gcd-blast-2026-11-02/gcd-prime-blocks/ODD_PRIME_BLOCKS_AND_EXCESS.md).
2. [Top-half prime-square classification](../pi-gcd-prime-adic-2026-11-02/gcd_blocks/PRIME_BLOCKS_AND_LIFTING.md).
3. [Exact all-index binary valuation of A_n](../pi-gcd-blast-2026-11-02/gcd-valuation-theory/two-adic-proof/TWO_ADIC_NUMERATOR_VALUATION.md).
4. [Binary saturation of the actual gcd](../pi-error-exponent-2026-11-02/exponent-audit/TWO_ADIC_GCD_SATURATION.md), with sole exceptional index n=3.
5. [Positive-integral error, factorial normalization and sharp asymptotics](../pi-joint-arithmetic-2026-11-02/joint-arctan-integral/DIRECT_INTEGRAL_GCD_THRESHOLD.md).
6. [Uniform Legendre coefficient bounds](../pi-error-exponent-2026-11-02/exponent-prime-budget/NET_PRIME_BUDGET_AND_ERROR_EXPONENT.md).

These paths already existed at the starting commit; relative links work in a full checkout. The new edge-height estimate needs the numerator's binary theorem but does not require gcd saturation. The quality-exclusion implication uses the latter as well.

## Primary and official external references

Tony D. Noe, *On the Divisibility of Generalized Central Trinomial Coefficients*, Journal of Integer Sequences 9 (2006), Article 06.2.7. [Publisher PDF](https://cs.uwaterloo.ca/journals/JIS/VOL9/Noe/noe35.pdf), Sections 6–8, provides the classical Schur/Holt framework for Legendre and scaled Legendre coefficients. The specialized coefficient reversal is written out in this checkpoint. The folded-degree and edge-capacity deductions are direct arguments for the specified fixed sequence; no uniform diagonal bound is attributed to that paper.

[NIST DLMF §27.12](https://dlmf.nist.gov/27.12), especially equations 27.12.4–27.12.5, records the classical prime number theorem. Its consequences θ(x)=x+o(x) and ψ(x)=x+o(x) justify cancellation of the explicit global prime weights. They do not control Legendre zero events at changing moduli. No prime number theorem is needed for the finite divisibility or edge-capacity certificates.

Both sources were directly checked. No secondary search excerpt is used as the premise of a proof.

## Scope distinctions

- The full finite scan covers indices n≤100000 in the specified family, not all rational denominators up to 100000.
- First-layer p∣A_(n−p) permits gcd valuation one or two. Only valuation two is a positive excess beyond n!, and both layers are recorded for dangerous clusters.
- The sparse abstract zero-set model in the independent review is a logical obstruction to an inference from local properties. It is not a counterexample for the actual recurrence.
- The uniform S(n) bound, the resulting exclusion of every sparse successful subsequence, and any global assertion about π remain unproved.
