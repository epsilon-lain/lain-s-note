# Independent finite reflection-factor certificate audit

2026-11-02. This audit does not rely on the supplied reviews or on a symbolic factorization routine. The input files were read without modification.

**Result:** all ten stored reflection factors, for `n=4,6,8,10,12` and both signs, are exactly the normalized determinants claimed in `REFLECTION_HEIGHT.md`. All ten are irreducible over the rationals. This is a finite result only.

## Direct verification of every stored coefficient

For `Z=2x`, multiply row `s` of the defining reflection matrix by `2^s s!`. Its entry in column `a` is the integer

\[
\sum_{k=0}^{\min(a,s)}
\binom ak x^{a-k}2^k\frac{s!}{(s-k)!}
\bigl(1+\varepsilon(-1)^{a+s-k}\bigr).
\]

An independent integer Bareiss implementation evaluates its determinant. Division by the known row-scaling product `2^(n(n-1)/2) A_n` recovers the defining `L_epsilon(2x)`.

Every column `a` has polynomial degree at most `a`, giving the unconditional determinant degree bound `D=n(n-1)/2`. The proposed normalized expression `sign * Z^r F_epsilon(Z)/B_n` has degree `r^2<=D`. Equality was checked at all `D+1` distinct arguments `Z=2,4,...,2(D+1)`. This certifies the entire polynomial identity, rather than only selected values or coefficients.

Across the ten fixtures there were 330 exact determinant evaluations. The computed `B_n=A_n/A_r^2`, leading coefficients, degrees, evenness, and signed constant terms all match the input. In particular the signs of the constants are:

| n | sign of F_+(0) | sign of F_-(0) |
|---|---|---|
| 4 | + | − |
| 6 | − | − |
| 8 | − | + |
| 10 | + | + |
| 12 | + | − |

Their magnitudes are `(n-1)B_n` and `B_n` respectively, as claimed.

## Independent irreducibility verification

For each of the 22 supplied prime certificates, the audit independently checks primality, squarefreeness of the reduced fixture, and the exact product of the supplied monic modular factors. Every one of the 71 modular factors was proved irreducible using the implemented Rabin criterion: for degree `d`, verify `X^(p^d)=X` modulo the factor, and verify that the factor is coprime to `X^(p^(d/l))-X` for every prime divisor `l` of `d`.

The independent subset-sum computation reproduces every listed surviving factor-degree set. The final intersection is empty for each fixture. A monic rational factor would be monic integral and preserve its degree after reduction at each prime, so this empty intersection rules it out. There is no circular reliance on the input field `Q_irreducible` or on SymPy's `is_irreducible`.

The certificate prime lists are:

| n | F_+ primes | F_- primes |
|---|---|---|
| 4 | 7 | 5 |
| 6 | 7, 11 | 7, 13 |
| 8 | 11, 17 | 11, 13 |
| 10 | 11 | 11, 13, 19, 37 |
| 12 | 13, 19 | 13, 17, 19, 23, 29 |

The supplied verifier's subset-degree method is mathematically valid. Its use of a factorization library is independently covered here by modular multiplication and Rabin tests. The finite evidence does not establish irreducibility for untested dimensions, or the height and degree of a particular irreducible root factor in general.

Reproduce with `python reviews/finite-height/independent_finite_certificates.py`. Machine-readable results are in `independent_finite_certificates.json` beside the script. The implementation uses only Python's standard library and exact integer/modular arithmetic.
