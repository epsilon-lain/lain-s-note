# Sources, dependencies, and scope

Checked 2026-11-02. Repository parent at the start of this checkpoint:
`873b5c24420869a4f1a06b0e7e8c5f64566838a8`.

## Exact previous proof inputs

1. [Direct positive integral, precise asymptotics, and actual reduced gcd threshold](../pi-joint-arithmetic-2026-11-02/joint-arctan-integral/DIRECT_INTEGRAL_GCD_THRESHOLD.md).
2. [Legendre–arctangent identities and exact all-index error bounds](../pi-joint-arithmetic-2026-11-02/joint-log-pade/LEGENDRE_ARCTAN_JOINT_COST.md).
3. [All-index 2-adic numerator valuation](../pi-gcd-blast-2026-11-02/gcd-valuation-theory/two-adic-proof/TWO_ADIC_NUMERATOR_VALUATION.md).
4. [Odd-prime Wronskian caps, local poles, and integrality](../pi-gcd-blast-2026-11-02/gcd-integrality-obstruction/LOCAL_POLES_AND_EXCESS_SIEVE.md) and [prime blocks](../pi-gcd-blast-2026-11-02/gcd-prime-blocks/ODD_PRIME_BLOCKS_AND_EXCESS.md).
5. [Supplied prime-square lift and partial combined scan](../pi-gcd-prime-adic-2026-11-02/README.md), already preserved before this checkpoint. No duplicate copy is made here. The newly proved saturation theorem strengthens its previously conjectural equality for the 2-part of the gcd.

These files exist at the stated parent commit. Relative links also remain usable in a full repository checkout. The analytic asymptotics and q_n→∞ are previous proof inputs; no numerical extrapolation replaces them here.

## Supplied artifact integrity

Source archive: `pi_gcd_prime_adic_20261102.zip`.

SHA-256: `c3e64443dbd2d77fa8eca77541097ab210f9dfb3eaf271ed2e0087c808082445`.

All 20 entries of its original `GCD_RESEARCH_MANIFEST.json` matched both byte counts and SHA-256. The arithmetic audit can independently repeat this with its optional `--source` argument. The original package is not needed to run the new mathematical verifiers. Its n≤50000 partial-prime scan was read and scoped, not rerun. In particular, absence of pairs with at least one tested prime p≤500 does not exclude pairs with both primes above 500.

## External framework and notation

Nadav Ben David, Guy Nimri, Uri Mendlovic, Yahel Manor, Carlos De la Cruz Mengual, and Ido Kaminer, *On the Connection Between Irrationality Measures and Polynomial Continued Fractions*, Arnold Mathematical Journal 10, 529–566 (2024), DOI [10.1007/s40598-024-00250-z](https://doi.org/10.1007/s40598-024-00250-z). The publisher full text was checked, especially Theorem 2 and Appendix D.

Their error convention uses exponent 1+δ, so our c=1+δ. Specializing their gcd-growth framework to leading coefficients A=2, B=1 gives the same exponential map c=2a/(a−τ), with a=ln(1+√2) and τ the growth rate of g_n/n!. This is an existing general framework. Factorial reduction allows negative τ; it does not ensure c=2. The exact saturation, critical prefactor, signed odd-prime budget, finite integer count certificate, and simultaneous-large-prime consequence are proved directly in this checkpoint rather than asserted from that framework.

The weak guaranteed exponent additionally uses the classical prime number theorem. The official [NIST DLMF §27.12](https://dlmf.nist.gov/27.12), equations 27.12.4–27.12.5, states the relevant prime-count asymptotic. Partial summation yields θ(x)~x; prime powers beyond the first contribute O(√x (ln x)²)=o(x), so ln lcm(1,…,n)=ψ(n)=n+o(n). The rest of the weak guarantee is the explicit denominator bound proved in the exponent note.

## No global approximation conclusion

All exponents and exclusions here concern the specified reduced sequence, including its subsequences. A bad-approximability conclusion requires all reduced rationals; non-bad-approximability requires a constructed unbounded sequence with q² error tending to zero. Neither follows merely from a fixed C/q² inequality, a finite scan, or a conditional growth-rate identity. No new record rational approximation or optimal family exponent is claimed.
