# Sources, dependencies, and exact claim boundary

2026-11-02. Input repository checkpoint: c4d1ddfc235c7f32ebb901de78349e54c6bd22a9. The supplied S6-period bridge is used through its displayed elliptic differential system and normalized Ramanujan special value. No global topological conclusion about S6 is assumed.

## Earlier repository inputs

- [Current construction direction](../CURRENT_DIRECTION.md).
- [Full-parity primitive kernel](../pi-parity-beta-redesign-2026-11-02/new-kernel-design/PARITY_FREE_ENDPOINT_FAMILY.md): exact endpoint kernel, odd-prime period image, no-log degree condition, 4/121 kernel bound.
- [The supplied bridge](../pi-s6-period-bridge-2026-11-02/README.md) and its [normalized period subsystem](../pi-s6-period-bridge-2026-11-02/periods/MODULAR_SUBSYSTEM.md).

The main Smith period quotient, Hilbert/Schur objective, charged integer nearest-plane lift, forced-index-D sublattice and rational certificate selection are direct new deductions with proofs in period-main-kernel/ and independent audits in bridge-review/. The main generator imports the earlier parity-free kernel through a repository-relative search; that earlier folder must remain present for generation. Saved-output verifiers need only standard-library arithmetic and their saved JSON.

## Primary references read

1. Jesus Guillera, *A method for proving Ramanujan series for 1/pi*, arXiv1807.07394v4 (2018), Formula2 and Legendre relation Formula6. https://arxiv.org/pdf/1807.07394 . This supplies the classical special value, not a denominator theorem for new approximants.
2. NIST DLMF16.12.2, Clausen's formula. https://dlmf.nist.gov/16.12 . Substitution a=b=1/4 gives the smaller second-order period whose square is the three-half-parameter series. The integer square coefficient normalization is proved directly in the new note.
3. Bing He, *On the divisibility properties concerning sums of binomial coefficients*, The Ramanujan Journal43 (2017),313-326; published online29March2016. DOI10.1007/s11139-016-9780-6. https://link.springer.com/article/10.1007/s11139-016-9780-6 . The official abstract explicitly states the all-n>1 divisor 2n*binom(2n,n) of the Ramanujan partial-sum integer. The subscription full WZ proof was not read; the theorem is an imported published result, with bounded checks of the specialization. The tail control does not use up this theorem's extra divisor gain; the note states that limitation explicitly.

The bridge's geometry manuscript and S6 topology are not additionally re-audited in this checkpoint; their validation is unnecessary for the standalone hypergeometric identities and the independent integer primitive construction used here.

## Exactness and interpretation

- All construction coordinates, Gram matrices, Schur complements, rounding and gcds are rational/integer.
- Pi/Machin intervals are used only for post-construction comparisons. `select_by_certificate.py` reads only the proved rational quality-squared upper expression, never those intervals.
- For a fixed period quotient coordinate, zero-period lifting changes representation norm and its certificate, but not the true fraction, gcd, integral or quality.
- The new congruence really restricts the two-dimensional period lattice; its index and transformed metric are charged. Primitive integer and multiplier contents are recorded in the certified benchmark.
- The elliptic shared-denominator control has an explicit differential compatibility defect. A common polynomial denominator is not by itself a full preservation of the second-order system or the physical 2pi*i normalization.
- Finite recovery of355/113 is an implementation and certificate benchmark. No new approximation record or universal quadratic-quality decay is asserted.
