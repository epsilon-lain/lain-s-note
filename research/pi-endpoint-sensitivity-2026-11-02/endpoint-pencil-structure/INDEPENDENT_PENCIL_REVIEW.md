# Independent review of the endpoint pencil

Reviewed 2026-11-02 by a separate research agent. Scope: `ENDPOINT_PENCIL_STRUCTURE.md`, including the finite certificate added as equation (19). Verdict: the stated identities and consequences are correct, conditional on the explicitly cited earlier period and Gram estimates. No unresolved mathematical blocker was found.

## Algebra and signs

Factoring `ET` on the left gives exactly

\[
C(Z)=T(Z)^{-1}E^{-1}T(Z)E,
\qquad G_n(Z,Y)=\det(YI-C(Z)).
\]

The inverse-conjugacy is in the asserted orientation:
`C(-Z)=T(Z) C(Z)^{-1} T(Z)^{-1}`. Because `det C=1`, the resulting sign is `(-1)^n`, and the factor in the reciprocal identity is `Y^n`. The constant coefficient is consequently `(-1)^n`, so there are no zero endpoint roots.

At imaginary `Z`, coefficientwise conjugation changes `Z` to `-Z`. This proves reciprocal-conjugate root pairing, including algebraic multiplicity. An odd degree forces a unit-circle root, but not the particular root `Y=1`. The derivative calculation gives precisely `Re(G_Y/G)=n/2`; it imposes no upper bound on the imaginary part.

The moment-pencil scaling is

\[
G_n(Z,Y)=Z^{n^2}A_n^{-2}\det(YV(Z)-U(Z)).
\]

Both endpoint moment determinants equal `A_n^2/Z^(n^2)` without an additional sign. Thus, for `M=(V-U)^{-1}V`, both `det M` and `det(M-I)` equal `1/P_n(omega)`. In particular the transformations between endpoint roots and eigenvalues never encounter an eigenvalue zero or one. The map `Y=1-1/mu` sends reciprocal-conjugate pairing to `mu -> 1-conjugate(mu)`, with multiplicities intact. Diagonalizability is unnecessary throughout.

## Spectral and counting estimates

The normalization is a similarity of the relative operator:

\[
(SHS)^{-1}(SVS)=S^{-1}MS.
\]

This preserves eigenvalues. The proof does not claim that this similarity preserves singular values. Every singular value of `SHS` is at most one, while their product exceeds `1/6`. Each therefore exceeds `1/6`, yielding inverse norm less than six. The resulting upper spectral-radius bound `rho(M)<6L_n` is valid. The determinant lower bound gives `rho(M)>=|P_n(omega)|^(-1/n)`.

Taking reciprocal eigenvalue moduli establishes the exact nearest-root identity `delta_n=1/rho(M)`, hence equation (17). The large-eigenvalue count uses the product of all moduli and an upper bound on each modulus; its threshold, denominator, asymptotic coefficient `(1-alpha)/(2-alpha)`, and multiplicity convention are correct.

The finite certificate (19) is also sufficient. If fewer than `k` eigenvalues exceed `T_n`, order eigenvalues by modulus: the largest `k-1` are bounded by `6L_n` and the rest by `T_n`. This bounds their complete product by the expression used in the note. The strict rational upper bound `2*pi<44/7` gives the strict contradiction claimed. The certified counts are conservative lower bounds, not numerical root estimates.

## Independent checks and limits

I reran `python verify_endpoint_pencil.py`; it passed, including the certificates `(n,k)=(64,1),(96,3),(128,5)`.

Separately, an independent exact SymPy computation tested sizes `n=2,...,6` at `Z=3/2, 2i, 7i/3`, with endpoint parameters `Y=3/2,i`. Together with the companion exterior audit this made 170 successful exact checks. Relevant checks here included the commutator inverse-conjugacy, raw/moment determinant scaling, reciprocal identity, endpoint translation, both relative determinant identities, and the real trace. No floating-point approximation to pi was used. These checks catch finite algebra and indexing errors; the all-size assertions follow from the proofs, not from finite experiments.

The conclusions locate roots in independent `Y` at fixed `Z=omega`. They do not locate roots of `P_n(Z)=G_n(Z,1)` or on `Y=exp(Z)`. No estimate for rational approximation to pi follows without the additional coupling argument stated in the note. This limitation is necessary and has been retained.
