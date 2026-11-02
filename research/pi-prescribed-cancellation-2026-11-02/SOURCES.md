# Sources and dependency boundary

2026-11-02. The user's active direction is construction, with the old prime-pool route paused. This checkpoint follows that direction.

## Repository inputs

Parent research checkpoint: db1290c827479ff85beaeca8606e8258d72f0c69.

- [Current direction](../CURRENT_DIRECTION.md): authoritative scope and candidate accounting requirements.
- [Constructive redesign](../pi-constructive-redesign-2026-11-02/README.md): fixed family, integer-antiderivative coordinates and candidate data.
- [Integer antiderivative family](../pi-constructive-redesign-2026-11-02/integer-kernel/INTEGER_ANTIDERIVATIVE_FAMILY.md): invertible rational primitive map, determinant, endpoint conditions and period image.
- [Rational integral certificate](../pi-constructive-redesign-2026-11-02/certificates/integral_certificate.py): centered polynomial expansion and rational remainder bound. The new full-truncation construction uses that expansion as a generator, not just as a posterior certificate.
- [Existing fixed family data](../pi-constructive-redesign-2026-11-02/integer-kernel/new_family_results.json): Git blob 3cd2fd50b288c470977f776c9c34e6b85c1dccde. The new comparison reads this exact source and generates no additional searched candidate.

The new identities, base periods, period independence, full moment-space parametrization, Fourier extraction, K=13m rule and joint bounds have proofs in this directory. Computational checks audit their algebra and implementation; they do not substitute for an asymptotic quality theorem.

## Primary references checked

1. NIST DLMF, section18.5(ii), Rodrigues formulas and Table18.5.1. https://dlmf.nist.gov/18.5 . This identifies the classical orthogonal-polynomial form. Our proof of the rational moment-space bijection is direct integration by parts and dimension counting; no statement that Rodrigues cancellation improves rational approximation quality is attributed to DLMF.
2. Doron Zeilberger and Wadim Zudilin, *The irrationality measure of pi is at most 7.103205334137...*, author manuscript revised7January2020. https://sites.math.rutgers.edu/~zeilberg/mamarim/mamarimPDF/pimeas.pdf . Checked as a construction-method reference already named by the user. Its irrationality-exponent bound is not a theorem about Q_m tending to zero and is not used to prove any new approximation claim here.

The rational constant22/7 used in the upper certificate is justified by the elementary positive integral identity integral_0^1 x^4(1-x)^4/(1+x^2) dx=22/7-pi. It does not serve as a target approximant in multiplier generation.

## Unproved objective

The improved norm gives the sufficient condition q_m^2||H_m||_1/|alpha_m|=o(16^m sqrt(m)). The explicit truncation generator has its own signed-integral joint upper bound, but neither condition has been proved to yield Q_m->0. Convergence of the generated rationals to pi and divergence of their reduced denominators do hold; no monotonicity of those denominators is asserted. No claim of a new record, pi bad approximability, or pi non-bad-approximability is made.
