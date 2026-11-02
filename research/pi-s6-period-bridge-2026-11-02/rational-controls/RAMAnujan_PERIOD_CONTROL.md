# Elliptic-period output: rational Ramanujan construction and denominator-aware modification

2026-11-02. This is a concrete trial of the elliptic-period machinery adjacent to the six-sphere construction. It does not require the global manifold to be a six-sphere, and it does not claim a new identity or a new approximation record.

## A period identity which really eliminates the extra constants

The external analytic input is the classical rational identity

\[
\frac{16}{\pi}=\sum_{n\ge0}\frac{(42n+5)\binom{2n}{n}^{3}}{4096^n}.
\]

A primary proof reference is Jesús Guillera, *A method for proving Ramanujan series for 1/π*, arXiv:1807.07394v4, formula (2). Sections 2–3 explain the mechanism: a modular transformation, Clausen's product identity and Legendre's relation eliminate elliptic periods and leave 1/π. This is an existing positive bridge, rather than an inference that integral monodromy alone controls rational approximations.

Primary URLs:
- https://arxiv.org/pdf/1807.07394
- https://www.davidhbailey.com/dhbpapers/pi-formulas.pdf, formula (36) (author's catalogue)

The output formula contains only rational coefficients. There is no residual square root, elliptic value, logarithm or connection constant to be silently ignored.

## A first integer family and its complete cost

Put
\[
a_n=\frac{(42n+5)\binom{2n}{n}^{3}}{4096^n},\quad
S_N=\sum_{n=0}^N a_n,\quad
C_N=4096^N S_N\in\mathbb Z.
\]

Thus C_0=5 and
\[
C_N=4096 C_{N-1}+(42N+5)\binom{2N}{N}^{3}.
\]

The rational approximation is r_N=16/S_N=16·4096^N/C_N, reduced completely. The exact gcd is
\[
\boxed{\gcd(16\cdot4096^N,C_N)=2^{3s_2(N)},}
\]
where s_2(N) is the binary digit sum. Hence
\[
\boxed{p_N=2^{12N+4-3s_2(N)},\qquad q_N=C_N/2^{3s_2(N)}.}
\]

Proof: v_2 binom(2n,n)=s_2(n), and 42n+5 is odd. In the defining sum for C_N, the final term has valuation 3s_2(N). Every previous j<N has valuation
12(N−j)+3s_2(j)>3s_2(N), because s_2(N)−s_2(j)≤N−j. The minimal valuation is unique. Since the un-reduced numerator is a power of two, this accounts for every common factor, including the absence of any odd cancellation.

This arithmetic calculation is useful for designing the next approximation: fast digit convergence of the period series is not being mistaken for a small reduced denominator.

## Fully rational error certificates

For all n≥0,
\[
\frac{a_{n+1}}{a_n}
=\frac{(42n+47)(2n+1)^3}{512(42n+5)(n+1)^3}<\frac1{32}.
\]

After clearing denominators, this follows from positivity of
\[
336n^4+1216n^3+1440n^2+588n+33.
\]

Therefore, with T_N=16/π−S_N,
\[
a_{N+1}<T_N<\frac{32}{31}a_{N+1}.
\]

The endpoint equality is excluded because all later terms are positive and the ratio estimate is strict. If T lies between those two rational bounds, the exact approximation error is
\[
r_N-\pi=\frac{16T}{S_N(S_N+T)}.
\]

The right side is increasing in T. Multiplication by q_N² gives a rational interval for the normalized quality. No decimal digits of π are used.

## A constructive modification: direct period Padé before inversion

Let
\[
H(t)=\sum_{k\ge0}(42k+5)\binom{2k}{k}^3t^k,
\qquad t_0=1/4096,\quad H(t_0)=16/\pi.
\]

Choose Q_n(0)=1 with degrees P_n,Q_n≤n and impose
\[
Q_n(t)H(t)-P_n(t)=O(t^{2n+1}).
\]

This exact rational linear system uses coefficients of H only; it does not use π or a target p/q. At t_0 define
\[
r_n=16 Q_n(t_0)/P_n(t_0),
\]
provided P_n(t_0),Q_n(t_0) are nonzero, and reduce completely. These nonzero and nonsingularity properties are checked in the stated finite range; they are not asserted for every order without proof.

If R_n=Q_nH−P_n, then
\[
r_n-\pi=\frac{\pi R_n(t_0)}{P_n(t_0)}.
\]

The exact reduced p,q and rational period-tail interval are used together in the verifier. Thus the arithmetic cost of the Padé coefficients is retained.

### A genuinely improved finite design at equal information budget

Both constructions below use exactly H_0,H_1,H_2 (three series coefficients).

- Direct truncation gives 33554432/10680675, with normalized quality between approximately 1.0731·10^9 and 1.1078·10^9.
- The [1/1] period Padé construction has
\[
Q_1=1-\frac{2403}{47}t,\qquad
P_1=5+\frac{5657}{47}t,
\]
and yields
\[
\boxed{r_1=3041744/968217.}
\]
Its normalized quality is about 567103.3362. Its reduced denominator is smaller by the exact factor 10680675/968217≈11.03, and its quality is rigorously improved by a factor greater than 1892, using the conservative rational tail certificates.

This is an effective finite modification, not a good Diophantine approximation: its quality is still much larger than one. Order [2/2] also improves the equal-coefficient truncation in this finite comparison. From order [3/3] onward in this small scan, coefficient growth erases the benefit; the code reports those cases too, rather than extrapolating a monotone improvement.

## Completed exact checks

`ramanujan_exact.py` completes truncations N=0,...,50, verifies the exact gcd formula at every such N, checks the positive ratio polynomial, and constructs all diagonal Padé orders 1,...,10 with every formal cancellation equation checked exactly. It uses a 220-term rational period-tail interval to certify each Padé quality.

`ramanujan_checks.json` contains the full finite data, same-coefficient comparisons and reduced p,q. No assertion about an infinite quality-zero sequence follows. The constructive lesson is specific: the period identity does remove the unwanted constants, and resummation can reduce denominator cost at low order; extending that gain requires a structured period recurrence/denominator design rather than a generic high-order Hankel solve.

Run from this directory:

```bash
python ramanujan_exact.py
```

SymPy supplies exact rational linear algebra. All error bounds use Fraction arithmetic and the cited identity, not floating-point π. Floating-point fields in the JSON are descriptive displays only; rational bounds are authoritative.
