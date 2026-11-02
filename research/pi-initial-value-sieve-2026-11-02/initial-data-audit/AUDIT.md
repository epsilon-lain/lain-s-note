# Independent audit: first-coordinate Wronskian filter

2026-11-02. Scope: the fixed recurrence and initial values at checkpoint
`05e6e140c6a41a5ed62f38d6d583487d0455bc7c`.

The original differential identity is correct. At the evaluation point,
however, it reduces exactly to the previously established discrete
Wronskian. It should not be advertised as an independent new identity.

The direct Wronskian proof gives the stronger threshold
\(k=\lfloor\log_p(n-1)\rfloor-v_p(n)\), rather than
\(\lfloor\log_p n\rfloor\). If
\(m=v_p(a_n)>k\), the first term of
\(n(a_nb_{n-1}-a_{n-1}b_n)=(-1)^n\)
has positive valuation. Its second term is therefore a unit, yielding
\[
v_p(b_n)=-v_p(n)-v_p(a_{n-1})\le0.
\]
This is valid for every odd prime and every \(n\ge2\), including
\(p\mid n\) and the negative-threshold case \(n=p^r\). All quantities
are rational; no valuation choice in an extension field is needed.

Complete removal of the affected prime power from the positive-surplus
multiple is legitimate. A positive surplus requires \(1\le m\le k\),
so the new multiple includes \(p^m\) precisely on that range. Other
primes contribute no positive factor. The exact negative surplus can
also certify a deficit; its prime power must not be double-counted with
another certificate for the same prime.

The exact rational-polynomial audit covers degrees 0 through 25. A
separate integer recurrence audit through index 1000 checks the
determinant, all odd-prime implications including \(p=n\), the exact
forced surplus, and divisibility of the gain-only multiple by the true
positive-surplus integer. It passes 225 forced nonpositive cases, 208
strict deficits, 184 prime-power boundary cases, and 32 strict
gain-multiple improvements. These are finite validation counts.

The height elimination is also valid:
\(2W-W_1\le\ln C_n\), combined with
\(W+W_1\ge(1+\tfrac12\ln2)n-o(n)\), implies
\[
W_1\ge\left(
\frac{2+\ln2-\ln(2+\sqrt2)}3-o(1)\right)n.
\]
The coefficient is approximately 0.4884000010868099. This only imposes
necessary simple-prime mass; it neither proves that such primes actually
produce positive surplus nor excludes all successful subsequences.
