# An exact early stopping search for the fixed arctangent recurrence

2026-11-02. This is an implementation and refinement of the existing
integer certificates, not a claimed new infinite-index theorem.

Let the integer sequences A,B and g_n be those of repository checkpoint
05e6e140c6a41a5ed62f38d6d583487d0455bc7c. Put
U_n=2^n A_n/n!, C_n=odd(U_n), d_n=v_2(U_n)-n.
The numerator U has only O(n) bits, whereas A_n has O(n log n) bits.
It can be streamed by an exact integer recurrence:

\[
(n+1)U_{n+1}=2(2n+1)U_n+4nU_{n-1},\qquad U_0=1,U_1=2.
\]

No full B_n is needed by the sieve. Let

\[
L_n=\operatorname{lcm}(1,\ldots,n-1),\quad
G_0=\operatorname{odd}\bigl(L_n/\gcd(L_n,n)\bigr),\quad
G_{\rm old}=\gcd(C_n,G_0).
\]

The known normalized Wronskian cap proves E_n divides G_0. Positive
odd-prime surplus also divides C_n, because it is bounded by the
valuation of the normalized A_n. Therefore E_n divides G_old.

The actual-initial-value identity improves this cap. For an odd prime
put k_p=floor(log_p(n-1))-v_p(n) and m_p=v_p(C_n). Positive excess at p
is impossible if m_p>k_p. Thus replace the old cap by

\[
G_n=\prod_{\substack{p\text{ odd}\\1\le m_p\le k_p}}p^{m_p},
\qquad E_n\mid G_n\mid G_{\rm old}.
\]

No factorization of C_n is required: factor only G_old, and when its
p-valuation equals k_p test C_n modulo p^(k_p+1). A zero deletes that
entire p^k_p factor. The proof uses the normalized Wronskian
n(ell_n w_(n-1)-ell_(n-1) w_n)=(-1)^n and the denominator bound for
w_(n-1); its complete statement is in the separately audited initial
identity note.

For an odd prime p<=n write n=Np^e+r with p^e<=n<p^(e+1), 1<=N<p.
If p does not divide C_n and does not divide B_N, the established
digit residue formula gives exactly

\[
v_p(g_n)-v_p(n!)=-e.
\]

For any selected subset of such certified primes, put
\(\widehat D=\prod p^e\). It divides the actual odd deficit D_n.
Since v_2(g_n/n!)<=d_n (equality except at n=3), one has rigorously

\[
\frac{g_n}{n!}\le2^{d_n}G/\widehat D,
\qquad
Q_n>\frac{2}{2n+1}2^{-2d_n}\frac{\widehat D^2}{G^2}.
\]

The last strict inequality is the exact positive Legendre error lower
bound, not an asymptotic approximation. Hence every prime may be tested
adaptively and testing may stop as soon as the right side exceeds a
specified rational target. Failure of a certificate does not label the
index good; it merely leaves it for more detailed exact evaluation.

For p>sqrt(n), N<=sqrt(n); for p<=sqrt(n), N<p<=sqrt(n).
Thus testing every odd prime still needs only B_1,...,B_floor(sqrt(n)),
and early stopping usually requires a much shorter prefix. Higher
prime-power generic deficits are retained, unlike a first-level-only
product of certified primes.

An alternative before computing C_n is to use G_0 directly. For a
certified prime whose positive cap is u_p=max(0,floor(log_p(n-1))-v_p(n)),
replace the provisional factor p^u_p by the actual p^-e. This divides
the upper gcd budget by p^(u_p+e). In the top-half range with p not
dividing n, that factor is p^2. This alternative gives a weaker but
still exact partial-state certificate.

## Implementation and validation

`adaptive_exact_sieve.py` builds U_n and the lcm, computes both caps,
then adds deficits in descending-prime order only when needed. Every
accept/reject decision compares rational integers. The optional --audit
mode independently builds the original full A_n,B_n and checks each
upper gcd budget against its exact value. This audit is separate from
the main sieve and does not change what input the sieve itself needs.

Examples: at n=4 the prime3 raises the exact lower quality certificate
to8. At n=7, primes7 and5 give the lower certificate245/24>1. At n=21,
the initial G=17 bound131072/12427 already exceeds1, so no prime tests
are needed. Such early stopping is deliberately aimed at rejecting a
candidate, not maximizing an already sufficient lower bound.

Run:

```bash
python adaptive_exact_sieve.py --N 10000 --audit --output audit_10000.json
python adaptive_exact_sieve.py --N 50000 --output sieve_50000.json
python adaptive_exact_sieve.py --N 50000 --gmp --output sieve_50000_gmp.json
```

The target defaults to 1; `--target 2718/1000` or any other positive
rational changes the rejection threshold. The JSON records exact finite
coverage, survivors, operation counts and trace certificates. Surviving
indices are not asserted successful; the sieve supplies no infinite
exclusion theorem.

The full independent audit through 10000 passed every integer check.
It certified Q_n>1 for every 4<=n<=10000 and left only indices 1,2,3.
The old and improved initial caps each rejected 9991 indices without
any deficit tests. The initial-value improvement nevertheless removed
121 prime powers at 120 indices; these were real improvements to the
upper gcd budget, although they did not change this threshold outcome.
At most 2 deficit-prime tests per index were needed; the only B value
requested by the sieve was B_1=1. The audit additionally constructs
the full B_n solely to independently check the certificate.

The main exact sieve, with no full B_n and only Python stdlib, completed
through 50000 in 150.15 seconds. It rejects all 49997 indices from 4
through 50000. The initial cap alone rejected 49991 indices; only 8
deficit-prime tests were performed over the entire run. The improved
initial-value cap deleted 300 complete prime-power factors at 298
indices. At the endpoint the exact certificate is

\[
Q_{50000}>\frac{2^{49989}}{268440656479721}>1.
\]

These are finite exact certificates. Neither their density nor the
observed rarity of large G proves that every larger index is rejected.
`--gmp` optionally accelerates gcd using the system GMP shared library
through stdlib ctypes. It does not change integer arithmetic or the
certificate. It was separately checked against Python's exact gcd on
506 pairs: zero operands, signed operands, large powers and 500
reproducibly generated large-integer pairs. The complete mathematical
output of the two backends was also compared through 10000. Both passed.
The wrapper supports conventional 64-bit Linux with 64-bit GMP limbs;
these ABI assumptions are checked at runtime. Omit it when GMP is
unavailable; the default version has no external dependency.

The optional GMP extension completed through 100000 in 246.68 seconds,
again leaving only 1,2,3. It removed 394 whole prime powers at 392
indices; the total deficit-prime test count remained 8. Its endpoint
certificate is Q_100000>2^99989/5339604497889. The stdlib 50000 run is
the main finite record; the GMP run is an additional exact experiment.

Files carrying the completed evidence are
`adaptive_sieve_full_audit_10000.json`, `adaptive_sieve_50000.json`,
`adaptive_sieve_100000_gmp.json`, and `gmp_backend_checks.json`.
Run `python verify_gmp_backend.py` to repeat the backend comparison.
