# Independent review of the supplied prime-adic package and prime-square lift

2026-11-02. The original supplied package declares baseline commit
`dff048f73da83e63978e1c01576d38ea71af1546`. This is original-input provenance,
not a statement of the current checkpoint's parent commit.

**Verdict:** the new top-half prime-square classification in
`gcd_blocks/PRIME_BLOCKS_AND_LIFTING.md` is correct. Its local hypotheses, transfer indices, and upper valuation bound are sufficient. It provides an exact extra-factor mechanism, not an aggregate gcd-growth theorem.

The package manifest was checked independently: all 20 listed files have the stated SHA-256 and byte length. This validates the supplied package's internal integrity; it does not by itself verify any remote repository state.

## 1. The prime-square classification

For p odd and n=p+r with 0<r<p, the singular coefficient p^2 in the recurrence at index p vanishes modulo p^2. Thus

\[
X_{p+1}\equiv(2p+1)X_p\pmod{p^2}.
\]

All subsequent steps before n=2p propagate this one-dimensional initial state. The shifted recurrence defining C_r(u) gives exactly

\[
X_{p+r}\equiv C_r(p)X_p\pmod{p^2}.
\]

This is valid for both A and B, with no assumption about X_(p-1). Since the coefficients of C_r are integers,

\[
C_r(p)\equiv C_r(0)+pC_r'(0)=A_r+pD_r\pmod{p^2}.
\]

Differentiation of the shifted recurrence gives the displayed recurrence for D_r, including the two inhomogeneous terms 2A_r+2rA_(r-1). The initial derivative D_1=2 is correct.

The earlier leading-pole identity gives B_p nonzero modulo p, while A_p is divisible by p. Hence the valuation of C_r(p) distinguishes whether the common divisor has exponent zero, one, or at least two. The last possibility is exactly two because the Wronskian at n has p valuation

\[
2v_p((n-1)!)=2,
\]

and p^k dividing both A_n and B_n forces p^k to divide that Wronskian. Here p<n<2p is essential. Consequently

\[
v_p(g_{p+r})=
\begin{cases}
0,&p\nmid A_r,\\
1,&p\mid A_r,\quad A_r/p+D_r\not\equiv0\pmod p,\\
2,&p\mid A_r,\quad A_r/p+D_r\equiv0\pmod p.
\end{cases}
\]

The examples (p,r)=(73,34) and (131,35) fall in the last branch; (17,4) falls in the middle branch. Since v_p(n!)=1 in this range, the last branch is a genuine positive excess beyond the factorial baseline.

## 2. Supporting identities and their limits

The constant-term Lucas proof is valid because the final factor of degree r<p has exponent support [-r,r]. Only its zero exponent can cancel an exponent divisible by p. The reflected-coefficient identity follows from reversing the degree p-1 polynomial D(z)^((p-1)/2); its sign uses (-1/p), as stated.

The leading-pole residue only determines the exact valuation when that residue is nonzero. If it vanishes, the package correctly claims a one-step improvement over the generic lower bound and requires further lifting. It does not identify every vanishing residue with positive surplus.

The transfer block beginning at index zero must be distinguished from the initial B data. The package does so: its nilpotent block is applied starting at index p, followed by the singular step at 2p, and does not assume B_1=B_0. The matrix product order and the conclusions about n>=2p are consistent with that choice.

No issue was found with the criterion's factor 4: it is a unit at the odd primes under discussion.

## 3. New theorem removing the remaining 2-adic conjecture

The supplied discussion leaves open whether the 2-adic part of g_n always reaches that of A_n. The accompanying
`TWO_ADIC_GCD_SATURATION.md` proves

\[
v_2(g_n)=v_2(A_n)\quad(n\ne3),\qquad
v_2(g_3)=v_2(A_3)-1.
\]

The proof needs only the already established A valuation theorem, the integer Wronskian, and r_4=160/51. It proves that all adjacent rational differences from index five onward have nonnegative 2-adic valuation, so every later reduced denominator is odd. An exact B valuation formula is unnecessary.

This removes the extra-loss variable from the exponent budget. It does not remove the need to control both positive and negative odd-prime contributions.

## 4. Scope of the finite data

The package explicitly distinguishes its full gcd verification through n=5000 from its separate partial-prime scan through n=50000. The latter covers p<=500 and then computes full gcds at its 115 candidates. Its exclusion concerns combined surplus with at least one tested prime; two simultaneous surplus primes both above 500 are outside that scan. The summary states this limitation correctly.

This review does not rerun that large partial scan. The new local lifting theorem and manifest are checked directly, and the companion verifier independently checks the exact prime-square classification in a bounded prime range. Finite scans, including the absence of combined gains in their stated coverage, do not establish a uniform asymptotic bound for the net odd-prime contribution.

The independent verifier passed all 20 manifest entries, 5,001 actual gcd-saturation checks through n=5000, 5,000 adjacent-difference valuation checks, 12,540 shifted transfer equalities modulo p^2, and 4,180 exact lift classifications for odd primes p<=199. The only negative adjacent-difference valuations are at indices three and four, as used in the all-index proof.

The verifier is repository-portable: its default invocation runs every mathematical check without requiring the original package directory. It reports the optional manifest check as not requested. To repeat the original-input integrity audit, pass `--source` with the exact directory containing `GCD_RESEARCH_MANIFEST.json`; use `--output` to preserve that provenance result separately. The checkpoint retains the original 20-file integrity result in `original_input_prime_lift_verification.json`, while `prime_lift_and_saturation_verification.json` records the portable default run.

The step from local lifting to an estimate of the form C/q^c still requires a global reduced-denominator budget. The new 2-adic saturation theorem makes that budget cleaner, but supplies no such global odd-prime estimate by itself.
