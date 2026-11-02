# Independent audit of the odd-prime and 2-adic gcd attack

2026-11-02. Reviewed:

- `../gcd-prime-blocks/ODD_PRIME_BLOCKS_AND_EXCESS.md`
- `../gcd-valuation-theory/two-adic-proof/TWO_ADIC_NUMERATOR_VALUATION.md`
- `../gcd-valuation-theory/ODD_EXCESS_BARRIER.md`

**Verdict:** the displayed general identities and inequalities are sound. In particular, the 2-adic formula is an all-index theorem with a complete proof, not an empirical recurrence pattern. The result forces exponentially large positive odd-prime excess for good quality in this fixed family, but no theorem excluding all such excess has been obtained.

## 1. Odd-prime valuation bounds

The generating-function identity G=F integral(F), with F=(1-2t-t^2)^(-1/2), gives the claimed convolution for w_n. Since ell_n is in Z[1/2], this proves the lower valuation bound for w_n at every odd prime.

For e_p(n)=min(v_p(ell_n),v_p(w_n)), the backward Wronskian has both terms of valuation at least e_p(n)-floor(log_p(n-1)); its nonzero right-hand side has valuation -v_p(n). This proves the stated bound. The forward Wronskian gives the second bound in exactly the same manner. The odd-part divisibility is equivalent to the backward inequality prime by prime; negative e_p values cause no problem because it compares integer valuations after multiplying by n and n!.

The stronger positive-excess assertion was derived independently during this audit and added to the companion note: if e_p(n)>0, then

\[
v_p(w_{n-1})=-v_p(n)-v_p(\ell_n),\qquad
v_p(w_{n+1})=-v_p(n+1)-v_p(\ell_n).
\]

The other term in either Wronskian has strictly positive valuation, so these exact equalities follow without a cancellation hypothesis. They bound the full valuation of ell_n whenever excess occurs.

The Lucas argument is self-contained: in odd characteristic the two square roots in
F(t)=D(t)^((p-1)/2)F(t^p), D=1-2t-t^2, have the same square and constant coefficient. Coefficient comparison proves the digit factorization. In the leading-pole formula, after multiplication by p^h, precisely the summands whose integration denominator is a multiple of p^h survive modulo p. Their digit factors give the claimed chi_p^h*ell_r*w_b. This establishes a leading residue, not higher p-adic valuations when that residue vanishes.

The transfer order C_s*S_m^k is correct for n-1=km+s because the coefficient matrices repeat with period m modulo m. The prime block sends the A initial vector to zero, fixing the column sum; the leading-pole formula and Wilson pairing fix the B image. This gives the stated rank-one nilpotent block. None of these arguments assume invertibility of the block modulo p.

The first odd-excess counterexample n=72,p=17 is reproduced independently. Thus replacing the excess bound by e_p(n)<=0 would be an actual error.

## 2. Exact 2-adic valuation theorem

The alternative coefficient formula

\[
\ell_n=\sum_k\frac{n!}{2^k(k!)^2(n-2k)!}
\]

follows by expanding (1-2z-z^2)^(-1/2) as
(1-z)^(-1)(1-2z^2/(1-z)^2)^(-1/2). Reversing this finite sum gives the author's even and odd prefactors exactly, including their powers of two.

For the even sum, the central-binomial valuation gives the tail estimates. When m is odd, R_2+R_3 is divisible by sixteen, the terms of index at least four are divisible by eight, and 1+m^2 is two modulo four. When m is even, the sum is odd. This proves the even-index formula with no unexamined cancellation.

The subtle odd-index step also checks. For c_j=(j!)^2/(2j+1)!, the coefficient ratio gives

\[
z(4-z)F'(z)+(2-z)F(z)=2.
\]

Writing H(z)=(2-z)F(z(4-z))-2F(z), direct differentiation gives the same homogeneous equation for H and H(0)=0. Over the rational formal power series this forces H=0. In evaluating at z=2 in the 2-adic sense, the expanded terms of c_j z^j(4-z)^j have valuations at least 2j-s_2(j). This tends to infinity uniformly in the finitely many terms at each j, justifying their regrouping. Thus F(2)=0 is valid in Q_2; the proof does not identify it with the real sum of the same rational terms.

After subtracting that zero series, the tail terms in T(m)/(m+1) have valuation at least v_2(j!)-floor(log_2 j), which is at least two for every j>=6. The finite H_5 polynomial has an odd denominator and yields H_5(2v)=2 modulo four. Hence T(m)/(m+1) is exactly two modulo four for odd m. The resulting valuations translate correctly to

\[
v_2(A_{4m})=2m,\quad v_2(A_{4m+1})=2m,
\]
\[
v_2(A_{4m+2})=2m+2,\quad
v_2(A_{4m+3})=2m+3+v_2(m+1).
\]

This supplies a proof for all indices, including the low-index cases.

## 3. The odd-excess barrier and the precise remaining gap

The inequalities v_2(g_n)<=v_2(A_n), v_2(n!)=n-s_2(n), and s_2(n)<=floor(log_2 n)+1 give the author's deliberately loose bound

\[
g_n/n!\le8n^2\,2^{-n/2}u_n,
\]

where u_n is the positive odd numerator of the reduced rational g_n/n!. Both versions of the quality bound check:

\[
Q_n\ge\frac{2^n}{64n^5u_n^2}
\quad\text{for all sufficiently large n},
\]
\[
Q_n>\frac{2^n}{32(2n+1)n^4u_n^2}
\ge\frac{2^n}{96n^5u_n^2}
\quad(n\ge1).
\]

The first uses the established asymptotic analytic multiplier tending to sqrt(2); the second uses the exact Legendre error floor and has no unspecified threshold. A subcritical exponential upper bound for u_n would close this particular family. Such an upper bound has not been proved, and the finite small-excess observations do not supply it.

## 4. Independent sieve and finite validation

The root investigation's U_n-only sieve was checked independently. The companion `LOCAL_POLES_AND_EXCESS_SIEVE.md` also proves a stronger whole-prime removal rule: if v_p(ell_n) exceeds either admissible neighboring cap, p contributes no positive excess. Gcd saturation implements this rule without factorization.

The independent verifier `verify_local_poles_sieve.py` passed 48 exact polynomial differential-identity checks and compared the refined bound against actual gcds at all 999 indices from n=2 through n=1000. Its lower bound alone certifies Q_n>1 at 991 of those indices. The survivors are exactly 2,3,4,5,6,7,11,15; survival means only that this sufficient bound did not exclude the index. Whole-prime removal strengthens the bound at 32 indices in this range.

The two-adic author verifier, the odd-excess counterexample verifier, and the odd-prime block verifier were rerun and passed. The last includes 7,000 digit checks, 7,000 leading-pole checks, 216 prime-power block checks and all eleven positive-excess neighbor identities in its range. Finite checks audit algebra, constants and counterexamples; the general proofs above are what justify the all-index statements.

No result here proves or disproves pi's bad approximability, and no all-index exclusion of quality-zero subsequences in this one family is claimed.
