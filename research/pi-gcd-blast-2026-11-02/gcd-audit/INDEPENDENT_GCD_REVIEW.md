# Independent audit of the gcd attack

2026-11-02. Reviewed the odd-prime block note, the direct scan and U-only sieve, `SIEVE_THEOREM_AND_FINITE_EXCLUSION.md`, their stored results, the all-index 2-adic numerator proof, and the odd-excess barrier including its final uniform bound. Verdict: **the claims within their stated scope pass; no unresolved mathematical blocker was found.** The full `N=20000` computation was code-reviewed and its reported output inspected; it was not redundantly rerun. A separate small exact implementation checks the mechanisms below.

The all-index statements in this review follow from proofs. The finite exclusion through 20000 concerns one specified recurrence family. Neither establishes a lower bound for every rational approximation to pi.

## 1. Generating functions, Wronskians, and odd-prime bounds

Let `ell_n=A_n/n!`, `w_n=B_n/n!`, and `g_n=gcd(A_n,4B_n)`. The dyadic binomial expression for `ell_n` proves its integrality at each odd prime. The recurrence and initial values independently give

\[
F(t)=(1-2t-t^2)^{-1/2},\qquad
G(t)=F(t)\int_0^tF(s)\,ds.
\]

Thus `w_n=sum_(j=1)^n ell_(j-1) ell_(n-j)/j`, and its p-adic valuation is at least `-floor(log_p n)`. This is a rational denominator bound; no analytic convergence or numerical logarithms are needed.

The signed Wronskian is

\[
\ell_nw_{n-1}-\ell_{n-1}w_n=(-1)^n/n.
\]

At an odd prime, `e_p(n)=v_p(g_n)-v_p(n!)=min(v_p(ell_n),v_p(w_n))` may be negative. Applying the denominator bound to the neighboring `w` values yields both caps

\[
e_p(n)\le\min\{\lfloor\log_p(n-1)\rfloor-v_p(n),
\lfloor\log_p(n+1)\rfloor-v_p(n+1)\},\quad n\ge2.
\]

Their signs and endpoints are correct. The corresponding odd-part divisibility `odd(n g_n) | odd(n! lcm(1,...,n-1))` is valid prime by prime. It does not say that `odd(g_n)` divides `odd(n!)`.

The stronger necessary condition for positive excess is also sound. If `e_p(n)>0` and `u_p=v_p(ell_n)`, one term in each neighboring Wronskian has strictly positive valuation, whereas the Wronskian itself has nonpositive valuation. The other term therefore has exactly that latter valuation. This proves

\[
v_p(w_{n-1})=-v_p(n)-u_p,
\qquad v_p(w_{n+1})=-v_p(n+1)-u_p.
\]

Consequently `u_p` itself must satisfy both caps. If it exceeds either cap, that prime cannot contribute positive excess. The current sieve uses the weaker backward cap; this is conservative, not erroneous.

## 2. Lucas congruence, leading poles, and transfer blocks

In odd characteristic, the identity

\[
F(t)=D(t)^{(p-1)/2}F(t^p),\qquad D(t)=1-2t-t^2,
\]

follows by squaring: both sides square to `1/D(t)` and have constant term one. Since `2` is invertible, this square root is unique. The polynomial prefactor has degree `p-1`, which proves the stated Lucas digit factorization without an unproved interchange of infinite sums.

For `n=b p^h+r`, multiplying the exact convolution for `w_n` by `p^h` removes all terms except `j=k p^h` modulo `p`. The digits of `k p^h-1` contribute `chi_p^h ell_(k-1)`, giving

\[
p^h w_n\equiv\chi_p^h\ell_r w_b\pmod p.
\]

If this residue is nonzero, the valuation is exactly `-h`, and so is `e_p(n)`. This conclusion remains correct at `h=0`. If the residue vanishes, it supplies no exact lifted valuation. In particular, digit zeros are not additive valuation contributions.

For the transfer matrix, the product order is correct:

\[
v_n\equiv C_s S_m^k v_1\pmod m,
\qquad n-1=km+s.
\]

The prime block sends the A initial vector to zero modulo `p`. Its first column is determined by the leading-pole formula for B and Wilson's identity. This gives precisely

\[
S_p=-\chi_p\begin{pmatrix}1&-1\\1&-1\end{pmatrix}\pmod p,
\qquad S_p^2=0\pmod p.
\]

No invertibility is assumed. Reduction modulo `p` alone does not determine the blocks or valuations modulo higher powers.

## 3. The counterexample is exact

Independent recurrence reconstruction gives

\[
v_{17}(A_{72})=v_{17}(B_{72})=5,
\qquad v_{17}(72!)=4.
\]

After removing `17^5`, the first nonzero residues of A and B are both `2 mod 17`. Thus `e_17(72)=1` is an exact counterexample to the proposed no-odd-excess divisibility. Checking all earlier indices confirms that 72 is the first odd-excess index through this range.

Independent modular recurrence checks also reproduce the base-17 repeated-digit examples at `n=4,72,1228,20880`, with `(v_17 ell_n,v_17 w_n,e_17(n))` respectively `(1,0,0),(1,1,1),(2,-1,-1),(2,0,0)`. These disprove that particular finite extrapolation; they do not establish a periodic or infinite pattern.

## 4. The U-only sieve is a valid necessary condition

Write `U_n=2^n ell_n`, `s_n=v_2(U_n)`, `M=lcm(1,...,n-1)`, and

\[
E_n=M/\gcd(M,n),\qquad
G_n=\gcd(\operatorname{odd}(U_n),\operatorname{odd}(E_n)).
\]

For each odd prime, positive excess is bounded both by `v_p(U_n)` and by the nonnegative part of the backward Wronskian cap. Therefore its positive contribution divides `G_n`. Negative odd valuations can be discarded for an upper bound. Since `g_n` divides `A_n`, its 2-adic gain over `n!` is at most `s_n-n`. Hence

\[
g_n/n!\le 2^{s_n-n}G_n.
\]

The earlier positive Legendre error integral gives the strict bound

\[
Q_n>\frac2{2n+1}(n!/g_n)^2
\ge\frac{2\,2^{2(n-s_n)}}{(2n+1)G_n^2}=:L_n.
\]

Thus `Q_n<=1/5` requires `L_n<1/5`. The code's strict comparison `5*lower_num<lower_den` has the correct boundary behavior: equality must be excluded. Its two branches for negative and nonnegative powers of two implement the rational quantity exactly.

The integer recurrence for U, the update from `M_(n-2)` to `M_(n-1)`, the index `n-1` in the lcm, and the initialization at `n=1` were all checked. The second stage uses the actual gcd and actual reduced fraction; passing the first sieve is not treated as evidence of good approximation.

The new theorem note's coefficient identity `U_n=[x^n](1+2x+2x^2)^n` is correct: choosing `k` quadratic terms, `n-2k` linear terms, and `k` constant terms gives `n! 2^(n-k)/(k!^2(n-2k)!)`, exactly `2^n` times the binomial expression for `ell_n`. This supplies an independent integrality proof.

The reported completed run has `N=20000`, excludes 19995 indices in the first stage, and leaves exactly `2,3,6,7,15`. Their strict actual lower quality bounds are respectively `1/10`, `9/14`, `50/13`, `245/24`, and `4058104050/31`. The upper bound at `n=2` is strictly `1/5`, so the other four candidates are excluded and `n=2` uniquely minimizes quality over the checked family. It represents `3/1`. These logical deductions from the completed output are valid.

The separate direct scan through 5000 uses the same integer recurrence and actual gcd. Its quantities `n(g_n/n!)^2`, quality floors, and reduced odd numerator are computed exactly. Its reported maxima and minima are finite results, not all-index hypotheses. The theorem note says the maxima are “attained at” the recorded indices, which correctly avoids asserting unverified uniqueness of those maxima. In contrast, uniqueness of the best quality at `n=2` follows from the strict quality bounds. The stronger `Q_n>9/14` floor is restricted to the direct-scan range, as required. No numerical pi appears in either scan.

## 5. The exact 2-adic numerator theorem

The full proof in `two-adic-proof/TWO_ADIC_NUMERATOR_VALUATION.md` was checked independently, including its convergence step.

For even indices, reversing the binomial expression gives a factor with valuation `m` and a remaining sum whose first two terms are `1+m^2`. The bound `v_2(R_j)>=v_2(j!)`, and the combined expression for `R_2+R_3`, give exactly the stated parity-dependent extra valuation.

For odd indices, the inverse-binomial series has coefficients `c_j=(j!)^2/(2j+1)!`, of valuation `-s_2(j)`. Its differential equation and formal duplication identity are correct. This series is not the generating function `F(t)` from the odd-prime section; the reused letter is local to the proof.

Substitution of `z=2` into the duplication identity is justified. Both `F(2)` and `F(4)` converge 2-adically. When the composed term `c_j z^j(4-z)^j` is expanded and then evaluated at 2, every summand has valuation at least `2j-s_2(j)`. That bound tends to infinity uniformly over the inner expansion, permitting regrouping. Thus the proof legitimately obtains `F(2)=0` in the 2-adic field. It does not identify this limit with the different real limit.

For `u=m+1` even, the identity for `P_j(u)-1` and the divisibility of `binom(u,r)` give

\[
v_2\left(t_j\frac{P_j(u)^2-1}{u}\right)
\ge v_2(j!)-\lfloor\log_2j\rfloor\ge2\quad(j\ge6).
\]

Consequently only five terms remain modulo four. Independent exact polynomial expansion reproduces all ten coefficients of the displayed `W(v)` and the denominator `155925`, which is odd. Therefore `T(m)/(m+1)` is `2 mod 4`, proving the exact odd-index valuation.

This establishes the four all-index formulas for `v_2(A_n)`, including the exceptional term at `n=4m+3`. Finite recurrence checks are supporting checks only; they are not the reason the result holds for every index.

## 6. The remaining odd-excess barrier and verification scope

The exact 2-adic theorem implies the loose upper bound used in `ODD_EXCESS_BARRIER.md`. If `u_n` is the odd part of the reduced numerator of `g_n/n!`, then

\[
g_n/n!\le8n^2 2^{-n/2}u_n.
\]

Combining with the previously proved joint asymptotic yields the stated eventual inequality

\[
Q_n\ge\frac{2^n}{64n^5u_n^2}.
\]

The qualifier “for all sufficiently large n” matters for this constant `64`, because this step uses an asymptotic lower bound.

The final appended all-index estimate is independently valid and needs no asymptotic threshold. Combining the exact strict quality floor with the same gain bound gives, for every `n>=1`,

\[
Q_n>\frac2{2n+1}\frac{2^n}{64n^4u_n^2}
=\frac{2^n}{32(2n+1)n^4u_n^2}
\ge\frac{2^n}{96n^5u_n^2}.
\]

The last comparison uses only `2n+1<=3n`, including `n=1`. Thus the constant `96` is uniform over all positive indices; `64` belongs to the separate eventual estimate. If `Q_n<=M`, the uniform version forces the strict lower bound `u_n>2^(n/2)/(sqrt(96M)n^(5/2))`. The strictness follows from the first strict quality inequality, even when the final elementary denominator comparison is an equality.

Both versions imply the stated necessary exponential odd-excess growth and conditional exclusion under a subcritical bound. No such global upper bound on `u_n` has been proved here.

The companion `verify_independent_gcd.py` uses only exact standard-library arithmetic. It checks a short full range through 300, 2800 Lucas identities, 2800 leading-pole identities, 2786 two-sided valuation caps, prime and prime-power transfer matrices, the counterexamples, the U-only bound against actual gcds, the formal inverse-binomial identities, and the independent five-term polynomial certificate. Its results are in `verify_independent_gcd.json`. It also checks the repeated-digit examples via modular recurrence up to 20880 without performing the full large scan.

The full 20000-index result remains an audited finite computation, not an all-index theorem. The missing task is simultaneous control of positive odd-prime excess and its lifts. Neither the counterexample nor the 2-adic loss decides that task, and even an exclusion of this entire family would not prove pi badly approximable.
