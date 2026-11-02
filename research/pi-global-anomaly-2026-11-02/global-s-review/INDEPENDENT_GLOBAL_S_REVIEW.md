# Independent review of the global first-layer reduction

2026-11-02. Reviewed the new budget, reflection, and structural-experiment files independently. **The stated reduction and boundary-removal results pass. The eventual estimate S(n)<=n/6 remains unproved.** This is an independent research-agent audit, not external peer review.

## 1. Exact signed budget

Use the specified integer recurrence, with g_n=gcd(A_n,4B_n). Let E_n/D_n be the odd part of g_n/n!. Let H_n be the product of odd primes n/2<p<n, P_n the product of the exceptional primes p|A_(n-p) in that interval, and L_n the odd part of lcm(1,...,n-1).

The exact divisibility

\[
E_n\operatorname{odd}(n)H_n^2\mid D_nP_n^2L_n
\]

is correct for every n>=2. At lower-half primes, including p=n/2 if present, its valuation inequality is exactly the existing Wronskian cap. In the upper half, v_p(n!)=1 and p does not divide n. The nonexceptional excess is -1 and an exceptional excess is at most one, giving the needed inequality e_p(n)<=-1+2·1_exception. At p>=n the same cap handles the endpoint p=n. No missing endpoint contributes a linear term.

If T_n is the product of upper-half primes for which p²|g_n, the exact local classification gives

\[
e_p(n)=-1+\mathbf1_{p\mid A_{n-p}}+\mathbf1_{p^2\mid g_n}.
\]

Hence the stronger divisibility with P_nT_n in place of P_n² also holds. Its logarithmic budget involves S(n)+ln T_n. Replacing T_n by its divisor upper bound P_n recovers the simpler first-layer target; no first-layer exception is silently counted as an actual positive excess.

Writing J_n=ln(E_n/D_n), the finite logarithmic consequence is

\[
J_n\le2S(n)+\psi_{\rm odd}(n-1)-2\ln H_n-\ln\operatorname{odd}(n).
\]

The classical prime number theorem gives psi_odd(n-1)=n+o(n) and ln H_n=n/2+o(n). The remaining integer log is O(ln n), so J_n<=2S(n)+o(n) follows along all indices, with no density assumption or subsequence exclusion.

## 2. Consequence for the actual approximation quality

Let Q_n=q_n²|pi-p_n/q_n|, where the fraction is fully reduced, and let d_n=v_2(g_n/n!). The known strict all-index error bound gives

\[
Q_n>\frac2{2n+1}\,2^{-2d_n}(D_n/E_n)^2.
\]

Thus the new integer divisibility supplies the rigorous finite lower bound

\[
Q_n>
\frac{2^{1-2d_n}\operatorname{odd}(n)^2H_n^4}
{(2n+1)P_n^4L_n^2}.
\]

Replacing P_n⁴ by P_n²T_n² strengthens it. Binary saturation and the exact numerator valuation give d_n=-n/2+O(ln n). Combining with PNT yields ln Q_n>=n ln2-4S(n)+o(n).

Consequently, the proposed eventual S(n)<=n/6 is indeed sufficient to force Q_n to infinity exponentially. The coefficient ln2-2/3 is positive. More generally, any proved limsup S(n)/n<ln2/4 would suffice. The equality case is not settled by this exponential comparison. The implication does not establish its S-hypothesis.

## 3. Reflection and folded degrees

Coefficient reversal of (1-2z-z²)^((p-1)/2) correctly yields

\[
\ell_{p-1-r}\equiv(-1)^{(p-1)/2+r}\ell_r\pmod p.
\]

Its zero equivalence transfers to A_r only for r<p, where r! is a p-unit; this condition holds throughout the upper-half application. Thus k(n,p)=min(n-p,2p-n-1) is a valid folded degree and satisfies 0<=k<=(n-1)/3. For fixed k, at most the two indicated prime values n-k and (n+1+k)/2 occur. Their coincidence at the transition must be counted once, as the note does.

Let U_k=2^k A_k/k! and C_k=odd(U_k). The existing exact binary theorem gives v_2(U_k)>=k/2. Combined with ell_k<=lambda^k, this proves the exact finite bound C_k<=(2+sqrt(2))^k, without any O(ln k) loss.

Distinct exceptional primes folded to k<=K have product dividing the product of C_k for k<=K. This remains correct when both distinct branch primes occur at one k, since each divides the same integer. Therefore

\[
S_{\rm edge}(n,K)\le\frac{\ln(2+\sqrt2)}2K(K+1).
\]

For K=o(sqrt(n)), the bound is uniformly o(n). The logarithmic empty strips and the reduction to n^(1/3)<k<=(n-1)/3 are valid. They remove no fixed positive fraction of all possible folded degrees and do not control the interior diagonal mass.

## 4. Finite experiments and their scope

The reflected scan is complete in its stated finite range: r=0,1 cannot be zeros; the recurrence tests every lower member of a reflection pair that can produce an index p+r<=N; every upper member in range is then included. The open interval condition is preserved. There is no sampled-prime omission.

Its comparisons against n/6 and the reference maximum use exact integer inequalities, not floating logarithms. In particular P^6<2^n<e^n suffices for most indices; the remaining small case uses e>8/3. The reference comparison is P^21<17^n, with equality only at n=21. Floating values appear only in labeled rankings or descriptions. The saved root-heavy samples describe zero degrees visible at the chosen index limit, not an asymptotic maximum per prime.

An independent implementation using the integral U recurrence without reflection reproduces the complete n<=10000 result: 957 anomaly events, 911 nonempty indices, maximum cluster size three, and unique maximum S(n)/n=ln17/21. All 10000 inequalities were certified by integers. The author extends its reflected finite scan to 100000; that extension is not a proof at unbounded indices and was not independently rerun in full here.

The six observed three-prime clusters were also checked independently by iterating both original A and B recurrences directly modulo p² for their 18 prime pairs. In all 18 cases, p divides both residues while at least one is nonzero modulo p²; hence v_p(g_n)=1=v_p(n!) and the positive excess is zero. This confirms that these particularly visible first-layer clusters do not themselves give extra-factor gains. It says nothing about unobserved larger indices or all two-prime clusters.

## 5. Why per-prime sparsity alone cannot close the uniform target

A useful logical distinction can be made precise. Consider an abstract family of symmetric zero sets, not the actual recurrence. Choose N_j=4^j so the prime windows (N_j/2,N_j) are disjoint. For each prime in one such window, prescribe the zero pair

\[
\{N_j-p,\ 2p-N_j-1\}.
\]

Omit the O(1) primes that would put an element at 0,1,p-2,p-1 or make the two zeros adjacent. Every remaining set has at most two zeros, is invariant under r↦p-1-r, and has no adjacent zeros. Nevertheless its moving-diagonal sum at N_j has weight N_j/2+o(N_j), by PNT. All other primes may have empty zero sets.

This is **not** a counterexample to the actual S-estimate: the abstract sets do not satisfy the full Legendre recurrence. It shows why even a uniform two-zero-per-prime bound, together with reflection and no-adjacent-zero restrictions, would not alone imply the desired all-index diagonal estimate. A further arithmetic relation across the changing moduli is needed.

## 6. Independent verification

`verify_independent_global_s.py` uses only the Python standard library and imports no author verifier. It passes both aggregate divisibilities and the corresponding finite quality comparisons for n=2,...,500, 11070 upper-half valuation checks, 4225 reflection values, 501 exact odd-height bounds, and 2951 folded edge-product checks. It additionally runs the unreflected n<=10000 scan and the 18 direct mod-p² cluster checks described above. Its JSON records the exact finite coverage, residues, and event-stream hash.

All infinite-index assertions accepted here follow from the displayed algebraic arguments and the explicitly stated PNT input. No average-case result, finite maximum, class/parity observation, or reflection symmetry is promoted into the missing eventual uniform bound. Neither global exclusion of this family nor a conclusion about pi's bad approximability has been established.
