# Independent review of the reduced-denominator exponent map

2026-11-02. Independently reviewed the final `REDUCED_DENOMINATOR_EXPONENT_MAP.md`, `NET_PRIME_BUDGET_AND_ERROR_EXPONENT.md`, `TWO_ADIC_GCD_SATURATION.md`, their verification code, and the substantive claims in the checkpoint README. **The stated results pass within their declared scope. No unresolved mathematical blocker was found.**

The target constants `C>0` and `c>0` are fixed independently of the index. The results concern the specified reduced Legendre–arctangent approximants. They do not identify the optimal exponent of this family, prove a global upper bound for pi's irrationality exponent, or decide bad approximability.

## 1. Exact threshold and its prefactor

Write `R_n=g_n/n!`, so the genuine reduced denominator is exactly `q_n=ell_n/R_n`. Positivity immediately gives

\[
\delta_n\le Cq_n^{-c}
\iff R_n\ge\ell_n(\delta_n/C)^{1/c}.
\]

This is an exact equivalence, including the non-strict boundary. It uses the reduced gcd, not the unreduced integer denominator.

Using the established expansions

\[
\ell_n=L n^{-1/2}\lambda^n(1+O(n^{-1})),
\qquad \delta_n=D_0\lambda^{-2n}(1+O(n^{-1})),
\]

with `lambda=1+sqrt(2)`, `L=sqrt(lambda/(2*pi*sqrt(2)))`, and `D_0=4*pi/lambda`, gives exactly

\[
T_n(C,c)=K(C,c)n^{-1/2}\lambda^{(1-2/c)n}(1+O(n^{-1})),
\]
\[
K(C,c)=\frac{\lambda^{1/2-1/c}(4\pi/C)^{1/c}}
{\sqrt{2\pi\sqrt2}}.
\]

The powers of lambda, pi, and two are correct. The errors are for each fixed `C,c`; uniformity in a changing exponent is not asserted. At `c=2`, cancellation gives `K(C,2)=2^(1/4)/sqrt(C)` and `D_0 L^2=sqrt(2)`.

The normalized logarithmic quantity `S_n` in the map satisfies `ln(delta_n q_n^c/C)=-c S_n+O(1/n)`. Its convergence statements are correct. In particular, a limit `S_n=0` only gives a normalized-error limit of one; it does not determine which side of the exact inequality holds eventually.

## 2. Polynomial boundaries and rate-free exponents

For a conditional profile `R_n=K_0 n^{-a} exp(tau n)(1+o(1))`, direct substitution gives the stated factor `n^{c(a-1/2)}` and exponential rate `c(alpha-tau)-2alpha`. Thus, at the critical exponential rate, the polynomial threshold is exactly `a=1/2`. At that polynomial equality, the constant prefactor must still be compared with the prescribed `C`.

For `c=2`, the equivalent condition for quality tending to zero along a specified diverging subsequence is `sqrt(n) R_n -> infinity`. Merely having the observed exponent tend to two does not imply this condition. Being unbounded above permits extraction of a further diverging subsequence; it is not the same as tending to infinity along an already specified subsequence. The final notes correctly distinguish these statements.

The limsup/liminf formulas need no existence of an exponential growth-rate limit. If `alpha=ln lambda`, `x_n=ln R_n/n`, and `y_n=ln q_n/n`, then

\[
y_n=\alpha-x_n+o(1),\qquad
(-\ln\delta_n)/n\longrightarrow2\alpha>0.
\]

Eventually `q_n>1`, so `y_n>0`. Taking reciprocals of a positive sequence gives precisely

\[
\limsup c_n=\frac{2\alpha}{\alpha-\limsup x_n},\qquad
\liminf c_n=\frac{2\alpha}{\alpha-\liminf x_n}.
\]

If the first denominator is zero, the limsup is positive infinity. If the second denominator is infinite, the liminf is zero. The condition `q_n -> infinity` does not require `y_n` to stay away from zero, and the proof does not assume that it does. These extended-value boundary cases have been handled correctly.

## 3. Exact 2-adic saturation of the gcd

The saturation proof correctly starts from the integer Wronskian and obtains

\[
v_2(r_n-r_{n-1})
=2+2v_2((n-1)!)-v_2(A_n)-v_2(A_{n-1}).
\]

Substitution of the already established all-index A-valuation formula gives all four displayed residue-class formulas. In the class `n=4m`, the identity `s_2(m-1)=s_2(m)-1+v_2(m)` accounts for the coefficient `-3v_2(m)`; its sign is correct.

The nonnegativity proof covers every `n>=5`. For the difficult first class, writing `m=2^t u` gives the asserted lower bound `u(4*2^t-2)-3(t+1)`. At `t=0` the allowed odd `u` is at least three; at `t>=1` the bound is at least `4*2^t-3t-5`, zero at `t=1` and increasing thereafter. The other classes follow from `s_2(m)<=m` and `v_2(m+1)<=m` where used.

The base value `r_4=160/51` is correct. Addition preserves the ring of rational numbers with odd denominator, so the induction proves saturation at every later index. The initial values account separately for `n=0,1,2,3`, with exactly one missing factor of two at `n=3`. Therefore

\[
v_2(g_n)=v_2(A_n)-\mathbf1_{\{n=3\}}
\]

is an all-index theorem, rather than a conjecture inferred from the table. No formula for `v_2(B_n)` is needed.

The two exponent notes use slightly different beta conventions: the error map incorporates the `n=3` correction into its actual binary exponent, whereas the budget note uses `beta_n=v_2(A_n)-v_2(n!)+n/2` and subtracts `t_n` in `J_n`. Both are exact, and the combined quantities agree at the exceptional index.

## 4. Signed odd budget and uniform finite-index bounds

The exact decomposition is `R_n=2^{d_n} E_n/D_n`, where `E_n` and `D_n` are coprime positive odd integers. The positive numerator and negative deficit enter with opposite signs. Removing `D_n` only gives a one-sided upper bound for the gcd gain; it cannot give a sufficient approximation criterion.

The coefficient proof of `lambda^n/(2n+1)<=ell_n<=lambda^n` is valid. After the indicated rescaling, the coefficients are symmetric and unimodal. The explicit convolution-difference formula preserves unimodality, including its center case. The central coefficient consequently lies between the average coefficient and their total.

Combining this with the strict positive-error bounds verifies the two inequalities (7) in the budget note. When `c<2`, the exponent on `ell_n` is negative; the displayed positive-part terms correctly reverse the coefficient bound. The necessary condition (8) is strict because the error lower bound is strict. The sufficient condition (9) can be non-strict because its error upper bound is strict. No integrality assumption on `c` is required.

At `c=2`, the polynomial terms in those two tests match and the remaining logarithmic gap is `(ln 2)/2`. For fixed `C=e`, only the constant terms change. The sharper asymptotic formula retains the explicit binary correction, `ln n` term, and prefactor; it does not pretend that an exponential-rate equality settles the fixed-constant question.

The conditional exclusion under `ln E_n=o(n)` is correct only for

\[
c>c_*:=\frac{2\ln(1+\sqrt2)}{\ln(1+\sqrt2)+(\ln2)/2}.
\]

It does not exclude smaller exponents, and at equality the polynomial factors and odd deficits remain decisive. The actual subexponential excess hypothesis is not proved or silently assumed.

## 5. An unconditional, weaker exponent guarantee

The exact convolution for `w_n`, together with dyadic denominators for `ell_n`, proves `v_p(w_n)>=-floor(log_p n)` for every odd prime. Therefore the odd deficit `D_n` divides the odd part of `lcm(1,...,n)`. The exact binary formula gives `d_n>=-n/2`, including the exceptional index three. It follows that

\[
R_n\ge2^{-n/2}/\operatorname{lcm}(1,\ldots,n).
\]

Invoking the prime number theorem in the standard form `ln lcm(1,...,n)=n+o(n)` yields

\[
\liminf c_n\ge
\frac{2\ln(1+\sqrt2)}{\ln(1+\sqrt2)+(\ln2)/2+1}
=0.7911979206\ldots.
\]

For every fixed `C>0` and every strictly smaller positive target exponent, the normalized error tends to zero. Thus the eventual `e/q^(3/4)` assertion follows. The equality case is properly left open, and no effective first valid index is asserted.

The elementary alternative also checks. In `lcm(1,...,2m)/lcm(1,...,m)`, each new prime appears with exponent one and its new prime power forces a carry in `binom(2m,m)`. For the odd endpoint, a possible new prime power equal to `2m+1` is supplied by the extra factor `2m+1`. This proves the two quotient divisibilities. Repeated halving gives `ln lcm(1,...,n)<=n ln4+O((ln n)^2)` and hence the weaker sufficient exponent quoted in the note, without the prime number theorem.

These are actual sufficient guarantees for this family, but they are weak bounds, not optimality or approximation-record claims.

## 6. Necessary density of simultaneous large-prime gains

For every positive odd-prime excess, the Wronskian cap implies `e_p(n) ln p<=ln(n-1)`. Thus `ln E_n<=N_+(n) ln(n-1)`. At a successful index, the fixed-constant necessary budget supplies `ln E_n>Theta_c n-O_{C,c}(ln n)`.

Small primes `p<=sqrt(n-1)` contribute at most `sqrt(n-1) ln(n-1)=o(n)` in total; this argument needs no prime number theorem. Every positive prime above that cutoff has integer excess exactly one. Hence, whenever `Theta_c>0`, success at unbounded indices requires

\[
N_{\rm large}(n)\ge(\Theta_c-o(1))\frac n{\ln n}
\]

along those successful indices. This is a necessary count of distinct primes at the same index. It is not a construction, not a count observed in the archive's partial-prime scan, and not sufficient after negative deficits are restored.

## 7. Exact necessary counts at index 50000

The count certificate was independently reconstructed without computing `A_50000` or `B_50000`. Here `s_2(50000)=6`, `d^A=-24994`, and the strict lower error bound, with `e<3`, forces

\[
3\cdot100001\cdot49999^{2N_+}>2^{49989}.
\]

Direct integer comparison confirms failure at `N_+=1600` and passage at `1601`.

Independently enumerating the 47 odd primes at most 223 and applying both neighboring Wronskian caps reproduces the 188-digit small-prime capacity product. With that product squared included, the corresponding inequality fails at `N_large=1560` and passes at `1561`.

Thus the two necessary counts in the note and README are correct: **at least 1601 positive-excess odd primes, including at least 1561 above 223**, would be needed for the stated target at index 50000. Neither comparison asserts that these primes actually occur. The supplied partial-prime search was not rerun and is not used to infer absence of all possible large-prime combinations.

## 8. Verification and final scope

The independent `verify_independent_exponent.py` passes 501 actual gcd-saturation checks, 500 rational-difference valuation checks, 500 signed factorizations and odd-deficit/lcm checks, 499 elementary lcm-halving checks, and 320 rational-exponent uniform-bound checks in `Q(sqrt(2))`. It also independently reconstructs both index-50000 count certificates. Results are stored in `verify_independent_exponent.json`.

The final author verifiers were inspected and rerun successfully: 3200 exact error classifications, 2400 reduced-threshold comparisons, and the budget verifier's 500 uniform exponent-bound checks. Their rational exponents are tested by raising positive quantities to integer powers. The pi and Euler-e intervals are rational certificates, not floating-point constants.

The checkpoint README's mathematical statements agree with these results. Its independent executable should be named `exponent-audit/verify_independent_exponent.py`. During review, the budget's subsequence wording was tightened from merely “unbounded above” to “tends to positive infinity” for a specified quality-zero subsequence.

The missing question remains the growth of the net odd ratio `E_n/D_n`, particularly at critical exponential and polynomial scales. The audit supplies no unproved rate assumption. All new all-index assertions follow from displayed arguments, while finite classifications retain their finite scope.
