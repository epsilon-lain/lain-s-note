# A pi-free necessary-condition sieve and a complete finite exclusion

2026-11-02. This note concerns only the fixed Legendre–arctangent family in the preceding checkpoint. It proves an inexpensive all-index quality lower bound, then applies it to a complete finite search. It does not exclude successful indices beyond the checked range and does not scan all rational numbers.

## 1. Fixed model and inherited exact error bounds

Let A_0=A_1=1, B_0=0,B_1=1, and let both integer sequences satisfy

\[
X_{n+1}=(2n+1)X_n+n^2X_{n-1},\qquad n\ge1.
\]

Set

\[
\ell_n=A_n/n!,\quad w_n=B_n/n!,\quad
g_n=\gcd(A_n,4B_n),\quad
p_n=4B_n/g_n,\quad q_n=A_n/g_n,
\quad Q_n=q_n^2|\pi-p_n/q_n|.
\]

The prior positive Legendre error identity gives the strict bounds, for every n≥1,

\[
\boxed{\frac{2}{2n+1}\left(\frac{n!}{g_n}\right)^2
<Q_n<\frac{4}{2n+1}\left(\frac{n!}{g_n}\right)^2.}
\tag{1}
\]

They are exact analytic inequalities. All finite comparisons below therefore require integer arithmetic only, with no decimal enclosure of pi.

## 2. A necessary sieve using U_n and an lcm only

Define the positive integer

\[
U_n=2^n\ell_n=[x^n](1+2x+2x^2)^n.
\]

The coefficient formula proves integrality. Equivalently,

\[
U_0=1,\quad U_1=2,\quad
nU_n=2(2n-1)U_{n-1}+4(n-1)U_{n-2}\quad(n\ge2).
\tag{2}
\]

For n≥2 let

\[
M_{n-1}=\operatorname{lcm}(1,\ldots,n-1),\quad
E_n=\frac{M_{n-1}}{\gcd(M_{n-1},n)},\quad
G_n=\gcd(\operatorname{odd}(U_n),\operatorname{odd}(E_n)),
\]

where odd(k) is the positive odd part of the integer k. Write s_n=v_2(U_n).

The [odd-prime Wronskian bound](../gcd-prime-blocks/ODD_PRIME_BLOCKS_AND_EXCESS.md) proves, for every odd prime p,

\[
e_p(n):=v_p(g_n)-v_p(n!)
\le\lfloor\log_p(n-1)\rfloor-v_p(n).
\]

If e_p(n)>0, then also e_p(n)≤v_p(ell_n)=v_p(U_n). Since

\[
v_p(E_n)=\max\{\lfloor\log_p(n-1)\rfloor-v_p(n),0\},
\]

the odd numerator of the reduced fraction g_n/n! divides G_n. Dropping its odd denominator can only enlarge that fraction. For the power of two, g_n divides A_n, so

\[
v_2(g_n/n!)\le v_2(\ell_n)=s_n-n.
\]

Consequently the following unconditional bound holds for n≥2:

\[
\boxed{\frac{g_n}{n!}\le2^{s_n-n}G_n,\qquad
Q_n>\frac{2\,2^{2(n-s_n)}}{(2n+1)G_n^2}.}
\tag{3}
\]

The power of two is a rational power when its exponent is negative, so implementations must not use truncating integer division. This bound uses no unproved formula for s_n and does not require computing B_n.

For n=1, use the explicit A_1=B_1=1 and (1), which already excludes quality at most 1/5. The code's E_1=G_1=1 convention agrees with that check.

The forward Wronskian and the stronger neighboring-pole condition permit further sieve improvements. The executable certificate deliberately uses only the backward bound (3); those optional refinements are unnecessary for its result.

## 3. Exact finite search to n=20000

For the threshold Q_n≤1/5, a necessary condition is that the lower bound in (3) be strictly below 1/5. If it is equal to 1/5, the strict inequality already excludes the index.

`sieve_gcd.py --limit 20000` checked every index from 1 through 20000. The only survivors of this necessary condition were

\[
2,\ 3,\ 6,\ 7,\ 15.
\]

Survival is not a good-approximation certificate. The program then constructed the actual integer sequences and gcds only through the last survivor, n=15. Applying (1) to these five indices leaves only n=2. At n=2,

\[
p_2/q_2=3/1,\qquad 1/10<Q_2<1/5.
\]

Every other index in the complete checked range has Q_n>1/5. Thus n=2 is the unique best-quality index within this finite family. This says nothing about the best approximation among all rational numbers, and the fraction 3/1 is not a new approximation record.

The certificate includes the U_n/lcm/G_n chain SHA-256, all five exact survivor checks, and selected bit-size summaries. It can be independently rerun with the Python standard library. If an alternative limit leaves a survivor beyond the optional exact-stage cap, the report explicitly marks that stage incomplete and makes no complete-exclusion claim.

## 4. Separate full-gcd cross-check through n=5000

`scan_gcd.py --limit 5000` directly constructs A_n,B_n,n!,g_n at every index, independently of the U/lcm sieve. It certifies

\[
\max_{1\le n\le5000} n(g_n/n!)^2=8\quad\text{attained at }n=2,
\]

\[
\max_{3\le n\le5000} n(g_n/n!)^2=4/3\quad\text{attained at }n=3.
\]

In that range, all n≠2 have the stronger floor Q_n>9/14. This floor is not claimed for the additional sieve-only range up to n=20000.

The direct scan also finds the failure of the no-odd-excess shortcut at n=72 and p=17. Through n=5000 the largest odd reduced numerator of g_n/n! found was 2677 at n=4762. Its apparent smallness is a finite observation, not a uniform bound.

## 5. What remains

The all-index theorem here is the necessary condition (3), not a finite search extrapolation. To exclude every sparse successful subsequence, a uniform quantitative bound for the relevant odd-prime product is still needed. The exact 2-adic theorem in the companion analysis makes the necessary odd growth exponential, but the present work has not proved that this odd growth is impossible.
