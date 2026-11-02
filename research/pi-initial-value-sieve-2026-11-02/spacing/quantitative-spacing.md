# Quantitative spacing from common large prime divisors

Let `U_n=2^n A_n/n!`. It is the generalized central trinomial sequence
`[x^n](1+2x+2x^2)^n`, with `U_0=1,U_1=2` and

\[
(n+1)U_{n+1}=2(2n+1)U_n+4nU_{n-1}.
\]

Let `t_n` be the lower cutoff in the existing multi-band theorem, so that
`t_n >= sqrt(n)` and `K_n=floor(n/t_n)=O(n^(1/3))`. Define

\[
\mathcal P_n=\{p:\ t_n<p\le n,\ p\text{ odd},\ p\mid U_n\},\qquad
W_n=\sum_{p\in\mathcal P_n}\log p.
\]

For all integers `m>n>=1`, putting `h=m-n`, one has the finite inequality

\[
W_n+W_m\le \vartheta(m)+h\log(5m)
 +\frac{\log6}{2}K_n(K_n+1).
\tag{1}
\]

The existing quality theorem implies `W_n >= c n-o(n)` whenever `Q_n` is
bounded, with `c=1/2+log(2)/4`. Consequently, for any fixed `M` and positive
epsilon, there is a threshold after which any two distinct indices with
`Q_n,Q_m<=M` satisfy

\[
m-n\ge \left(\frac{\log2}{2}-\epsilon\right)\frac n{\log n}.
\tag{2}
\]

It follows that their counting function up to `N` is at most

\[
\left(\frac1{\log2}+o(1)\right)(\log N)^2.
\tag{3}
\]

This does not prove finitely many successful indices and does not settle pi's
bad approximability. It is stronger than the previous fixed-gap exclusion.

## Proof of (1)

Define integer polynomials `F_h(n),G_h(n)` through

\[
F_0=1,\ G_0=0,\quad F_1=2(2n+1),\ G_1=4n,
\]
\[
V_{h+1}=2(2n+2h+1)V_h+4(n+h)^2V_{h-1}\quad(V=F,G).
\]

Induction gives the exact integer identity

\[
\prod_{j=1}^{h}(n+j)\,U_{n+h}
 =F_h(n)U_n+G_h(n)U_{n-1}.
\tag{4}
\]

Thus a common prime divisor of `U_n,U_m` either divides `U_(n-1)` or
divides `G_h(n)`. This requires no inversion of `(n+1)...(n+h)` modulo `p`.
The latter polynomial is positive for `n,h>=1`. For all `k<=h` its recurrence
coefficients are bounded above by `4m` and `4m^2`, respectively, so induction
using `4/5+4/25<1` gives

\[
0<G_h(n)\le(5m)^h.
\]

For the adjacent-zero alternative, write `n=bp+r`, `0<=r<p`. As `p>sqrt(n)`,
we have `b<p`, and the Lucas congruence gives `U_n=U_b U_r (mod p)`.
If `r=0`, then `p|U_n` already implies `p|U_b`. If `r>=1` and `p` does
not divide `U_b`, adjacent zeros would force `U_r=U_(r-1)=0 (mod p)`.
The three-term recurrence, run backwards while indices remain below `p`,
would force `U_0=0`, a contradiction. Therefore every bad adjacent prime in
`P_n` divides

\[
L_n=\prod_{b=1}^{K_n}U_b.
\]

The recurrence also proves `U_b<=6^b`: after substituting the induction
bounds, the normalized factor is `(28b+12)/(36b+36)<=1`. Hence
`log L_n <= log(6) K_n(K_n+1)/2`.

Since the common prime product of `P_n,P_m` divides `G_h(n)L_n`, its log
is at most the last two terms of (1). The union prime mass is at most
`theta(m)`, proving (1).

## Proof of (2) and (3)

If (2) fails along any infinite sequence, then `h=O(n/log n)=o(n)`.
The left side of (1) is at least `c(n+m)-o(n)`, while
`theta(m)=m+o(n)`. It follows that

\[
\frac{\log2}{2}n\le h\log n+o(n),
\]

a contradiction. For (3), consecutive successful indices therefore increase
`(log n)^2` by at least `log(2)-o(1)`, and sum these increments.

## Relevant primary literature and an alternative target

Tony D. Noe, *On the Divisibility of Generalized Central Trinomial
Coefficients*, Journal of Integer Sequences 9 (2006), Article 06.2.7:
https://cs.uwaterloo.ca/journals/JIS/VOL9/Noe/noe35.pdf
contains this exact sequence as `(a,b,c,d)=(1,2,2,-4)`, OEIS A006139, and
proves generalized Schur/Lucas and Holt reflection congruences.

Nadav Kohen, *Density and Symmetry in the Generalized Motzkin Numbers Modulo
p*, Integers 26 (2026), A22:
https://math.colgate.edu/~integers/aa22/aa22.pdf
Proposition 4 gives Lucas congruences for any constant term of a Laurent
polynomial of degree one. Its density results are fixed-prime statements;
they do not bound the moving-prime mass `W_n`.

Define the active primes
`Z={p odd: exists 0<=r<p with p|U_r}`. Lucas implies `W_n<=theta_Z(n)`.
Therefore a uniform upper density bound
`theta_Z(x) <= (0.6732867951-eta)x` would eliminate all bounded-quality
subsequences. Exact root scans through 5000 give
`theta_Z(5000)/theta(5000) approx 0.41455`; this is evidence for choosing a
research target, not a proof of its asymptotic density. Fixed-prime density
one of divisibility, when a zero digit exists, must not be mistaken for
rarity of moving-prime divisibility.
