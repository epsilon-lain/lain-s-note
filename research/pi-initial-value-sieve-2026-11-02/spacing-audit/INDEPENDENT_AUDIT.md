# Independent audit: quantitative spacing

The transfer and prime-overlap argument is valid. It upgrades the earlier
qualitative separation result but does not rule out an infinite sparse
successful subsequence.

Let `U_0=1`, `U_1=2` and

\[
(n+1)U_{n+1}=2(2n+1)U_n+4nU_{n-1}.
\]

Write `B_n=floor(n^(1/3))`. Select odd primes satisfying

\[
p\le n,\quad p^2>n,\quad p(B_n+1)>n,\quad p\mid U_n.
\]

Their squarefree product is `P_n`; thus `log P_n=W(n)` for the selected
multi-band set. Put

\[
K_n=\prod_{b=1}^{B_n}U_b.
\]

For `h>=1`, the exact identity

\[
\left(\prod_{j=1}^h(n+j)\right)U_{n+h}
=F_h(n)U_n+G_h(n)U_{n-1}
\]

has `F_0=1`, `G_0=0`, `F_1=2(2n+1)`, `G_1=4n`, and

\[
V_{h+1}=2(2n+2h+1)V_h+4(n+h)^2V_{h-1}.
\]

No division modulo a prime is used. If a prime divides `U_n` and
`U_{n+h}` but does not divide `U_{n-1}`, it divides `G_h(n)`.

The remaining primes divide `K_n`: write `n=bp+r`. The Lucas congruence
gives `U_n=U_b U_r (mod p)`. If `r=0`, central divisibility already forces
`p|U_b`. If `r>0` and `p` does not divide `U_b`, consecutive zeros
`U_r=U_(r-1)=0 (mod p)` follow. The recurrence can be run backwards
because `1<=r<p` and `p` is odd; it eventually contradicts `U_0=1`.
Finally `p(B_n+1)>n` implies `b<=B_n`.

Consequently the following finite integer certificate holds, for `m=n+h`:

\[
\boxed{P_nP_m\ \mid\
\left(\prod_{\substack{p\le m\\p\text{ odd}}}p\right)K_nG_h(n).}
\]

For `m=n+h`, direct induction gives

\[
0<G_h(n)\le(5m)^h.
\]

Indeed the recurrence coefficients for intermediate indices are at most
`4m` and `4m^2`, and `4/5+4/25<1`. Similarly `U_b<=6^b`, hence

\[
\log K_n\le\tfrac12(\log6)B_n(B_n+1)=O(n^{2/3}).
\]

Taking logarithms gives the completely finite overlap bound

\[
W(n)+W(m)\le\vartheta(m)+h\log(5m)+O(n^{2/3}).
\]

Assume the prior multi-band theorem, namely every fixed-quality event
`Q_k<=M` satisfies

\[
W(k)\ge c k-o(k),\qquad c=(2+\log2)/4.
\]

Then for every fixed `M>0` and `epsilon>0`, there is `N(M,epsilon)` such
that whenever `m>n>=N` and both `Q_n,Q_m<=M`,

\[
\boxed{m-n\ge
\left(\tfrac12\log2-\epsilon\right)\frac n{\log n}.}
\]

For a rigorous quantifier argument, a sequence violating this estimate
would have `h=O(n/log n)=o(n)`. On that sequence the prime number theorem
and the finite overlap bound force

\[
h\log(5m)\ge (2c-1)n-o(n)
=\tfrac12(\log2)n-o(n),
\]

which contradicts the assumed violation. Large relative gaps need no
separate asymptotic approximation.

The number of fixed-quality indices up to `N` is therefore at most

\[
\left(\frac1{\log2}+o(1)\right)(\log N)^2.
\]

This follows by summing increments of `(log n)^2`; successive successful
indices have an increment at least `log2-o(1)`.

The refinement that distinguishes exact first valuations does not by
itself improve the constant: its prime weights and the shared-prime
capacity can both double. An additional distribution theorem is still
required for an improved universal constant or for finiteness.

Independent exact verification checks 12,000 `(n,h)` pairs with
`1<=n<=400`, `1<=h<=30`; it verifies the transfer identity, the integer
divisibility certificate and its height bound. These tests support the
algebra and are not a replacement for the asymptotic argument.
