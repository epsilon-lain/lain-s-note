> 历史资料恢复说明：本文件是 2026-09-06 项目原稿从 Library 分段读取后重建的文本；分页/换行及字节不保证与原始文件完全相同。原始文件：[Library 原稿](https://chatgpt.com/library?file_id=libfile_ba1b267839d08191a17f48795005b687&account_id=personal)。其证明、有限实验和引用尚未在本次逐条独立核验；日期和状态均为历史记录。

# P2848 — PROCESS-DNA / BBP / RADIX-ORDER ROUND
## Consolidated theorems, derivations, stop boundaries, and active direction

**Date:** 2026-09-06  
**Status:** `P2848 NOT PROVED`

This note consolidates the full recent round beginning with the two-time process quantity, its collapse to ordinary critical approximation, the base-\(b\) generalization, the BBP/base-16 attack, the carry obstruction, the smooth-factor stop boundary, and the surviving continued-fraction × radix-order direction.

The user's process object \(p\) remains the complete synchronized digit-generation / left-shift process. Finite quantities \(P_t^{(b)}=\lfloor b^t\pi\rfloor\) are coordinate shadows only.

---

# 0. Target

A sufficient form of P2848 is:

\[
\exists\ \text{reduced }a_j/q_j,\qquad
q_j\to\infty,\qquad
q_j^2\left|\pi-\frac{a_j}{q_j}\right|\to0.
\]

Indeed,

\[
|q_j\pi-a_j|
=
q_j\left|\pi-\frac{a_j}{q_j}\right|,
\]

so if \(a_j/q_j\to\pi\),

\[
a_j|q_j\pi-a_j|
=
\frac{a_j}{q_j}
q_j^2\left|\pi-\frac{a_j}{q_j}\right|
\to0.
\]

Since

\[
|\sin a_j|
\le
|a_j-\pi q_j|,
\]

we obtain

\[
a_j|\sin a_j|\to0.
\]

Thus every route below is judged by whether it gives an ordinary reduced rational sequence with

\[
q^2|\pi-a/q|\to0.
\]

---

# 1. Two-time decimal process: first exact collapse

Define

\[
P_t=\lfloor10^t\pi\rfloor,
\qquad
E_t=10^t\pi-P_t.
\]

Then

\[
0<E_t<1.
\]

For \(s\ge0,m\ge1\), define the true \(m\)-digit block following position \(s\):

\[
B_{s,m}
:=
P_{s+m}-10^mP_s.
\]

Since

\[
P_{s+m}
=
\lfloor10^m(P_s+E_s)\rfloor
=
10^mP_s+\lfloor10^mE_s\rfloor,
\]

we have the exact identity

\[
\boxed{
B_{s,m}=\lfloor10^mE_s\rfloor.
}
\tag{1.1}
\]

Also,

\[
E_{s+m}
=
10^{s+m}\pi-P_{s+m}
=
10^mE_s-B_{s,m},
\]

hence

\[
\boxed{
E_{s+m}-E_s
=
(10^m-1)E_s-B_{s,m}.
}
\tag{1.2}
\]

Therefore

\[
\boxed{
\left|
E_s-\frac{B_{s,m}}{10^m-1}
\right|
=
\frac{|E_{s+m}-E_s|}{10^m-1}.
}
\tag{1.3}
\]

So the pair

\[
\gcd(B_{s,m},10^m-1)
\]

and

\[
|E_{s+m}-E_s|
\]

are respectively the reduction factor and approximation error of the single rational

\[
\frac{B_{s,m}}{10^m-1}
\]

approximating \(E_s\).

---

# 2. The old \(\Xi_{s,m}\) Boss collapses to a critical rational score

Let

\[
M:=10^m-1,
\qquad
d:=\gcd(B_{s,m},M).
\]

Write

\[
B_{s,m}=db,
\qquad
M=dq,
\qquad
\gcd(b,q)=1.
\]

From (1.2),

\[
E_{s+m}-E_s
=
d(qE_s-b).
\]

The previous Boss was

\[
\Xi_{s,m}
=
\frac{10^sM|E_{s+m}-E_s|}{d^2}.
\]

Substitution gives

\[
\boxed{
\Xi_{s,m}
=
10^sq\,|qE_s-b|.
}
\tag{2.1}
\]

Using

\[
E_s=10^s\pi-P_s,
\]

we get

\[
qE_s-b
=
10^sq\,\pi-(qP_s+b).
\]

Define

\[
n:=10^sq,
\qquad
k:=qP_s+b\in\mathbb Z.
\]

Then

\[
\boxed{
\Xi_{s,m}
=
n|n\pi-k|
=
n^2\left|\pi-\frac{k}{n}\right|.
}
\tag{2.2}
\]

If \(k/n\) is not fully reduced and \(r=\gcd(k,n)\), then after reduction

\[
\frac{k}{n}=\frac{a}{q'},
\qquad
q'=n/r,
\]

and

\[
\boxed{
(q')^2\left|\pi-\frac{a}{q'}\right|
=
\frac{\Xi_{s,m}}{r^2}
\le
\Xi_{s,m}.
}
\tag{2.3}
\]

Thus \(\Xi\to0\) is sufficient, but the numerical quantity itself is not a new invariant: it is a process-generated critical rational approximation score.

If \(d\ge2\), then

\[
|qE_s-b|
=
\frac{|E_{s+m}-E_s|}{d}
<
\frac12,
\]

so \(b\) is the nearest integer to \(qE_s\), and

\[
\boxed{
\Xi_{s,m}
=
n\|n\pi\|.
}
\tag{2.4}
\]

### Stop boundary

The new content cannot be the identity \(\Xi=n\|n\pi\|\) itself.  
Any gain must come from special arithmetic restrictions on the process-generated integers \(n\).

---

# 3. Exact full reduction of the two-time periodic shadow

Keep

\[
M=10^m-1,
\qquad
B=B_{s,m},
\qquad
d=\gcd(B,M),
\]

and write

\[
B=db,\qquad M=dC,\qquad (b,C)=1.
\]

The eventually periodic rational shadow is

\[
R_{s,m}
=
\frac{P_{s+m}-P_s}{10^s(10^m-1)}.
\]

Since

\[
P_{s+m}-P_s
=
MP_s+B
=
d(CP_s+b),
\]

we obtain

\[
R_{s,m}
=
\frac{CP_s+b}{10^sC}.
\]

Define

\[
K:=CP_s+b.
\]

Because

\[
(K,C)=1,
\]

the only remaining common factor comes from \(10^s\). Define

\[
h_{s,m}:=\gcd(K,10^s).
\]

Then the fully reduced denominator is

\[
\boxed{
q_{s,m}
=
\frac{10^sC}{h_{s,m}}.
}
\tag{3.1}
\]

Moreover,

\[
\boxed{
\gcd(P_{s+m}-P_s,\ 10^s(10^m-1))
=
d\,h_{s,m}.
}
\tag{3.2}
\]

Let

\[
\Delta_{s,m}
:=
|E_{s+m}-E_s|.
\]

The raw secant error is

\[
\left|\pi-R_{s,m}\right|
=
\frac{\Delta_{s,m}}{10^s(10^m-1)}.
\]

Hence the fully reduced critical score is

\[
\boxed{
\Theta_{s,m}
:=
q_{s,m}^2|\pi-R_{s,m}|
=
\frac{
10^sC^2\Delta_{s,m}
}{
(10^m-1)h_{s,m}^2
}.
}
\tag{3.3}
\]

This is exact.

A sufficient condition is therefore

\[
\boxed{
\Theta_{s_j,m_j}\to0.
}
\tag{3.4}
\]

Define the logarithmic surplus

\[
\boxed{
\mathcal S_{s,m}
=
\log_{10}(10^m-1)
-s
-2\log_{10}C
+2\log_{10}h
-\log_{10}\Delta.
}
\tag{3.5}
\]

Then

\[
\Theta_{s,m}=10^{-\mathcal S_{s,m}},
\]

up to the exact use of \(\log_{10}(10^m-1)\), so

\[
\boxed{
\mathcal S_{s_j,m_j}\to+\infty
\Longrightarrow
P2848.
}
\tag{3.6}
\]

### Interpretation

The previous apparent gains

- cyclotomic reduction \(d\),
- \(2,5\)-adic prefix reduction \(h\),
- Archimedean tail return \(\Delta\),

are representation-dependent components of one invariant reduced rational score.

---

# 4. General base-\(b\) version

The process theory is not restricted to base 10.

For any fixed integer base \(b\ge2\), define

\[
P_t^{(b)}:=\lfloor b^t\pi\rfloor,
\qquad
E_t^{(b)}:=b^t\pi-P_t^{(b)}.
\]

For \(s,m\),

\[
M:=b^m-1,
\]

\[
B_{s,m}^{(b)}
:=
P_{s+m}^{(b)}-b^mP_s^{(b)}.
\]

Let

\[
d=\gcd(B,M),
\qquad
B=d\beta,
\qquad
M=dC,
\qquad
(\beta,C)=1,
\]

and

\[
K:=CP_s^{(b)}+\beta,
\qquad
h:=\gcd(K,b^s).
\]

Then the fully reduced two-time shadow has denominator

\[
\boxed{
q_{s,m}^{(b)}
=
\frac{b^sC}{h}.
}
\tag{4.1}
\]

Let

\[
\Delta_{s,m}^{(b)}
=
|E_{s+m}^{(b)}-E_s^{(b)}|.
\]

Then

\[
\boxed{
\Theta_{s,m}^{(b)}
=
\left(q_{s,m}^{(b)}\right)^2
\left|
\pi-R_{s,m}^{(b)}
\right|
=
\frac{
b^sC^2\Delta_{s,m}^{(b)}
}{
(b^m-1)h^2
}.
}
\tag{4.2}
\]

Thus

\[
\boxed{
\Theta_{s_j,m_j}^{(b)}\to0
\Longrightarrow
P2848.
}
\tag{4.3}
\]

This allows us to choose a radix adapted to a special \(\pi\)-formula.

---

# 5. Why base 16 was selected: BBP

The Bailey--Borwein--Plouffe formula is

\[
\boxed{
\pi
=
\sum_{k=0}^{\infty}
\frac1{16^k}
\left(
\frac4{8k+1}
-\frac2{8k+4}
-\frac1{8k+5}
-\frac1{8k+6}
\right).
}
\tag{5.1}
\]

Define

\[
c_k
:=
\frac4{8k+1}
-\frac2{8k+4}
-\frac1{8k+5}
-\frac1{8k+6}.
\]

Combining the four fractions gives

\[
\boxed{
c_k
=
\frac{
120k^2+151k+47
}{
(2k+1)(4k+3)(8k+1)(8k+5)
}.
}
\tag{5.2}
\]

In particular \(c_k>0\).

A crude bound sufficient for the later tail estimates is, for \(k\ge1\),

\[
\boxed{
\frac1{21k^2}<c_k<\frac1{k^2}.
}
\tag{5.3}
\]

For the lower bound, after cross multiplication the difference is

\[
2008k^4+2147k^3+275k^2-194k-15>0
\]

for \(k\ge1\).

---

# 6. Direct BBP truncation is a dead route

Let

\[
S_N:=\sum_{k=0}^{N}16^{-k}c_k.
\]

The tail

\[
T_N:=\pi-S_N
\]

satisfies, using only the first omitted term,

\[
T_N
>
\frac{
16^{-(N+1)}
}{
21(N+1)^2
}.
\tag{6.1}
\]

Now inspect the odd part of the reduced denominator of \(S_N\).

For every prime

\[
4N<p\le8N,
\qquad
p\equiv1\text{ or }5\pmod8,
\]

the prime \(p\) occurs in exactly one BBP denominator \(8k+1\) or \(8k+5\) with \(k\le N\), and no other term has a \(p\)-adic denominator.

The relevant numerator is nonzero modulo \(p\), so \(p\) survives in the reduced denominator.

By the prime number theorem in arithmetic progressions,

\[
\sum_{\substack{4N<p\le8N\\p\equiv1,5\ (8)}}\log p
=
2N+o(N).
\]

Therefore the reduced denominator \(q_N\) satisfies

\[
\boxed{
q_N\ge\exp((2+o(1))N).
}
\tag{6.2}
\]

Combining with (6.1),

\[
q_N^2|\pi-S_N|
\ge
\frac{
\exp((4-\log16-o(1))N)
}{
O(N^2)
}.
\]

Since

\[
4-\log16
=
1.227411\ldots>0,
\]

we get

\[
\boxed{
q_N^2|\pi-S_N|
\to+\infty.
}
\tag{6.3}
\]

### Stop boundary

\[
\boxed{
\text{Direct denominator-cleared BBP truncation is DEAD.}
}
\]

Its critical score diverges exponentially in the wrong direction.

---

# 7. BBP modular dynamics survives, but ordinary orbit tracking is too weak

Define a finite BBP head at radix position \(s\):

\[
F_s
:=
\sum_{k=0}^{s}16^{s-k}c_k.
\]

Write

\[
16^s\pi
=
F_s+T_s,
\]

where

\[
T_s
=
\sum_{k=s+1}^{\infty}16^{s-k}c_k.
\]

The tail obeys

\[
\boxed{
T_s>\frac1{336(s+1)^2}.
}
\tag{7.1}
\]

Indeed the first omitted term contributes

\[
16^{-1}c_{s+1}
>
\frac1{336(s+1)^2}
\]

up to an immaterial index shift in the crude bound.

The BBP modular algorithm can compute the fractional part of \(16^s\pi\) through modular powers of \(16\) modulo the linear denominators.

This gives a genuine \(\pi\)-specific perturbed base-16 dynamical system.

However, polynomial tracking such as \(O(s^{-2})\) is astronomically weaker than the exponential shrinking scale required by P2848.

Thus ordinary distribution / normality / equidistribution statements do not close the target.

---

# 8. A single exact modular return does NOT force a process-block divisor

A tempting false implication was

\[
r\mid16^m-1
\quad\Longrightarrow\quad
r\mid B_{s,m}^{(16)}.
\]

This is false.

Take

\[
r=17=8\cdot2+1.
\]

Since

\[
16^2\equiv1\pmod{17},
\]

the \(r=17\) BBP component has exact modular return after \(m=2\).

But the hexadecimal expansion begins

\[
\pi=(3.243F6A88\ldots)_{16}.
\]

For \(s=2,m=2\), the true two-digit block is

\[
B_{2,2}^{(16)}
=
(3F)_{16}
=
63.
\]

Yet

\[
63\not\equiv0\pmod{17}.
\]

Therefore

\[
\boxed{
\operatorname{ord}_r(16)\mid m
\centernot\Longrightarrow
r\mid B_{s,m}^{(16)}.
}
\tag{8.1}
\]

The missing mechanism is carry / phase.

---

# 9. General divisor-carry theorem

This theorem is independent of BBP and is one of the cleanest results of this round.

Fix a base \(b\ge2\).

Let

\[
M=b^m-1,
\qquad
E_s=\{b^s\pi\},
\qquad
B_{s,m}=\lfloor b^mE_s\rfloor.
\]

Take any divisor

\[
D\mid M,
\qquad
L:=M/D.
\]

Write

\[
LE_s=N+u,
\qquad
N\in\mathbb Z,
\quad
0\le u<1.
\]

Then

\[
B_{s,m}
=
\lfloor(M+1)E_s\rfloor
=
DN+\lfloor Du+E_s\rfloor.
\]

Therefore

\[
\boxed{
B_{s,m}
\equiv
\lfloor D\{LE_s\}+E_s\rfloor
\pmod D.
}
\tag{9.1}
\]

Since

\[
0\le D\{LE_s\}+E_s<D+1,
\]

we obtain the exact criterion

\[
\boxed{
D\mid B_{s,m}
}
\]

if and only if either

\[
\boxed{
\{LE_s\}
<
\frac{1-E_s}{D}
}
\tag{9.2a}
\]

or

\[
\boxed{
\{LE_s\}
\ge
1-\frac{E_s}{D}.
}
\tag{9.2b}
\]

The two target intervals have total length exactly

\[
\boxed{
\frac1D.
}
\tag{9.3}
\]

Since

\[
E_s=b^s\pi-P_s
\]

and \(LP_s\in\mathbb Z\),

\[
\{LE_s\}
=
\{b^sL\pi\}.
\]

Thus the divisor condition is equivalent to a shrinking-target condition for the ordinary \(\pi\)-orbit:

\[
\boxed{
D\mid B_{s,m}
\iff
\{b^sL\pi\}
\text{ lies in a total window of width }1/D.
}
\tag{9.4}
\]

If this occurs, then

\[
\|b^sL\pi\|<\frac1D,
\]

and hence

\[
\boxed{
b^sL\,\|b^sL\pi\|
<
\frac{b^sL}{D}
=
\frac{b^sL^2}{b^m-1}.
}
\tag{9.5}
\]

### Consequence

A huge smooth / cyclotomic factor of \(b^m-1\) does not itself force a large process gcd.  
Factorization abundance and phase forcing are separate problems.

---

# 10. Full BBP-head synchronization and the carry integer

Write the finite BBP head in reduced form

\[
F_s=\frac{A_s}{R_s},
\qquad
(A_s,R_s)=1.
\]

The denominator \(R_s\) is odd, because after combining \(c_k\), all denominators in (5.2) are odd and the powers of \(16\) lie in the numerator of \(F_s\).

Assume

\[
R_s\mid16^m-1.
\]

Define

\[
M:=16^m-1,
\qquad
G:=\frac{M}{R_s}.
\]

Using

\[
E_s
=
F_s+T_s-P_s,
\]

we have

\[
(M+1)E_s
=
M(F_s-P_s)+MT_s+E_s.
\]

The first term is an integer multiple of \(G\). Therefore

\[
\boxed{
B_{s,m}^{(16)}
\equiv
C_{s,m}
\pmod G,
}
\tag{10.1}
\]

where

\[
\boxed{
C_{s,m}
:=
\lfloor MT_s+E_s\rfloor.
}
\tag{10.2}
\]

Hence

\[
\boxed{
\gcd(B_{s,m}^{(16)},G)
=
\gcd(C_{s,m},G).
}
\tag{10.3}
\]

This is the exact BBP carry reduction.

---

# 11. The full-head carry condition itself is another shrinking target

Since

\[
R_sT_s
=
16^sR_s\pi-A_s,
\]

we have

\[
\boxed{
\{R_sT_s\}
=
\{16^sR_s\pi\}.
}
\tag{11.1}
\]

Let

\[
x_s:=\{16^sR_s\pi\}.
\]

Because

\[
M=R_sG,
\]

write

\[
R_sT_s=N+x_s.
\]

Then

\[
C_{s,m}
=
GN+\lfloor Gx_s+E_s\rfloor.
\]

Therefore

\[
\boxed{
G\mid B_{s,m}^{(16)}
}
\]

if and only if

\[
\boxed{
x_s<
\frac{1-E_s}{G}
}
\tag{11.2a}
\]

or

\[
\boxed{
x_s\ge
1-\frac{E_s}{G}.
}
\tag{11.2b}
\]

Again the total window width is

\[
1/G.
\]

Thus

\[
\boxed{
\text{full BBP modular synchronization does not force the gcd;}
}
\]

it only produces the modulus \(G\). A separate phase theorem at scale \(1/G\) is still required.

---

# 12. Strong lower bound for the BBP-head denominator

The finite BBP head denominator \(R_s\) is much larger than \(16^s\) on the exponential scale.

Let

\[
N(k)=120k^2+151k+47.
\]

For a prime \(p\) appearing in one of the large linear factors, the following reductions hold:

\[
p=8k+1
\quad\Longrightarrow\quad
64N(k)\equiv1920\pmod p,
\tag{12.1}
\]

\[
p=8k+5
\quad\Longrightarrow\quad
64N(k)\equiv-32\pmod p,
\tag{12.2}
\]

\[
p=4k+3
\quad\Longrightarrow\quad
16N(k)\equiv20\pmod p.
\tag{12.3}
\]

Thus sufficiently large such primes do not cancel against the numerator.

Now consider:

### Range I

\[
\frac{8}{3}s<p\le4s.
\]

Every odd prime in this interval appears uniquely in one of

\[
8k+1,\quad8k+5,\quad4k+3
\]

for \(k\le s\), and no larger odd multiple can occur in another BBP denominator factor.

Its total logarithmic prime weight is

\[
\frac43s+o(s).
\]

### Range II

\[
4s<p\le8s,
\qquad
p\equiv1,5\pmod8.
\]

These primes appear uniquely in \(8k+1\) or \(8k+5\), with total logarithmic weight

\[
2s+o(s).
\]

Therefore

\[
\boxed{
\log R_s
\ge
\left(\frac{10}{3}+o(1)\right)s.
}
\tag{12.4}
\]

Equivalently,

\[
\boxed{
R_s
\ge
16^{(1.2022458\ldots+o(1))s}.
}
\tag{12.5}
\]

---

# 13. Consequences for full-head synchronization

If

\[
R_s\mid16^m-1,
\]

then

\[
16^m-1\ge R_s.
\]

Hence

\[
\boxed{
m
\ge
(1.2022458\ldots+o(1))s.
}
\tag{13.1}
\]

So the exact BBP modular clock closes more slowly than the digit-position clock.

Also from the positive BBP tail,

\[
T_s>\frac1{336(s+1)^2}.
\]

Therefore

\[
(16^m-1)T_s
\ge
R_sT_s
>
\frac{
\exp((10/3+o(1))s)
}{
336(s+1)^2
}
\to\infty.
\tag{13.2}
\]

Thus the carry

\[
C_{s,m}
=
\lfloor(16^m-1)T_s+E_s\rfloor
\]

is exponentially large along full-head synchronization.

Hence

\[
\boxed{
\text{full-head synchronization + carry-free mechanism is DEAD.}
}
\tag{13.3}
\]

---

# 14. Carry-free synchronization cannot be a hidden shortcut

Suppose, hypothetically, that carry-free synchronization gave

\[
C_{s,m}=0.
\]

Then

\[
(16^m-1)T_s<1.
\]

If the resulting full \(G\)-saving were enough for a critical hit, we would need

\[
\frac{
16^sR_s^2
}{
16^m-1
}
\to0.
\]

Combining with

\[
16^m-1<1/T_s
\]

would imply

\[
\boxed{
16^sR_s^2T_s\to0.
}
\tag{14.1}
\]

But the direct BBP rational

\[
\frac{A_s}{16^sR_s}
\]

has error

\[
\left|
\pi-\frac{A_s}{16^sR_s}
\right|
=
\frac{T_s}{16^s},
\]

so its raw critical score is exactly

\[
\boxed{
(16^sR_s)^2
\frac{T_s}{16^s}
=
16^sR_s^2T_s.
}
\tag{14.2}
\]

Therefore if carry-free synchronization were strong enough to solve P2848, the BBP head itself would already solve it.

Thus carry-free synchronization is not an independent hidden gain.

---

# 15. Smooth-factor / lane-alignment route: stop boundary

A large divisor

\[
D\mid16^m-1
\]

made from small primes or BBP denominator lanes gives a potentially large algebraic modulus.

However Theorem 9.1 shows that

\[
D\mid B_{s,m}
\]

still requires a phase hit in a target of total width

\[
1/D.
\]

Therefore:

\[
\boxed{
\text{large smooth part of }16^m-1
\text{ does not by itself feed the process gcd.}
}
\]

External results on the \(Km\)-smooth part of \(a^m-1\) further indicate that exponentially large small-prime support is exceptional for most \(m\), but even exceptional factorization does not remove the phase condition.

Thus:

\[
\boxed{
\text{smoothness + BBP lane alignment alone is DEAD as an independent feeder.}
}
\]

---

# 16. Full BBP carry gcd also collapses back to an ordinary P2848 witness

Let

\[
G=\frac{16^m-1}{R_s},
\qquad
d:=\gcd(C_{s,m},G),
\]

and write

\[
C_{s,m}=dc,
\qquad
G=dq.
\]

The exact relation is

\[
E_{s+m}-E_s
=
(16^m-1)T_s-C_{s,m}.
\]

Because

\[
16^m-1=R_sG,
\]

we get

\[
E_{s+m}-E_s
=
d(qR_sT_s-c).
\]

Using

\[
R_sT_s
=
16^sR_s\pi-A_s,
\]

we obtain

\[
E_{s+m}-E_s
=
d\left(
16^sR_sq\,\pi-(qA_s+c)
\right).
\]

Define

\[
\boxed{
n:=16^sR_sq
=
\frac{
16^s(16^m-1)
}{d}.
}
\tag{16.1}
\]

Then for an integer \(k=qA_s+c\),

\[
\boxed{
E_{s+m}-E_s
=
d(n\pi-k).
}
\tag{16.2}
\]

If \(d\ge2\), then \(|E_{s+m}-E_s|<1\), so \(k\) is the nearest integer and

\[
\boxed{
n\|n\pi\|
=
\frac{
16^s(16^m-1)|E_{s+m}-E_s|
}{
d^2
}.
}
\tag{16.3}
\]

Thus a successful BBP carry gcd is already an ordinary P2848 witness in disguise.

### Consequence

A carry-gcd attack remains meaningful only if one proves an independent \(\pi\)-specific divisibility theorem for \(d\). Merely renaming the gcd does not create new arithmetic input.

---

# 17. Legendre gate: every successful process shadow eventually becomes a convergent

Let \(a/q\) be a reduced process-generated rational shadow.

If

\[
q^2\left|\pi-\frac aq\right|<\frac12,
\]

then

\[
\left|\pi-\frac aq\right|<\frac1{2q^2}.
\]

By Legendre's theorem,

\[
\boxed{
a/q
\text{ is a continued-fraction convergent of }\pi.
}
\tag{17.1}
\]

Therefore any process-shadow sequence satisfying

\[
q_j^2|\pi-a_j/q_j|\to0
\]

must eventually lie among the true convergents of \(\pi\).

So a successful process theory cannot avoid continued fractions at the critical scale.

---

# 18. Exact critical score on a convergent

Let

\[
\frac{p_k}{q_k}
\]

be a convergent of \(\pi\), and let

\[
\alpha_{k+1}
=
[a_{k+1};a_{k+2},\ldots]
\]

be the continued-fraction tail.

The standard exact formula is

\[
\boxed{
q_k^2
\left|
\pi-\frac{p_k}{q_k}
\right|
=
\frac1{
\alpha_{k+1}+q_{k-1}/q_k
}.
}
\tag{18.1}
\]

Hence

\[
q_k^2
\left|
\pi-\frac{p_k}{q_k}
\right|
\to0
\]

along a subsequence if and only if

\[
a_{k+1}\to\infty
\]

along a subsequence.

Thus a successful process mechanism must eventually manufacture convergents followed by very large partial quotients.

---

# 19. Canonical base-16 complexity of a rational

Let

\[
\frac aq
\]

be a reduced rational.

Write

\[
q=2^uC,
\qquad
C\text{ odd}.
\]

The minimal hexadecimal preperiod length is

\[
\boxed{
\ell
=
\left\lceil\frac{u}{4}\right\rceil.
}
\tag{19.1}
\]

The odd denominator \(C\) is coprime to \(16\), and the minimal repeating period is

\[
\boxed{
r
=
\operatorname{ord}_C(16).
}
\tag{19.2}
\]

Thus

\[
\ell+r
\]

is the canonical "one preperiod + one full period" hexadecimal complexity of the rational.

---

# 20. Critical decimal / hexadecimal compression theorem

Call a rational \(R=a/q\) a **self-generated base-16 process shadow** of \(\pi\) if its canonical hexadecimal expansion

- has preperiod length \(\ell\),
- period length \(r\),

and its first \(\ell+r\) fractional hexadecimal digits agree with the first \(\ell+r\) hexadecimal digits of \(\pi\).

Then

\[
\boxed{
\left|\pi-\frac aq\right|
<
16^{-(\ell+r)}.
}
\tag{20.1}
\]

Define the compression surplus

\[
\boxed{
\Omega(R)
:=
\ell+r-2\log_{16}q.
}
\tag{20.2}
\]

Then

\[
q^2
\left|\pi-\frac aq\right|
<
16^{-\Omega(R)}.
\tag{20.3}
\]

Therefore

\[
\boxed{
\Omega(R_j)\to+\infty
\Longrightarrow
P2848.
}
\tag{20.4}
\]

This is an invariant, reduced-height formulation. It no longer counts fake gcd amplification from nonminimal representations.

---

# 21. Compression implies a large next partial quotient

Assume the self-generated shadow \(R=a/q\) is also a convergent

\[
R=p_k/q_k.
\]

For convergents,

\[
\left|
\pi-\frac{p_k}{q_k}
\right|
>
\frac1{
(a_{k+1}+2)q_k^2
}.
\]

Combining with (20.1),

\[
\frac1{
(a_{k+1}+2)q_k^2
}
<
16^{-(\ell+r)}.
\]

Therefore

\[
\boxed{
a_{k+1}+2
>
\frac{16^{\ell+r}}{q_k^2}
=
16^{\Omega_k},
}
\tag{21.1}
\]

where

\[
\boxed{
\Omega_k
=
\ell_k+r_k-2\log_{16}q_k.
}
\tag{21.2}
\]

Thus:

\[
\boxed{
\text{strong self-generated radix compression}
\Longrightarrow
\text{large continued-fraction partial quotient}.
}
\tag{21.3}
\]

---

# 22. IMPORTANT CORRECTION: the converse requires a digit-boundary margin

A previous exploratory statement said that

\[
a_{k+1}>16^{\Omega_k}
\]

would by itself force the convergent to share the first \(\ell_k+r_k\) hexadecimal digits with \(\pi\).

That implication is **not valid without an additional boundary condition**.

Two numbers can be closer than \(16^{-N}\) while lying on opposite sides of a base-16 cylinder boundary.

The correct sufficient criterion is as follows.

Let

\[
N:=\ell_k+r_k.
\]

Let \(I_N(\pi)\) be the unique length-\(16^{-N}\) base-16 cylinder containing \(\pi\), and define the boundary margin

\[
\eta_N(\pi)
:=
\operatorname{dist}
\big(
\pi,\partial I_N(\pi)
\big).
\]

If

\[
\boxed{
\left|
\pi-\frac{p_k}{q_k}
\right|
<
\eta_N(\pi),
}
\tag{22.1}
\]

then \(p_k/q_k\) has the same first \(N\) hexadecimal digits as \(\pi\).

Using the convergent upper bound

\[
\left|
\pi-\frac{p_k}{q_k}
\right|
<
\frac1{a_{k+1}q_k^2},
\]

a sufficient condition is

\[
\boxed{
a_{k+1}
>
\frac1{
q_k^2\,\eta_{\ell_k+r_k}(\pi)
}.
}
\tag{22.2}
\]

Thus the surviving CF × radix-order direction must control not only \(r_k\), but also the base-16 cylinder margin of \(\pi\).

---

# 23. Exact arithmetic clues

These are navigation clues, not proofs of an infinite mechanism.

## 23.1 The convergent \(355/113\)

We have

\[
113=8\cdot14+1,
\]

so \(113\) lies exactly in one BBP denominator lane.

Also,

\[
\boxed{
16^7\equiv1\pmod{113}.
}
\]

Indeed

\[
16^7-1
=
113\cdot2375535.
\]

Thus

\[
\operatorname{ord}_{113}(16)=7.
\]

Its hexadecimal expansion is

\[
\boxed{
\frac{355}{113}
=
(3.\overline{243F6F0})_{16}.
}
\tag{23.1}
\]

Meanwhile

\[
\pi
=
(3.243F6A8885A3\ldots)_{16}.
\]

They share the first five hexadecimal fractional digits

\[
243F6.
\]

The rational has \(\ell=0,r=7\), so

\[
\Omega
=
7-2\log_{16}(113)
\approx
3.5899105.
\]

Hence

\[
16^\Omega
\approx21022.43.
\]

The actual next partial quotient after \(355/113\) is

\[
292.
\]

So \(355/113\) is a striking three-way alignment

\[
\text{CF convergent}
+
\text{BBP lane}
+
\text{short base-16 order},
\]

but it is still far from the full critical compression threshold.

---

## 23.2 The convergent \(104348/33215\)

We have

\[
33215=5\cdot7\cdot13\cdot73
\]

and

\[
\boxed{
33215\mid16^9-1.
}
\]

In fact

\[
\operatorname{ord}_{33215}(16)=9.
\]

Its hexadecimal expansion begins

\[
\frac{104348}{33215}
=
(3.\overline{243F6A89F})_{16}.
\]

It agrees with \(\pi\) for the first seven hexadecimal fractional digits before diverging.

Here

\[
\ell=0,
\qquad
r=9,
\]

and

\[
\Omega
=
9-2\log_{16}(33215)
\approx1.49022635.
\]

Thus

\[
16^\Omega
\approx62.289.
\]

But the next partial quotient is only

\[
1.
\]

Again the arithmetic order is unusually short, but the CF gain does not synchronize with it.

---

## 23.3 A genuine self-generated semiconvergent

Inside the large \(292\) continued-fraction corridor,

\[
\boxed{
\frac{42223}{13440}
=
\frac{118\cdot355+333}{118\cdot113+106}.
}
\]

So it is an intermediate / semiconvergent.

Its hexadecimal expansion is

\[
\boxed{
\frac{42223}{13440}
=
(3.24\overline{3F6})_{16}.
}
\tag{23.2}
\]

Since

\[
\pi
=
(3.243F6\ldots)_{16},
\]

this rational is genuinely self-generated by the first five hexadecimal fractional digits:

\[
\ell=2,\qquad r=3.
\]

However,

\[
\boxed{
13440^2
\left|
\pi-\frac{42223}{13440}
\right|
\approx70.7514773,
}
\]

so it is nowhere near a P2848 critical hit.

### Interpretation

Self-generated eventually-periodic shadows really do intersect the CF semiconvergent structure. The missing issue is quantitative strength, not existence of the objects.

---

# 24. Current stop boundaries from the BBP round

The following routes are now frozen unless genuinely new input appears.

### STOP 1 — Raw \(\Xi\) language

\[
\Xi_{s,m}
\]

is an ordinary critical rational score in disguise.

### STOP 2 — Fake higher gcd gain

Nonminimal periodic representations move factors among \(d,h,\Delta\); only the fully reduced rational score matters.

### STOP 3 — Direct BBP truncation

Its critical score diverges exponentially.

### STOP 4 — Single BBP modular return

\[
\operatorname{ord}_r(16)\mid m
\]

does not imply

\[
r\mid B_{s,m}.
\]

The \(r=17\) example is an exact counterexample.

### STOP 5 — Full BBP-head synchronization

The head denominator obeys

\[
\log R_s\ge(10/3+o(1))s,
\]

so synchronization requires \(m\gtrsim1.202s\), and the tail carry becomes exponentially large.

### STOP 6 — Carry-free synchronization

If strong enough, it would imply that the direct BBP head already solves P2848.

### STOP 7 — Smooth/cyclotomic abundance alone

A divisor \(D\mid16^m-1\) helps only if a separate phase condition of width \(1/D\) is met.

### STOP 8 — Carry gcd without a new theorem

A successful carry gcd is exactly another \(n\|n\pi\|\) witness. Without an independent divisibility theorem it is only a reformulation.

### STOP 9 — Ordinary equidistribution / normality

Even perfect fixed-scale equidistribution does not supply the required \(1/q\)-scale shrinking target.

---

# 25. CURRENT ACTIVE BOSS

The most promising surviving formulation is now:

\[
\boxed{
\textbf{continued fractions}
\times
\textbf{base-16 multiplicative order}
\times
\textbf{digit-cylinder depth}.
}
\]

Let

\[
\frac{p_k}{q_k}
\]

be the convergents of \(\pi\).

Write

\[
q_k=2^{u_k}C_k,
\qquad
C_k\text{ odd}.
\]

Define

\[
\boxed{
\ell_k
=
\left\lceil\frac{u_k}{4}\right\rceil,
}
\]

\[
\boxed{
r_k
=
\operatorname{ord}_{C_k}(16),
}
\]

and the radix complexity

\[
\boxed{
K_k:=\ell_k+r_k.
}
\]

Define the CF approximation depth

\[
\boxed{
N_k
:=
-\log_{16}
\left|
\pi-\frac{p_k}{q_k}
\right|.
}
\]

If

\[
N_k\ge K_k,
\]

and the digit-boundary issue is satisfied, then the convergent's complete canonical preperiod + one period is visible inside the true hexadecimal process of \(\pi\).

Define the compression surplus

\[
\boxed{
\Omega_k
=
K_k-2\log_{16}q_k.
}
\tag{25.1}
\]

For a genuinely self-generated convergent,

\[
\boxed{
a_{k+1}+2>16^{\Omega_k}.
}
\tag{25.2}
\]

And if a sequence of genuinely self-generated shadows satisfies

\[
\boxed{
\Omega_k\to+\infty,
}
\tag{25.3}
\]

then

\[
q_k^2
\left|
\pi-\frac{p_k}{q_k}
\right|
\to0
\]

and P2848 follows.

---

# 26. The next exact research questions

The next round should NOT return to ordinary prime ideals or generic BBP modular returns.

The sharp questions are:

## A. Radix-order question

Can one prove for infinitely many \(\pi\)-convergent denominators \(q_k\) that

\[
\operatorname{ord}_{C_k}(16)
\]

is unusually small relative to

\[
\log q_k?
\]

In particular, can one obtain an infinite subsequence with

\[
\ell_k+r_k
-
2\log_{16}q_k
\to+\infty
\]

or any nontrivial positive surplus?

## B. Synchronization question

Can short multiplicative order be synchronized with a large next partial quotient

\[
a_{k+1}?
\]

The examples \(355/113\) and \(104348/33215\) show that each phenomenon can appear strongly without the other being strong enough.

## C. Cylinder-margin question

For

\[
N=\ell_k+r_k,
\]

can one control

\[
\eta_N(\pi)
=
\operatorname{dist}
(\pi,\partial I_N(\pi))
\]

from the BBP / process structure strongly enough that a CF error bound forces actual digit matching?

This is required for the converse direction from "very close rational" to "self-generated periodic shadow".

## D. Semiconvergent question

The exact self-generated example

\[
42223/13440=(3.24\overline{3F6})_{16}
\]

lies inside the \(292\)-corridor.

Can large partial-quotient corridors systematically generate low-radix-complexity semiconvergents?

If so, can one obtain a quantitative critical score improvement rather than merely a structured rational?

## E. Cross-base question

Base \(16\) was chosen because BBP is \(\pi\)-specific.

Could another base or another digit-extraction / hypergeometric identity give a stronger relation between

\[
\operatorname{ord}_{C}(b)
\]

and \(\pi\)'s convergent denominators?

Any such change must use a real \(\pi\)-specific theorem, not merely re-encode the same arbitrary rational approximation.

---

# 27. One-line handoff

```text
P2848 CURRENT ACTIVE DIRECTION

Do not attack ordinary primality again.
Do not treat BBP modular return as a digit-block divisor.
Do not count nonminimal cyclotomic gcd factors as independent gain.

The surviving object is the canonical reduced radix complexity of
continued-fraction convergents (and possibly semiconvergents).

For a pi convergent p_k/q_k, write

    q_k = 2^{u_k} C_k,  C_k odd,

    ell_k = ceil(u_k/4),

    r_k = ord_{C_k}(16),

    Omega_k = ell_k + r_k - 2 log_16(q_k).

A genuinely self-generated base-16 periodic shadow satisfies

    q_k^2 |pi-p_k/q_k| < 16^{-Omega_k},

and if it is a convergent then

    a_{k+1}+2 > 16^{Omega_k}.

Thus an infinite sequence of self-generated shadows with

    Omega_k -> +infinity

would prove P2848.

The active research problem is to derive, from pi-specific structure,
an infinite synchronization among:

    (1) unusually small base-16 multiplicative order r_k,
    (2) large continued-fraction next coefficient a_{k+1},
    (3) sufficient hexadecimal cylinder margin / digit agreement.

Known exact clues:
    355/113: q=113 is a BBP 8k+1 denominator, ord_113(16)=7,
             and it shares the hex prefix 243F6 with pi.

    104348/33215: ord_33215(16)=9 and it shares seven hex digits with pi.

    42223/13440 = (3.24 overline{3F6})_16 is an exact self-generated
                   semiconvergent, but its critical score is about 70.75.

P2848 remains unproved.
```

---

# 28. Current verdict

\[
\boxed{\text{P2848 NOT PROVED}.}
\]

The ordinary "micro-prime" idea did not survive as ordinary prime-element theory.

What survived after the entire reduction chain is more precise:

\[
\boxed{
\text{process DNA}
\to
\text{canonical radix compression}
\to
\text{continued-fraction gate}
\to
\text{large partial quotient}.
}
\]

The BBP/base-16 attack did not supply the final theorem, but it exposed a concrete \(\pi\)-specific arithmetic interface:

\[
\boxed{
\text{CF denominator}
\longleftrightarrow
\operatorname{ord}_{C}(16)
\longleftrightarrow
\text{hex digit process}.
}
\]

The next attack should be concentrated there.
