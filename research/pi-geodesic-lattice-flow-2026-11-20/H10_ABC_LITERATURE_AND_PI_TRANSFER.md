# Hilbert's tenth problem, the abc conjecture, and possible tools for the pi lattice-flow project

Date: 2026-11-20  
Scope: literature review and research transfer assessment. This note distinguishes established theorems, fresh preprints, conditional results, and speculative links to our bad-approximability project.

## 1. Executive conclusion

The H10 literature is genuinely useful for our elliptic-curve/formal-group branch, but it does **not** directly imply that \(\pi\) is badly approximable.

The most transferable pattern is:
\[
\text{local valuation control}
\longrightarrow
\text{many-prime divisibility}
\longrightarrow
\text{global height inequality}
\longrightarrow
\text{contradiction}.
\]
Pasten–Vidaux (2026) implements this pattern using an elliptic formal-group approximation lemma and an \(abc\)/Vojta height inequality. Our current construction has the local formal-group ingredient, but lacks the global bridge that converts a very small real error into sufficiently many distinct prime divisors or a controlled radical.

A second, more directly relevant tool is the theory of **linear forms in elliptic logarithms**. Our real recurrence \(V=[s]P\to O\) is exactly a small-divisor problem for the elliptic logarithm of \(P\). David–Hirata-Kohno give effective lower bounds for such linear forms. This is a natural way to quantify the height cost of making the formal parameter small.

## 2. What is known about H10 as of 20 November 2026

It is essential to distinguish three rings/fields.

### 2.1 The integers \(\mathbb Z\)

The Davis–Putnam–Robinson–Matiyasevich theorem gives a negative solution: no algorithm decides solvability of arbitrary polynomial equations over \(\mathbb Z\).

### 2.2 The rational field \(\mathbb Q\)

A major new preprint is:

- Ashay Burungale, *A Diophantine definition of \(\mathbb Z\) in \(\mathbb Q\)*, arXiv:2611.16130 (12 Nov 2026).

It claims an unconditional Diophantine definition of \(\mathbb Z\) in \(\mathbb Q\), with a polynomial of total degree at most four. The abstract describes quadratic twists of the congruent-number elliptic curve, analytic rank one via the Tian–Yuan–Zhang criterion (in Smith's formulation), a height estimate on a quotient of a Shimura curve, and Green–Tao–Ziegler prime patterns. It explicitly says it does not use the 2-converse or the 2-adic Fontaine–Mazur input in the OpenAI preprint.

Pasten–Vidaux, in their 19 Nov 2026 preprint, state that in their opinion Burungale's version is complete and that this yields a negative answer to H10 over \(\mathbb Q\). This is strong contemporary expert endorsement, but the result is still a very recent preprint, not a long-settled peer-reviewed theorem. The careful status is: **a claimed unconditional solution with substantial independent expert confidence, still deserving normal proof scrutiny**.

Related items:
- OpenAI, *Hilbert's tenth problem over the rational numbers* (6 Oct 2026): claims a negative answer, but Pasten–Vidaux say its special 2-adic inputs were yet to be verified.
- Alexander Kozachinskiy, Hector Pasten, Xavier Vidaux, *Why should \(\mathbb Z\) be Diophantine in \(\mathbb Q\)?*, arXiv:2611.17124 (12 Nov 2026): an accessible conditional route assuming (X) finiteness of the 2-primary parts of Sha groups, (Y) abc, and (Z) Dickson's conjecture.

### 2.3 A general number field \(K\)

Hector Pasten and Xavier Vidaux, *On Hilbert's tenth problem for number fields*, arXiv:2611.24748 (18/19 Nov 2026), prove:
\[
\left.
\begin{array}{l}
abc\text{ over }K,\\
0<\operatorname{rank}E(\mathbb Q)
=\operatorname{rank}E(K)
\text{ for some }E/\mathbb Q
\end{array}
\right\}
\Longrightarrow
\mathbb Z\text{ is Diophantine in }K.
\]
For quadratic number fields, the required rank-stable elliptic curve is supplied by Zywina's work, so the remaining hypothesis is \(abc\) over \(K\).

Do not confuse this with H10 for the ring of integers \(\mathcal O_K\). The latter is a different problem. Alpöge–Bhargava–Ho–Shnidman and collaborators proved rank stability for abelian varieties in quadratic extensions, yielding an unconditional negative answer for H10 over \(\mathcal O_K\) for every number field. That does not by itself solve H10 over the whole field \(K\).

## 3. The most relevant mechanism in Pasten–Vidaux

Their proof has four layers.

### 3.1 The Pheidas predicate

They use a Diophantine set \(\Phi\subset K\) with:
- all rational numbers with square denominator belong to \(\Phi\);
- outside a finite set of primes, elements of \(\Phi\) cannot have negative odd valuation.

This turns parity of local valuations into a Diophantine condition.

### 3.2 Elliptic formal-group local approximation

Their Lemma 4.1 says: for integers \(k\) and \(M\ge2\), there is a rational elliptic-curve quantity \(\beta\) such that
\[
\beta-k\in M\mathbb Z_p\qquad(p\mid M).
\]
The proof chooses a non-torsion point \(P\), multiplies it into sufficiently small formal-group subgroups at all \(p\mid M\), and uses
\[
[k](z)=kz+z^2F_k(z),\qquad F_k\in\mathbb Z[[z]].
\]
This is essentially the same local formal-group mechanism we have been studying. Crucially, the lemma imposes **local congruences**, not a real approximation estimate and not a denominator-height bound.

### 3.3 Odd valuations propagate to conjugate differences

The Diophantine predicate forces certain rational witnesses to be close to an algebraic element at every prime where a chosen polynomial has positive odd valuation. Applying this to all primes over a rational prime makes those primes divide a conjugate difference
\[
b-\sigma(b).
\]
Thus a local parity condition becomes divisibility by many distinct rational primes.

### 3.4 \(abc\)/Vojta converts height into a prime-counting bound

The number-field \(abc\) conjecture is formulated through a Vojta-type inequality on \(\mathbb P^1\) with truncated counting functions. The paper only needs the corresponding Vojta inequality for one fixed squarefree polynomial \(G\) of degree at least seven. It gives a lower bound for the number of primes at which \(G(b)\) has positive odd valuation. The conjugate-difference argument gives an upper bound for that count in terms of the height of \(b-\sigma(b)\). Comparing the two bounds forces \(b=\sigma(b)\), hence rationality.

This is the precise local-to-global pattern worth borrowing.

## 4. What abc says, and its current status

For coprime positive integers \(a,b,c\) with \(a+b=c\), define
\[
\operatorname{rad}(abc)=\prod_{p\mid abc}p.
\]
The abc conjecture says that for every \(\varepsilon>0\), there is a constant \(C_\varepsilon\) such that
\[
c\le C_\varepsilon\operatorname{rad}(abc)^{1+\varepsilon}.
\]
Informally: an integer cannot be enormously larger than the product of the distinct primes dividing it, except by a small power.

Over number fields there is an analogous height-versus-truncated-prime-counting statement. In the Pasten–Vidaux paper, this is identified with the Vojta conjecture on \(\mathbb P^1\) for the divisor \([0]+[1]+[\infty]\), and Belyi maps transfer it to other divisors.

**Status:** abc is not currently an accepted theorem. The IUT papers of Mochizuki claim a proof, but there is no broad mathematical consensus that the proof has been completed. Project LANA's official July 2026 interim report says it is still suspending judgment, with specific unresolved formalization/interpretation issues around IUT-III, Theorem 3.11 to Corollary 3.12. Thus we must treat abc as a conjectural input, not a proved tool.

## 5. Papers most worth reading for our project

### Priority A — directly transferable

1. **Pasten–Vidaux (2026), arXiv:2611.24748.**  
   Read Lemma 4.1 and Sections 2, 5, 6. This is the closest template for combining formal groups, local valuations, and global heights.

2. **Kühn–Müller (2015), *A height inequality for rational points on elliptic curves implied by the abc-conjecture*.**  
   Under uniform abc, they derive a strong inequality relating the coordinates and denominator of a rational elliptic point to the radical of its denominator:
   \[
   \max\{\tfrac12\log|a_P|,\log d_P\}
   \le(1+\varepsilon)\log\operatorname{rad}(d_P)+O_\varepsilon(1).
   \]
   This is the most directly relevant abc-to-elliptic-height result in the list. It controls the support of denominators, not the real closeness of \(\xi(P)\) to zero.

3. **David–Hirata-Kohno (2009), *Linear forms in elliptic logarithms*.**  
   They prove effective lower bounds for linear forms in elliptic logarithms and solve a conjecture of Lang in the elliptic case. This directly addresses the real recurrence \( [s]P\to O \) used in our strengthened formal-group construction.

4. **Stroeker–Tzanakis (1994), *Solving elliptic Diophantine equations by estimating linear forms in elliptic logarithms*.**  
   A practical demonstration of the method: use canonical heights for an upper bound on an elliptic-logarithm linear form, David's theorem for a lower bound, then bound the integer coefficients in the Mordell–Weil group.

### Priority B — understand the H10 architecture

5. **Burungale (2026), arXiv:2611.16130.**  
   Important new claimed solution for \(\mathbb Q\). Read its height estimate on the Shimura-curve quotient and the construction of the Diophantine model.

6. **Kozachinskiy–Pasten–Vidaux (2026), arXiv:2611.17124.**  
   A useful map of the conditional route through abc, Dickson, and Sha finiteness.

7. **Poonen (2002), *Using elliptic curves of rank one towards the undecidability of Hilbert's Tenth Problem over rings of algebraic integers*.**  
   The classical rank-stability/Diophantine-definition mechanism underlying the later work.

8. **Anscombe et al. (2023), *A Survey of Local–Global Methods for Hilbert's Tenth Problem*.**  
   Background on valuations, quadratic forms, central simple algebras, and definability.

## 6. What can actually help our \(\pi\) project?

We need to separate three possible transfers.

### A. Transfer the local-to-global proof architecture — promising

Our formal-group construction already controls local congruences at selected primes. The missing ingredient is a global arithmetic object whose local behavior forces many distinct primes to divide a single nonzero integer or rational numerator.

A candidate research question is:

> Can the high-order formal-group combination \(B_m(V)-1\), together with its prescribed \(p\)-adic congruences, force a lower bound for the radical of its reduced numerator or denominator?

If yes, an abc/Vojta inequality might turn that radical information into a height bound. At present there is no proof that our construction forces enough new prime divisors.

### B. Use abc to control elliptic-point denominators — relevant but insufficient alone

Kühn–Müller show that uniform abc would control the height of an elliptic rational point in terms of the radical of its denominator. But our high-order cancellation uses a **sum of rational functions at several multiples** \([j]V\). A height bound for each individual point does not automatically control the denominator after the weighted sum and cancellation.

We would need a theorem controlling the radical and common factors of the resulting linear combination, not just the coordinates of each \([j]V\).

### C. Quantify the real recurrence using elliptic logarithms — the best immediate experiment

On the identity component,
\[
E(\mathbb R)^0\simeq\mathbb R/\mathbb Z.
\]
If \(u=\phi(P)\), then the formal parameter of \([s]P\) near \(O\) is comparable to the distance
\[
\|su\|=\min_{m\in\mathbb Z}|su-m|.
\]
Thus our choice of \(s\) is an ordinary one-dimensional Diophantine approximation problem for the elliptic logarithm \(u\).

David–Hirata-Kohno provide effective lower bounds for the corresponding elliptic logarithm linear forms. On the other hand, Dirichlet's pigeonhole principle gives recurrent multiples with \(\|su\|\ll 1/N\) for some \(s\le N\).

This gives a quantitative framework for the tension:
- recurrence gives upper bounds on \(|\xi([s]P)|\);
- elliptic-logarithm estimates give lower bounds;
- canonical height grows quadratically in \(s\);
- our order-\(m\) cancellation raises the power of the small formal parameter, but its integer weights have total size \(2^{m+1}-1\), and the crude output-height bound grows like \(O(m^3\widehat h(V)+m^2)\).

The next experiment should compute these quantities for one explicit rank-one curve and a generator \(P\), then compare
\[
-\log|\xi([s]P)|,\quad
\widehat h([s]P),\quad
h(B_m([s]P)),\quad
-\log|B_m([s]P)-1|.
\]
This will tell us whether high-order cancellation can beat its height cost in any meaningful regime.

## 7. The crucial obstruction: abc does not directly see \(\pi\)

The classical abc conjecture concerns algebraic integers and rational/algebraic heights. The number \(\pi\) is transcendental, so it cannot be substituted as a point of a number field in the H10/abc arguments above.

A proof of
\[
\inf_{q\ge1}q\|q\pi\|>0
\]
would need a new bridge from rational approximants \(p/q\) to an algebraic object on which abc or Vojta can act. The bridge must preserve the scale \(q|p-\pi q|\) and must not introduce uncontrolled heights or denominators.

Therefore:
- H10 and abc supply powerful techniques for algebraic local-global arithmetic.
- They do not presently give a lower bound for \(q\|q\pi\|\).
- The most promising concrete import is the elliptic-logarithm height analysis, followed by a search for a radical/divisibility statement for our high-order formal-group combinations.

## 8. Source list

- Pasten–Vidaux, arXiv:2611.24748: https://arxiv.org/abs/2611.24748
- Burungale, arXiv:2611.16130: https://arxiv.org/abs/2611.16130
- Kozachinskiy–Pasten–Vidaux, arXiv:2611.17124: https://arxiv.org/abs/2611.17124
- OpenAI, *Hilbert's tenth problem over the rational numbers* (6 Oct 2026): https://github.com/openai/math/tree/main/preprints/Hilberts-tenth-problem-over-the-rational-numbers-October-6-2026
- Poonen, *Using elliptic curves of rank one...*: https://math.mit.edu/~poonen/papers/ants5.pdf
- Kühn–Müller, arXiv:1210.6543: https://arxiv.org/abs/1210.6543
- David–Hirata-Kohno, *Linear forms in elliptic logarithms*: https://doi.org/10.1515/CRELLE.2009.018
- Stroeker–Tzanakis, *Solving elliptic Diophantine equations...*: https://users.math.uoc.gr/~tzanakis/Papers/Acta94-v2.pdf
- Anscombe et al., arXiv:2309.14987: https://arxiv.org/abs/2309.14987
- Project LANA official July 2026 status: https://zen.ac.jp/news/zmcpostevent0717e
