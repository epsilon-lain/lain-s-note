# Prime-pool capacity: primary-source scope check

2026-11-02. The fixed sequence here is

\[
U_n=[x^n](1+2x+2x^2)^n=\operatorname{ct}[(x^{-1}+2+2x)^n].
\]

The target concerns the varying-prime set
\(\mathcal Z=\{p\text{ odd}:\exists r<p,\ U_r\equiv0\pmod p\}\), not the density of zero indices after fixing one prime. No theorem establishing the needed uniform prime-pool upper bound was found in the primary sources checked below. This is a report of the checked sources, not a claim that all literature has been searched.

## 1. Kohen 2026 already identifies a heuristic of almost exactly the observed scale

Nadav Kohen, *Density and Symmetry in the Generalized Motzkin Numbers Modulo p*, Integers 26 (2026), #A22. https://math.colgate.edu/~integers/aa22/aa22.pdf . DOI 10.5281/zenodo.18305088.

Remark 1, page 7, discusses the **usual** central trinomial sequence. The observed inactive-prime fraction below one million is about 0.6075. Reflection halves the distinct zero tests. Pretending these tests are independent with zero probability 1/p gives heuristic inactive proportion exp(-1/2), hence active proportion 1-exp(-1/2), about 0.39346934. This is not an independence theorem or density bound.

The current U sequence's observed 0.399696 is compatible with that heuristic. **This is an inference, not an identification or density proof.**

Theorem 1 proves reflection. Proposition 4 proves digit multiplication for every integral Laurent polynomial of degree one, covering x^-1+2+2x. The main density theorem fixes p and counts Motzkin indices; it does not bound Z.

Parameter caution: symmetric T_n^{a,b} would require a^2=2,b=2 here, not integer pair (2,2). The general Laurent formulation avoids this ambiguity.

## 2. Noe 2006 gives exactly the local input, without a prime-density theorem

Tony D. Noe, *On the Divisibility of Generalized Central Trinomial Coefficients*, Journal of Integer Sequences 9 (2006), Article 06.2.7. https://cs.uwaterloo.ca/journals/JIS/VOL9/Noe/noe35.pdf .

His parameters (a,b,c)=(1,2,2), discriminant d=b^2-4ac=-4, give our sequence. Formula (13) and Theorems 8.5–8.6 supply digitwise divisibility and reflection. Theorems 8.7–8.9 evaluate a few endpoint terms modulo p. These are fixed-prime statements. They certify finite classification of every p, but do not estimate the number or weight of active primes as p varies.

## 3. General Lucas formulation removes parameter ambiguity

Joel A. Henningsen and Armin Straub, *Generalized Lucas congruences and linear p-schemes*, author manuscript dated 2021-11-16, published JCTA 2022. https://arminstraub.com/downloads/pub/generalized-lucascongruences.pdf .

Corollary 3.1 states Lucas congruences for ct[(a x^-1+b+c x)^n] with integer a,b,c. Lemma 3.2 and its discussion recall the discriminant endpoint and reflection congruences. Thus the implication p|U_n => p in Z is directly covered for (a,b,c)=(1,2,2). These results do not assert a cross-prime density.

## 4. Fixed-prime recurrence classifications have the same limitation

Nadav Kohen, *Uniform Recurrence in the Motzkin Numbers and Related Sequences mod p*, Electronic Journal of Combinatorics 32(2) (2025), #P2.44. https://www.combinatorics.org/ojs/index.php/eljc/article/download/v32i2p44/pdf/ . DOI 10.37236/13089.

For each fixed p, the generalized central trinomial sequence has digit multiplication. If one zero digit exists, the set of zero indices has density one; otherwise the sequence is uniformly recurrent. Theorem 11 and Corollaries 13–14 classify this dichotomy. The dichotomy answers what an active prime does after it is chosen, not how frequently active primes occur among primes.

Narad Rampersad and Jeffrey Shallit, *Congruence Properties of Combinatorial Sequences via Walnut and the Rowland–Yassawi–Zeilberger Automaton*, EJC 29(3) (2022), #P3.36. https://www.combinatorics.org/ojs/index.php/eljc/article/download/v29i3p36/pdf/ . DOI 10.37236/10799. Its automata and decidability machinery are likewise for fixed prime powers; they do not turn finitely many prime classifications into an eventual weighted capacity bound.

## 5. Other tempting results do not answer this axis of variation

John Cullinan and Farshid Hajir, *On the Galois groups of Legendre polynomials*, author PDF https://faculty.bard.edu/cullinan/papers/legendre.pdf . The main problem concerns factorization and Galois groups of a fixed-degree polynomial while its argument varies. The present pool fixes the argument i and permits all degrees r<p to vary with p. No uniform theorem converting the former problem into the latter is supplied there; some broad Galois assertions in the paper are conjectural.

The 2026 paper by Alekseyev, Amdeberhan, Shallit and Vukusic, *On the p-adic valuations of values of Legendre polynomials*, DOI 10.1007/s40993-026-00743-0, concerns P_n evaluated **at the prime p** according to its official abstract. That is a different specialization from P_n(i). Full publisher text was inaccessible in the retrieval used here, so no theorem from it has been relied on.

Kohen's dissertation was located via the author's publications page, but its ProQuest URL returned HTTP 403. It has not been treated as read.

## 6. What would actually close the cited route

Reflection alone gives pairs of non-midpoint zeros; it does not make different zero events independent. A rigorous first-moment or sieve estimate must control those events uniformly while the degree range grows proportionally to p. A theorem for each fixed degree or fixed prime does not justify exchanging those ranges with the prime limit.

Therefore a proof of theta_Z(x) <= (0.6732867951...-eta)x would be additional mathematical input beyond the checked digit and symmetry theorems. The observed near-0.4 scale has a published heuristic explanation, but its rigorous upper bound remains the specific missing step.
