# Sources, scope, and dependency record

Parent checkpoint: 9ac9c9065a038a87a6ab587b8d87c65f689f393b. The new work fixes the same U_n and approximation family; it neither modifies the meaning of a successful subsequence nor derives any global bad-approximability claim for pi.

## Existing repository inputs

- [Initial-value sieve and active pool](../pi-initial-value-sieve-2026-11-02/README.md): the new capacity target supplied by the user, exact finite search scope, and W(n)<=theta_Z(n).
- [Active-prime capacity](../pi-initial-value-sieve-2026-11-02/active-primes/ACTIVE_PRIME_CAPACITY.md) and its existing CSV: local classification and the finite table to 20000. CSV Git blob: f9c06a0b04955a0328d5c5a25616b07135d04629. The new summary does not generate additional primes.
- [Multi-band deficit certificate](../pi-cross-prime-budget-2026-11-02/multi-band-budget/MULTI_BAND_DEFICIT_CERTIFICATE.md): ln Q_n >= (2+ln2)n-4W(n)+o(n).
- [Earlier midpoint interpretation](../pi-cross-prime-budget-2026-11-02/midband-arithmetic/EDGE_CAPACITY_AND_MIDPOINT.md): midpoint coefficient, trace congruence and the weaker necessary residue restriction p=7 mod8. The new deduction is the non-CM certificate and logarithmic zero-density removal.

## Primary analytic input

J.-P. Serre, *Quelques applications du theoreme de densite de Chebotarev*, Publications Mathematiques de l'IHES 54 (1981), pp.123-201, Section8.2, Theorem20 and Corollary2, printed p.189. https://www.numdam.org/item/PMIHES_1981__54__123_0.pdf . For a non-CM elliptic curve over Q, the trace-zero primes have count o(x/ln x), unconditionally. This implication is used only for the fixed midpoint curve E. It supplies no bound on z_p at other degrees.

## Non-CM certificate support

- Andrew Sutherland, MIT18.783, Spring2019 Lecture22, reduction injection, equation(3). https://math.mit.edu/classes/18.783/2019/LectureNotes22.pdf . Good reduction injects endomorphisms; pass to a finite field of definition of the geometric endomorphisms before specializing.
- Andrew Sutherland, MIT18.783, Fall2025 Lecture13, ordinary/supersingular classification, Theorem13.6 and Corollary13.20. https://math.mit.edu/classes/18.783/2025/LectureNotes13.pdf . For an ordinary elliptic curve over a finite field, the geometric endomorphism algebra is the imaginary quadratic Frobenius field. This is why the two distinct fields at p=3,5 certify no geometric CM.
- Denis Xavier Charles, *Complex Multiplication Tests for Elliptic Curves*, 2004, Sections5.1 and6.2. https://arxiv.org/pdf/math/0409501 . Theorem5.8 is stated with order/conductor conditions; our proof uses the good-reduction injection on rational endomorphism algebras, avoiding an unwarranted assertion that the full integral orders agree at every prime. Corollary6.2 also sketches the non-CM trace-zero density argument.

## Local congruences and heuristic boundary

- Tony D. Noe, *On the Divisibility of Generalized Central Trinomial Coefficients*, JIS9 (2006),06.2.7. https://cs.uwaterloo.ca/journals/JIS/VOL9/Noe/noe35.pdf . Exact parameters(a,b,c)=(1,2,2), discriminant-4; Lucas and Holt congruences only.
- Nadav Kohen, *Density and Symmetry in the Generalized Motzkin Numbers Modulo p*, Integers26(2026),A22. https://math.colgate.edu/~integers/aa22/aa22.pdf . Remark1, p.7, gives a heuristic inactive-prime proportion exp(-1/2) for the usual central trinomial sequence. The inference of a comparable scale for our U_n is explicitly heuristic. Proposition4 covers our integral Laurent polynomial; the main fixed-p density theorem does not settle the varying-prime pool.
- Joel A. Henningsen and Armin Straub, *Generalized Lucas congruences and linear p-schemes*, Corollary3.1. https://arminstraub.com/downloads/pub/generalized-lucascongruences.pdf . Direct integer Laurent-polynomial formulation of the local digit rule.

Further source checks and access limitations are recorded in prime-pool-literature/PRIMARY_SOURCE_SCOPE.md. References involving varying parameters, fixed primes, or fixed polynomial degrees are not substituted for a bound on our fixed-parameter, varying-prime pool. The sufficient first-moment estimate is not supplied by any theorem used here.
