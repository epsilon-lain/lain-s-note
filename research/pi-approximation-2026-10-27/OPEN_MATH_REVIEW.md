# Independent mathematical review task

Prepared 2026-10-27 for later use with Open Math or another research system.
This is a self-contained review prompt, not a platform-specific API schema.
The package has not been submitted to Open Math and contains no Lean-checked proofs.

## Scope and priority

Audit the assertions below before attempting extensions. Report a statement as
proved, false with a counterexample, or incomplete with the exact missing step.
Do not assume that pi is badly approximable. Do not infer small true approximation
errors from the failure of a lower-bound certificate. Do not infer an infinite
statement from any finite scan in the data files.

### Task A: construction and transferred lower bounds

For integers m>=2 let T_m=floor(m^2*pi^2). Set N_m=T_m unless T_m is a
perfect square, in which case N_m=T_m-1. Define

    F_m=sqrt(N_m)/m,
    g_m=gcd(N_m,m^2), A_m=m^2/g_m, B_m=N_m/g_m,
    c_m=1/[A_m(4+F_m)], eta_m=pi-F_m,
    H(m,q)=max(0,c_m/q^2-eta_m), q a positive integer.

Check all of the following, including boundary cases and unreduced fractions:

1. N_m is positive and nonsquare; 3<F_m<pi<4; eta_m<1/(3m^2).
2. For every integer a and positive integer q,
   |F_m-a/q|>=c_m/q^2 and |pi-a/q|>=H(m,q).
3. With delta_m=A_m*pi^2-B_m>0,

       eta_m=delta_m/[A_m(pi+F_m)],
       c_m/eta_m=(pi+F_m)/[(4+F_m)delta_m].

4. H(m,q)>=C/q^2 is equivalent, for C>0, to

       1/(4+F_m)-q^2 delta_m/(pi+F_m)>=A_m C.

The detailed proofs are in PURE_ROOT_FUNCTIONS.md (updated section 6) and
ROOT_ENVELOPE_OBSTRUCTION.md, sections 1-2.

### Task B: infinite-envelope obstruction (highest priority)

Let x be any real number with 3<x<4 and x^2 irrational. Let R be the set
of coprime positive integer pairs (A,B) such that 3<sqrt(B/A)<4 and
sqrt(B/A) is irrational. For alpha=sqrt(B/A), put

    c(A,B)=1/[A(4+alpha)], eta(A,B)=|x-alpha|,
    U_x(q)=sup_{(A,B) in R} max(0,c(A,B)-q^2 eta(A,B)).

Proposed theorem: U_x(q) is finite for every positive integer q, and
U_x(q)->0 as q->infinity.

Proposed argument: choose K with 1/[7(K+1)]<epsilon. Terms with A>K are
bounded by this number. There are only finitely many pairs with A<=K
because 9A<B<16A; none has alpha=x because x^2 is irrational. Their
positive eta values make every one of their contributions vanish for all
sufficiently large q. Check especially the supremum, the quantifier order,
and whether any hidden uniform convergence assumption is being used.

Apply the result to x=pi, using only the classical fact that pi^2 is
irrational. With G_infinity(q)=sup_m H(m,q), conclude

    q^2 G_infinity(q)->0.

This should rule out only the current polynomial constant followed by
triangle-inequality transfer, even when m may depend arbitrarily on q.
It should not rule out improved constants for the reference roots.

As a diagnostic, verify that x=2+sqrt(2) satisfies this same obstruction
while also satisfying |x-a/q|>=1/(4q^2) for all a,q. This is intended to
prevent an incorrect inference about actual approximation errors.

Detailed proof: ROOT_ENVELOPE_OBSTRUCTION.md, sections 3-4.

### Task C: independent exact-arithmetic verification

Run `python3 reproduce.py` with assertions enabled (do not use `python -O`).
Inspect the interval orientation, root enclosure, flooring, reduction by gcd,
and all claimed q thresholds. Recompute, rather than trust, the JSON files.
In particular verify

    |pi-sqrt(440933025)/6684|<1.3e-13,
    c_6684>9e-9,

and the finite-height certificate at K=200:

    q>=32 => q^2 G_infinity(q)<=1/1407.

The latter concerns an upper bound on the strength of a lower-bound
certificate. It is NOT an upper bound on |pi-a/q|.

### Task D: digit-wave statements, separately

For irrational x in (0,1), audit the equivalence between

    exists q_j->infinity,n_j in Z: q_j^2|x-n_j/q_j|->0

and the following variable-base condition: there are integer bases b_j->infinity
whose first two fractional digits d1,d2 agree and |d3-d1|/b_j->0.
Here b=q+1 and n=d1. Audit the corresponding period-two statement for
b=q-1 in WAVE_RESULTS.md. Do not replace relative third-digit deviation
by a fixed absolute deviation or by equality of all three digits.

Also audit the equivalence between being Liouville and having unbounded
initial identical-digit prefix lengths across integer bases. The Mahler
bound used for pi's length cap is an external theorem, not a numerical
observation. It is independent of Tasks A-B.

### Task E: focused research extension

Only after A-D, investigate whether a stronger rigorously computable
reference constant can evade Task B. For example, for an irrational
alpha consider kappa(alpha)=inf_{q>=1} q*||q alpha||. The polynomial
constant c(A,B) is a lower bound for kappa, not generally its exact value.
Do not assume that kappa tends to zero merely because c tends to zero.

A useful deliverable would be either a new proved bound, a counterexample
to a proposed uniformity statement, or a precise conditional lemma
identifying a new pi-specific input. Ordinary numerical closeness of
roots to pi, or rephrasing bounded partial quotients, is insufficient.

## Requested output

- A theorem-by-theorem verdict with proof gaps or counterexamples.
- Independent exact-computation results and reproducible commands.
- For any formalization, identify theorem dependencies and clearly state
  whether Lean actually ran; do not use `sorry` or an unverified axiom as
  evidence that the target theorem was checked.
- If continuing research, state exactly which new estimate is needed and
  which of the current bounds or definitions it changes.

## Follow-up: weighted slopes, tight prefix certificates, oscillatory control

Audit WEIGHTED_WAVE_RETURNS.md and WAVE_CERTIFICATES_AND_CONTROL.md.
For x=pi-3 in base p, check (p-1)x=d1+sum_j(d[j+1]-d[j])/p^j.
For a prefix A=floor(p^L*x), check that the exact infimum/supremum
of E=(p-1)*||(p-1)x|| over the prefix cell are obtained from the
integer and half-integer breakpoints in the scaled cell closure.
Verify the resolution bound (p-1)^2/p^L and rational-threshold termination
for fixed irrational x when arbitrarily many digits can be certified.

Check the explicit control y=46/99+sum_{j>=1}10^(-j!): all adjacent
decimal digits differ by at least 2; y is Liouville; its mean squared
adjacent difference tends to 4; in bases 99*10^(k!)+1 it has k identical
initial digits for k>=4. This refutes fixed-base amplitude/energy criteria,
not the full signed weighted phase or an all-base criterion.

Run `python3 wave_certificates.py`. Check the finite pi witness E<1/20000
from its three-digit certificate, including the exact denominator and
large digit values in wave_certificate_results.json. It is neither a
new numerical approximation record nor evidence that the liminf is zero.

## Follow-up: real-base lines and integer-endpoint needle areas

Independently audit REAL_BASE_NEEDLE_MODEL.md. For irrational x in (0,1),
let beta_n=1+n/x, and check the greedy repeated digit n and the identity
q^2*|x-n/q|=x*q*|(q+1)-beta_n|. Check all quantifiers and the proof that
liminf_n n*||beta_n|| equals liminf_q q*||qx||.

Define the needle as the rigid segment from (0,0) to (q,a), with q>=1
and a an integer; rotate about the origin to slope pi. Audit the sector
formula and prove that its area is asymptotic to q^2*|pi-a/q|/2 when
a/q approaches pi. Check the liminf statement for nearest numerators and
the equivalence of the area-to-zero property for arbitrary numerators.

For target slope 2+sqrt(2), verify the stated strictly positive global
area bound 1/85, including negative and outside-[3,4] slopes. Recompute
the exact interval showing the pi needle (113,355) has area <1/85.
Do not infer pi is well or badly approximable from that finite comparison.
Run `python3 real_base_needle.py` with assertions enabled.

## Follow-up: infinite rays, tube overlaps, and a Kakeya-inspired bridge

Audit INFINITE_RAYS_KAKEYA_BRIDGE.md independently. Rays have a shared
origin, slope pi or a/q, and are truncated by 0<=t<=q. Check that the
continuous rotation sweeps a triangle of area q^2*|pi-a/q|/2; check that
the full infinite-ray sweep has infinite area for distinct slopes.

For vertical half-width 1/q, verify that each tube has area 2 and their
overlap fraction is Phi(E), E=q^2*|pi-a/q|, where Phi(E)=1-E/4 if
0<=E<=2 and Phi(E)=1/E if E>=2. Check the area of the union of all
intermediate-slope tubes, and note these tubes shear rather than undergo
rigid Euclidean rotation. Verify the pi (113,355) certified interval and
the global bound Phi(E)<13/14 for target 2+sqrt(2).

Audit the dyadic max M_Q and its equivalence to liminf q*||q*pi||=0.
Prove the stated union-bound estimate for the exceptional sets, and
verify it implies only an almost-everywhere shrinking gap, not a fixed
positive gap or a pointwise conclusion for pi. Audit the Gram kernel
and its exact Fourier/Plancherel formula.

Compare—not conflate—the unit delta-tube scaling with section 9.2 of
the OpenAI preprint. P2848-relevant pairs have direction gaps o(delta),
whereas that section's Lemma 9.4 assumes delta-separated directions.
Do not claim the preprint's PDE argument was transferred, or that an
Lp or dimension result alone establishes a pi-specific assertion.
Run `python3 ray_tube_overlap.py` with assertions enabled.

## Follow-up: real-base wave space and constant-digit cells

Audit BASE_SPACE_MODEL.md, including all endpoint inequalities. For real
base beta>1 and target y with fractional part x, define d_j by the greedy
recursion; check that F(beta,y;v,u) is bounded and that the digits plus
the separately recorded integer part reconstruct y.

For n>=1 and beta>n+1, verify that the first K greedy digits equal n
if and only if L_K(beta)<=x<U_K(beta), including the left endpoint and
excluding the right endpoint. Audit the local interval [beta_minus,
beta_plus) and its asymptotics as K tends to infinity for fixed x,n.

For fixed x and K, independently verify width ~ x^(K-2)*n^(1-K) as
n tends to infinity, paying attention to the different limits. Check the
almost-every-offset finite-hit result for K>=3 by the convergent sum of
interval lengths, and verify that it says nothing about the specific
integer-grid offset zero. Check the Liouville example where offset zero
is exceptional.

Review the polynomial cell inequalities and the finite-prefix stability
at the all-equal wave. Verify certified digits in bases 8 and 114 and the
rational bisection intervals. Run `python3 base_space.py` with assertions
enabled. Do not infer a new pi irrationality measure or P2848 result from
a generic offset statement or finite prefixes.
