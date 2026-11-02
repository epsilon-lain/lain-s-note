# Zuuv v2 static audit — finite rational CF versus certified pi prefix
2026-11-02

Scope
Repository: sinandredemption/zuuv
Pinned commit: 496dede87362b2e8bb7a7acec8ecc11301599bc6 (2021-09-11)
Files read in full: source/cruncher.cpp, continuant.cpp, continuant.h, params.h,
disk_mpz.cpp, disk_mpz.h, io.cpp, io.h, mul.cpp, utils.cpp, cruncher.h, README.md.
No 653.52-billion-term recomputation was performed, and no binary/source identity
attestation for the 2025 author's actual executable was established here.

Overall result
The central mathematical mechanism is a standard, valid integer half-GCD /
continued-fraction computation. Continuant and full-size correction formulas are
integer identities; the computation is not based on floating-point quotients.
However, the pinned disk-mode source does not emit a rigorous certificate that
each chosen prefix is shared by an interval enclosing pi. It also relies on a
fixed ten-term backoff and checks only residual sign, not the complete regular-CF
tail-domain predicate. Therefore this source alone does not justify claiming the
2025 record has been independently rigorously verified. This is a certification
gap, not proof that the reported record is wrong.

1. What is exact and mathematically sound

cruncher.cpp 16-78: the base case computes Euclidean quotients with integer MPIR
division/subtraction. The optional most-significant-limb fast path accepts a
quotient only when all four coarse endpoint ratios agree; otherwise it falls
back to full long division. This is intended as an interval/corner test for that
single quotient (subject to ordinary C++ unsigned arithmetic edge cases).

continuant.cpp: recursive continuants use the standard split identity; block
combination multiplies the two 2x2 convergent matrices. Entries are MPIR/YMP
integer products and sums. There is no transcendental approximation in this step.

cruncher.cpp 147-167 and 443-466 compute the inverse-Mobius correction. For an
input rational N/D and a candidate prefix matrix [p,p_prev;q,q_prev], the residual
tail is
  t=(D*p_prev-N*q_prev)/(N*q-D*p).
It exactly reconstructs N/D as (p*t+p_prev)/(q*t+q_prev). Three thousand
independent Fraction examples verify this identity in verify_static_formulas.py.

disk_mpz.cpp 187-219 implements these products/differences over all stored integer
blocks, with carries. The intended correction arithmetic is full precision,
even though the candidate coefficients were selected using leading blocks.

2. The disk-mode certification gap is concrete

cruncher.cpp 420-435:
- read only the top two stored blocks;
- trim them to a working leading size;
- compute a candidate half expansion;
- pop Params::TruncateConvergents (ten) terms;
- recompute its integer continuants.

params.h line20 fixes TruncateConvergents=10 and labels this a safety margin,
suggesting four might suffice. No quantitative bound here links that fixed
term count to an enclosing input interval or to a verified matrix-growth bound.
The threshold/reduction factors at params.h 9-18 are chosen numerical cutoffs
(2.0,2.8,3.0), not an exposed certificate of common CF digits.

After full correction, disk mode at cruncher.cpp 469-477 only normalizes when
both residual numbers are negative and throws if signs differ. It does not
check correction_num >= correction_den. At line488 its own comment states the
reasoning in probabilistic language (almost-certain positivity).

In contrast, the RAM/internal routine at cruncher.cpp 175 includes the predicate
correction_num < correction_den and invokes perform_corrections, and RAM main
has the corresponding tail>=1 assertion at line347.

Why sign is insufficient: input N/D=7/5, guessed prefix [1;1], gives correction
(-2)/(-3), normalized to 2/3. Both numbers are positive, so the disk predicate
accepts, but a regular-CF residual tail must be at least1. The actual regular CF
is [1;2,2]. This is an explicit counterexample to the sign-only predicate being
a certificate. It is NOT a demonstrated end-to-end counterexample surviving
this program's preceding ten-term trim; a proof that trimming always prevents
such states would address this issue, but none was found in these source files.

3. Exact CF of a finite digit rational is not automatically CF of pi

disk_mpz.cpp 10-77 explicitly interprets a finite hexadecimal digit file, ignores
the integer part, and constructs the exact rational fractional part. For pi this
changes the first coefficient from3 to0 but leaves subsequent coefficients and
their positions unchanged when the initial coefficient is counted consistently.
The disk reader is hexadecimal; pushback_digits at 303-309 uses base16.

The source stores only one numerator/denominator. It has no representation of
the upper endpoint obtained by adding one unit in the last input digit, no
propagated pi tail bound, and no comparison of the two endpoint CF prefixes.
The requested number of CF terms is supplied by the user. README's conversion
from number of digits to number of CF terms is explicitly a rule of thumb.
Thus the source does not by itself prove that the retained terms are insulated
from truncation of pi's finite input digits.

An appropriate certifier would start with rigorous L<pi<U and verify either:
(a) transformed interval endpoints share each integer quotient; or
(b) the entire interval lies in the CF cylinder of each stored matrix block.
For a regular prefix through p_n/q_n, that cylinder has endpoints p_n/q_n and
(p_n+p_{n-1})/(q_n+q_{n-1}), ordered by parity. The interval must lie strictly
inside it to avoid rational-boundary ambiguity.

4. Peripheral source issues must not be conflated with this record

RAM read_from_file at io.cpp 23-26 constructs a power-of-two denominator from
binary length even when the UI selected decimal digits. For decimal input
3.14159 this represents314159/131072, not314159/100000. The hex/pi disk-mode
record does not use this decimal RAM path, so this is not evidence that its
reported coefficients are wrong.

Very short-list and maximum-limb edge cases deserve separate implementation
testing. They were not needed to infer the certification gap and have not been
used to reject the record. No unexercised corner case is claimed to have occurred
in the author's run.

5. What may safely be concluded for the user

The method searches already-known pi information efficiently: convert a long
digit prefix into continued fractions, and identify unusually large partial
quotients. It does not construct an infinite sequence with quality tending0,
nor prove a final positive minimum. The integer identities underlying the
strategy check out. The exact 2025 execution remains a reported computation
until input bounds, safe-prefix certification and sufficiently independent
recomputation/certificates are inspected.

Canonical pinned sources
https://github.com/sinandredemption/zuuv/blob/496dede87362b2e8bb7a7acec8ecc11301599bc6/source/cruncher.cpp
https://github.com/sinandredemption/zuuv/blob/496dede87362b2e8bb7a7acec8ecc11301599bc6/source/continuant.cpp
https://github.com/sinandredemption/zuuv/blob/496dede87362b2e8bb7a7acec8ecc11301599bc6/source/params.h
https://github.com/sinandredemption/zuuv/blob/496dede87362b2e8bb7a7acec8ecc11301599bc6/source/disk_mpz.cpp
https://github.com/sinandredemption/zuuv/blob/496dede87362b2e8bb7a7acec8ecc11301599bc6/source/io.cpp

Upstream source is referenced by pinned links; this directory does not bundle it.
