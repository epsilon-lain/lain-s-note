# External inputs and provenance

The elementary derivations in this directory are working proofs, not a claim
of bibliographic novelty. Computation is performed locally using rational
intervals. No external source is taken as independent validation of our proofs.

- Kurt Mahler, *On the approximation of π*, Indagationes Mathematicae 15
  (1953), 30–42. Theorem 1, original p.33, reprint p.561: the explicit
  bound |π−p/q|>q^(-42), for positive p and integers q>=2.
  [Publisher](https://ems.press/books/dms/252/4986),
  [original-paper reprint PDF](https://content.ems.press/assets/public/full-texts/books/252/chapters/online-pdf/252-chapter-4986.pdf).
  This is used only in the pi repeated-prefix cap, not the root-envelope
  obstruction. The theorem itself was rechecked in the original reprint.
- Classical fact: π is transcendental, hence π² is irrational. Only the
  latter consequence is needed in ROOT_ENVELOPE_OBSTRUCTION.md. We do
  not supply a new proof of transcendence. The general theorem is stated
  separately under the explicit assumption x² irrational.
- Standard continued-fraction facts used in QUADRATIC_BARRIER_RESULTS.md:
  the convergent error estimates and Legendre's criterion. An exposition
  with proofs is [Bruce Ikenaga, Approximation by Rational Numbers](https://sites.millersville.edu/bikenaga/number-theory/approximation-by-rationals/approximation-by-rationals.html).
- [Archive of Formal Proofs: Liouville Numbers](https://isa-afp.org/entries/Liouville_Numbers.html)
  is a reference for the definition. It does not formally verify this
  directory's variable-base equivalences.
- [Open Math Harness](https://github.com/openmathmodel/open-math-harness)
  is the intended future research environment, not a mathematical input.
  Its README, inspected 2026-10-27, says the full runtime relies on the
  Open Math OS deployment environment and does not yet provide a
  standalone local setup for interactive agents or proof search.
  OPEN_MATH_REVIEW.md is deliberately a portable text task, not an
  invented integration schema.

The earlier historical report mentions an OpenAI preprint about the
irrationality exponent of pi. That result is not used as an input to any
current root-envelope theorem or certificate here. The supplied local-smoothing
wave-equation paper has not been turned into a proof of any Diophantine claim
in this package.

No third-party full papers or private account/application data are included.
