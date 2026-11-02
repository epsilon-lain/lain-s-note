# Independent review of the finite same-column period certificates

2026-11-02. Reviewer: target-spatial-filtration agent. Independently read `PHASE_PERIOD_NONVANISHING.md` and `verify_phase_period.py`, ran the verifier, and checked the scaled matrix formula by a second direct coefficient calculation. The finite nonvanishing certificates and the stated limitations pass review.

## 1. The finite-field calculation is a valid one-sided certificate

The modulus 1,000,033 is prime; the verifier tests all possible integer divisors up to its square root. The selected image of i is 649,529, and its square is -1 modulo that prime. Since 7 is invertible, the evaluation map from Z[i,1/7] to the field is well-defined. A nonzero determinant residue therefore proves that the characteristic-zero Gaussian rational determinant is nonzero. No injectivity of the reduction map is needed.

Multiplication of each row by t! times the product of the transverse factorials clears coefficient-normalization denominators. These scalars are nonzero in characteristic zero. The displayed falling-factorial entry formula agrees with direct coefficient extraction; I independently checked 2,430 entries with exact rational arithmetic at the real formal evaluation Z=2 before reducing modulo the same prime. This comparison is separate from the main modular determinant algorithm.

The finite-field square determinant is recomputed after choosing the source pivot columns. The target uses precisely those same raw columns. Its nonzero residue is not inferred from source rank, from another target basis, or from support matching.

## 2. The actual-period inference is exact

For every recorded fixed row and raw-column selection, the unscaled determinant is a polynomial in Z with rational coefficients. Nonvanishing at 44i/7 proves that polynomial is not identically zero. The transcendence of 2πi then proves nonvanishing there. The columns are fixed after selection and before treating Z as a formal variable, so there is no hidden Z-dependent basis substitution.

This inference proves only nonzero. It provides neither a lower bound on the determinant magnitude nor proximity of any algebraic root to the actual period. The text keeps this distinction explicit.

## 3. The complete source, target and budgets are retained

The source loops use the full strict weighted cutoff, including every admissible spatial multi-index and every admissible nonnegative time index. The ceiling endpoints implement that strict cutoff correctly, including when a weight ratio is integral. The raw-column loops implement the non-strict weighted budget correctly for the rational model weights.

Both parameter families satisfy the displayed half-density normalization exactly. Their finite source/column counts are separately enumerated; the asymptotic half-density is not used as an exact finite identity. The time-budget parameter inequality is checked explicitly.

The target construction is the stated dilation-two, binary-phase construction with the low-order identity correction. It retains distinct labels, nonnegative times, conserved total order, coordinatewise upward movement and the target transverse budget. The selected-column coordinate surpluses and all intrinsic time cuts are checked on the same selected raw family.

The exact Pascal calculation uses the full lower-set bucket. Integer dilation two makes each Newton diagonal factor 2^|β|; the verifier independently checks its determinant with integer Bareiss elimination. This is stronger than an entrywise matching certificate.

## 4. Verification and scope

The executed verifier passed all 16 listed complete packets, with row counts 32 through 240. Every packet has two nonzero residues on a common raw column set. The finite transport and time-cut checks also passed.

The proof and code establish those 16 finite cases. They do not establish an all-heights theorem, an infinite-dimensional limit, a logarithmic-denominator weight realization, a bound on the complete aggregate, or a Diophantine consequence for π. Both the main note and output JSON state these limitations.

One textual correction was requested and incorporated in the comparison with earlier zero-center rigidity: strict weighted descent there invokes v0*T_i>v_i, not merely >=. This affects only the account of the prior argument, not any present certificate. The ordinary shear raises the old weighted order when v_i>v0, so its nontrivial zero-center movements are not contradicted by that earlier truncated-tail argument.
