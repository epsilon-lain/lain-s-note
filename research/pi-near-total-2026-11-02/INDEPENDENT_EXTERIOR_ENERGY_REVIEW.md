# Independent review of exterior column energy

2026-11-02. Reviewer: near-total-cost-coupling agent. Reviewed `EXTERIOR_COLUMN_ENERGY.md` and ran `verify_exterior_energy.py` independently.

**Result: passed; no mathematical correction required.**

- Equations (2) and (3) follow directly from the nonnegative Cauchy--Binet expansion. The direction of the factor comparison is correct for 0<q<1. The note correctly avoids promoting the determinant inequality to an unsupported Loewner inequality.
- The top-fraction energy formula (4) is the tail first moment of the Beta(d,1) cost distribution. Its half-density expansion (5) has coefficient 1-ln(2), with the stated sign.
- The full local-jet coordinate change has determinant one: its linear shear and its one-variable formal substitution each induce determinant-one maps on the complete truncated jet space. This justifies preserving both finite and infinite Gram determinants.
- The Pascal block minors, the exact geometric Gram determinant, and the composition identities for E_0 and T agree. The optimized old cutoff has q_*=(n+d)/(3n+2d-1), and its normalized positive loss is exactly (3 ln(3)-2 ln(2))/(2(d+1)).
- The repaired estimate uses the valid inequality E_max<=E_0+dM. At q=1/(n+1), its normalized logarithm is bounded by [d ln(n+1)+1]/(2n), so the loss disappears at the leading MH scale.
- The example has theta=1 and equal source weights. The note accurately excludes a conclusion about the critical separated-weight construction or pi.

The verifier passed 12 complete packet models, 1,215 exact local-coordinate entries, 18 exact geometric Gram blocks, and 18 finite weighted-simplex selections. The asymptotic assertions rely on the accompanying integral and exact-formula proofs, rather than on those finite checks.

The outstanding issue is unchanged: obtaining a sufficiently strong actual-period estimate for the near-total-escape rows in the critical parameter family, with the same columns' arithmetic normalization.
