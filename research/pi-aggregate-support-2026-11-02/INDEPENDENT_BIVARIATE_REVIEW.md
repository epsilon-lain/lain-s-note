# Independent review of the exact bivariate-block appendix

2026-11-02. Reviewed `EXACT_BIVARIATE_BLOCKS.md` and reran `verify_bivariate_blocks.py` independently of their author. No correction is required.

The condition v_0>=v is used in the correct direction. A nonzero coefficient of z^s u^l in (e^z-1)^r(z+u)^a has a>=l and r<=s-a+l, so its direct-coordinate cost is at most v_0s+vl-(v_0-v)(a-l). Every required direct row therefore belongs to the same packet. Ordering by time index makes the row change triangular, with diagonal coefficient one. The shifted-Y column change also preserves the full budget by lowering the Y exponent.

The resulting time-exponent blocks are ordinary polynomial Hermite interpolation at distinct points jZ. Their ranks are exactly min(K n_h,d_h), including the possibility of an empty row or column block. This proves the finite rank formula for every nonzero Z. It does not assert that the same transformation preserves a packet when v_0<v. The stated Gram factorization correctly uses the shifted columns, without identifying their metric with the raw-column metric.

For the 15-row counterexample, (T-Y+1)^4 has permitted column degree and vanishes to weighted order at least eight after substituting Y=e^z and T=z+u. The exact nonzero 14-row minor independently supplies the complementary rank lower bound. Thus rank 14, rather than merely a nontrivial kernel, is certified.

The verifier passed all 36 exact models and 14,592 coordinate-identity comparisons, the shifted-column Gram identities, the rank-14 example with determinant 1/640, and eight invisible-polynomial-family cases. These calculations support the written arguments and do not extrapolate the shear-invariant formula to the critical small-time-weight regime.
