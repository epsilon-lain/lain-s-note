# Independent review of the escape-count Gram envelope

2026-11-01. Reviewed `ESCAPE_GRAM_ENVELOPE.md`, sections 1--8. No mathematical error found under the stated normalization, rank, and fixed-parameter hypotheses. This is a review of the algebraic reductions, not an estimate of an actual-period spectrum.

1. With X=C B^(-1), the identity C=X B gives B^*B+u C^*C=B^*(I+u X^*X)B. The relative period minor after replacing source index set I by outside index set S is, up to row-order sign, det X_(S,I). This requires B nonsingular exactly where the note says so.

2. Expanding det(I+eta R X) by its principal minors and applying Cauchy--Binet gives the stated sum over matching I,S. The squared sums of the transport minors and period minors are the respective elementary symmetric polynomials. Thus the fixed-count Cauchy--Schwarz bound, and its comparison to the weaker scalar-fugacity bound, have the correct indices and powers.

3. The source-singular-safe Gram pencil is correct. With G=B^*B+C^*C positive definite, whitening gives Q and I-Q. The multiplicity of zero in Q is M-rank B; the multiplicity of one is M-rank C. The determinant coefficients are a Poisson-binomial distribution times det G, with no probabilistic assumption about the period. If the stack is rank deficient, every full minor is zero.

4. For the unrestricted period target sets, the largest coordinate entry of the kth exterior power is at most the product of the k largest singular values. Conversely the squared Frobenius norm of that exterior power is at least the square of its operator norm and is a sum of N_k coordinate squared minors. These observations give both sides of the claimed maximum-versus-spectrum inequality. The logarithmic discrepancy is uniformly O(M), and hence O(1/H) after division by MH in the fixed-weight packet. The note correctly does not transfer the lower bound or equality to the smaller transport-compatible class.

5. The scalar convex transform follows from log(1+exp y)-max(y,0) lying between zero and log 2. The infimum formula for the sum of the largest k normalized log singular values also handles k=0 and k=rank X as stated; zero singular values remain at minus infinity and contribute zero to the positive-part sum.

6. The generic aligned/disjoint examples have identical separate Gram spectra but different products R X. The scalar optimum occurs at u=r/s and gives the displayed bound. The single-exit example indeed changes a source-relative logarithmic determinant by order MH at escape fraction 1/M. Both are clearly labeled generic linear algebra, not actual-period counterexamples.

7. The E,D refinement is a valid similarity of R X. Its first determinant expands into terms |det R_(I,S)|^2 times product E_I divided by product D_S; its second has the reciprocal weight monomial and |det X_(S,I)|^2. Thus each logarithm is log-sum-exp of affine forms in the diagonal logarithms. The claimed convexity and disjoint-support limiting recovery follow directly. This still discards phases and does not prove exactness of the variational infimum in general.

The newly proved small-escape Hardy bound has a different normalization: it controls absolute period minors using the complete Hardy source Gram and its explicit Hermite right inverse. It does not control C B^(-1) for an arbitrarily small selected polynomial source minor. Therefore the note's source-relative spike warning and the actual-period Hardy continuity result are compatible.
