# Independent review of the primitive reflection factors

2026-11-02. Reviewed `PRIMITIVE_REFLECTION_FACTORS.md` independently, including the leading coefficients, zero orders, arithmetic normalization, and the selected-factor comparison. No substantive mathematical error was found. The factor argument remains a structural supplement; its direct comparison has the stated quadratic logarithmic loss and is weaker than the separate normalized-Gram comparison developed in this round.

## Orders and signs

At Z=0, the plus entry has multiplier 1+(-1)^s and the minus entry has multiplier 1-(-1)^s after the coefficient extraction is evaluated. Thus each determinant has r entire rows divisible by Z. The already established exact order 2r of their product forces both individual orders to be exactly r.

The parity calculation has no missing column sign. On replacing Z by -Z and then the coefficient variable z by -z, the sign (-1)^a from the spatial power cancels the corresponding interchange of the two exponentials. The remaining row sign is (-1)^s for the plus factor and (-1)^(s+1) for the minus factor. Since n=2r, both determinant signs are (-1)^r. Consequently division by Z^r leaves an even polynomial.

For the leading term, grouping even columns before odd columns has sign (-1)^(r(r-1)/2). The two binomial evaluation determinants each contribute 2^(r(r-1)/2). The powers of Z/2 contribute 2^(-r^2), leaving 2^(-r) before the coefficient-vector determinant.

The confluent Vandermonde coefficient matrix W at nodes +1/2,-1/2 has determinant (-1)^(r^2) A_r^2/A_n. Passing to the sum/difference columns multiplies it by (-2)^r, giving the positive factor 2^r A_r^2/A_n. This proves the displayed leading coefficient of L_+. Swapping the two r-column groups contributes (-1)^(r^2)=(-1)^r and gives the displayed coefficient of L_-. These signs also reproduce the known product leading coefficient (-1)^r/d_n.

As a separate arithmetic check, direct exact matrix determinants verified the two binomial formulas and the confluent matrix formula for r=1,...,6: 18 checks passed. The general argument above does not depend on those finite checks.

## Integer factors and rational values

The normalizations make F_+ and F_- monic polynomials over Q. Their product times Z^n is the already proved monic integer polynomial (-1)^r d_nP_n. Gauss's lemma therefore applies to each monic rational factor separately and gives integer coefficients. No assumption of separate entrywise denominator cancellation is hidden in this step. Their leading coefficient one implies primitivity.

Their degree Delta=r(r-1) is even. At ia/b in lowest terms, multiplying the value by b^Delta gives an integer whose residue at each prime dividing b is i^Delta a^Delta, hence nonzero. The reduced denominator is therefore exactly b^Delta, and the value is nonzero for b>=2. Passing to a polynomial in Z^2 does not improve this denominator: it halves the degree while squaring the argument's denominator.

## Derivative and selection of the smaller factor

The reflection entry bounds are correct: on |Z|<8 its magnitude is at most 2*5^a, while its derivative is at most a*5^(a-1). The factor 1/2 from differentiating Z/2 cancels the factor 2 from summing the two exponential contributions. The determinant product rule consequently gives

\[
|L_\pm'|\le n!2^{n-1}S_n5^{S_n-1},\qquad
S_n=n(n-1)/2.
\]

On the comparison segment, |Z|>1. Differentiating Z^(-r)L_\pm then gives the bound B_n n!2^n5^(S_n)(r+S_n/10), which is at most the displayed M_n. There is no extra missing factor n: the product rule has already summed the spatial degrees into S_n.

The actual-period product identity guarantees that at least one factor has magnitude at most V_n. It is legitimate to choose that factor: both signs have the same degree, rational-value floor, and derivative bound. The sign need not be identified uniformly in n. With V_n<=1/(2b^Delta), the triangle inequality yields a segment-length lower bound 1/(2b^Delta M_n), and converting segment length to |pi-p/q| introduces the further factor 2. Thus the constant 4 in the final bound is correct. The treatment of segments longer than one removes the initial proximity condition.

## Asymptotics and scope

The expansion log(V_n)/Delta=-log n+2C+o(1) retains the same threshold kappa=exp(2C) as the full-polynomial comparison. Choosing n=cb+O(1), c>kappa, makes the true-period sufficient condition hold eventually. The two leading costs in the resulting error bound are Delta log b and log B_n; together they give

\[
-\tfrac12c^2b^2\log b+O_c(b^2).
\]

The use of 44/7 in a finite rational certificate is sufficient, not equivalent to the exact pi-dependent threshold. Nothing here proves general irreducibility, excludes additional factors, or establishes a quadratic Diophantine lower bound. The note preserves all three limitations.
