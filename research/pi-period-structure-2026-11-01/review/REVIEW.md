# Independent review: actual-period cancellation and grouped estimates

2026-11-01. This review checks the four notes on diagonal order, grouped Hardy norms, period stability, and the Laurent gauge. It also reruns their exact verification scripts. It does not independently certify the manuscript's Theorem 2.1 or establish any result on the bad approximability of π.

**Verdict:** the four principal arguments are valid under their stated hypotheses. Two qualifications were tightened during review: a product-gain bound needed a positive part when logarithmic gains can be negative, and the grouped Gram inequality must not be described as preserving every phase cancellation in the original determinant. Both corrections have been incorporated.

## 1. Exact vanishing at the actual period

Let $D$ be a nonzero polynomial over an algebraic number field. In coordinates

\[
z=r_1,\qquad \delta_i=r_i-r_1\quad(i\ge2),
\]

write $D=\sum_\gamma F_\gamma(z)\delta^\gamma$, and let $k$ be the least $|\gamma|$ with $F_\gamma\ne0$. Known transcendence of $\omega=2\pi i$ implies $F_\gamma(\omega)\ne0$ for every nonzero coefficient polynomial. Thus substitution of $z=\omega+\varepsilon_1$ and $\delta_i=\varepsilon_i-\varepsilon_1$ produces no terms below degree $k$, and its degree-$k$ part is a nonzero polynomial in the independent normal coordinates. This verifies the equality between the local Taylor order at the actual period and the diagonal-ideal order.

In particular, vanishing of $D$ at $(\omega,\ldots,\omega)$ is equivalent to the polynomial identity $D(Z,\ldots,Z)=0$. The conclusion uses transcendence and is qualitative: it does not bound a small nonzero value there. Along a special direction of the error vector, even the nonzero leading homogeneous polynomial may vanish, causing additional cancellation.

The weighted-degree bound also checks. For fixed columns with transverse support budgets $b_\ell$, each nonzero matrix entry has weighted center degree at most $b_\ell-w\beta$. Thus

\[
\deg_wD\le E:=\sum_\ell b_\ell-\sum_\rho w\beta_\rho.
\]

When $w_1\le w_i$, the change $r_i=z+\delta_i$ does not increase weighted degree if $z$ has weight $w_1$ and $\delta_i$ weight $w_i$. Hence $k\min_{i\ge2}w_i\le E$. No preservation of the weighted source packet under a change of transverse jet coordinates has been assumed.

Both complete four-row examples were checked by the exact verifier: the raw monomial minor is $r_1(r_1-r_2)$, and the mixed-column minor is $(r_2-r_1)^4$. In the latter case the arithmetic lower bound $(2/\operatorname{lcm}(q_1,q_2))^4$ and the error upper bound give precisely the ordinary rational-spacing inequality. The fourth power does not by itself improve that comparison.

For an explicit product of difference factors, weighted degree charges the edge $(i,j)$ by $\max(w_i,w_j)$. Under $|e_j|\le |e_i|/2$, the factor gain is at most $\alpha_i=-\ln|e_i|$ per power. Replacing the total charged degree by its upper bound $E$ is valid only when the final multiplier is nonnegative. The revised formula uses the positive part of the maximum gain-to-weight ratio; this fixes the negative-gain case without imposing an unnecessary extra approximation hypothesis. The residual polynomial factor and sums of different normal monomials remain uncontrolled by this particular estimate.

## 2. The whole-shift Hardy factorization

The factorization $A=ZC$ is an exact actual-period identity. The test rows of $Z$ are square summable because their evaluation points lie strictly inside the Cauchy disk, and the coefficient columns of $C$ are square summable by Parseval. Finite truncations followed by Cauchy–Binet and Cauchy–Schwarz give

\[
|\det A|^2\le\det(ZZ^*)\det(C^*C).
\]

Taking the truncation limit is legitimate because the finite products and the two finite Gram matrices converge. If each entire function has disk-bound $D_H$, the separate-type direct sum gives column squared norm at most $N_A D_H^2$. Hadamard then yields the displayed relaxed coefficient-Gram bound. The factor $N_A$ is required; its logarithmic cost is lower order at fixed dimension.

The reproducing-kernel formula has the correct conjugation: the second actual-period node is $-k\omega$, while the coefficients of $\ln(1+t)$ are real. The finite transverse factor is the exact inner product over the type label. This keeps the shift dependence and row correlations represented by $Z$; the exterior Cauchy–Schwarz inequality can still lose relative phase information between the test and coefficient exterior vectors. The final wording now preserves this distinction.

At zero shift, types split into independent Gram blocks. For $n$ consecutive time jets at each of $K$ centers, the triangular change from $\ln(1+t)$ jets to normalized ordinary derivatives has diagonal one. Rescaling to the unit disk contributes

\[
R^{-2K\sum_{s=0}^{n-1}s}=R^{-Kn(n-1)}.
\]

The confluent Cauchy determinant supplies the remaining quotient of node-difference and kernel factors. There are no extra factorials because normalized derivatives are used.

The resulting period exponent is consistent term by term. In half the logarithm of a block determinant, the coefficient of $n^2\ln R$ is $-K^2/2$; the remaining linear term is $(Kn/2)\ln R$. Summing over types yields exactly the stated $\Psi_{K,R}$ formula. Since

\[
\frac{\sum_\beta n_\beta^2}{MH}\longrightarrow
\frac{2}{K v_0(m+2)},
\]

the normalized saving is $2\Psi_{K,R}/[K v_0(m+2)]$. Complete consecutive jets allow the separate $\ln2/v_0$ coefficient-extraction charge to disappear here. That refinement is not automatically available for arbitrary sparse selections.

The large-$K$ expression follows from the weighted logarithmic sum of node differences and the uniformly convergent kernel series. The integrable logarithmic endpoint requires the indicated cutoff argument. At $R=14K$, the rational bounds in the note establish $\phi(14)>11/10$; this is a limit coefficient for the zero-error full packet, not a uniform estimate for all nonzero error vectors.

The finite nonzero-shift example uses polynomial columns without a $Y$ factor, so its rational imaginary base node does not falsely invoke exponential periodicity. Its grouped bound is less than $3/2$ times the actual determinant modulus and more than a factor $100$ better than the stated rowwise Hardy relaxation. This compares two specific finite bounds and does not certify the endpoint asymptotics.

## 3. Stability and forced relative spectral growth

For nonsingular $B$, the relation

\[
A=B(I+E),\qquad E=B^{-1}(A-B)
\]

gives the stated determinant stability bounds whenever a dominating norm of $E$ is below one. A fixed norm bound below one suffices for negligible change at normalization $MH$; convergence of that norm to zero is not required.

The Gaussian-rational preconditioner construction is rigorous. Its rational coefficient majorant dominates the complex row-sum norm throughout the polydisc. If the majorant is less than one, every matrix in the polydisc is invertible, with the displayed determinant bounds. The two finite disks containing the actual period are certified by exact Machin enclosures, not floating-point evaluation.

For the asymptotic comparison, the new $q_iL_i$ scaled and mixed columns give arithmetic lower exponent $b_{\kappa,H}$. The period upper exponent pays $\kappa_{\max}$ once. Their difference is therefore

\[
\Gamma=S_{\rm per}-R_{\rm ar}-E_{\rm per}.
\]

If the determinant ratio is at least $\exp(MH(\Gamma-o(1)))$, at least one eigenvalue of $B^{-1}A$ has modulus at least $\exp(H(\Gamma-o(1)))$. Subtracting one gives the lower bound on the spectral radius of $E$. This does not require diagonalizability. A further common invertible column-basis change conjugates $E$, so this spectral conclusion is basis invariant.

This is not a condition-number lower bound by itself: such a bound would also require quantitative information about the absolute perturbation. Nor does it apply when the chosen period minor is singular. The common mixed-column construction explicitly ensures that both square matrices are nonsingular before this argument is used.

For the full period packet, the identity-type second moment is $2K/[v_0(m+2)]$, giving the reported high-dimensional saving. The exact margins at dimensions $32$ and $40$ retain the unknown interpolation thresholds in addition to the explicit minimum-weight requirement. They imply the stated exponential spectral-growth bounds only for sufficiently large admissible $H$; no explicit such $H$ is claimed.

The separate fixed-$355/113$ staircase family is independent of that interpolation theorem. Its counts satisfy

\[
M_n=20n(n+1),\qquad
W_n=\frac{10}{3}n(n+1)(2n+1),\qquad
W_n/M_n=(2n+1)/6.
\]

The staircase columns obey the polynomial budget for $H=20n$. The exact determinant formula, rather than a materialized huge matrix, justifies the large finite example. Its positive binomial-sum certificate gives the reported lower bound greater than $100$. A very small fixed center displacement therefore does not force uniformly small relative matrix perturbation as the packet grows.

## 4. The Laurent gauge and its joint budget

Multiplication of every column by $Y^{-d}$ acts on source rows by convolution with $(1+t)^{-d}$. The complete packet is closed under lowering the time index at fixed center and transverse type, so the row operation is unitriangular. Every full-row minor on its fixed columns is unchanged.

Negative integer powers of $1+t$ have integer binomial coefficients. Thus the same $q_iL_i$ row-and-column normalization remains Gaussian integral; no additional denominator charge is needed. At the actual period, integer $d$ preserves the required exponential periodicity.

For a scaled column under the single budget $w_0h+w\alpha\le H$, its logarithmic norm is bounded by

\[
R|h-d|+(H-w_0h)b_R.
\]

For $0\le d\le H/w_0$, convexity in $h$ reduces this to the two endpoint values. Balancing those endpoints and rounding $d$ down proves

\[
\ln D_H\le H\max\{b_R,(R/w_0+b_R)/2\}+R.
\]

The extra $R/H$ is explicitly controlled for fixed radius, so allowing $d=d_H$ to grow with $H$ creates no hidden uniformity gap. The arithmetic scaling is already included in $b_R$ and must not be added again.

The sparse-rank qualification is essential and correctly included: an arbitrary selected sparse test minor is not invariant under this row convolution. Instead, one applies the gauge to the complete auxiliary one-center packet, preserves its independence, and then takes the required subset. The period-polynomial nonvanishing and common-column steps are subsequently repeated.

Finally, the radius optimization of this new joint envelope was checked independently. If $\kappa=\max\kappa_i$, then

\[
G_R-\kappa\ge \tfrac12\max(0,R/w_0-\kappa).
\]

With $A=7K/w_0$ and $d_0=1/(2v_0)$, the upper objective is
$d_0\ln u-\tfrac12(Au-\kappa)_+$ for $u>1$. When $\kappa>A$ and $\kappa>2d_0$, it increases up to $u=\kappa/A$ and decreases thereafter. This gives the stated upper bound $\ln(\kappa/A)/(2v_0)$. Along the specified shell family it is $O(m/v_0)\to0$. Thus the earlier low-score witness survives this particular improved uniform envelope and per-dimension radius optimization. This remains a limitation of that bound, not a lower bound on an actual test determinant or a prohibition on sharper grouped estimates.

## 5. Explicit coefficient height and the direct normal-order bound

The final coefficient-height addition was checked independently. For a raw monomial column, replace the coefficients of the formal logarithm by their absolute values and evaluate this majorant at $|t|=1/2$. Its sum is $\ln2<1$. In normal-center coordinates the polynomial $jr_i$ has coefficient norm at most $2j$. Thus a source entry has normal-polynomial coefficient norm at most

\[
2^s(3/2)^h\binom\alpha\beta(2j+\ln2)^{|\alpha|-|\beta|}
\le2^s(3/2)^h(4K)^{|\alpha|}.
\]

This is a coefficient-majorant argument, not an inference from a small function supremum to a coefficient one-norm. Submultiplicativity of that one-norm and the $M!$ determinant permutations give exactly the stated bound for $L=\|D\|_1$ in normal coordinates. The row time budget bounds $\sum s/(MH)$ by $1/v_0$, while the single column budget bounds the combined $h$ and transverse powers by the maximum of their normalized costs.

Because the normal weighted degree is at most $MH$, its ordinary degree is at most $MH/w_*$. For $0<\eta\le1/2$, the normal errors have size at most $2\eta\le1$, and the common center has size below eight. A normal order $k$ therefore yields the stated upper bound containing

\[
-\frac{k}{MH}\ln\frac1{2\eta}
+\frac{\ln2}{v_0}
+\max\{\ln(3/2)/w_0,\ln(4K)/w_*\}
+\frac{\ln8}{w_*}+\frac{\ln M}{H}.
\]

The condition $\eta>0$ is needed when this estimate is written with logarithms. At zero error the preceding polynomial inequality can instead be interpreted directly.

For arbitrary fixed polynomial columns, their actual coefficient one-norms must replace the raw-column factors. In particular, the arithmetic pre-scaling cannot be omitted from the analytic height charge. The note states this cost explicitly.

Finally, if $L_\eta=\ln(1/(2\eta))$ is at most the least normal weight, the degree bound implies $kL_\eta/(MH)\le1-\bar b_H$. This prevents closure of the particular **uniform** arithmetic comparison whose required saving is $1-\bar b_H$ plus its displayed positive errors. It does not rule out a sharper column-specific arithmetic threshold, a directional error estimate, a smaller actual coefficient norm, or additional cancellation. The final statement correctly retains this limited scope.

## 6. Verification and scope

All four accompanying verifiers were run successfully during review. Their checks cover the actual-diagonal polynomial identities, the confluent and shifted Gram formulas, actual-period rational preconditioner certificates, exact staircase and high-dimensional margins, and the Laurent row identities and envelope arithmetic. The general statements rest on the proofs; finite verification does not replace them.

The order of choices remains important: dimension, weights, approximation tuple and Cauchy radius are fixed before each $H$ limit. The Laurent gauge is the explicitly controlled exception allowed to vary with $H$. The diagonal classification uses known transcendence of π, and high-dimensional rank existence uses the stated manuscript theorem. No proof of π's bad approximability, non-bad approximability, or a new rational-approximation record has been obtained.
