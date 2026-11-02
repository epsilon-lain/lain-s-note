# Independent review of the period logarithmic derivative

2026-11-02. This review independently rederives the central calculation in
`PERIOD_LOG_DERIVATIVE.md` and checks its consequence for the largest complex
disk at the fixed threshold `2 U_n`. It uses the previously proved Jacobi
variance and cosine lower bound, and the previously proved explicit disk
bound. It does not assume a central limit theorem or a limiting value for the
Jacobi characteristic function.

## Verdict and scope

The endpoint-insertion argument proves the stated factorial-scale asymptotic
for the actual derivative of the formal polynomial. Combined with the prior
disk estimate and Cauchy's inequality, it proves

\[
\log\mathcal R_n=-2n\log n+O(n)
\]

for the largest complex disk centered at `2 pi i` on which this polynomial
stays below `2 U_n`. This does not determine the distance to a root, the
behavior on a selected ray, or a rational-approximation exponent for pi.

Two drafting points were corrected during review: the logarithm used below
is the sum of the individual factor logarithms, not necessarily the principal
logarithm of their product; and polynomial parity follows from reflection of
the entire moment functional, not parity of each individual moment.

## Exact insertion and phase

Put `t=2 pi`, `omega=i t`, `m=n-1`, and let `Delta_m` be the Hilbert determinant.
The endpoint functional satisfies, for every polynomial `f`,

\[
\mathcal V_\omega(f)
=-\frac{i}{t}\int_0^\infty e^{-v}f(1+iv/t)\,dv.
\]

Expanding a monomial verifies the finite moment formula term by term; there
is no contour approximation. Differentiating the determinant's
Vandermonde-square formula inserts this functional once, with coefficient
`n/n! = 1/m!`. Hence, for `D_n(Y)=det[mu_(i+j)(omega,Y)]`,

\[
D_n'(1)=-\frac{i\Delta_m}{t}
\int_0^\infty e^{-v}\mathbb E_m\left[
e^{it\sum u_j}\prod_{j=1}^m(1+iv/t-u_j)^2\right]dv.
\]

Factor `iv/t` out of every linear factor. The resulting phase `(-1)^m`
cancels the centering phase `exp(itm/2)=(-1)^m`. With independent
`V ~ Gamma(2m+1,1)`, this proves exactly

\[
D_n'(1)=-i\Delta_m\frac{(2m)!}{t^{2m+1}}
\mathbb E[e^{itX_m}C_m],\qquad
C_m=\prod_{j=1}^m\left(1-\frac{it(1-u_j)}V\right)^2.
\]

Thus the centered prefactor is `-i`, independent of parity. The parity in
the final answer enters through
`D_n(1)=(-1)^n Delta_n chi_n(t)`.

## The unbounded Gamma tail is controlled

On `V >= m`, for `m >= 2t`, define the additive logarithm

\[
L_m=2\sum_{j=1}^m\operatorname{Log}
\left(1-\frac{it(1-u_j)}V\right),\qquad C_m=e^{L_m},
\]

where each logarithm is given by its convergent power series. The quadratic
remainder estimate gives

\[
L_m=-it\frac mV+2it\frac{X_m}V+R,
\qquad |R|\le\frac{2t^2}{m}.
\]

The first two terms are purely imaginary. Since `exp(-it/2)=-1`,
`E X_m^2 <= 1/15`, and `E(V-2m)^2=2m+2`,

\[
\mathbb E[|C_m+1|;V\ge m]
\le \frac{2t^2e^{2t^2/m}}m
+\frac{t\sqrt{2m+2}}{2m}
+\frac{2t}{\sqrt{15}\,m}.
\]

Concentration by itself would not control the remaining event: `C_m` is
unbounded at zero. Instead its singularity cancels the Gamma density:

\[
\mathbb E[|C_m|;V<m]
\le\frac1{(2m)!}\int_0^m e^{-v}(v+t)^{2m}dv
\le e^t\Pr(V\le m+t).
\]

The elementary Laplace-transform bound
`Pr(V <= a) <= exp(a) 2^(-(2m+1))` therefore implies

\[
\mathbb E[|C_m+1|;V<m]
\le \frac{e^{2t}+1}{2}(e/4)^m.
\]

The second summand in the prefactor accounts for `Pr(V<m)` and must not be
omitted. These estimates prove the claimed `O(m^(-1/2))` bound in `L1`.
Multiplication by `exp(it X_m)` preserves it, so

\[
\mathbb E[e^{itX_m}C_m]
=-\chi_m(t)+O(m^{-1/2}).
\]

The previously proved lower bound `chi_m(t)>1/6` converts this absolute
error into relative error for `m>=2`. In particular, `chi_1(2pi)=0` must
not be included in that bound. Only sufficiently large `n` are needed for
the asymptotic theorem.

## Factorial ratio, chain rule, and linear coefficient

The Hilbert determinant quotient is

\[
\frac{\Delta_m}{\Delta_{m+1}}
=\frac{(2m)!(2m+1)!}{m!^4}.
\]

Therefore

\[
\operatorname{Tr}(H_n^{-1}V_n)
=i(-1)^nK_n\frac{\chi_{n-1}(t)}{\chi_n(t)}
\bigl(1+O(n^{-1/2})\bigr),
\qquad
K_n=\frac{(2m)!^2(2m+1)!}{m!^4t^{2m+1}}.
\]

At the actual period, `F_n(Z)=G_n(Z,exp Z)` and
`P_n(Z)=G_n(Z,1)` give

\[
\frac{P_n'}{P_n}
=\frac{F_n'}{F_n}-\operatorname{Tr}(H_n^{-1}V_n).
\]

The omitted chain-rule term would be the dominant term, so it is essential
that the endpoint variable remain independent until this step. Independently,
differentiating the Jacobi integral shows `F_n'/F_n=O(n^2)`; the bound
`|chi_n'/chi_n| <= 6/sqrt(15)` also gives the more precise displayed formula
in the main note. This polynomial-size term is negligible relative to `K_n`.

Stirling's formula yields

\[
\log K_n
=2m\log m+(6\log2-2-2\log t)m+O(\log m).
\]

Using `t=2pi`, replacing `m` by `n-1`, and using
`1/6 < chi_(n-1)/chi_n < 6`, gives exactly

\[
\log\left|\frac{P_n'(\omega)}{P_n(\omega)}\right|
=2n\log n+[\log(16/\pi^2)-2]n+O(\log n).
\]

## Parity and the explicit eventual sign

Individual formal moments do not have pure parity: for example,
`mu_2(Z,1)=1/Z-2/Z^2`. The valid identity is the functional reflection

\[
\mathcal L_{-Z}(f(u))=\mathcal L_Z(f(1-u)),
\qquad \mathcal L_Z(u^k)=\mu_k(Z,1).
\]

One proof uses
`L_Z(f'+Zf)=f(1)-f(0)`, which uniquely determines the polynomial functional.
Replacing `f(u)` by `f(1-u)` proves the reflection relation. The reflected
monomial basis changes the moment matrix by `R M R^T`; `det(R)^2=1`.
Thus its determinant is unchanged, while the prefactor `Z^(n^2)` changes
by `(-1)^n`. Consequently `P_n(-Z)=(-1)^n P_n(Z)`, and the logarithmic
derivative at a nonzero imaginary argument is purely imaginary.

The explicit error certificate for `m>=10000` is valid. At that threshold
the four rational majorants in the main note sum to less than `1/12`;
each decreases with `m`. Together with `chi_m>1/6` this gives the stated
sign of the dominant endpoint term. The elementary bound

\[
K_n\ge (2m)!/7^{2m+1}
\ge (m/49)^m/7
>200^m/7>2(m+1)^2+24
\]

for `m>=10000` dominates the correction from `F_n'/F_n`. Thus the eventual
sign threshold `n>=10001` follows; it is not asserted to be sharp.

## Maximal disk conclusion

For every radius qualifying in the definition of `mathcal R_n`, Cauchy's
estimate and the established period lower bound give

\[
R\le\frac{2U_n}{|P_n'(\omega)|}
<\frac{12}{|P_n'(\omega)/P_n(\omega)|}.
\]

This applies to every qualifying radius, so it also bounds their supremum.
The old explicit disk has
`log r_n=-2n log n+O(n)` and satisfies `mathcal R_n>=r_n`. The two bounds
prove the claimed leading scale of `mathcal R_n`.

This upper bound concerns uniform smallness on the entire complex disk.
A large derivative alone does not establish a zero in that disk, or an
obstruction for other shapes, column families, or approximation methods.

## Executed checks

The reviewer read and reran `verify_period_log_derivative.py` after the two
drafting corrections above. It passed 27 exact endpoint insertion checks,
75 endpoint monomial integrations, 81 Gamma product normalization checks,
297 formal reflection moment checks, and 29 factorial-scale checks. The
explicit rational error majorant is `59983/930000 < 1/12`. These finite
checks validate formulas and constants; the uniform asymptotic is justified
by the analytic bounds reviewed above.


A separate [root-distance review](../nearest-root-scale/INDEPENDENT_NEAREST_ROOT_REVIEW.md) checks the additional finite-degree and bounded-disk arguments. The present review scopes its no-root conclusion to the derivative and Cauchy estimate by themselves.
