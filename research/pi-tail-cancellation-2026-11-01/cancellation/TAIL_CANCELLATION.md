# Exact cancellation of logarithmic tails in full-row determinants

2026-11-01. This note continues the sparse-shell checkpoint. Its main finite-jet identity is unconditional algebra. Using it in the original nonzero-determinant construction still requires the manuscript's interpolation theorem. No claim about π's bad approximability follows.

The new fact is stronger than preservation of rank: **the entire full-row determinant is unchanged when the logarithmic tails are removed at the weighted cutoff.** Consequently the tail-bearing summands of the period translation cancel exactly in aggregate. This improves both translation errors and arithmetic denominator costs, but leaves the sparse constant-translation terms intact.

## 1. The precise row packet and determinant identity

Fix positive rational weights `v₀,w₁,…,wₘ`, a positive `θ`, a cutoff `H>0`, and finitely many centers `cⱼ=(cⱼ₁,…,cⱼₘ)` over a field of characteristic zero. Write `wβ=Σᵢwᵢβᵢ` for a multi-index `β∈ℤ≥0ᵐ`. At every center include **all** rows

\[
\mathcal R_{H,j}=\{(j,s,\beta):s\in\mathbb Z_{\ge0},\quad
                 v_0s+w\beta/\theta<H\}.
\]

Let `L(t)=ln(1+t)` be its formal Taylor series. Let `τᵢ(t)` be any formal series with

\[
\operatorname{ord}_t\tau_i\ge T_i\ge1,
\qquad v_0T_i\ge w_i/\theta.                 \tag{1}
\]

Choose any fixed finite family of polynomial columns `P`. A full-row square minor uses every row in the union of the packets and exactly that many columns. Its row entries, with a formal parameter `z`, are

\[
A(z)_{(j,s,\beta),P}
 =[t^su^\beta]P(1+t,c_{j1}+L(t)+u_1+z\tau_1(t),\ldots,
                    c_{jm}+L(t)+u_m+z\tau_m(t)).
\]

Then

\[
\boxed{A(z)=U(z)A(0),\qquad \det U(z)=1.}     \tag{2}
\]

Here `A` may be rectangular. In particular **every full-row square minor, on its same chosen columns, is independent of z**. The same statement holds after arbitrary fixed column combinations or column scalings.

**Proof.** Expanding the transverse shift gives

\[
A(z)_{j,s,\beta}
 =\sum_{a\ge\beta}\binom a\beta z^{|a-\beta|}
   \sum_{\ell=0}^{s}
   [t^{s-\ell}]\prod_i\tau_i(t)^{a_i-\beta_i}
   A(0)_{j,\ell,a}.                          \tag{3}
\]

A nonzero coefficient with `d=a−β` obeys

\[
s-\ell\ge\sum_iT_id_i,
\quad
v_0\ell+wa/\theta
\le v_0s+w\beta/\theta
 -\sum_i(v_0T_i-w_i/\theta)d_i<H.             \tag{4}
\]

Thus every row on the right belongs to the **same original packet**. Equality in (1) causes no problem, because its row cutoff is strict. If `a≠β`, then `ℓ<s`; if `a=β`, the only contribution is `ℓ=s` with coefficient one. Order each center's rows by increasing `s`. The square row-change matrix `U(z)` is triangular with diagonal entries all one, proving (2).

One may alternatively take the exterior product of all rows. In a nonidentity choice from (3), at least one time order decreases and none increases. A bijection of the complete row set would preserve the sum of time orders, so this choice cannot be a bijection; some row repeats and its exterior product vanishes. This explains the determinant identity without estimating the size of the tails. ∎

The identity also holds with independent formal markers `zⱼᵢ` for every center and coordinate. Every coefficient of positive total tail degree, or any specified nonzero multidegree in these markers, is zero.

## 2. What cancels in the manuscript's exact period expansion

Set `cⱼᵢ=j rᵢ` and `ω=2πi`. For each chosen column define

\[
f_{a,P}(z)=[u^a]P(e^z,z+u_1,\ldots,z+u_m).
\]

At the **actual** period, `e^{hjω}=1` for integer monomial degree `h`. Put `εⱼᵢ=j(rᵢ−ω)`. Expanding `A(z)` about the common period center gives exactly the manuscript's row identity (3.13), except that every tail multiplicity vector `d` carries the extra marker `z^{|d|}`. This is an identity between finite coefficient polynomials; each original column has finite transverse degree.

By (2), the determinant is constant in the marker. Therefore the sum of all determinant summands with any fixed positive total tail degree is exactly zero. At `z=1`, the whole determinant equals its degree-zero part:

\[
\boxed{
\Delta_H
=\det\left(
\sum_{a\ge\beta}\binom a\beta
       \prod_i[j(r_i-\omega)]^{a_i-\beta_i}
       ([t^s]f_{a,P}(j\omega+\ln(1+t)))_P
\right)_{(j,s,\beta)}.}                     \tag{5}
\]

Only types supported by the chosen columns contribute. For monomials in the original domain `w₀h+wα≤H`, they satisfy `wa≤H`.

Thus one can perform the analytic translation using **only `d=k=0` and `ℓ=s`**. This is an exact regrouping of the full determinant, not an approximation dropping small terms. Individual tail-bearing summands after constant translation need not vanish separately; their sums in each marked degree vanish. The earlier sparse-shell construction already uses `d=k=0`, so it survives this cancellation.

When a formal period variable is useful, first use the real period identity above and then algebraize the resulting formula. One must not assume `e^{hjW}=1` at arbitrary formal `W`. The shift identity itself can also be written directly with `P(1+t,…)`, which contains no such exponential assertion.

## 3. Reduced translation error

Under the manuscript's approximation hypothesis

\[
|\pi-p_i/q_i|\le q_i^{-\nu},\qquad
w_i=\lceil\ln q_i\rceil,\qquad r_i=2ip_i/q_i,
\]

one has `|j(rᵢ−ω)|≤2K e^ν e^{-νwᵢ}`. In (5), using
`binom(a,β)≤2^{|a|}` and `|a|≤H/w_*`, where `w_*=minᵢwᵢ`, gives

\[
|\xi^0_{j,s,\beta,a}|
\le\exp\{-\nu w(a-\beta)+HE_{\rm tr}^{0}\},
\qquad
\boxed{E_{\rm tr}^{0}
=\frac{\ln2+\ln(2K)+\nu}{w_*}.}             \tag{6}
\]

Compared with the original translation bound, the `ν/F₀` term disappears, the scalar `ln2/v₀` term disappears, and the scalar binomial constant improves from `ln4` to `ln2`. There are at most `(⌊H⌋+1)^m` type choices per row when the transverse weights are positive integers, so the normalized logarithm of the number of choices still tends to zero for fixed parameters.

The collision estimate continues to have its own coefficient-extraction cost `ln2/v₀`; (6) does not remove that separate term. Combining (6) with the unchanged collision estimate gives

\[
E_{\rm an}^{0}
=\frac{100K}{w_0}+\frac{\ln2}{v_0}
 +\frac{\ln(800K^2)+\nu}{w_*}.               \tag{7}
\]

This reduces error terms. It does not increase the main scalar-and-collision score of a constant-translation selection.

## 4. Minimal truncation and improved arithmetic normalization

By (1), choose the cutoff without any auxiliary large `F₀`:

\[
T_i=\left\lceil\frac{w_i}{\theta v_0}\right\rceil,
\qquad
G_i(t)=\sum_{1\le k<T_i}\frac{(-1)^{k+1}t^k}{k},
\quad
L_i=\operatorname{lcm}(1,\ldots,T_i-1),       \tag{8}
\]

with the empty least common multiple equal to one. Every full-row determinant formed using `Gᵢ` equals the same-column determinant formed using `ln(1+t)` exactly.

For a monomial column `Y^hX^α`, scale that column by
`Qα=∏ᵢ(qᵢLᵢ)^{αᵢ}` and row `(j,s,β)` by `Qβ^{-1}`. The resulting entry is a Gaussian integer: each polynomial

\[
q_iL_i(jr_i+G_i(t))\in\mathbb Z[i][t],
\]

and the coefficient formula contains only products, binomial coefficients, and `(1+t)^h`. Hence a nonzero full-row minor satisfies the exact lower bound

\[
\ln|\Delta_H|
\ge-\sum_{\text{columns }\alpha}\sum_i\alpha_i\ln(q_iL_i)
   +\sum_{\text{rows }\beta}\sum_i\beta_i\ln(q_iL_i).      \tag{9}
\]

This simultaneous row-and-column normalization avoids multiplying every entry by a common product of all denominators.

Let `Λ=4 ln2` and `b̄_H=(MH)^{-1}Σ_rows wβ`. Using the manuscript's elementary bound `ln lcm(1,…,n)≤Λn` and

\[
T_i-1<\frac{w_i}{\theta v_0},\qquad
\ln(q_iL_i)\le w_i\left(1+\frac\Lambda{\theta v_0}\right),
\]

while the positive row term is at least `MH(b̄_H−θ/w_*)`, yields

\[
\boxed{
\frac{\ln|\Delta_H|}{MH}
\ge-(1-\bar b_H)-E_{\rm ar}^{0},
\qquad E_{\rm ar}^{0}
=\frac\Lambda{\theta v_0}+\frac\theta{w_*}.}               \tag{10}
\]

The improved error has **no factor `m`** and no `Σᵢ1/wᵢ` ceiling loss. Its proof uses the same arithmetic data as the manuscript and the exact identity, not a lower bound for π's approximation errors.

### Compatibility with the common mixed-column extension

For clarity, this uses the new `qᵢLᵢ` column scaling, not the earlier `qᵢ`-only normalization. Let `A` be the rational-center matrix, `B` a full-rank test matrix, and `Q` the diagonal monomial-column matrix above. The common bounded mixing lemma supplies an integer matrix `C` with entries in `{0,1,2}` such that both `AQC` and `BQC` are nonsingular. The scaled matrix `diag(Qβ^{-1})AQC` is Gaussian integral, so

\[
\frac{\ln|\det(AQC)|}{MH}\ge b_{\kappa,H},
\quad
b_{\kappa,H}=\frac1{MH}\sum_{\rm rows}\sum_i\beta_i\ln(q_iL_i).
\]

Define `κ_max=maxᵢ ln(qᵢLᵢ)/wᵢ`. The mixed-column analytic norm costs at most `κ_max+o(1)` at scale `MH`: each of the `N` raw columns is scaled by at most `exp(κ_max H)`, and the bounded sum contributes only `ln(2N)/H=o(1)`. Therefore the main arithmetic-versus-analytic threshold becomes

\[
\kappa_{\max}-b_{\kappa,H}
\le 1-\bar b_H+\frac\Lambda{\theta v_0}+\frac\theta{w_*}.  \tag{11}
\]

So the common mixed-column extension retains the same improved comparison bound. The exact tail identity applies to these fixed mixed columns as well. Rank existence and the common mixing lemma remain separate inputs.

## 5. Threshold and completeness are necessary hypotheses

A smallest below-threshold counterexample has `m=K=1`, `v₀=w₀=w₁=1`, `θ=1/2`, `H=2`. The full row packet consists of `(s,β)=(0,0),(1,0)`. Take columns `1,X` and any center `c`. Replacing `ln(1+t)` by `ln(1+t)+zt` changes its determinant from `1` to `1+z`. In particular the too-short truncation `G=0` changes it to zero. Here the start order `T=1` violates the required `T≥2`.

Even at valid threshold `T=2`, arbitrary sparse subpackets need not have the identity. At `H=3`, retain only the row `(s,β)=(2,0)` and column `X`. Its entry is `−1/2` with the logarithm, but zero with `G=t`. The full packet contains the missing lower-time transverse row that makes the complete determinant invariant.

## 6. Reproduction and remaining gap

Run `python3 verify_tail_cancellation.py`. The standard-library verifier checks four examples with 4, 12, 22, and 8 rows. It verifies the **entire formal matrix identity** coefficient by coefficient, reconstructs the constant-plus-tail period-reduced expansion independently, and checks the determinant polynomial at more exact points than its degree. The selected minor values are respectively `1`, `59049`, `105468750`, and `−625`. The weak threshold equality is exercised. It also checks both counterexamples and 552 supported entries after Gaussian arithmetic scaling, including 126 entries with nonzero imaginary part.

These finite examples certify the stated computations; the general identity rests on the proof above. They do not verify the unknown successive weight thresholds of the source interpolation theorem.

The resulting analytic comparison may replace `E_ar+E_an` by `E_ar^0+E_an^0` and sum only constant-translation types. However the previous sparse-shell family is already inside this surviving class. Its endpoint score tending to `0.3` remains below the comparison threshold tending to `0.5`. Thus the exact cancellation is a useful simplification and removes avoidable losses, but **does not close the endpoint argument or determine π's bad approximability**. Further work must address actual constant-offset magnitudes, their remaining coefficient sums, or stronger function estimates.
