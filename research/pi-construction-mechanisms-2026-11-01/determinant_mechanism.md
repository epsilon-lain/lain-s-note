# π construction laboratory: arithmetic height, exact periods, and Taylor collisions

2026-11-01. This is an independent finite model and a conditional audit of the determinant stage of `pi017`, pinned at `4eaee87548250b2682affda7d05d947c05554bae`. It does not establish the separated-weight interpolation theorem, the claimed value of μ(π), or a badly approximable endpoint. When the original source cache or pinned submodule is available, `period_lab.py` verifies every Git blob SHA against `pinned_source_manifest.json`. Without local originals it explicitly reports the byte audit unavailable and still runs the independent finite model. Its complete exact output is `results_period_lab.json`.

Reading scope for this note: the full bodies of `build/main.tex`, `build/sections/determinant.tex`, and `build/sections/conclusion.tex`; the statement of the interpolation theorem and its conventions in `build/sections/interpolation.tex`; and the previous finite-rank scripts and rank note. No full-reading ledger was changed.

## 1. What makes this a construction for π

The construction puts two structures on one determinant. Its entries at the approximate centers belong to Q(i), so denominators give a lower bound whenever the determinant is nonzero. After an exact row expansion, the same determinant is expressed using entire functions tested at the exact logarithmic periods `jω`, where `ω=2πi`. The identity `exp(hjω)=1`, for integer `h,j`, makes these entire-function tests match the original `Y=1+t` coordinate.

The approximate centers are `(Y,X₁,…,Xₘ)=(1,jr₁,…,jrₘ)`, with `rᵢ=2ipᵢ/qᵢ`. They are not points `(exp(jr),jr,…,jr)` on an exponential curve. Replacing ω by a near period inside an exponential would introduce a phase `exp(hjr)` which is absent from the actual algebraic matrix. Exact translation followed by transverse displacement is what avoids that error.

The columns share one weighted budget:

\[
P=Y^hX^\alpha,\qquad w_0h+w\alpha\le H.
\]

The entire functions attached to a transverse index `a` are

\[
f_{a,P}(z)=[u^a]P(e^z,z+u_1,\ldots,z+u_m)
=\binom\alpha a e^{hz}z^{|\alpha|-|a|},
\]

with zero entry if `a` is not componentwise bounded by `α`. The exponential contributes many polynomial-column directions, while its growth is controlled by `h≤H/w₀`. The transverse index identifies which functions are being tested; grouping rows by that index is essential.

## 2. The arithmetic lower bound, derived on the actual minor

Fix a nonzero square minor using all `M` rows. This nonvanishing must be supplied separately; denominator clearing cannot produce it. For its columns and rows, multiply column `α` by `∏qᵢ^{αᵢ}` and row `(j,s,β)` by `∏qᵢ^{-βᵢ}`. An entry becomes

\[
\prod_i\binom{\alpha_i}{\beta_i}
[t^s](1+t)^h\prod_i\bigl(2ijp_i+q_iG_i(t)\bigr)^{\alpha_i-\beta_i}.
\]

This is zero unless `α≥β`. Put `Lᵢ=lcm(1,…,Tᵢ−1)` and

\[
D_H=\prod_iL_i^{\lfloor H/w_i\rfloor}.
\]

Every scaled entry has denominator dividing `D_H`. Thus the determinant of the matrix obtained by multiplying each entry by `D_H` is a nonzero Gaussian integer. Its squared modulus is a positive integer, so its modulus is at least 1. Undoing all scalings gives the exact inequality

\[
\log|\Delta_H|\ge-M\log D_H
-\sum_{\rm columns}\alpha\cdot\log q
+\sum_{\rm rows}\beta\cdot\log q.\tag{1}
\]

The signs in (1) matter: transverse row derivatives recover denominator factors that were charged to the columns. With `wᵢ=ceil(log qᵢ)`, each selected column costs at most `H`, whereas

\[
\sum_{\rm rows}\beta\cdot\log q
\ge MH\bar b-MH\theta/w_*.
\]

The elementary bound `log lcm(1,…,n)≤(4 log 2)n` gives

\[
\frac{\log D_H}{H}
\le\frac{(4\log2)F_0m}{v_0}
+(4\log2)\sum_i\frac1{w_i}.
\]

Consequently

\[
\frac{\log|\Delta_H|}{MH}
\ge -(1-\bar b)-E_{\rm ar}.
\]

This mechanism needs neither an upper bound on the integer numerators nor any random-matrix heuristic. It needs a nonzero determinant and an integral structure after controlled scaling. A huge numerator can only increase the modulus of a nonzero Gaussian integer.

The fixed truncations help control that scaling. For a fixed jet packet the coefficients of the full logarithm are already rational, but letting their denominator complexity grow freely with `H` would not yield the stated small uniform arithmetic error. The fixed `Tᵢ` keep the denominator cost explicit. Their omitted tails are recovered on the analytic side, with a separate budget.

## 3. A new exact 6×6 minor, including the exponential column

Take the six columns

\[
(1,Y,X_2,X_2^2,X_2^3,X_1)
\]

and the six rows `(j,s,β)`, with `j=0,1`, `s=0,1,2`, and `β=(0,0)`. For all these jets `Gᵢ=t−t²/2` is sufficient. The actual coefficient matrix is

\[
A(r_1,r_2)=
\begin{pmatrix}
1&1&0&0&0&0\\
0&1&1&0&0&1\\
0&0&-1/2&1&0&-1/2\\
1&1&r_2&r_2^2&r_2^3&r_1\\
0&1&1&2r_2&3r_2^2&1\\
0&0&-1/2&1-r_2&-3r_2^2/2+3r_2&-1/2
\end{pmatrix}.
\]

Exact elimination, a separate 720-term Leibniz expansion, and independent differentiation of the logarithmic expressions give

\[
\boxed{\det A=-3r_2^2(r_1-r_2).}\tag{2}
\]

The factor `r₁−r₂` has a transparent source: subtracting the `X₂` column from the `X₁` column leaves `(r₁−r₂)(0,0,0,1,0,0)^T`. Thus the original matrix has rank 6 when `r₂≠0` and `r₁≠r₂`.

This is an actual whole-row-packet minor of the source's coefficient matrix for finite data

\[
m=2,\quad K=2,\quad H=15,\quad w_0=v_0=5,
\quad(w_1,w_2)=(2,5),\quad\theta=1/10.
\]

The row inequality excludes every nonzero `β` and gives exactly `s=0,1,2`. Every selected column lies in the weighted column simplex. Taking `F₀=21` gives `T=(9,21)`, whose higher terms do not change these jets. No assertion is made that these small weights satisfy the source's unknown successive separation thresholds; explicit nonvanishing makes that unnecessary for this finite example.

For a concrete arithmetic instance choose `r₁=44i/7`, `r₂=710i/113`. Then

\[
\det A=\frac{3024600i}{10100279}.
\]

The column multipliers are `(1,1,113,113²,113³,7)`. Multiplying every entry of this scaled matrix by 2 already makes it Gaussian integral, with determinant

\[
279307921036800i\ne0.
\]

The source's universal factors are `L=(840,232792560)` and `D_H=840⁷·232792560³`; they also clear all entries, and are much less sharp for this small packet. The example illustrates the lower-bound mechanism; its approximation quality and the asymptotic interpolation theorem are not certified by the computation.

### Three ranks that must stay separate

On these same six columns:

| Matrix or term | Exact rank | Meaning |
|---|---:|---|
| Original `A(r₁,r₂)` over `Q(r₁,r₂)` | 6 | Algebraic coefficient matrix at approximate centers |
| Exact-period translated choice `a=0` in all rows | 5 | One term of the row expansion, with `X₁=X₂` along the diagonal |
| A choice with `a=(1,0)` in row `(j,s)=(1,0)` and `a=0` elsewhere | 6 | One shifted term needed in the exact reconstruction |

For the shifted choice the determinant is `−3Z²`, where `Z` is a formal placeholder for ω after the exponential-period identity is applied. Its scalar is `e₁=r₁−ω`. The corresponding `a=(0,1)` choice has determinant `+3Z²` and scalar `e₂=r₂−ω`.

The original rank is not inferred from a translated term's rank, nor does the original rank force every translated term to have full rank. Equation (2) and exact row translation are separate statements about the same chosen columns.

## 4. Exact translation: errors, logarithmic tails, and cancellation

Write `εᵢ=j(rᵢ−ω)` and `τᵢ=Gᵢ−log(1+t)`. The exact coefficient-row identity is

\[
\begin{aligned}
\operatorname{row}_{j,s,\beta}
={}&\sum_{a\ge\beta}\binom a\beta
\sum_{0\le d\le a-\beta}\binom{a-\beta}{d}
\epsilon^{a-\beta-d}
\sum_{k=0}^s[t^k]\tau^d\\
&\hspace{12mm}\cdot
\bigl([t^{s-k}]f_{a,P}(j\omega+\log(1+t))\bigr)_P.
\end{aligned}\tag{3}
\]

The coefficients outside the last row vector do not depend on the column. That is why determinant multilinearity applies to (3). Expanding each entry with a different choice depending on its column would not justify a rowwise determinant estimate.

For the 6×6 model the tails vanish to the tested order. The three rows at `j=0` have zero displacement and need only `a=0`. For each of the three rows at `j=1`, the possible nonzero transverse indices are

\[
(0,0),(0,1),(0,2),(0,3),(1,0).
\]

The script enumerates all `5³=125` choices on the same minor. There are 16 nonzero vector determinants. Grouping their full contributions by total transverse degree gives:

| Total degree | Sum of contributions |
|---:|---|
| 1 | `−3Z²(e₁−e₂)` |
| 2 | `−6Ze₂(e₁−e₂)` |
| 3 | `−3e₂²(e₁−e₂)` |
| 4 | `0` |

Their exact sum is

\[
-3(Z+e_2)^2(e_1-e_2),
\]

which reconstructs (2) after `rᵢ=Z+eᵢ`. This gives two different cancellations. The all-zero-transverse choice vanishes because its vector rows have rank 5. There are also nonzero individual degree-4 terms whose contributions cancel after summation. If `e₁=e₂`, every degree sum cancels, although some individual translated terms are nonzero. A surviving term alone therefore cannot prove nonvanishing of the original determinant.

The source obtains its upper bound by taking absolute values before summing translated choices, so it does not depend on these additional cancellations. Certifying zero vector determinants can nevertheless improve a future upper bound because those terms can be removed before taking absolute values.

The lab also verifies (3) in 70 cases with nonzero `β` and different genuine truncation tails `T=(2,3)`. It exercises 268 nonzero tail summands. For example `τ₁=t²/2−t³/3+…`, and `τ₂=−t³/3+…`. A further 300 entries are independently checked by differentiation of the entire-function tests, and all 36 actual 6×6 entries are independently checked by differentiation of the true logarithm.

### Why a high index pays an approximation-error cost

The bad-approximation assumption gives

\[
|\epsilon_i|\le2Kq_i^{-\nu}
\le2Ke^\nu e^{-\nu w_i}.
\]

Thus increasing `a−β` normally gives the gain `−νw(a−β)` in the scalar exponent. A tail factor can replace a displacement factor, but it vanishes to order at least `Tᵢ`, so any nonzero tail term satisfies

\[
\sum_i d_iw_i\le(v_0/F_0)\sum_i d_iT_i
\le v_0s/F_0<H/F_0.
\]

It can erase at most the budget `νH/F₀`. The coefficient estimate on `|t|=1/2`, binomial factors, and rounding `log qᵢ` to `wᵢ` produce the stated remaining translation errors. High transverse indices are therefore not a free way to escape Taylor collisions.

## 5. Why equal transverse indices force distinct Taylor degrees

Fix the columns once. On `|z|≤R`, write the row of functions for index `a` as

\[
f_{a,P}(z)=\sum_{d\ge0}c_{a,d,P}(z/R)^d.
\]

A translated row is a linear functional `Lρ` applied to these functions:

\[
L_\rho(f_{a,P})=\sum_dL_\rho((z/R)^d)c_{a,d,P}.
\]

In a determinant summand, if two rows choose the same pair `(a,d)`, their coefficient vectors `(c_{a,d,P})_P` are identical; the remaining functional values are only row scalars. The summand is zero. This is an exact alternating-multilinear cancellation, independent of any numerical closeness of centers.

Within a group of `nₐ` rows sharing `a`, every nonzero Taylor summand therefore uses `nₐ` distinct nonnegative integers. Their least possible sum is

\[
0+1+\cdots+(n_a-1)=\frac{n_a(n_a-1)}2.\tag{4}
\]

The quadratic term comes from the sum of these forced distinct degrees. It is not a dimension count for the original algebraic matrix and is not a statement that the original rows themselves coincide.

### Exact 5×5 collision model

Use `m=1`, columns `(1,X,X²,X³,Y)`, and tests

\[
(j,\ell)=(0,0),(0,1),(1,0),(1,1),(2,0),\qquad a=0.
\]

At exact periods the matrix is

\[
B(Z)=\begin{pmatrix}
1&0&0&0&1\\
0&1&0&0&1\\
1&Z&Z^2&Z^3&1\\
0&1&2Z&3Z^2&1\\
1&2Z&4Z^2&8Z^3&1
\end{pmatrix},\qquad
\det B=-6Z^5.
\]

It has rank 5 over `Q(Z)`. Before applying the functionals, its Taylor coefficient vector at degree `d` is

\[
c_d=(\mathbf1_{d=0},\mathbf1_{d=1},\mathbf1_{d=2},\mathbf1_{d=3},1/d!).
\]

Repeating any degree repeats a whole vector and kills its coefficient determinant. The first permitted degree set is `(0,1,2,3,4)`, with sum 10 and coefficient determinant `1/24`. Its moment determinant for these tests is `4Z⁸`, so this one Taylor selection contributes `Z⁸/6`.

This does **not** say that the exact-period determinant must have a factor `Z¹⁰`, or even `Z⁸`. The collision estimate is a bound using geometric decay at a fixed radius. It is not a vanishing-order assertion in a formal period variable. The full infinite sum, followed by the identity `exp(Z)=1` at nonzero exact periods, gives the displayed polynomial `−6Z⁵`.

One can see the distinction exactly. At general centers `j x`, keep `E=exp(x)` independent before specialization. The determinant is

\[
-x^4(-E^2+4Ex-4E+2x+5).
\]

Setting `E=1` gives the exact-period polynomial `−6x⁵`. Following the analytic path `E=exp(x)` near `x=0` instead gives

\[
x^8/6+2x^9/15+O(x^{10}).
\]

Those are different operations. The degree-10 Taylor selection, combined with two first-derivative tests, accounts for the latter order 8. Applying exact-period algebra while varying the period continuously near zero would be invalid.

### Equal degrees in different transverse groups need not cancel

Take columns `(X₁,X₂,YX₁,YX₂)` and coefficient labels

\[
((1,0),0),\ ((0,1),0),\ ((1,0),1),\ ((0,1),1).
\]

Their coefficient matrix is

\[
\begin{pmatrix}
1&0&1&0\\
0&1&0&1\\
0&0&1&0\\
0&0&0&1
\end{pmatrix},\qquad\det=1.
\]

Degrees 0 and 1 are both repeated globally. They are distinct within each of the two groups. The correct least total degree is `binom(2,2)+binom(2,2)=2`; requiring four globally distinct degrees would incorrectly give 6. The collision object is `(a,d)`, not `d` alone.

## 6. From collision cancellation to an exponential upper bound

The source sets `R=100K`. Cauchy's estimate on the entire functions gives

\[
|c_{a,d,P}|\le\mathcal D_H
=\exp\!\left(H\left(\frac R{w_0}+\frac{\log(2R)}{w_*}\right)\right).
\]

For `|t|=1/2`, every exact period test lies inside `|z|≤R/2`. Hence

\[
|[t^\ell]((j\omega+\log(1+t))/R)^d|\le2^\ell2^{-d}.
\]

Split `2^{-Σd}` into two halves. Equation (4) bounds the first half by

\[
2^{-\frac12\sum_a\binom{n_a}{2}}
=\exp\!\left(-\frac{\log2}{4}\sum_an_a^2+\frac{\log2}{4}M\right).
\]

The other half is summable independently in every row. The coefficient determinant is bounded by `M!𝒟_H^M`, and `Σℓ≤MH/v₀`. After normalization by `MH`, `log M!/MH` and the geometric-sum constants tend to zero because all parameters are fixed and `M` is polynomial in `H`. Thus

\[
\log|\det T|\le-c\sum_an_a^2+MH(E_{\rm hol}+o(1)),
\quad c=(\log2)/4.
\]

The calculation uses only equal-vector cancellation, a radius gap, and controlled growth. It is uniform over the originally selected columns and every translated row choice. No numerical rank estimate is used.

## 7. Why the two alternatives cover every translated choice

Choose a transverse cutoff `AH` and a row fraction `η`. Every outer determinant-expansion term falls into one of two cases.

If at least `ηM` indices have `wa≤AH`, those indices occupy at most

\[
D_A(H)\sim\frac{(AH)^m}{m!\prod_iw_i}
\]

groups. Cauchy–Schwarz forces

\[
\sum_an_a^2\ge\frac{(\eta M)^2}{D_A(H)}
\sim MH\frac{\eta^2K\theta^m}{(m+1)v_0A^m}.
\]

The high collision saving pays the arithmetic cost. In this case the approximation-error saving is nonnegative and may simply be discarded.

If fewer than `ηM` indices lie below the cutoff, more than `(1−η)M` lie above it. Their total transverse weight is at least `MHA(1−η)`, so the product of translation scalars supplies the normalized saving

\[
\nu(A(1-\eta)-\bar b).
\]

In this case the collision saving is nonnegative and may be discarded. This is a packing-versus-displacement principle: using few small groups causes many collisions; avoiding those groups requires large transverse shifts and therefore powers of the rational-approximation errors.

The row count has dimension `m+1`, while the transverse-group count has dimension `m`. Their ratio is what creates a saving of order `MH`, rather than only order `M`. Keeping the full column simplex is also what bounds the arithmetic column cost by exactly `H` per column. Replacing either simplex by an unrelated box requires redoing both comparisons.

The normalized upper bound is consequently

\[
\frac{\log|\Delta_H|}{MH}
\le E_{\rm an}+o(1)+
\max\{-cL_H,-\nu(A(1-\eta)-\bar b)\}.
\]

The large-index gap over the arithmetic main cost is

\[
\nu(A(1-\eta)-\bar b)-(1-\bar b)
\ge\nu(A(1-\eta)-\theta)-(1-\theta)=g.
\]

Thus the source's three determinant-comparison conditions are sufficient: `g>0`, `E_ar+E_an<g`, and `cL>1+E_ar+E_an`. This implication is conditional on the same nonzero full-row minor and all fixed-parameter hypotheses. The estimates do not independently establish those hypotheses.

## 8. Removing pieces: what the finite experiments establish

Collapsing the exponential factor to 1 on the 6×6 minor replaces its `Y` column by the `1` column. The original rank falls from 6 to 5; its all-zero-transverse rank falls from 5 to 4. This shows the exponential contributes an independent direction in this specific construction. It does not prove that no other polynomial-only construction could use a different set of columns.

Retaining only `a=0` in the row translation leaves the rank-5 exact-period term and loses the actual nonzero determinant. The shifted family restores rank-6 terms. Conversely, allowing several transverse groups reduces collision pressure: the 4×4 coefficient example has repeated degrees across groups and determinant 1. The analytic dichotomy is needed precisely because the transverse family both restores possible rank and provides an escape from large single-group collisions.

These two effects should not be conflated. Original nonvanishing concerns `A(r₁,…,rₘ)`. Collision cancellation concerns coefficient determinants inside a fixed translated choice. Exact translation connects them, and possible cancellation among different choices remains a third issue.

## 9. Audit and limits

The denominator-clearing signs, the binomial factors in exact row translation, the tail-budget inequality, the groupwise `nₐ(nₐ−1)/2` saving, the constants arising from splitting geometric decay, and the two-alternative comparison were checked independently. No local formula error was found in these determinant-stage implications. The source's parameter-limit algebra is compatible with these formulas when its interpolation input is assumed.

This audit does not verify the geometric interpolation proof or its independence of centers. The finite examples do not supply a uniform lower bound for asymptotic ranks, a uniform minimum surviving displacement cost, or a new irrationality bound. Rejecting the zero all-transverse term in one minor does not authorize rejecting it in every asymptotic minor. For arbitrary polynomial minors, the general rule transferring nonvanishing in `Q[Z]` to nonvanishing at the actual `2πi` uses its independent transcendence. The displayed monomial determinants need only `Z≠0`. The lab keeps `Z` formal and does not replace it by a numerical approximation.

Even a valid proof of `μ(π)=2` would give an eventual bound `|π−p/q|≥q^{−(2+ε)}` separately for each `ε>0`. It would not give a single `c>0` with `|π−p/q|≥c/q²`, bounded partial quotients, or bad approximability. Nothing here proves any of those endpoints.

Reproduction: `python research/pi-construction-mechanisms-2026-11-01/period_lab.py`. Only exact integer, rational, Gaussian-rational, polynomial, symbolic derivative, and symbolic series calculations are used.
