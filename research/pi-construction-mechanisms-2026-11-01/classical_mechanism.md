# Constructing rational obstructions to approximating pi

Date: 2026-11-01. This note reconstructs finite mathematical mechanisms and compares them with OpenAI017. It does not certify the theorem `mu(pi)=2`. Exact symbolic checks below establish only the displayed finite identities. Reading coverage and access gaps are recorded in `classical_coverage.json`.

## What a construction has to accomplish

An accurate decimal for pi supplies one rational approximation. An upper irrationality bound must instead obstruct **every** sufficiently large rational approximation that is too accurate. The usual strategy builds quantities that have two incompatible properties: arithmetic discreteness prevents a nonzero quantity from becoming arbitrarily small, while cancellation or integration makes it very small analytically. Construction, nonvanishing, and size are three separate obligations.

Pi fits naturally into this strategy because the exponential map has an algebraic value at an imaginary multiple of pi:

\[
 e^{i\pi/2}=i,\qquad e^{2\pi i}=1.
\]

The first relation lets one specialize integer polynomials at a Gaussian integer while retaining pi as a logarithm. The second lets one evaluate the *same entire functions* at many different lifts of the same algebraic point. The logarithm near 1 supplies local coordinates; its branches are encoded by adding periods to the lifted variable. No assertion that all these periods are values of one principal branch is needed.

## Classical source mechanisms, with exact locations

**Mahler (1953).** Original pp.30–34, §§1–6, start from integer polynomial coefficients `A_hk(x)` whose matrix is invertible and whose combinations with powers of `log x` are small. Specialize at `x=i`, replace `log i` by `ip/(2q)`, and select a nonzero combination using invertibility. Its denominator divides `(2q)^m`, so its modulus is at least `(2q)^(-m)`. A difference quotient converts this into a bound on `|pi-p/q|`; the degree is chosen in proportion to `log q`. The paper proves exponent 42 for all `q>=2`, and states exponent 30 eventually without giving that optimization. Its initial approximation lemma is imported from Mahler's separate logarithms paper; this note has not read that paper. Original pp.34–42 extend the arithmetic step using conjugates and apply the result to sine denominators. [Full primary scan, original p.30](https://ems.press/content/book-chapter-files/27418?nt=1#page=2), [p.34](https://ems.press/content/book-chapter-files/27418?nt=1#page=6), [p.41](https://ems.press/content/book-chapter-files/27418?nt=1#page=13).

**Mignotte (1974).** Original pp.122–123, I.1–I.2, explicitly construct the remainders as residues of

\[
 \frac{P}{2\pi i}\oint_C
 \frac{z^h x^{z+n}}{(z(z+1)\cdots(z+n))^{m+1}}\,dz,
 \quad
 P=m!\operatorname{lcm}(1,\ldots,n)^m(n!)^{m+1}.
\]

Pole multiplicities produce powers of `log x`; decay at infinity forces a high zero at `x=1`. A degree calculation proves `det(A_hk)=c(x-1)^((m+1)n)` with `c!=0`. Original pp.123–125 repeat the logarithm substitution, now with `m=5` and sharper remainder/coefficient estimates. Theorem 1 actually gives exponent **20.6** for every `q>=2`, then exponent **20** for `q>=exp(e^110)`. OpenAI017's “21” is a valid weakening, not the exact theorem. The finite-denominator completion uses an external computation of 21,000 continued-fraction convergents, while explicit prime estimates enter the denominator bounds. Appendices pp.129–131 analyze contours and coefficients. [Construction](https://www.numdam.org/item/10.24033/msmf.139.pdf#page=3), [specialization](https://www.numdam.org/item/10.24033/msmf.139.pdf#page=5), [Theorem 1](https://www.numdam.org/item/10.24033/msmf.139.pdf#page=6).

**Zudilin (2004), arXiv v2.** This is an author's exposition of constructions attributed to Hata, Rukhadze, and Rhin, not a substitute for checking those original proofs. Pages 2–5 isolate the rate comparison for integer linear forms and arithmetic improvements from factorial valuations. Pages 6–9, §2, integrate one rational function to different endpoints: its `z^(-1)` coefficient is the same for every endpoint, so the logarithmic forms share a coefficient. For endpoints `2` and `1+i`, their logarithms combine to yield pi. Saddle points control both remainders and coefficient growth; prime divisibility reduces the denominator multiplier. Theorem 2 gives `8.01604539...` for nonzero rational combinations of `log 2` and pi. Pages 10–13 generalize the numerator to integer polynomials, making analytic optimization and denominator compatibility explicit. Pages 14–16 discuss hypergeometric directions and references. [Rate criteria](https://arxiv.org/pdf/math/0404523v2#page=2), [simultaneous construction](https://arxiv.org/pdf/math/0404523v2#page=6), [arithmetic](https://arxiv.org/pdf/math/0404523v2#page=8).

**Chudnovsky (1982), access boundary.** Only original pp.299–300 were available through the publisher's two-page preview. The introduction describes Hermite's exponential approximation system, inverse Laplace transforms, and explicit asymptotics. Pages 301–322 were not read. Thus this note does not attribute a detailed construction, theorem proof, or numerical bound to that chapter. [Publisher](https://link.springer.com/chapter/10.1007/BFb0093516), [two-page preview](https://page-one.springer.com/pdf/preview/10.1007/BFb0093516).

## Own finite model 1: manufacture a high-order zero

Take the four-dimensional space of exponential polynomials

\[
 F(z)=a+bz+(c+dz)e^z.
\]

Demand that its coefficients of `1,z,z^2` vanish. These are three rational linear conditions:

\[
 a+c=0,\qquad b+c+d=0,\qquad c/2+d=0.
\]

One integer solution is

\[
 F(z)=(z-2)e^z+z+2
     =z^3/6+z^4/12+z^5/40+O(z^6).
\]

This is a first-type Hermite–Padé construction: a polynomial-coefficient linear combination of the functions `1,e^z` has a prescribed high-order zero. The cancellation is solved exactly before any numerical evaluation. More generally, coefficients of polynomials multiplying `1,e^z,...,e^(kz)` become unknowns in a rational linear system obtained by matching Taylor coefficients.

A dimension count finds a nonzero coefficient vector when there are more unknowns than homogeneous conditions. It does not by itself control coefficient heights, prove a useful remainder estimate at a distant point, or ensure an evaluation is nonzero. These are the remaining arithmetic and analytic tasks.

Changing variables to `x=e^z` turns this example into

\[
 R_0(x)=(x+1)\log x+2(1-x).
\]

The companion

\[
 R_1(x)=x-1-\log x
\]

comes from the same residue denominator in the case `m=n=1` of the displayed contour family. We independently computed its two residues at `z=0,-1`; no high-degree theorem is used in this calculation. The coefficient matrix is

\[
 A(x)=\begin{pmatrix}2(1-x)&x+1\\x-1&-1\end{pmatrix},
 \qquad \det A(x)=-(x-1)^2.
\]

At `x=i`, the determinant is `2i`. Therefore, for **every** complex `y`, the two numbers

\[
 U_0(y)=2(1-i)+(i+1)y,\qquad U_1(y)=i-1-y
\]

cannot both vanish: otherwise the invertible matrix would kill the nonzero vector `(1,y)`. At `y=ip/(2q)`, at least one `2q U_h(y)` is a nonzero Gaussian integer, giving

\[
 |U_h(y)|\ge\frac1{2q}.
\]

This is the precise role of the determinant: it makes the lattice lower bound applicable. The tiny example does not have a remainder small enough at `x=i` to prove a useful irrationality bound. Raising the approximation order improves analytic cancellation but also increases coefficient and denominator costs.

## Own finite model 2: why simultaneous approximation needs a common denominator

Let

\[
 f_+(x)=\log(1+x),\qquad f_-(x)=\log(1-x).
\]

Seek a common normalized denominator `q(x)=1+q_1x+q_2x^2`, and degree-two numerators `p_+,p_-`, with

\[
 q(x)f_\pm(x)-p_\pm(x)=O(x^4).
\]

Coefficients of degree zero through two define each numerator. The two degree-three conditions are

\[
 1/3-q_1/2+q_2=0,\qquad
 -1/3-q_1/2-q_2=0.
\]

They force `q_1=0`, `q_2=-1/3`. Clearing denominators gives the integer polynomials

\[
 Q(x)=6-2x^2,\qquad P_\pm(x)=\pm6x-3x^2,
\]

and the exactly checked remainders

\[
 R_\pm(x)=Q(x)\log(1\pm x)-P_\pm(x)
 =-x^4/2\pm8x^5/15-x^6/2+O(x^7).
\]

These are second-type simultaneous Padé conditions: two approximation requirements share one denominator. Evaluating the **exact definitions** at `x=i` and using the principal logarithms gives

\[
 \log(1+i)-\log(1-i)=i\pi/2,
 \qquad R_+(i)-R_-(i)=4i(\pi-3).
\]

Thus cancellation between two logarithms extracts a linear form in pi, even though each logarithm also contains `log 2`. This finite example yields the rough rational approximant 3. The `O(x^4)` statement is local and supplies no small numerical bound at `x=i` on its own; the exact identity remains valid there by the analytic definitions. Higher-dimensional constructions require an actual remainder analysis rather than substituting a boundary point into a local `O` estimate.

## Own rate calculation: denominators are part of the approximation

Suppose a raw rational-coefficient construction has remainder size `exp(-alpha n+o(n))`, coefficient size at most `exp(beta n+o(n))`, and requires a clearing multiplier `exp(delta n+o(n))`. After clearing, the useful integer form has rates

\[
 C_0=\alpha-\delta>0,\qquad C_1=\beta+\delta.
\]

The multiplier both weakens decay and enlarges height. A usual rate criterion leads to

\[
 \mu\le1+\frac{\beta+\delta}{\alpha-\delta},
 \qquad
 \frac{d}{d\delta}\left(1+\frac{\beta+\delta}{\alpha-\delta}\right)
 =\frac{\alpha+\beta}{(\alpha-\delta)^2}>0.
\]

Here is the elementary reason behind the criterion for ordinary integer forms. Write `L_n=b_n gamma-a_n`, and take `n` just large enough that `q |L_n|<=1/2`, so asymptotically `n~log q/C_0`. If `q a_n-p b_n!=0`, integrality yields

\[
 |b_n|\,|\gamma-p/q|\ge1/q-|L_n|\ge1/(2q).
\]

If that integer is zero, then `|gamma-p/q|=|L_n|/|b_n|`; an exact exponential decay rate also supplies a lower bound for `|L_n|`. Allowing an arbitrarily small error in both exponential rates yields exponent `1+C_1/C_0` in either case. A mere upper bound on the remainder, without nonvanishing or a lower rate, does not justify this second case.

This calculation explains why optimizing the maximum of an integrand alone is insufficient: a numerator producing better decay may introduce worse denominator arithmetic. If `alpha<=delta`, the cleared forms do not even tend exponentially to zero.

## Own finite model 3: repeated Taylor degrees kill determinant terms

For three evaluation rows, take

\[
 M(\varepsilon)=
 \begin{pmatrix}
 1&0&1\\
 1&\varepsilon&e^{\varepsilon}\\
 1&2\varepsilon&e^{2\varepsilon}
 \end{pmatrix}.
\]

Direct symbolic calculation gives

\[
 \det M(\varepsilon)
 =\varepsilon(e^{\varepsilon}-1)^2
 =\varepsilon^3+O(\varepsilon^4).
\]

The exponent 3 is not obtained by multiplying independent estimates of the nine entries. Expand the vector `(1,z,e^z)` in Taylor coefficient vectors. Any determinant summand choosing the same Taylor degree in two rows has two proportional rows and vanishes. Surviving degrees must be at least `0,1,2`, whose sum is 3.

For `n` rows drawn from one coefficient-vector family, the corresponding minimum sum is `n(n-1)/2`. On a smaller disk inside a Cauchy contour, every degree costs a fixed geometric factor, turning this sum into a quadratic saving. If rows belong to different families indexed by `a`, the repeated pair must be `(a,degree)`; matching degrees in different families need not cause cancellation. Consequently the saving is proportional to `sum_a n_a^2`, not automatically to the square of the total number of rows.

## What OpenAI017 adds to the mechanisms above

The following describes the argument's organization, conditional on its interpolation theorem. It is a comparison with the sources actually read, not a priority claim over all literature.

| Obligation | Classical finite models above | OpenAI017 implementation |
|---|---|---|
| Encode pi arithmetically | `log i=i pi/2`, or combine two algebraic logarithms | Many rational centers `c_ji=2ij p_i/q_i` close to common exact periods `2pi ij` |
| Produce nonvanishing | Explicit coefficient determinant | Surjectivity of a weighted jet map gives a nonzero square minor using every row |
| Preserve rational entries | Integer coefficients and residue denominators | Truncate `log(1+t)` by rational polynomials; preserve the required weighted coefficient space |
| Arithmetic lower bound | A nonzero Gaussian integer has modulus at least 1 | Column clearing costs `q_i^(alpha_i)` and row scaling returns `q_i^(beta_i)`; main normalized cost is `1-bbar` |
| Analytic upper bound | A high-order zero or a small integral | Translate to exact periods; each row becomes a sum of tests on entire functions `e^(hz) z^d` with an extra transverse index |
| Combine many constraints | Shared coefficients/denominators | A dichotomy: many low-index rows force Taylor collisions; many high-index rows contribute large powers of approximation errors |
| Choose parameters | Balance approximation order against `log q` | Fix dimension and geometric constants, then separated denominator weights and centers, and finally send degree `H` to infinity |

The new burden in this comparison is the strong, **center-independent separated-weight interpolation statement**, together with the high-dimensional balance that lets its two analytic alternatives beat the arithmetic cost for every fixed `nu>2`. It is not the use of exponentials, logarithms, determinants, or arithmetic discreteness individually.

In the determinant section the two alternatives are especially concrete. If at least a fraction `eta` of `M` rows have transverse weight at most `AH`, there are only `D_A(H)` possible transverse indices. Cauchy–Schwarz gives `sum n_a^2 >= (eta M)^2/D_A(H)`. Otherwise, most chosen indices exceed `AH`, and the assumed approximation quality supplies a factor of roughly `exp[-nu MH(A(1-eta)-bbar)]`. The relevant main arithmetic cost is `exp[-MH(1-bbar)]`. Keeping the row rebate `bbar` matters in the comparison.

The paper's parameter choice uses `A^2<theta` alongside `nu(A-theta)>1-theta`. With `K~C^m`, `w_0=B^(-m)`, and `v_0=2K theta^m w_0`, it arranges errors tending to zero while the normalized collision term grows as `(B/A)^m/(m+1)`. Then one **finite** sufficiently large dimension is fixed; the degree limit is taken afterwards. Interchanging those orders would require additional uniform estimates.

The exact rank claim cannot be replaced by coefficient counting. Likewise, the finite checks in this note do not establish the geometric interpolation theorem, its center uniformity, or the passage from its hypotheses to the exponent-two conclusion. Those remain the mathematically decisive parts of a full audit.

## Verification and boundaries

`check_classical_models.py` uses exact symbolic arithmetic to verify: the first-type cancellation order; both `m=n=1` residues; the determinant and its value at `i`; the common denominator obtained from the simultaneous linear conditions; the two remainder expansions; the three-row collision identity; and the derivative of the denominator penalty. The recorded result is `classical_model_checks.json`.

Mahler, Mignotte, and Zudilin were read through their complete extracted bodies. Selected Mignotte formulas were checked against rendered pages because the scan's text layer corrupts superscripts. Reading a published proof is distinguished here from independently proving its estimates or reproducing its external computations. No historical numerical bound was independently recomputed in full. Niven's primary PDF endpoints returned access errors in this subtask; no Niven coverage is claimed here. Chudnovsky remains limited to the two-page preview. The precise coverage manifest records these boundaries and the local comparison sections actually read.
