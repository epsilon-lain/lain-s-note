# How the separated-weight interpolation construction works

**Scope and status.** This is a reconstruction and audit of the interpolation section of the pinned OpenAI017 source, rather than an acceptance of its downstream arithmetic conclusion. The source pin is `4eaee87548250b2682affda7d05d947c05554bae`; `build/sections/interpolation.tex` has Git blob SHA `69dd30eac65dd0817226970299e4100e15e9ec0b`, and `build/main.tex` has blob SHA `fe28ac589090c44a8877aa326b60e753c863fa24`, as recorded in `source_context/pi017/closure.json`. Exact elementary calculations below are independently checked. The asymptotic conclusion additionally uses the algebraic-geometric inputs identified below. This note does not audit the determinant estimates or infer their arithmetic conclusion.

Source links below all point to that pin:

- [Interpolation statement, lines 10–57](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/interpolation.tex#L10-L57).
- [Local comparison, lines 72–175](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/interpolation.tex#L72-L175).
- [Curve argument, lines 180–425](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/interpolation.tex#L180-L425).
- [Passage to packets, lines 437–584](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/sections/interpolation.tex#L437-L584).

## 1. The types, and the actual rank problem

The domain is the finite-dimensional polynomial space

\[
\mathcal P_W(H)=\operatorname{span}\{Y^hX^\alpha:w_0h+\sum_iw_i\alpha_i\le H\}.
\]

The target at each center is a quotient of a *formal local ring*, retaining precisely the monomials

\[
t^su^\beta,\qquad v_0s+\sum_i(w_i/\theta)\beta_i<H.
\]

The map is linear over \(\mathbb C\); it sends a polynomial to coefficients after
\(Y=1+t,\ X_i=c_{ji}+u_i+\log(1+t)\). It is not evaluation on an algebraic logarithm variety. The logarithm is a formal series at the fixed value \(Y=1\). Different centers receive independent packets, even though the same polynomial produces them.

The leading ratio of target dimension to domain dimension is
\(K(w_0/v_0)\theta^m\). The additional inequality \(K\theta^m<1\) is the analogous ratio *inside the fiber \(Y=1\)*. A positive coefficient surplus allows a homogeneous system of vanishing conditions to have a nonzero solution. It does not imply that the coefficient map is onto. That missing implication is the entire content of the geometric argument.

In particular, the proof legitimately uses dimension counting to manufacture a nonzero auxiliary polynomial with many zero coefficients. It does **not** use the same counting argument to manufacture an arbitrary prescribed packet.

## 2. Coordinate derivatives and logarithmic directions differ

Set
\[
D_0=Y\partial_Y+\sum_i\partial_{X_i},\qquad D_i=\partial_{X_i}.
\]

The chain rule gives, exactly,

\[
D_0=(1+t)\partial_t,\qquad D_i=\partial_{u_i},\qquad
\partial_Y=\partial_t-(1+t)^{-1}\sum_i\partial_{u_i}.
\]

Thus a coordinate normal basis involving \(\partial_Y\) is not automatically the same normal basis as one involving \(D_0\). The proof exploits this difference. The positive-index fields agree with their coordinate directions, while the zeroth field mixes the \(Y\) direction with every \(X_i\) direction.

With \(z=\log(1+t)\), the logarithmic flow is \((Y,X)=(e^z,c+u+z)\), and \(D_0=\partial_z\). In particular, in one transverse variable,

\[
[t^su^b]P(1+t,c+u+\log(1+t))
=\frac{1}{s!b!}\left[\prod_{r=0}^{s-1}(D_0-r)\right]D_1^bP(1,c).
\]

The falling factorial in \(D_0\), rather than \(D_0^s\), comes from using \(t=e^z-1\). The accompanying script verifies this identity directly by formal expansion. Since the change from \(t\) to \(z\) is triangular with invertible diagonal, it preserves finite jet rank when all the relevant one-variable orders are retained.

Each \(D_a\) preserves the polynomial degree budget \(W\): \(Y\partial_Y\) preserves monomial degree, and every \(\partial_{X_i}\) lowers it. In the local series it loses at most its assigned jet cost \(v_a\). These two properties are different and both are needed.

## 3. Two exact ablations: surplus is not rank

All ranks here use exact rational arithmetic in `interpolation_lab.py`; its results are in `interpolation_lab_results.json`.

### A. Ambient volume alone misses a fiber obstruction

Take \(m=1,K=3,W=(1,1),V=(2,2),\theta=1/2,H=3\), and centers \((1,0),(1,1),(1,2)\). The target consists of value and first \(t,u\) derivatives at each center. There are **10 columns, 9 rows, rank 7**.

This rank can also be derived without software. First logarithmic jets are an invertible change of ordinary first jets at each point. Write a cubic uniquely to first order along \(Y=1\) as

\[
P(Y,X)=f(X)+(Y-1)g(X)\pmod{(Y-1)^2},
\quad \deg f\le3,\ \deg g\le2.
\]

At three collinear points the image therefore has dimension \(4+3=7\). Equivalently, the kernel is exactly \((Y-1)^2\mathcal P_{(1,1)}(1)\), of dimension 3. The ambient volume ratio is \(3/4<1\), while the fiber ratio is \(3/2>1\). This explains why the source has a second volume hypothesis.

### B. Both volume inequalities still do not replace separation

Take \(m=K=1,w_0=1,v_0=1/2,\theta=1/4,w_1=1/10\), and center \((1,0)\). Thus \(V=(1/2,2/5)\). Both source volume inequalities hold, with ratios \(1/4\) and \(1/2\), but the coordinate weights are not successively separated.

| H | Domain columns | Target rows | Exact rank |
|---:|---:|---:|---:|
| 1 | 12 | 5 | 4 |
| 2 | 33 | 14 | 10 |
| 3 | 64 | 29 | 21 |
| 4 | 105 | 48 | 35 |
| 5 | 156 | 73 | 53 |

The rank loss persists asymptotically; it is not a small-degree artifact. On the algebraic curve \(C:X=0\), put \(t=Y-1\), so \(u=-\log(1+t)\). The weighted degree is \(\deg_WC=1\), while contact at the center is

\[
h=\min\{1/(1/2),1/(2/5)\}=2.
\]

Thus even the unstrengthened inequality \(\deg_WC\ge h\) fails. More explicitly, a monomial of target weight at least \(H\) has total degree at least \(2H\), so restriction to this curve induces a surjection from the target packet to \(\mathbb C[[t]]/(t^{2H})\). Every domain polynomial restricts to a polynomial in \(Y\) of degree at most \(H\), an image of dimension at most \(H+1\). Hence surjectivity is impossible for every integer \(H>1\). These are counterexamples to deleting a hypothesis, not counterexamples to the source theorem.

### A genuine separated example also needs a degree threshold

For the same \(m,K,w_0,v_0,\theta\), choose \(\sigma=1/10\). The source's three slack inequalities hold, and its explicit separation condition is \(w_1>30\). Choose \(w_1=31\), so \(v_1=124\). For \(H<124\), only the \(u^0\) rows occur. The domain functions on the logarithmic curve are \(e^{hz}z^a\), with \(h+31a\le H\).

At \(H=50\) there are 100 rows and only 71 columns, so surjectivity fails. At \(H=100\) there are 200 rows and 218 columns, and rank is 200. The script certifies both ranks modulo the prime 1,000,003. A nonzero minor modulo a prime proves a nonzero rational minor. Entries are the integers

\[
\left.\frac{d^n}{dz^n}(e^{hz}z^a)\right|_{z=0}
=\begin{cases}n!h^{n-a}/(n-a)!,&n\ge a,\\0,&n<a.\end{cases}
\]

This is also the confluent Vandermonde mechanism: for multiplicities \(r_h\), its square determinant is, up to sign,
\(\prod_h\prod_{a<r_h}a!\,\prod_{h<k}(k-h)^{r_hr_k}\ne0\).
The example explains both why logarithmic restriction can create substantial rank and why “sufficiently large \(H\)” cannot be discarded. It does not test general transverse packets or multiple centers.

## 4. How the proof converts a rank obstruction into an impossible curve

The source proves the stronger curve inequality

\[
\deg_W C\ge(1+\sigma)\sum_{P\mapsto a_j}h_P.
\]

Here degree counts the largest weighted coordinate pole at each point of the complete curve, and \(h_P\) is the minimum weighted vanishing order of \(t,u_1,\ldots,u_m\) on a branch. This is designed so that:

- a polynomial of degree at most \(N\) has total pole order at most \(N\deg_WC\);
- a local series of weighted order at least \(M\) has branch order at least \(Mh_P\).

Assume a curve violates the inequality. Coefficient surplus produces \(F_N\ne0\) of degree at most \(N\) with local order at least \((1+3\sigma)N\) at every center. Derivatives of total cost at most \(\sigma N\) retain order at least \((1+2\sigma)N\). On the violating curve their total zero order exceeds their total possible pole order. The divisor identity on a complete curve forces every such derivative to vanish identically on the curve.

This is stronger than saying that one auxiliary polynomial contains the curve. There is a whole derivative ideal containing it. Form nested common-zero sets using cost budgets \(r\sigma N/(m+2)\). A nested choice of irreducible components containing the curve has too many levels for its possible dimensions. At adjacent equal levels, a component \(Z\) persists under an extra positive derivative budget. That is the input for the local multiplicity lemma.

### The local inequality compares two independently chosen normal bases

For a codimension-\(k\) persistent component, choose a coordinate normal basis \(A\) and a commuting-field normal basis \(B\). On a transverse slice, derivative persistence forces all monomials of \(B\)-cost below \(\varepsilon N\) to vanish in the equations. Thus local quotient length is at least

\[
\frac{(\varepsilon N)^k}{2k!\prod_{b\in B}v_b}.
\]

Generic linear combinations select \(k\) equations with an isolated zero. Weighted Bézout bounds their local intersection length above by \(N^k/\prod_{a\in A}w_a\). Comparing gives

\[
\prod_{a\in A}w_a\le C_*\prod_{b\in B}v_b.
\]

Crucially, this is uniform in \(N\), the coefficients, the component, and the analytic flow map. The lower bound counts monomials, so no radius or analytic norm appears. The upper bound is obtained by finite power substitutions reducing weighted degree to ordinary degree. The statement would not be usable here with a threshold depending on the moving component \(Z=Z_N\).

### Why successively separated scales force a logarithmic differential identity

The weights are chosen so that the reverse strict inequality holds whenever the **largest positive index where \(A,B\) differ belongs to \(A\)**. There are finitely many comparisons. Higher common indices contribute \(w_j/v_j=\theta\), so choosing \(w_i\) only requires earlier weights and fixed data.

Now let \(i\) be the largest positive coordinate direction nontangent to \(Z\). Some coordinate normal basis contains \(i\). If a frame normal basis omitted \(i\), the local inequality would contradict the separated product inequality. Consequently every frame normal basis contains \(D_i\).

Deleting \(D_i\) leaves precisely the hyperplane \(\ker(dX_i-dY/Y)\). Its failure to span the normal quotient says that the tangent space of \(Z\) lies in this hyperplane. Hence

\[
dX_i=dY/Y\quad\text{on }Z.
\]

This is where logarithms matter. If \(Y\) were nonconstant, restrict to an algebraic curve in \(Z\) on which it stays nonconstant and pass to a smooth complete model. At a zero or pole of \(Y\), \(dY/Y\) has the nonzero residue \(\operatorname{ord}_P Y\). Every exact meromorphic differential \(dX_i\) has zero residue, because differentiating a Laurent series cannot create a \(z^{-1}dz\) term. Contradiction. Thus \(Y\) is constant. If all positive coordinate directions were tangent initially, properness already forces \(dY=0\).

The source then repeats the derivative argument inside \(Y=1\), now using the second volume inequality. Coordinate and derivative frames coincide there. Separation forces a unique normal-basis index set, hence at least one coordinate \(X_i\) is constant on a persistent component. Coordinatewise distinctness of the centers makes it meet at most one center. On a nonconstant coordinate \(X_l\) of the curve, zero–pole comparison gives the sharper final bound
\(\deg_WC\ge\theta^{-1}\sum h_P>(1+\sigma)\sum h_P\), a contradiction.

The coordinatewise condition also has an elementary obstruction behind it. If \(m=2,K=4,\theta=2/5,w_0=v_0\), both volume ratios are \(16/25<1\). Put four distinct centers on the line \(Y=1,X_2=0\), varying only \(X_1\). Whatever the separation of the weights, this line has degree \(1/w_1\) and total contact \(4\theta/w_1=8/(5w_1)\), so the desired curve inequality fails. Distinct points alone would therefore not support the source theorem. The restriction to that line also gives an asymptotic rank obstruction by ordinary one-variable Hermite interpolation.

## 5. Why a curve inequality yields all coefficient packets

The final passage is an actual global generation argument, not another coefficient count.

Choose a divisible \(R\) and compactify affine space by every monomial of degree at most \(R\). With \(L=\mathcal O_X(1)\), pure powers ensure \(\deg(L|_C)=R\deg_WC\). Truncate each logarithm beyond the relevant local weight. The resulting algebraic coordinates \(z_0=t,z_i=u'_i\) differ from the formal coordinates by a filtration-preserving automorphism, so every packet and every branch contact is unchanged.

At each center set \(I_j=(z_0^{R/v_0},\ldots,z_m^{R/v_m})\), and blow up the resulting finite-support ideal \(I\). For \(I\mathcal O_{X'}=\mathcal O_{X'}(-E)\) and \(A=p^*L\), the curve inequality gives \(A-(1+\sigma)E\) nef. Exceptional curves have positive degree for \(-E\); curves away from the centers need no contact correction. Relative ampleness gives some ample \(aA-E\), and a positive linear combination makes \(A-E\) ample.

For sufficiently large \(n\), ordinary Rees-algebra recovery and relative vanishing identify \(p_*\mathcal O_{X'}(-nE)=I^n\) and kill its higher direct images. Serre vanishing for the ample \(A-E\) then gives
\(H^1(X,I^n\otimes L^n)=0\). The exact sequence therefore supplies arbitrary local data modulo \(I^n\).

The inclusion is in the important direction

\[
I^n\subseteq\{\text{series of weight at least }nR\}.
\]

The quotient modulo \(I^n\) is the larger packet, so surjectivity onto it implies surjectivity onto the desired smaller packet. Equality is unnecessary: for example \((t^2,u^2)\) omits \(tu\), even though \(tu\) has weight 2. Finally, eventual surjectivity from ambient homogeneous polynomials to \(H^0(X,L^n)\) turns the sections into actual polynomials of weighted degree at most \(nR\).

## 6. Audit boundary and parameter order

I found no specific defect in this interpolation argument under the stated hypotheses. The elementary substitutions, exact-rank ablations, the zero–pole implication, the normal-basis linear algebra, and the residue obstruction have been reconstructed rather than merely repeated. The finite experiments establish only their stated finite ranks; they cannot establish the asymptotic theorem.

The remaining theorem depends on standard but substantive inputs: isolated-local-multiplicity Bézout, existence of complete smooth curve models, blowup and relative ampleness properties, numerical nefness and nef-plus-ample, eventual recovery of ordinary Rees powers, higher-direct-image vanishing, and Serre vanishing. The source explicitly uses *ordinary powers* only eventually; no assumption of projective normality, smooth compactification, or equality with integral closures is needed in this reconstruction. I have checked the roles and internal hypotheses of these invocations, not independently reproved all those external theorems.

The quantifier order is essential:

\[
(m,K,w_0,v_0,\theta)\ \longrightarrow\ \sigma,C_*\
\longrightarrow\ w_1\ \longrightarrow\cdots\longrightarrow\ w_m\
\longrightarrow\ (c_{ji})\ \longrightarrow\ H\to\infty.
\]

The first choices involve finite numerical comparisons and cannot depend on the centers. The centers may affect the later ampleness/vanishing threshold for \(H\). Reversing these dependencies would be circular in an application where center locations are built from approximation denominators whose sizes also determine the coordinate weights; see [the main-file explanation, lines 229–239](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/main.tex#L229-L239). The construction is useful precisely because it controls this order, in addition to producing rank.
