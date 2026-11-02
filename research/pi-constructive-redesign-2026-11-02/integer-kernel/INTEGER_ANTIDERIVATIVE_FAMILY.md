# A constructive integer-antiderivative family

2026-11-02. This changes construction coordinates before integration. It does not continue the prime-pool exclusion route and does not claim a new pi approximation record.

## 1. Integer data with a dyadic boundary

Fix an integer c>=1. Choose integers r_0,...,r_(c-1),alpha and put

\[
R(x)=x\sum_{j=0}^{c-1}r_jx^{2j},
\]
\[
N(x)=(1+x^2)R'(x)-2(c-1)xR(x)+\alpha(1+x^2)^{c-1}.
\]

Then N is an even integer polynomial of degree at most 2c, and direct differentiation gives

\[
\frac{N(x)}{(1+x^2)^c}
=\frac{d}{dx}\frac{R(x)}{(1+x^2)^{c-1}}+\frac\alpha{1+x^2}.
\]

Consequently the following is an exact identity for every allowed integer input:

\[
\boxed{I:=\int_0^1\frac{N(x)}{(1+x^2)^c}\,dx
=\frac{R(1)}{2^{c-1}}+\frac\alpha4\pi.}
\]

The rational boundary now has only a power of two as denominator. This is an explicit design choice, not a claim that the odd arithmetic cost disappears.

For c>=3 and alpha!=0, set

\[
g=\gcd\bigl(R(1),2^{c-3}\alpha\bigr),\qquad
r=-\frac{R(1)}{2^{c-3}\alpha},\qquad
q=\frac{2^{c-3}|\alpha|}{g}.
\]

Here q is the fully reduced denominator of r. The exact quality identity is

\[
\boxed{q^2|\pi-r|=\frac{4q^2}{|\alpha|}|I|.}
\]

Thus an integral should be selected jointly with the integer period coefficient alpha, boundary value R(1), and actual gcd g. Making I small alone is insufficient.

## 2. Endpoint zeros are integer linear constraints

Write t=x^2 and N(x)=sum_(j=0)^c n_j t^j. With r_j=0 outside 0<=j<c,

\[
n_j=(2j+1)r_j+(2j-2c+1)r_{j-1}
+\alpha\binom{c-1}{j},
\]

where the binomial is zero for j=c. For fixed nonnegative integers a,b with a+b<=c,

\[
N(x)=x^{2a}(1-x^2)^bH(x^2),\qquad H\in\mathbb Z[t],
\]

is equivalent to the integer conditions

\[
n_j=0\quad(0\le j<a),\qquad
\sum_{j=k}^c\binom jk n_j=0\quad(0\le k<b).
\]

The quotient H is integral because division by t^a and the monic polynomial (t-1)^b preserves integer coefficients when exact. Its degree is at most d=c-a-b.

These a+b linear constraints are independent over Q: on the full polynomial coefficient space they are exactly divisibility by t^a(t-1)^b. The resulting rational solution space has dimension d+1. A Smith decomposition with its unimodular right transform supplies the full saturated integer kernel, rather than only separately cleared rational nullspace vectors.

## 3. The exact arithmetic cost of the coordinate choice

Let M_c map (r_0,...,r_(c-1),alpha) to the coefficient vector (n_0,...,n_c). It is nonsingular, and

\[
\boxed{|\det M_c|=2^{c-1}(c-1)!.}
\]

Therefore M_c Z^(c+1) is an explicitly known full-rank sublattice of index 2^(c-1)(c-1)! in all even integer numerators of degree at most 2c.

Proof: replace the R basis by R_j=x(1+x^2)^j for 0<=j<c, and the N basis by (1+x^2)^j for 0<=j<=c. Both basis changes are integral triangular matrices with diagonal one. For the differential part L_c,

\[
L_c R_j=2(c-j-1)(1+x^2)^j
+(2j-2c+3)(1+x^2)^{j+1}.
\]

The last R column is (1+x^2)^c, while the alpha column is (1+x^2)^(c-1). Removing those two columns and rows leaves a triangular block with diagonal 2(c-1),2(c-2),...,2. The determinant formula follows, including c=1 where the empty product is one.

This proves that pre-clearing the boundary is not a free denominator gain. Odd factors can reappear in alpha or in the size of the required integer numerator. The construction engine handles that cost explicitly instead of concealing it in an integration lcm.

## 4. A non-circular finite selection rule

For fixed c,a,b:

1. Build the exact integer endpoint constraint matrix.
2. Obtain its full integer kernel by Smith decomposition.
3. Map a candidate to values N(x)/(1+x^2)^c at the rational sample points x=1/10,...,9/10. Add explicit cost coordinates lambda*alpha and lambda*R(1), where lambda is a rational power of two.
4. Clear the feature denominators and apply exact integer LLL; retain the original integer coordinates in the feature vector so the result can be reconstructed and checked.
5. Discard alpha=0 and evaluate only nondegenerate fractions using the exact quality identity.

No value of pi, continued fraction, target rational or target numerator/denominator appears in steps 1-4. Pi is used only afterwards, through a rigorous Machin alternating-series enclosure, to certify quality. The sampled values are selection heuristics; they are not integral-error upper bounds.

The total-derivative subspace must not be counted as success. In particular alpha=0,R(1)=0 gives I=0 but no approximation to pi. The code rejects all alpha=0 rows, including such zero-integral rows. Allowing arbitrary unbounded integer coordinates would merely include a period lattice and does not itself supply a quality-zero theorem; the meaningful research question is whether endpoint constraints and a controlled coefficient budget produce a growing family with provably improving quality.

## 5. Completed evidence and limits

The initial rational-nullspace sublattice search through c=20 examined 1994 LLL rows, discarded 973 with alpha=0, and produced 448 distinct positive fractions. It recovered 355/113 with

\[
c=7,\quad a=1,\quad b=3,\quad\alpha=-339,
\]
\[
(r_0,...,r_6)=(339,1921,4473,5455,3560,1200,92),
\]

so R(1)=17040 and r=17040/(16*339)=355/113 after division by gcd 48. Its quality is rigorously between 0.00340631193013 and 0.00340631193015. This is a known baseline, not a new record.

The rank-one control a=b=n,c=2n gives primitive alpha values
2,33,135,52395,281295,... and demonstrates explicitly that odd factors have moved into the coefficient bill. Full tables, saturated searches, comparisons at equal c, and exact algebraic checks accompany this note. Their finite outputs do not prove any infinite-subsequence conclusion.

The full-kernel script includes a local correction to SymPy 1.14's LLL nearest-integer floor. In this runtime math.floor on a very large PythonMPQ rounds through a float, whereas numerator//denominator is exact. The correction leaves SymPy's integer LLL algorithm and its final exact size/Lovasz assertions intact. This numerical implementation detail is separately tested; it is not part of the mathematical construction.

## 6. A proved design rule: asymmetric endpoint orders reduce the odd period bill

This is a constructive parameter-selection result, rather than an exclusion theorem for the old recurrence.

Let a,b>=1 and c>=a+b+1. Let Lambda_(c,a,b) be the image of the full integer endpoint kernel under

\[
(r_0,\ldots,r_{c-1},\alpha)\longmapsto(\alpha,R(1)).
\]

For every odd prime p, its localization is exactly

\[
\boxed{\Lambda_{c,a,b}\otimes\mathbb Z_p
=p^{\lfloor\log_p M\rfloor}\mathbb Z_p\times\mathbb Z_p,
\qquad M=\max(2a-1,b).}
\]

Equivalently, the odd part of the covolume of this rank-two period image is exactly

\[
\boxed{L_{a,b}=\operatorname{odd}\operatorname{lcm}(1,\ldots,\max(2a-1,b)).}
\]

### Necessity at zero

Write F=R/(1+x^2)^(c-1). The endpoint condition gives

\[
F(x)=-\alpha\arctan x+O(x^{2a+1}).
\]

Multiplication by (1+x²)^(c-1) and its inverse are integral triangular transformations at zero. Thus integrality of the first a odd coefficients is equivalent to alpha/(2j+1) being p-integral for 0<=j<a. For odd p this requires exactly v_p(alpha)>=floor(log_p(2a-1)).

### Necessity at one

Put u=x-1 and z=R(1). The local condition is

\[
F(1+u)=\frac z{2^{c-1}}
-\alpha\bigl(\arctan(1+u)-\pi/4\bigr)+O(u^{b+1}).
\]

Write

\[
\frac1{2+2u+u^2}=\sum_{j\ge0}h_j u^j.
\]

The recurrence 2h_j+2h_(j-1)+h_(j-2)=0 shows that every h_j is zero or a signed power of two. In particular h_(k-1) is a nonzero dyadic unit when k is odd. The coefficient of u^k in the arctangent difference is h_(k-1)/k. Hence all these coefficients through degree b become p-integral exactly when v_p(alpha)>=floor(log_p b): sufficiency follows from the denominator k, and necessity from k=p^e. The constant term z/2^(c-1) is p-integral for every z in Z_p.

### Sufficiency with the degree bound

The two local requirements just described specify R modulo x^(2a+1) at zero and modulo (x-1)^(b+1) at one, with the odd parity built into R=xS(x²).

At zero they prescribe S(t) modulo t^a. At one the substitution t-1=2u+u² is a power-series isomorphism over Z_p because p is odd; division by x is also a unit. Thus they prescribe S(t) modulo (t-1)^(b+1). The two monic ideals are comaximal over Z_p. Hermite/Chinese-remainder interpolation therefore produces S in Z_p[t] with degree at most a+b. Since c>=a+b+1, this degree is at most c-1, exactly the permitted degree. No uncharged enlargement of R is needed.

This proves the full local image and not only a necessary divisibility condition. Localization commutes with this integer kernel and linear image, so the conclusion applies to the saturated integer period lattice.

### Fair finite parameter comparison

Keep the same pole order c=16 and the same total endpoint budget a+b=14. Exact Smith/Hermite normal forms give

\[
\Lambda_{16,7,7}=
45045\mathbb Z\times1024\mathbb Z,
\]
\[
\Lambda_{16,5,9}=
315\mathbb Z\times32768\mathbb Z.
\]

The odd period bill falls by 45045/315=143. The dyadic boundary bill grows by 32. The complete period-image covolume falls by the exact factor 143/32, from 46,126,080 to 10,321,920. This is a certified arithmetic improvement at equal c and a+b. It is not yet a factor-143/32 improvement in Q: the changed endpoint shape also changes the integral and coefficient norm.

In general, fixing s=a+b makes M=max(2a-1,b) smallest near b=2a. The mandatory odd period bill consequently has logarithm (2s/3)+o(s), compared with s+o(s) for balanced endpoints, by the prime number theorem. This identifies a parameter direction to pursue before searching for a good rational approximation. The full dyadic and analytic bills must still be included. The observed dyadic pattern for some b,c is not asserted as a universal formula.

## 7. A constructive upper certificate, available before evaluating pi

For s=a+b>0, c>=s, and H(t)=sum h_j t^j, put h=sum|h_j|. On 0<=t<=1 one has |H(t)|<=h. A direct one-variable maximum calculation gives

\[
\max_{0\le y\le1}\frac{y^a(1-y)^b}{(1+y)^s}
=\frac{a^a b^b}{2^a s^s},
\qquad y_{\max}=\frac a{a+2b}.
\]

It follows that every constructed nondegenerate integer vector satisfies the fully rational, unconditional upper certificate

\[
\boxed{Q\le\frac{4q^2h}{|\alpha|}
\frac{a^a b^b}{2^a(a+b)^{a+b}}.}
\]

This inequality can be used as a selection objective without knowing pi. It is an upper bound on achieved quality, unlike the old exclusion route's lower bounds. For signed H it may be quite weak because it discards cancellation; it should not be substituted for the actual quality interval.

For example, the known 355/113 vector above has H(t)=-9t+67t²-92t³, h=168. Its simple certificate is 21357/16, although its actual quality is about 0.00340631. Thus this certificate has not itself verified the good baseline. The separately supplied ../certificates/integral_certificate.py now gives an arbitrarily refinable signed-integral certificate; proving a growing sequence still requires coefficient/period control. The file constructive_upper_certificates.json records this distinction.

## 8. Completed equal-c construction controls

The saturated search completed **4<=c<=12**, d=1,2,3, with endpoint allocations near a:b=1:1 or 1:2. It examined 698 LLL output rows, rejected 279 alpha-zero rows, and evaluated 90 distinct positive rational fractions. No LLL assertion failed after the exact-floor correction. The larger experimental c<=32 run was stopped without a completion claim; only the stated completed range is evidence.

The following compares a search with free H against the best *constant-H* kernel over every positive endpoint allocation at the same c. It does not hold the individual endpoint orders fixed.

| c | Best constant-H control | Control Q (display) | Free-H candidate | Candidate Q (display) |
|---:|---:|---:|---:|---:|
| 4 | 47/15 | 1.85834706 | 22/7 | 0.06195997 |
| 6 | 424/135 | 15.52611167 | 333/106 | 0.93505573 |
| 7 | 4288/1365 | 373.97198484 | 355/113 | 0.00340631 |
| 12 | 458752/146025 | 136206.21335 | 100798/32085 | 9.31839894 |

These are actual improvements in the bounded design experiment, with known classical good fractions serving as baselines. They do **not** exhibit monotone improvement as c grows, and the larger denominator at c=12 has Q>1. In particular, no growing quality-zero sequence has been obtained. The selection rule, all full-kernel output statistics, and the per-c results are in antiderivative_small_saturated.py and small_saturated_results.json; the constant-H controls are in compare_single_kernels.py and single_kernel_control.json.

## 9. A fixed growing family, actually tested through m=6

The following family was specified before running its search:

\[
a=m,\qquad b=2m,\qquad c=4m+1,\qquad\deg H\le m+1,
\qquad m\ge1.
\]

It follows the proved odd-cost rule exactly: M=max(2a-1,b)=2m. Its odd period bill is odd lcm(1,...,2m), instead of the roughly 3m cutoff for balanced endpoints with the same total endpoint order 3m. The analytic shape bound is also explicit:

\[
\boxed{Q_m\le\frac{4q_m^2\|H_m\|_1}{|\alpha_m|}
\left(\frac2{27}\right)^m.}
\]

Producing integer vectors in this family for which the right-hand side tends to zero would be sufficient. No such coefficient-growth estimate has been proved; the bound is a usable construction objective, not an achieved asymptotic conclusion.

`new_family_search.py` uses the full saturated endpoint kernel and only five rational feature-cost weights, 2^(-k) for k in {0,8m,12m,16m,20m}. Each m has a hard 30-second runtime limit. All m=1,...,6 completed, taking about 8.53 seconds in total. There were 30 distinct-positive-fraction candidates across the six cases; all original integer R, alpha, N and H coordinates are preserved in new_family_results.json.

| m | Best observed candidate | Reduced q | Q (display only) |
|---:|---:|---:|---:|
| 1 | 22/7 | 7 | 0.0619599741 |
| 2 | 355/113 | 113 | 0.00340631193 |
| 3 | 355/113 | 113 | 0.00340631193 |
| 4 | 104348/33215 | 33215 | 0.3658638091 |
| 5 | 1459097/464445 | 464445 | 1.6202067974 |
| 6 | 85563208/27235615 | 27235615 | 0.6380658635 |

The values are supported by exact Machin intervals after candidate generation; the generation itself has no pi value, continued fraction or target rational as input. The root workflow can alternatively certify I directly from its integer integrand. These are concrete growing-denominator outputs of the redesigned family, not new best approximations of pi. They do not show monotone quality improvement, and m=5 fails Q<1.

The crude H-coefficient-norm bounds are much larger than the actual quality, e.g. about 2.46e8 at m=6. Thus the family has already generated usable finite approximants, but an efficient signed-integral bound and a provable joint coefficient/denominator estimate remain necessary for an infinite quality-zero construction.
