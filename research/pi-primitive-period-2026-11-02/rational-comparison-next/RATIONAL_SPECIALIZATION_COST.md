# What primitive smallness gives at a rational center

2026-11-02. The improved coefficient clearing makes the two-frequency actual-period determinant tend to zero. This note keeps the denominator introduced by rational evaluation and derives a concrete comparison bound. It explains the remaining gap to a quadratic Diophantine bound. It does not rule out a sharper use of the same family.

## 1. The even family is a primitive integer polynomial

Let n>=4 be even, write n=2r, and put

\[
A_j=\prod_{a=0}^{j-1}a!,\quad
d_n=\frac{A_n^2}{A_r^4},\quad N_n=\frac{n^2}{2},\quad D_n=n(n-1).
\]

The endpoint-content theorem proves that d_n is the exact least positive integer coefficient denominator of the genuine two-frequency source polynomial

\[
P_n(Z)=\det\bigl([z^s]e^{hz}(jZ+z)^a\bigr),
\quad j,h\in\{0,1\},\quad 0\le s,a<n.
\]

Define R_n(Z)=d_n P_n(Z). It has integer coefficients, degree N_n, leading coefficient (-1)^r, and only even powers. The last assertion follows from the commutator formula P_n(Z)=det(ET(Z)-T(Z)E): since T(-Z)=T(Z)^(-1) and det T(Z)=1, one has P_n(-Z)=(-1)^nP_n(Z). The integer polynomial R_n is primitive. Let omega=2*pi*i and

\[
U_n=(2\pi)^{n^2}\frac{A_n^2}{A_{2n}}.
\]

The preceding Fourier--Hankel theorem gives U_n/6<|P_n(omega)|<=U_n for every n>=2. The present note uses its upper bound and the new even-degree and denominator formulas. Degree N_n is distinct from the safe pre-cancellation degree surplus D_n.

## 2. Rational evaluation introduces its full degree denominator

Take a reduced rational approximation p/q, where p is an integer, q>=1, and gcd(p,q)=1. Write

\[
\frac{2p}{q}=\frac ab,\qquad
g=\gcd(2,q),\quad a=2p/g,\quad b=q/g.
\]

Then gcd(a,b)=1. Assume b>=2. Since R_n is even, R_n(ia/b) is a real rational number. Multiplying by b^N_n gives an integer whose residue modulo each prime divisor of b is its leading term

\[
(-1)^r i^{N_n}a^{N_n}.
\]

This residue is nonzero: N_n is even and a is coprime to b. Thus the resulting integer is coprime to b, and

\[
\boxed{\operatorname{den}_{\rm reduced} R_n(2ip/q)=b^{N_n},
\qquad |R_n(2ip/q)|\ge b^{-N_n}>0.}
\tag{1}
\]

This is stronger than merely assuming the rational-center determinant is nonzero. No zero at a nonintegral rational imaginary center is possible for this even, unit-leading-coefficient polynomial. It does not say that its value is an integer before multiplication by b^N_n.

The large factor b^N_n remains even though every coefficient of R_n was already cleared optimally. For example, when n=10 and p/q=355/113, the reduced value denominator is exactly 113^50. This is a value-denominator statement, not a claim about the approximation quality of 355/113.

## 3. An explicit derivative bound on the comparison segment

For |Z|<=8, each normalized matrix entry in the h=0 or h=1 column with spatial exponent a has absolute value at most 9^a. For j=1 it is a partial sum of

\[
\sum_{k=0}^{\min(s,a)}\binom ak\frac{h^{s-k}}{(s-k)!}Z^{a-k};
\]

the coefficients h^(s-k)/(s-k)! have magnitude at most one. Summing the corresponding full binomial expansion gives (|Z|+1)^a. The j=0 entries are constant in Z and satisfy the same bound. Differentiating a column gives the bound a*9^(a-1), with zero for a=0.

There are (2n)! determinant permutations and the total spatial column degree is D_n. The product rule therefore gives

\[
|P_n'(Z)|\le(2n)!D_n9^{D_n-1}
\le(2n)!10^{D_n},
\]

where the last inequality is one term of the binomial expansion of (9+1)^D_n. Hence

\[
\boxed{|R_n'(Z)|\le M_n:=d_n(2n)!10^{D_n}\quad(|Z|\le8).}
\tag{2}
\]

This is a deliberately explicit upper estimate. It is not asserted to be optimal near omega.

## 4. The resulting rigorous comparison

Choose any even n>=4 satisfying

\[
\boxed{d_nU_n\le\frac1{2b^{N_n}}.}
\tag{3}
\]

If |2ip/q-omega|<=1, the joining segment lies in |Z|<8, since 2*pi+1<8. Equations (1)--(3), the triangle inequality, and integration of the polynomial derivative along that segment imply

\[
b^{-N_n}\le |R_n(2ip/q)|
\le |R_n(\omega)|+M_n|2ip/q-\omega|.
\]

Consequently

\[
\boxed{\left|\pi-\frac pq\right|
\ge\frac1{4b^{N_n}d_n(2n)!10^{D_n}}.}
\tag{4}
\]

If the segment length is greater than one, the left side exceeds 1/2 and (4) still holds. Thus (4) applies to every reduced p/q with b>=2 whenever (3) holds, without a prior proximity assumption. The small-denominator cases b=1 can be treated directly and are irrelevant to the asymptotic comparison.

Condition (3) can be certified without evaluating pi numerically: replace 2*pi in U_n by 44/7 and compare positive rational numbers. Every fixed b admits a sufficiently large even n, by the new denominator asymptotic. A search using this exact sufficient condition is finite; the existence argument is analytic, not inferred from a search cutoff.

## 5. Height needed by this specific comparison

The exact denominator formula and the two-sided period theorem give

\[
\ln(d_n|P_n(\omega)|)
=-\tfrac12 n^2\ln n+
\left[\tfrac34+\ln(\pi/\sqrt2)\right]n^2+O(n\ln n).
\]

The same expansion holds with |P_n(omega)| replaced by U_n. Define the positive constant

\[
\kappa=\frac{e^{3/2}\pi^2}{2}.
\]

After division by N_n=n^2/2,

\[
\frac{\ln(d_nU_n)}{N_n}=-\ln n+\ln\kappa+o(1).
\tag{5}
\]

For any fixed c>kappa, choose n to be the least even integer at least cb. As b tends to infinity, condition (3) then holds, because

\[
\frac{\ln(b^{N_n}d_nU_n)}{N_n}\longrightarrow\ln(\kappa/c)<0.
\]

If 0<c<kappa, the analogous actual-value condition |R_n(omega)|<=1/(2b^N_n) fails eventually. This is the threshold for this **particular baseline-value comparison**, not a necessary dimension for all possible arguments.

With c>kappa, the right side of (4) has logarithm

\[
\boxed{-c^2 b^2\ln b+O_c(b^2).}
\tag{6}
\]

Indeed ln d_n=N_n ln n+O(n^2), ln((2n)!)=O(n ln n), and D_n ln10=O(n^2). The two dominant logarithms are N_n ln b and N_n ln n. Formula (6) describes the explicit bound produced by (2)--(4); it is not an asymptotic claim about the true approximation error.

Thus this direct comparison gives a positive bound of exp(-O(q^2 ln q)) scale. A bad-approximability conclusion would require a fixed positive constant times q^(-2), which is much stronger. Primitive integer polynomials becoming small at pi do not by themselves supply that improvement. Stronger local derivative control, a different rational-value comparison, or a construction whose effective degree grows more slowly could change the result; those possibilities have not been excluded.

## Verification

The companion verifier uses the previous exact n=4,6,8,10 coefficient certificates to check primitive coefficients, exact degree, rational-value denominators, and the derivative bound. It also finds exact sufficient heights for a finite set of reduced center denominators b by comparing the rational upper bound using 44/7. Large comparisons are summarized by sizes and certified inequalities rather than printing enormous integers. The general comparison theorem and asymptotic scale follow from the proofs above.
