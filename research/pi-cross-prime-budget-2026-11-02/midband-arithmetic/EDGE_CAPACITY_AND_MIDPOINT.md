# Edge-capacity audit and an exact midpoint interpretation

2026-11-02. This is an auxiliary audit of the fixed recurrence

\[
A_0=A_1=1,\qquad A_{r+1}=(2r+1)A_r+r^2 A_{r-1},
\quad \ell_r=A_r/r!.
\]

It does **not** prove the proposed eventual bound \(S(n)\le n/6\), nor a new bound on simultaneous anomalies across different primes.

## 1. Removing sublinear edge windows is elementary

Recall the already established reflection reduction

\[
k(n,p)=\min(n-p,2p-n-1),\qquad
p\mid A_{n-p}\iff p\mid A_{k(n,p)}
\]

for odd primes \(n/2<p<n\). Define \(S_{\rm edge}(n,K)\) by restricting \(S(n)\) to \(k(n,p)\le K\), with integer \(K\ge0\).

For every fixed \(k\), its primes belong to the two-element set

\[
\{n-k,\ (n+1+k)/2\}.
\]

Consequently, even if **every** candidate prime were anomalous,

\[
\boxed{S_{\rm edge}(n,K)\le 2(K+1)\log n.}
\tag{1}
\]

Thus \(K=o(n/\log n)\) already suffices for \(S_{\rm edge}=o(n)\), with no coefficient-height estimate or prime-distribution input.

Ordinary PNT improves the range further. Put \(\vartheta(x)=\sum_{p\le x}\log p\). The condition \(k(n,p)\le K\) places \(p\) in the union

\[
[n-K,n)\ \cup\ (n/2,(n+1+K)/2].
\]

Therefore, for \(0\le K<n/3\),

\[
S_{\rm edge}(n,K)\le
\vartheta(n-1)-\vartheta(n-K-1)
+\vartheta((n+1+K)/2)-\vartheta(n/2).
\tag{2}
\]

Using \(\vartheta(x)=x+o(x)\) separately at these four endpoints shows

\[
\boxed{K=o(n)\ \Longrightarrow\ S_{\rm edge}(n,K)=o(n).}
\tag{3}
\]

There is no short-interval asymptotic claim here: the errors are only required to be \(o(n)\), not \(o(K)\). This deduction is valid for each integer-valued function \(K(n)=o(n)\).

The previous height estimate

\[
S_{\rm edge}(n,K)\le \tfrac12\log(2+\sqrt2)K(K+1)
\]

remains valid. Its consequence for \(K=o(\sqrt n)\) is weaker than (1) and (3). The separate height-based *empty* logarithmic-width strips remain genuine arithmetic exclusions. Merely discarding a sublinear-width edge window does not provide new control over concentration among primes in the remaining interior.

## 2. The exact midpoint is a Frobenius-trace zero

Set \(h=(p-1)/2\) and

\[
U_r=2^r\ell_r=[x^r](x^2+2x+2)^r\in\mathbb Z.
\]

Consider the single fixed elliptic curve

\[
E:\quad y^2=x(x^2+2x+2).
\]

Its discriminant is \(-256\), so every odd prime has good reduction. Write \(a_p=p+1-\#E(\mathbb F_p)\). Euler's criterion and finite-field power sums give

\[
a_p\equiv [x^{p-1}]\{x(x^2+2x+2)\}^h=U_h\pmod p.
\tag{4}
\]

For completeness, \(a_p=-\sum_x\chi(x(x^2+2x+2))\). Replace \(\chi(t)\) by \(t^h\) modulo \(p\). In its polynomial expansion the only positive exponent divisible by \(p-1\) is \(p-1\), since the degree is \(3h<2(p-1)\). Its power sum is \(-1\); all other power sums vanish. This proves (4), including its sign.

Because \(2^h/h!\) is a \(p\)-unit and Hasse's bound gives \(|a_p|\le2\sqrt p<p\) for \(p\ge5\),

\[
\boxed{p\mid A_{(p-1)/2}\iff a_p=0\qquad(p\ge5).}
\tag{5}
\]

This is a genuine interpretation by one arithmetic object shared across primes, but it describes only the exact midpoint. In \(S(n)\), the midpoint condition requires

\[
n=(3p-1)/2,\qquad p=(2n+1)/3.
\]

It therefore contributes at most one prime and weight at most \(\log n\) at any fixed \(n\).

## 3. One elementary necessary midpoint congruence

Equation (5) implies the rigorously proved restriction

\[
\boxed{p\mid A_{(p-1)/2}\ \Longrightarrow\ p\equiv7\pmod8
\qquad(p\ge5).}
\tag{6}
\]

Indeed, if \(p\equiv1\pmod4\), then \(x^2+2x+2\) splits into distinct roots in \(\mathbb F_p\). The full subgroup of four points of order dividing two is rational, so \(4\mid\#E(\mathbb F_p)\). But \(a_p=0\) would give \(\#E=p+1\equiv2\pmod4\), a contradiction.

Thus \(p\equiv3\pmod4\), and \((0,0)\) is the unique nonzero rational two-torsion point. Since \(\#E=p+1\) is divisible by four, the finite abelian group \(E(\mathbb F_p)\) contains a point \(Q\) of order four with \(2Q=(0,0)\): a finite abelian 2-group with unique nonzero two-torsion is cyclic. The duplication formula is

\[
x(2Q)=\frac{(x(Q)^2-2)^2}{4y(Q)^2}.
\]

As \(y(Q)\ne0\), this forces \(x(Q)^2=2\). The supplementary law \((2/p)=(-1)^{(p^2-1)/8}\), combined with \(p\equiv3\pmod4\), gives \(p\equiv7\pmod8\).

The stronger observed class \(31\pmod{32}\) is **not proved here** and is not assumed. A midpoint congruence cannot be extrapolated to all degrees: the existing exact examples show anomalies in all four odd residue classes modulo eight.

## 4. Scope and verification

The remaining obstacle is still a uniform bound for the moving diagonal \(p\mid A_{n-p}\) across the primes in a macroscopic interior interval. The midpoint curve supplies no such bound. A shrinking central strip around \(p=(2n+1)/3\), just like the edge strips, has total prime weight \(o(n)\) by ordinary PNT regardless of anomalies. Removing such strips is an auxiliary simplification, not a solution of the global problem.

`verify_edge_midpoint.py` verifies the affine fibers and interval containment for all \(n\le2000\), and independently compares recurrence values with direct elliptic-curve point counts for all odd primes \(p\le2000\). It checks (4), (5), (6), and the exact duplication formula for all order-four candidates at these primes. These finite checks supplement the proofs and are not evidence for an unproved global \(S\) bound.

Sources consulted:

- Andrew Sutherland, MIT 18.783, Spring 2021, Problem Set 4, Problem 1: the Hasse-invariant coefficient and Frobenius-trace congruence. https://math.mit.edu/classes/18.783/2021/ProblemSet4.pdf
- Andrew Sutherland, MIT 18.783, Spring 2021, Lecture 7, Theorem 7.3: Hasse's bound. https://math.mit.edu/classes/18.783/2021/LectureNotes7.pdf
- NIST DLMF, Section 27.12: the prime number theorem; its equivalent form \(\vartheta(x)\sim x\) also follows by partial summation from \(\pi(x)\sim x/\log x\). https://dlmf.nist.gov/27.12

All specialized consequences above are direct deductions for the stated curve and recurrence; the cited sources are not claimed to prove the moving-diagonal estimate.
