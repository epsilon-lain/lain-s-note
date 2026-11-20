# High-order formal-group cancellation with integral binomial weights

Date: 2026-11-20  
Status: **proved formal-group theorem**. It strengthens the real–\(p\)-adic simultaneous approximation construction to arbitrary prescribed order, but does not yet give a useful height bound for approximation to \(\pi\).

## 1. Formal multiplication polynomials

Let \(F(X,Y)\in R[[X,Y]]\) be a one-dimensional commutative formal group over a torsion-free characteristic-zero ring \(R\), with
\[
F(X,Y)=X+Y+\text{terms of total degree at least }2.
\]
Write its multiplication-by-\(k\) series as
\[
[k]_F(T)=kT+c_2(k)T^2+c_3(k)T^3+\cdots.
\]

**Lemma.** For every \(r\ge2\), \(c_r(k)\) is a polynomial in \(k\) of degree at most \(r\), and
\[
c_r(0)=c_r(1)=0.
\]

**Proof.** The identities \([0]_F(T)=0\) and \([1]_F(T)=T\) give the two vanishing values. For the degree bound, use
\[
[k+1]_F(T)=F([k]_F(T),T).
\]
Assume inductively that the coefficient of \(T^s\) in \([k]_F(T)\) is polynomial in \(k\) of degree at most \(s\) for all \(s<r\). In the coefficient of \(T^r\), every nonlinear monomial of \(F([k](T),T)\) contains at least one explicit factor \(T\). Its remaining coefficient has degree at most \(r-1\) in \(k\). Hence
\[
c_r(k+1)-c_r(k)
\]
is a polynomial of degree at most \(r-1\). A discrete antiderivative has degree at most \(r\), proving the claim. \(\square\)

## 2. An integer interpolation identity

For \(m\ge1\), define integer weights
\[
\boxed{
\lambda_{m,j}=(-1)^j\binom{m+1}{j-1},
\qquad 2\le j\le m+2.
}
\tag{1}
\]
For every polynomial \(f\) of degree at most \(m\),
\[
\boxed{
\sum_{j=2}^{m+2}\lambda_{m,j}f(j)=f(1).
}
\tag{2}
\]
Indeed, this is the \((m+1)\)-st forward-difference identity applied at \(1\), with the \(j=1\) term isolated. In particular,
\[
\sum_{j=2}^{m+2}\lambda_{m,j}j=1.
\tag{3}
\]

Apply (2) coefficientwise to \([j]_F(T)\). The linear coefficient \(j\) is sent to \(1\) by (3). For every degree \(2\le r\le m\), the coefficient polynomial \(c_r(j)\) has degree at most \(r\le m\) and vanishes at \(j=1\), so its weighted sum is zero. Therefore
\[
\boxed{
\sum_{j=2}^{m+2}\lambda_{m,j}[j]_F(T)
=T+O(T^{m+1}).
}
\tag{4}
\]
After division by \(T\),
\[
\boxed{
\sum_{j=2}^{m+2}\lambda_{m,j}
\frac{[j]_F(T)}{T}
=1+O(T^m).
}
\tag{5}
\]
This is an arbitrary-order cancellation identity with **integer** coefficients. No denominators are introduced by the interpolation weights.

## 3. Apply it to the elliptic formal group

Return to the elliptic curve \(E/\mathbb Q\) and formal parameter \(\xi\) from the Pasten–Vidaux construction. Let \(V\in E(\mathbb Q)\) be a non-torsion point in the chosen formal neighbourhood, and put \(z=\xi(V)\ne0\). Define
\[
\boxed{
B_m(V)=
\sum_{j=2}^{m+2}
(-1)^j\binom{m+1}{j-1}
\frac{\xi([j]V)}{\xi(V)}.
}
\tag{6}
\]
Every summand is rational, so \(B_m(V)\in\mathbb Q\). Equation (5) gives
\[
\boxed{
B_m(V)=1+O_{E,m}(z^m)
\qquad (z\to0\text{ in }\mathbb R).
}
\tag{7}
\]
The implied constant depends on \(E\) and \(m\), and no uniform-in-\(m\) bound is claimed.

Now suppose \(V\) was chosen so that
\[
z=\xi(V)\in M\mathbb Z_p
\qquad\text{for every }p\mid M.
\]
For each fixed positive integer \(j\), the formal-group expansion gives
\[
\frac{\xi([j]V)}{\xi(V)}
=j+zF_j(z),\qquad F_j(T)\in\mathbb Z_p[[T]].
\]
Hence every summand is congruent to \(j\pmod{M\mathbb Z_p}\). Using (3),
\[
\boxed{
B_m(V)-1\in M\mathbb Z_p
\qquad\text{for every }p\mid M.
}
\tag{8}
\]
Thus the same rational number is arbitrarily close to \(1\) in the real metric to order \(m\), while satisfying all prescribed local congruences.

For any rational target \(r=a/b\) with \(p\nmid b\) for all \(p\mid M\), the rational
\[
\beta_{m,r}(V)=rB_m(V)
\]
satisfies
\[
|\beta_{m,r}(V)-r|=O_{E,m,r}(|z|^m),
\qquad
\beta_{m,r}(V)-r\in M\mathbb Z_p\quad(p\mid M).
\tag{9}
\]
This is a real–\(p\)-adic simultaneous approximation statement for a prescribed rational target.

## 4. Exact check in the multiplicative formal group

For the multiplicative formal group
\[
F(X,Y)=X+Y+XY,
\qquad [j]_F(T)=(1+T)^j-1,
\]
the identity can be checked explicitly:
\[
B_m(T)=
\sum_{j=2}^{m+2}
(-1)^j\binom{m+1}{j-1}
\frac{(1+T)^j-1}{T}.
\]
The forward-difference identity gives
\[
\boxed{
B_m(T)=1+(-1)^m T^m+O(T^{m+1}).
}
\tag{10}
\]
So the order \(m\) is sharp for this model formal group.

The accompanying symbolic script verifies (2), (3), and (10) for a range of \(m\) using exact integer arithmetic.

## 5. Why this still does not solve the \(\pi\)-height problem

The construction supplies high-order smallness in the real formal parameter, but the rational height of \(B_m(V)\) is not controlled by (7)–(8).

There are two sources of growth:
1. the multiples \([j]V\), whose canonical heights grow quadratically in \(j\);
2. the interpolation weights, whose total absolute size is
   \[
   \sum_{j=2}^{m+2}|\lambda_{m,j}|=2^{m+1}-1.
   \]

A crude height inequality for sums of rationals gives a bound involving the sum of the heights of all \(\xi([j]V)\), which can grow at least like a polynomial in \(m\) times \(\widehat h(V)\). The recurrence step used to make \(V\) real-close to the identity gives no sufficiently sharp denominator bound in terms of \(|z|\). Thus the gain \(|z|^m\) cannot yet be compared with the denominator cost.

In particular, (9) approximates a **rational** target \(r\), not \(\pi\) itself. To transfer it to \(\pi\), one would have to combine it with rational approximants \(r=p/q\), while controlling the new denominator and the \(q^2\)-weighted error. No such quantitative estimate is proved here.

## 6. Research conclusion

The previous first-order formal-group construction has now been upgraded to arbitrary-order real cancellation with integral coefficients, preserving the same \(p\)-adic congruences. The next genuine bottleneck is no longer formal cancellation; it is a quantitative height theorem for the resulting rational combination.

This result is a new algebraic tool for the project, not a proof of bad approximability of \(\pi\), and not a claim that the elliptic-curve paper itself states the higher-order lemma.


## 7. A crude height budget

Let \(h\) be the logarithmic Weil height on \(\mathbb P^1(\mathbb Q)\), and \(\widehat h_E\) the canonical height on \(E\). Since \(\xi\) is a fixed rational function on \(E\), standard height comparison gives
\[
h(\xi(Q))\le C_E\widehat h_E(Q)+O_E(1).
\]
For \(Q=[j]V\),
\[
\widehat h_E([j]V)=j^2\widehat h_E(V).
\]
Using \(h(x+y)\le h(x)+h(y)+O(1)\) repeatedly in (6), one obtains the crude estimate
\[
\boxed{
h(B_m(V))
\le C_E\left(\sum_{j=2}^{m+2}j^2\right)\widehat h_E(V)
+O_E(m^2)
=O_E(m^3\widehat h_E(V)+m^2).
}
\tag{11}
\]
The \(O(m^2)\) term accounts for the logarithmic heights of the binomial weights and repeated rational additions. This is only an upper bound; it does not show that the reduced denominator actually has this size.

If \(V=[s]P\), then \(\widehat h_E(V)=s^2\widehat h_E(P)\). The compact-recurrence argument supplies some \(s\) for which \(|\xi(V)|\) is small, but the present construction does not give a sufficiently sharp relation between this \(s\), the real size of \(\xi(V)\), and the reduced denominator of \(B_m(V)\). Thus (11) does not convert the real error \(O(|\xi(V)|^m)\) into a useful power saving in the rational height.

For the weights themselves,
\[
\sum_{j=2}^{m+2}|\lambda_{m,j}|=2^{m+1}-1.
\]
The exact integer checks through \(m=30\) are stored in `elliptic_high_order_weight_checks.json` and can be regenerated with `verify_elliptic_high_order_weights.py`.


## 8. Minimal support and the binomial uniqueness barrier

The moment identity itself forces at least \(m+1\) distinct multiplication indices. If one insists that every index be at most \(m+2\), the support is forced to be exactly \(\{2,\ldots,m+2\}\), and the unique weights are the binomial finite-difference weights. Thus any different integer cancellation scheme must introduce some multiplier \(j\ge m+3\). This is proved in `ELLIPTIC_HEIGHT_AWARE_WEIGHTS_ACTUAL_TEST.md`.

The exact elliptic experiment there shows why this matters: the tested alternatives reduce \(\sum|\lambda_j|j^2\) but introduce larger multipliers, and on \(E:y^2=x^3-2\) their actual rational heights and errors are both worse in all nine tested cases.


## 9. Curve-specific improvement using CM symmetry

For \(E:y^2=x^3-2\), the sixth-order automorphism forces \([j]_F(z)\) to contain only powers \(z^{6k+1}\). The generic polynomial-moment filter therefore cancels many coefficients that are identically zero. A sparse three-point filter
\[
(\lambda_2,\lambda_3,\lambda_4)=(-364,291,-36)
\]
satisfies only the two necessary nonzero-coefficient constraints and gives
\[
\sum_{j=2}^4\lambda_j\frac{[j]_F(z)}z
=1-1\,031\,184\,000\,z^{12}+O(z^{18}).
\]
See CM_SPARSE_ORDER12_CANCELLATION.md for the derivation, minimum-multiplier argument, EDS height consequence, and exact elliptic experiment. This illustrates that the generic binomial scheme is not optimal once the automorphism-forced sparsity of a particular formal group is used.
