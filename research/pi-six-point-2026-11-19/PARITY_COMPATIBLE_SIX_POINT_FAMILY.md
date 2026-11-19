# Parity-compatible shifted six-point family

Date: 2026-11-19  
Status: **proved infinite auxiliary construction**. This repairs a gap in the earlier shifted-six-point note: the original construction controlled tail magnitudes but did not ensure that all six pairs had opposite parity. The present version does.

## 1. A parity-controlled rational approximant to the tail slope

Let
\[
\lambda=\frac{\log(2/5)}{\log(3/7)}.
\]
Let \(p_k/q_k\) be its continued-fraction convergents. Consecutive vectors
\((p_{k-1},q_{k-1})\), \((p_k,q_k)\) form a basis of \(\mathbb F_2^2\), because their determinant is \(\pm1\). Therefore the three nonzero residue vectors
\[
(p_{k-1},q_{k-1}),\quad(p_k,q_k),\quad
(p_{k-1}+p_k,q_{k-1}+q_k)
\]
represent all three nonzero elements of \(\mathbb F_2^2\).

For every \(k\), choose one of these three vectors whose residue is
\[
(p,q)\equiv(0,1)\pmod2.
\]
The resulting integers are coprime, since they are either a convergent pair or the sum of two consecutive convergent pairs. Also \(q\to\infty\), and
\[
\left|\frac pq-\lambda\right|\le\frac{2}{q^2}
\tag{1}
\]
for all sufficiently large \(k\). Indeed each of the two convergent errors
\(|p_i-\lambda q_i|\) is \(<1/q_i\), and the chosen denominator is at least the larger convergent denominator; summing gives (1), with a harmless fixed constant.

Thus we have infinitely many coprime approximants with \(p\) even and \(q\) odd.

## 2. Six points with opposite parity

Let \(c_0\) be an integer nearest to \(qd_*\), where
\[
d_*=\frac{\log C}{2\log(3/7)},\qquad
C=\lambda\frac{(2/5)}{(3/7)}
\frac{1+(3/7)^2}{1+(2/5)^2}.
\]
Take three offsets
\[
c_-=c_0-1,\quad c_0,\quad c_+=c_0+1.
\]
For each \(c\), solve
\[
p u_c\equiv-c\pmod q,\qquad 0\le u_c<q.
\]
Since \((p,q)=1\), there is a unique residue \(u_c\). Choose \(x_c\) to be one of \(u_c+q,u_c+2q\) so that
\[
x_c\not\equiv c\pmod2.
\tag{2}
\]
This is possible because \(q\) is odd, so the two representatives have opposite parity.

Put
\[
n_c=\frac{px_c+c}{q}.
\]
Because \(p\) is even and \(q\) is odd,
\[
n_c\equiv c\pmod2.
\]
Thus \(m=x_c\) and \(n=n_c\) have opposite parity.

Now take the paired point
\[
(x_c+2q,\ n_c+2p).
\]
Both coordinates retain their parity, so all six points have opposite-parity coordinates. They also retain the same integer offset:
\[
q(n_c+2p)-p(x_c+2q)=c.
\]

The six points are therefore
\[
(x_c,n_c),\quad(x_c+2q,n_c+2p),
\qquad c\in\{c_-,c_0,c_+\}.
\tag{3}
\]
All coordinates are positive and \(O(q)\).

For every point in (3),
\[
n-\lambda m
=\left(\frac pq-\lambda\right)m+\frac cq
=d_*+O(q^{-1}),
\tag{4}
\]
uniformly, by (1), \(m=O(q)\), and
\(|c/q-d_*|\le3/(2q)\).

The tail asymptotic from the preceding note then yields
\[
\boxed{
\max_{\text{six points}}
\left|
\frac{r_m(2/5)}{r_n(3/7)}-1
\right|=O(q^{-1}).
}
\tag{5}
\]
Since every pair has opposite parity, the signed error
\[
f_{m,n}=\frac\pi4-T_m(2/5)-T_n(3/7)
\]
is exactly, up to sign,
\[
f_{m,n}=\pm\bigl(r_m(2/5)-r_n(3/7)\bigr).
\tag{6}
\]

This is the parity-compatible infinite family needed to insert the signed tail error into the last determinant column.

## 3. Exact nonzero quadratic determinant

Let
\[
\kappa=\frac{2x_{c_0}-x_{c_-}-x_{c_+}}q.
\]
It is an integer because \(x_c\equiv-p^{-1}c\pmod q\). If \(\kappa=0\), replace \(x_{c_+}\) by \(x_{c_+}+2q\); this preserves (2) and changes \(\kappa\) to \(-2\). Thus \(\kappa\ne0\), and \(|\kappa|\) is bounded by an absolute constant.

For the six rows
\[
(1,m,n,m^2,mn,n^2),
\]
use \(c=qn-pm\) as the second coordinate. The transformed feature determinant is \(q^4\), while the paired \(m\)-coordinates differ by \(2q\). Direct expansion gives
\[
\boxed{\det A=\pm32\kappa.}
\tag{7}
\]
Hence this six-point family has an exact, bounded, nonzero integer determinant for every sufficiently large chosen denominator \(q\).

## 4. Exact tail-column insertion

For each of the six points define
\[
f_i=\frac\pi4-T_{m_i}(2/5)-T_{n_i}(3/7).
\]
Then every \(f_i\) is a signed difference of two positive tails by (6), and (5) gives the uniform \(O(q^{-1})\) relative cancellation.

Replacing the \(n^2\) column of the quadratic evaluation matrix by the column \(f_i\) produces a rational determinant after the common \(\pi/4\) term is removed using the all-ones column. Clearing the exact common denominator of the six Taylor truncations gives an integer determinant.

The remaining two tasks are still independent:

1. Prove this new tail-column determinant is nonzero (for example, by a modular certificate).
2. Prove its absolute value is smaller than the reciprocal of its **actual reduced denominator**.

The old quadratic determinant (7) does not by itself prove the tail-column determinant is nonzero; the column replacement changes the determinant. Likewise, a common-denominator upper bound is not a lower bound for the reduced determinant denominator.

## 5. Status

The parity gap in the earlier construction is repaired. We now have an infinite family with:

- all six pairs of opposite parity;
- six tail ratios tending to one at rate \(O(q^{-1})\);
- signed \(\pi/4\)-linear forms equal to differences of positive tails;
- a bounded nonzero quadratic integer determinant.

No claim is made that the tail-column determinant is nonzero for all members, that the denominator budget closes, or that \(\pi\) is badly approximable.

## 5. Finite parity-compatible certificate

A concrete parity-compatible choice is
\[
p/q=332/307,\qquad c_0=-5,
\]
with offsets \(-6,-5,-4\). The six points are
\[
(823,890),\ (1437,1554),\
(430,465),\ (1044,1129),\
(651,704),\ (1265,1368).
\]
All six pairs have opposite parity, and their offsets \(qn-pm\) are
\(-6,-6,-5,-5,-4,-4\).

Using 12 terms of the alternating rational tail series gives strict enclosures (entries are integers divided by \(10^{12}\)):
\[
\begin{array}{c|c}
(m,n)&r_m(2/5)/r_n(3/7)\\ \hline
(823,890)&(1005082313506,1005082315561)/10^{12}\\
(1437,1554)&(1011709117395,1011709119475)/10^{12}\\
(430,465)&(1006362522816,1006362524848)/10^{12}\\
(1044,1129)&(1013046226429,1013046228506)/10^{12}\\
(651,704)&(1014368986904,1014368988971)/10^{12}\\
(1265,1368)&(1021066594109,1021066596206)/10^{12}.
\end{array}
\]
For these points, the rational determinant with final column
\(T_m(2/5)+T_n(3/7)\), after clearing its exact common denominator, is nonzero and congruent to \(4\pmod5\). The accompanying script performs all operations with exact fractions and integer arithmetic:
[verify_parity_tail_column.py](verify_parity_tail_column.py).

This is a finite nonzero certificate for the tail-column determinant, not an infinite-family proof.

## 6. Status

The parity gap is repaired. We have an infinite family with all six pairs of opposite parity, tail ratios converging to one at rate \(O(q^{-1})\), signed \(\pi/4\)-linear forms equal to differences of positive tails, and a bounded nonzero quadratic integer determinant. For one explicit member, the tail-column determinant also has an exact modular nonzero certificate.

No claim is made that the tail-column determinant is nonzero for all members, that the denominator budget closes, or that \(\pi\) is badly approximable.
