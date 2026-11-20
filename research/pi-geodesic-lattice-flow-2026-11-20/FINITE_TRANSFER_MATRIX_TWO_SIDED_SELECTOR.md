# A finite transfer-matrix formula for the two-sided lattice selector

Date: 2026-11-20  
Status: **proved exact reformulation**. It removes the need to construct intermediate numerator-denominator pairs explicitly. The past remains a finite prefix statistic; the future requires only the explicit lookahead window proved in the preceding note.

## 1. Normalize away the scale

Let \(\alpha\) be irrational, with convergents \(p_n/q_n\), errors
\[
E_n=|q_n\alpha-p_n|,\qquad E_{n+1}<E_n.
\]
For a fixed candidate \(n\), define
\[
P_n^*=\max\left(1,\max_{i<n}
\frac{q_n^2-q_i^2}{E_i^2-E_n^2}\right),
\]
and
\[
F_n=\min_{j>n}
\frac{q_j^2-q_n^2}{E_n^2-E_j^2}.
\]
The candidate is active exactly when \(P_n^*<F_n\).

Set the dimensionless scale
\[
\rho_n=\frac{q_n}{E_n},
\]
and normalized ratios
\[
R_i=\frac{q_i}{q_n},\quad Z_i=\frac{E_i}{E_n}\quad(i<n),
\]
\[
\widehat R_k=\frac{q_{n+k}}{q_n},\quad
\widehat Z_k=\frac{E_{n+k}}{E_n}\quad(k\ge1).
\]
Then
\[
\boxed{
\frac{P_n^*}{\rho_n^2}
=
\max\left(\rho_n^{-2},
\max_{i<n}\frac{1-R_i^2}{Z_i^2-1}\right)
}
\tag{1}
\]
and
\[
\boxed{
\frac{F_n}{\rho_n^2}
=
\min_{1\le k\le K_n}
\frac{\widehat R_k^2-1}{1-\widehat Z_k^2},
}
\tag{2}
\]
where \(K_n\) is any valid future-lookahead bound from `EXPLICIT_LOOKAHEAD_FROM_PARTIAL_QUOTIENTS.md`. The common scale \(\rho_n^2\) cancels. Thus activity is decided by comparing the two dimensionless quantities in (1) and (2).

## 2. Future ratios from a two-by-two transfer matrix

Define
\[
M(a)=\begin{pmatrix}a&1\\1&0\end{pmatrix},
\qquad
N(a)=\begin{pmatrix}-a&1\\1&0\end{pmatrix}.
\]
The denominator recurrence gives
\[
\binom{q_{m+1}}{q_m}
=M(a_{m+1})\binom{q_m}{q_{m-1}}.
\]
Let
\[
r_n=\frac{q_{n-1}}{q_n},
\qquad z_n=\frac{E_{n+1}}{E_n}
=\frac1{\alpha_{n+2}}.
\]
Then for \(k\ge1\),
\[
\boxed{
\widehat R_k
=
(1,0)\,
M(a_{n+k})\cdots M(a_{n+1})
\binom{1}{r_n}.
}
\tag{3}
\]
For \(k=1\), this gives \(a_{n+1}+r_n\), as required.

The error recurrence
\[
E_{m+1}=E_{m-1}-a_{m+1}E_m
\]
gives
\[
\binom{E_{m+1}}{E_m}
=N(a_{m+1})\binom{E_m}{E_{m-1}}.
\]
Consequently
\[
\boxed{
\widehat Z_1=z_n,\qquad
\widehat Z_k
=
(1,0)\,
N(a_{n+k})\cdots N(a_{n+2})
\binom{z_n}{1}
\quad(k\ge2).
}
\tag{4}
\]
All matrices are explicit and have determinant \(-1\). Equations (3)–(4) compute the needed future ratios from the finite partial-quotient block and the single complete-quotient parameter \(z_n\).

## 3. Past ratios from reversed transfer matrices

For \(i<n\), the denominator recurrence run backward is
\[
q_{m-1}=q_{m+1}-a_{m+1}q_m.
\]
Hence
\[
\binom{q_{m-1}}{q_m}
=N(a_{m+1})\binom{q_m}{q_{m+1}}.
\]
Starting from
\[
\binom{q_{n-1}/q_n}{q_n/q_n}
=\binom{r_n}{1},
\]
we obtain, for \(h\ge1\),
\[
\boxed{R_{n-1}=r_n,}
\]
and for \(h\ge2\),
\[
\boxed{
R_{n-h}
=
(1,0)\,
N(a_{n-h+2})N(a_{n-h+3})\cdots N(a_n)
\binom{r_n}{1}.
}
\tag{5}
\]

For errors, use
\[
E_{m-1}=a_{m+1}E_m+E_{m+1}.
\]
Since
\[
\frac{E_{n-1}}{E_n}=a_{n+1}+z_n,
\]
we get
\[
\boxed{
Z_{n-1}=a_{n+1}+z_n.
}
\]
For \(h\ge2\),
\[
\boxed{
Z_{n-h}
=
(1,0)\,
M(a_{n-h+2})M(a_{n-h+3})\cdots M(a_n)
\binom{a_{n+1}+z_n}{1}.
}
\tag{6}
\]
Thus the entire past statistic in (1) is a finite maximum of explicit matrix expressions in the reversed finite continued-fraction prefix and the same scalar \(z_n\).

## 4. The two-sided selector is now a finite matrix expression

Combining (1)–(6), define
\[
\mathcal P_n(z_n)=
\max\left(
\rho_n^{-2},
\max_{1\le h\le n}
\frac{1-R_{n-h}^2}{Z_{n-h}^2-1}
\right),
\]
and
\[
\mathcal F_n(z_n)=
\min_{1\le k\le K_n}
\frac{\widehat R_k^2-1}{1-\widehat Z_k^2}.
\]
Then
\[
\boxed{
q_n\text{ is active}
\iff
\mathcal P_n(z_n)<\mathcal F_n(z_n).
}
\tag{7}
\]
For fixed \(n\), this is a finite collection of two-by-two matrix products and finite max/min operations. It does not require constructing all intermediate pairs \((p_j,q_j)\), nor scanning arbitrary rational denominators.

The only nonlocal input is \(z_n=1/\alpha_{n+2}\), the forward complete quotient. It can be enclosed by a continued-fraction cylinder interval from a finite block of future partial quotients. Substituting that interval into (7), with outward-rounded interval arithmetic, can certify activity or inactivity without knowing the infinite tail exactly.

## 5. What this does and does not achieve

**Achieved:**
- The past and future thresholds share the same scale \(\rho_n^2\), which cancels.
- All denominator ratios are transfer-matrix products of partial quotients.
- All error ratios are transfer-matrix products plus one complete-quotient parameter \(z_n\).
- The future matrix products have a proved finite length \(K_n\), controlled by \(a_{n+2},a_{n+3}\).
- The resulting activity test is a finite matrix/max/min formula for each candidate.

**Not achieved:**
- The past maximum still ranges over the entire finite prefix \(1,\ldots,n-1\); no uniform compression of that history is proved.
- The complete quotient \(z_n\) is not determined by only two partial quotients; a finite cylinder gives an interval, not its exact value.
- This does not prove that \(\pi\) is badly approximable or give a positive lower bound for \(q\|q\pi\|\).

The useful next step is to use interval arithmetic on the complete-quotient cylinder and determine when the sign of
\[
\mathcal F_n(z)-\mathcal P_n(z)
\]
is constant throughout that interval. If it is, a finite block of partial quotients certifies the candidate's status. If not, refine the cylinder by reading more future digits.

## 6. Connection with the user's original one-dimensional idea

The original clock used one scalar to enumerate all lattice points. The current construction uses one scalar flow parameter and a finite transfer-matrix selector. It does not visit all rational pairs: it selects only the candidates that can support the lower envelope. The past and future are represented by opposite matrix products, while the complete quotient couples them.

This is a precise two-sided arithmetic boundary-value formulation. It is not physical retrocausality; it is a global optimization criterion expressed through finite continued-fraction data.


## 7. Finite continued-fraction cylinders for the tail parameter

The only infinite-tail input in the selector is
\[
z_n=\frac1{\alpha_{n+2}}.
\]
A finite block of future digits gives a rigorous interval for it. Suppose the known prefix is
\[
\alpha_{n+2}=[b_0;b_1,\ldots,b_m,\tau],
\qquad \tau>1,
\]
where \(b_j=a_{n+2+j}\). Let \(P_m/Q_m\) and \(P_{m-1}/Q_{m-1}\) be the convergents of the finite word \([b_0;\ldots,b_m]\). Then
\[
\alpha_{n+2}
=\frac{P_m\tau+P_{m-1}}{Q_m\tau+Q_{m-1}}.
\]
As \(\tau\) ranges over \((1,\infty)\), this Möbius function is monotone, so its range is the interval between
\[
\frac{P_m}{Q_m}
\quad\text{and}\quad
\frac{P_m+P_{m-1}}{Q_m+Q_{m-1}}.
\]
Taking reciprocals gives an explicit rational interval \(I_m\) containing \(z_n\).

Substitute \(z\in I_m\) into the finite matrix formulas (1)–(7), using outward-rounded interval arithmetic. If the resulting interval for
\[
\mathcal F_n(z)-\mathcal P_n(z)
\]
lies strictly above zero, the candidate is certified active; if it lies strictly below zero, it is certified inactive. If the interval contains zero, reveal another future partial quotient and refine the cylinder.

This gives an adaptive finite-digit certificate. It is not guaranteed to decide every candidate from a uniformly bounded number of digits: a candidate whose exact selector margin is zero (a hull tie), or arbitrarily close to zero, can require more precision. For strict activity/inactivity the nested cylinders shrink to the true irrational tail, so interval refinement eventually decides once the selector has a nonzero margin.

## 8. Precise interpretation of “direct”

The matrix selector avoids constructing intermediate numerator-denominator pairs and avoids scanning arbitrary rational numbers. It still performs:
- a finite scan of the already known past prefix (or an equivalent stored past summary);
- a bounded, partial-quotient-controlled future matrix calculation;
- interval refinement of one complete-quotient parameter when needed.

Thus it is a finite transfer-matrix formula with adaptive tail certification, not a formula depending on only two partial quotients and not a constant-time local rule.
