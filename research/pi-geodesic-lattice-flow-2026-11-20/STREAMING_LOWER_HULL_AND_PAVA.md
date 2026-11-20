# Streaming lower-hull algorithm for the continuous lattice flow

Date: 2026-11-20  
Status: **proved exact finite-prefix algorithm**. It turns the geometric filtering rule into a one-pass stack process over continued-fraction convergents.

## 1. The key observation

After the active-vector theorem, only convergents need to be considered. Write their points as
\[
P_n=(X_n,Y_n)=(q_n^2,E_n^2),
\qquad E_n=|q_n\alpha-p_n|.
\]
The lattice-flow objective is the family of linear functionals
\[
\ell_n(x)=X_n+xY_n,\qquad x=e^{4t}\ge1.
\]
Since \(X_n\) strictly increases and \(Y_n\) strictly decreases, define the positive drop slope between two points \(P_i,P_j\), \(i<j\), by
\[
D(i,j)=\frac{Y_i-Y_j}{X_j-X_i}>0.
\]
For three consecutive points \(A,B,C\), the middle point \(B\) is strictly below the chord \(AC\) exactly when
\[
D(A,B)>D(B,C).
\]
If instead
\[
D(A,B)\le D(B,C),
\tag{1}
\]
then \(B\) lies on or above the chord \(AC\) and can never be a strict lower-hull vertex.

## 2. Exact streaming stack algorithm

Process the points in increasing denominator order. Start with the vertical lattice vector
\[
P_{-1}=(0,1),
\]
then append convergent points \(P_0,P_1,\ldots\).

Maintain a stack \(H\). When a new point \(C\) arrives:

1. While the stack has at least two points \(A,B\) at its end and \(D(A,B)\le D(B,C)\), remove \(B\).
2. Append \(C\).

The invariant is that the stack is exactly the lower convex hull of all points processed so far (with collinear interior points removed). This follows by induction: a new rightmost point can only invalidate a suffix of the old lower hull; condition (1) identifies the last vertex that has become hidden above a chord. Once a point is removed, it is above a chord joining two already processed points, so no later point can make it a vertex of the lower hull of the enlarged set.

Each point is pushed once and popped at most once. Thus for \(K\) convergents:
\[
\boxed{\text{time }O(K),\qquad \text{space }O(K).}
\]
This is a streaming algorithm, but not a fixed finite-memory rule: the stack can grow with \(K\).

## 3. Exact connection with weighted pool-adjacent-violators

For adjacent convergents define
\[
w_n=X_{n+1}-X_n=q_{n+1}^2-q_n^2>0,
\qquad
d_n=D(n,n+1).
\]
For a block \(i\le n<j\), telescoping gives
\[
D(i,j)=
\frac{\sum_{n=i}^{j-1}w_nd_n}
{\sum_{n=i}^{j-1}w_n}.
\tag{2}
\]
Therefore when (1) fails, the two adjacent blocks are pooled and their new drop slope is their weighted average:
\[
d_{\mathrm{merged}}=
\frac{W_1d_1+W_2d_2}{W_1+W_2},
\qquad W_r=\sum_{n\in\text{block }r}w_n.
\]
This is exactly the pool-adjacent-violators algorithm for projecting the sequence \(d_n\) onto nonincreasing sequences, with weights \(w_n\). The surviving hull vertices are the boundaries between the final pooled blocks.

This gives a second interpretation of the curve: it continuously reads the lower convex envelope, while the discrete convergent stream is compressed by weighted isotonic regression.

## 4. Recovering the continuous switching times

Let the final finite-prefix hull be
\[
H_0,H_1,\ldots,H_m
\]
with drop slopes
\[
\Delta_i=D(H_i,H_{i+1}).
\]
The hull property gives
\[
\Delta_0>\Delta_1>\cdots>\Delta_{m-1}>0.
\]
The switch from \(H_i\) to \(H_{i+1}\) occurs at
\[
\boxed{x_i=\frac1{\Delta_i},\qquad
t_i=\frac14\log\frac1{\Delta_i}.}
\tag{3}
\]
For \(x\ge1\), discard transitions with \(x_i<1\). A hull vertex is active precisely on the interval between its entering and leaving switch times, intersected with \(x\ge1\).

For an infinite stream, the last vertex of a finite prefix is only provisionally active after the final computed switch; future convergents may truncate that interval or remove that last vertex. Earlier vertices already certified by the processed prefix can still be affected only if they are in the current suffix, exactly as the stack update records.

## 5. Finite experiment for pi

Using 100 decimal digits and all convergent denominators \(q\le10^6\), the streaming stack gives the final finite-prefix hull:
\[
\begin{aligned}
&(0,1),\ (1,3),\ (7,22),\ (113,355),\\
&(33215,104348),\ (66317,208341),\\
&(99532,312689),\ (364913,1146408).
\end{aligned}
\]
It permanently removes the convergents
\[
(106,333),\quad(33102,103993),\quad(265381,833719)
\]
from this finite-prefix hull.

The switch times for the nonvertical vertices are approximately
\[
0.0050630444,\ 1.9461796353,\ 4.7263240240,\ 
10.4459621432,\ 11.3826074100,\ 11.5025193028,\ 
12.7701412190.
\]
The last vector is active only relative to the cutoff \(q\le10^6\); larger convergents may later change the final interval.

The experiment script generates convergents, applies the stack, and emits a JSON report. It is a finite numerical regression, not an infinite description of pi's active sequence.

## 6. Why this matters for the original continuous-velocity idea

The original triangular clock visited every integer pair. The lattice flow instead moves the whole lattice continuously, and the stack algorithm supplies a low-redundancy arithmetic description of the points that can control the motion.

The filtering is not achieved by manually deleting fractions from an enumeration. It is the lower-hull operation forced by the geometry. The stack gives a one-pass way to discover that hull as convergents arrive.

The limitation is also explicit: the stack is not bounded in size by this argument, and the algorithm still needs the continued-fraction data of alpha. For alpha=pi, the arithmetic challenge remains to prove a uniform lower bound for q||q pi||. The streaming hull is a structural reformulation and an efficient experiment, not a proof of bad approximability.

## 7. Reproducibility

Run pi_lattice_stack_experiment.py (requires mpmath). It writes pi_lattice_stack_experiment.json. The script uses the continued-fraction theorem to restrict candidates to convergents, then applies the exact stack inequalities to high-precision approximations. Exact rational interval certification of individual switch times is a separate optional step.
