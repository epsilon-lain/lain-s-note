# Exact extraction of the integer modular subquotient

2026-11-02. This is an explicit bridge calculation from the matrices in the introductory setup of <https://alpo.ge/s6.pdf>. It uses neither the asserted topology of the completed threefold nor any numerical value of pi.

## 1. Input and primitive invariant flag

In the ordered basis `(gamma,u,w,delta)`, take

\[
T_1=\begin{pmatrix}
1&0&-6&2\\0&-1&1&1\\0&-1&0&1\\0&0&0&1
\end{pmatrix},\qquad
T_2=\begin{pmatrix}
1&6&0&-3\\0&0&-1&1\\0&1&0&0\\0&0&0&1
\end{pmatrix}.
\]

They have determinant one and satisfy `T1^3=I`, `T2^4=I`. The integral direct-summand flag

\[
\mathbb Z\gamma\subset
W=\mathbb Z\gamma+\mathbb Zu+\mathbb Zw\subset V=\mathbb Z^4
\]

is invariant. Thus the middle subquotient `W/Z gamma` is a free rank-two integer local system, with no saturation index or clearing denominator. In the basis `(u,w)` its matrices are

\[
A_1=\begin{pmatrix}-1&1\\-1&0\end{pmatrix},\quad
A_2=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\quad
A_0=(A_1A_2)^{-1}=\begin{pmatrix}1&-1\\0&1\end{pmatrix}.
\]

Their fractional-linear actions on tau are exactly

\[
\tau\mapsto(\tau-1)/\tau,\qquad
\tau\mapsto-1/\tau,\qquad
\tau\mapsto\tau-1.
\]

These are the closed tau transformations displayed by the source's period system. The projective image is the usual modular group: `A2^2=-I` and `A1 A2` is a unit upper shear. This identifies a classical rank-two modular subsystem inside the rank-four construction. It does not assert that a generic complex two-torus of the family has an elliptic quotient.

## 2. The local unipotent mechanism is integral and linear in winding number

Direct computation gives

\[
T_0=(T_1T_2)^{-1}=
\begin{pmatrix}
1&0&0&1\\0&1&-1&0\\0&0&1&0\\0&0&0&1
\end{pmatrix},\qquad
N=T_0-I.
\]

`N^2=0`, and its Smith invariants are `(1,1,0,0)`. Indeed `N delta=gamma`, `N w=-u`, and `N gamma=N u=0`, so its image and kernel are the same primitive rank-two sublattice `Z gamma+Z u`. For every integer k, including negative k,

\[
\boxed{T_0^k=I+kN.}
\]

Consequently local branch differences remove a single logarithmic layer without exponential coefficient growth. In a suitable regular-singular flat-solution description, a logarithmic component has the form

\[
g(z)=h(z)+\frac{\log z}{2\pi i}f(z),
\]

and one positive circuit gives `g -> g+f`. The first branch difference is `f`; the second branch difference is zero. This is a precise cancellation operation supplied by the integer local system. To use it for rational approximation, the actual flat connection, normalization, and arithmetic of f and h must still be supplied.

For example `g(z)=log(z)/(2 pi i)+C` has the same integral monodromy for every complex C. More generally, multiplying both a holomorphic and logarithmic solution by an arbitrary nonzero holomorphic function preserves this local monodromy while changing its Taylor coefficients. Hence integer monodromy alone neither fixes connection constants nor bounds Taylor denominators.

## 3. What has actually been eliminated

The source writes the two-torus period matrix as

\[
\Pi=\begin{pmatrix}6\mu&\tau&1&0\\\beta&\mu&0&1\end{pmatrix}.
\]

The matrix subquotient and its tau transformation law are independent of mu and beta. There are two important qualifications.

First, this representation is not a rational direct sum containing the modular rank-two module. The exact intertwining equations

\[
T_iX=XA_i\quad(i=1,2),\qquad X\in M_{4\times2}(\mathbb Q)
\]

have only `X=0`. Equivalently, the invariant flag is a genuine extension rather than a free choice of two invariant columns.

An elementary check: the delta row of an intertwiner must vanish, since the modular two-dimensional representation has no fixed covector. Its u,w block must be `lambda I`, since A1 and A2 have only scalar common centralizers. Writing the gamma row as `(a,b)`, equivariance forces

\[
2a+b=0,\quad a+b=0,\quad b-a=6\lambda,
\]

hence `a=b=lambda=0`.

Second, if a constant integer vector v is required to make the *whole column* `Pi v` independent of both mu and beta, coefficient comparison forces `v_gamma=v_u=0`. Only the two constant columns remain. One can read tau from an individual matrix entry, but constructing its correctly normalized flat period basis requires the analytic differential equation; it is not obtained merely by integer column cancellation.

## 4. Constructive handoff

The useful handoff to the pi construction is therefore explicit:

1. Discard the extension variables at the level of the primitive modular subquotient, without paying an integer saturation denominator.
2. Use the source's modular uniformization or its differential equation to calculate the normalized rank-two periods and their Taylor recurrences.
3. Use an independently explicit period identity, such as an appropriate Wronskian or a relation between complete and quasi-periods, to eliminate the remaining period values and produce a rational linear form in pi or 1/pi.
4. Calculate the fully reduced rational denominator and the analytic remainder together. Neither a topological gluing identity nor the integer monodromy is a substitute for this joint estimate.

This calculation makes the proposed transfer concrete. It recovers a classical modular source of periods, so any resulting standard modular formula must be labelled a known baseline rather than a new consequence of the six-sphere topology.

## Reproduction

`python verify_monodromy.py` checks the matrices, the primitive flag, all cusp powers from -20 through 20, Smith invariants, the zero intertwiner space, and symbolic elimination in the full period matrix. The statements for all integer powers and the nonsplitting proof are given above; the finite checks are supplementary.
