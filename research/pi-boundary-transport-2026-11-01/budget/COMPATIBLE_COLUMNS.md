# Why compatibility of transverse components must be preserved

2026-11-01. This note accompanies the exact Gram displacement counterexample and boundary identity. It makes no assertion about pi's approximation constant.

## 1. A determinant with the same two transverse types

Take two centers 0 and r, with n consecutive time jets at each and transverse derivative index beta=0. Use the 2n columns

\[
Y^h,\quad XY^h,\qquad 0\le h<n,
\]

and the source substitution Y=1+t, X=jr+log(1+t), j=0,1. With rows grouped by center and columns by the displayed families, the matrix is

\[
A_n(r)=\begin{pmatrix}B&L\\B&L+rB\end{pmatrix},
\quad B_{s,h}=\binom hs,\quad
L_{s,h}=[t^s](1+t)^h\log(1+t).
\]

Subtracting the first block of rows from the second gives

\[
\boxed{\det A_n(r)=r^n,}                         \tag{1}
\]

because B is triangular with diagonal one. This is an exact polynomial identity, independent of any period or irrationality theorem. Thus for omega nonzero and |delta/omega|<=1/2,

\[
\frac{|\det A_n(\omega+\delta)|}{|\det A_n(\omega)|}
=|1+\delta/\omega|^n,
\quad
\left|\log\frac{|\det A_n(\omega+\delta)|}
{|\det A_n(\omega)|}\right|
\le 2n|\delta/\omega|.                           \tag{2}
\]

In contrast, in the grouped Hardy coefficient space with the same transverse types 0 and 1, the preceding obstruction note proves

\[
\frac{\det G_n(\delta)}{\det G_n(0)}
\ge |\delta|^{2n}
\left(\frac{R^2-|\omega|^2}{|\omega|^2}\right)^{n^2}.
\]

For R>sqrt(2)|omega| this can grow exponentially in n^2. Both statements can hold: the kernel allows arbitrary independent functions for the two transverse types, whereas the actual columns have

\[
(f_0,f_1)_{Y^h}=(e^{hz},0),\qquad
(f_0,f_1)_{XY^h}=(ze^{hz},e^{hz}).                \tag{3}
\]

The two components of an XY^h column are linked by multiplication by z. The amplified independent kernel directions need not be independently available among these columns. This comparison is a model calculation, not a counterexample to the source rank theorem. As n grows this row packet alone is not the full fixed-weight source simplex.

## 2. The exact missing factor in a Gram bound

For A=ZC, assume the two Gram matrices G=ZZ* and T=C*C are positive definite, and put

\[
U=G^{-1/2}Z,\qquad V=CT^{-1/2}.
\]

Then UU*=I and V*V=I, and

\[
\boxed{|\det A|^2=\det G\det T\,|\det(UV)|^2.}   \tag{4}
\]

Every singular value of UV is at most one. The usual Gram bound discards the last factor. Formula (4) is an identity identifying the discarded information; it is not itself an estimate of that factor. The row and column spaces are finite dimensional subspaces of the same Hardy coefficient Hilbert space, so the argument remains valid for the countably infinite intermediate basis.

## 3. General compatibility equations

For a polynomial P in the weighted column space, put

\[
F_P(z,u)=P(e^z,z+u_1,\ldots,z+u_m),\qquad
\mathcal D=\partial_z-\sum_i\partial_{u_i},
\quad d_0=\lfloor H/w_0\rfloor.
\]

On a monomial component Y^h X^alpha,
\(\mathcal D F_P=hF_P\). Therefore every allowed column satisfies

\[
\boxed{\prod_{h=0}^{d_0}(\mathcal D-h)F_P=0.}     \tag{5}
\]

The transverse polynomial support and the full weighted budget give additional restrictions. For a fixed h component, if f_a=[u^a]F_P, then

\[
f_a'(z)-hf_a(z)=\sum_i(a_i+1)f_{a+e_i}(z).       \tag{6}
\]

Equations (5)--(6) follow directly by differentiation. For mixed h columns, (6) must be applied after decomposing into h components; a single h cannot be assigned to such a column. Equation (5) is necessary, but it alone does not encode the weighted simplex budget.

The boundary transport identity keeps these compatible columns intact while canceling all internal row transport. Its remaining mixed minors, rather than an unconstrained comparison of G(delta) with G(0), are a concrete next target. Neither (4) nor (5) presently supplies the high-dimensional determinant estimate required to prove bad approximability.
