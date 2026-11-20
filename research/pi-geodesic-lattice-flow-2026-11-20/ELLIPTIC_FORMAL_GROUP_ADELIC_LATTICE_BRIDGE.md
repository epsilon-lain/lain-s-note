# Elliptic formal groups meet the continuous lattice flow: a two-sided local/global bridge

Date: 2026-11-20  
Status: **proved strengthening of the paper's local approximation lemma, plus an exact obstruction to the naive congruence-lattice approach**. No new lower bound for \(q\|q\pi\|\) is obtained.

## 1. The paper's mechanism, in the form relevant to our project

Pasten–Vidaux, arXiv:2611.24748, Lemma 4.1, fixes an elliptic curve \(E/\mathbb Q\) with a non-torsion point in a suitable rank-stable subgroup \(\Gamma\). For each prime \(p\mid M\), a finite-index open subgroup \(H_p\subset E(\mathbb Q_p)\) is chosen so that its formal parameter \(\xi=-x/y\) satisfies
\[
\xi(H_p)\subseteq M\mathbb Z_p.
\]
If \(L\) is divisible by all indices \([E(\mathbb Q_p):H_p]\), then \(V=[L]P\) lies in every \(H_p\). For a nonzero integer \(k\), put \(U=[k]V\). The formal group expansion
\[
\xi([k]V)=k\xi(V)+\xi(V)^2F_k(\xi(V)),
\qquad F_k\in\mathbb Z[[T]],
\]
gives
\[
\beta:=\frac{\xi(U)}{\xi(V)}\in\mathbb Q,
\qquad
\beta-k\in M\mathbb Z_p\quad(p\mid M).
\tag{1}
\]
This is a simultaneous local approximation at every prime dividing \(M\).

## 2. Strengthening: impose real closeness at the same time

The local construction can be synchronized with the real topology.

**Lemma (real–\(p\)-adic simultaneous approximation).**  
Assume \(E(\mathbb Q)\) contains a non-torsion point \(P\), and let \(\Gamma\) be a subgroup containing \(P\) and stable under integer multiplication. Fix \(M\ge2\), a nonzero integer \(k\), and \(\varepsilon>0\). Then there exist \(U,V\in\Gamma^\circ\) such that
\[
\beta=\frac{\xi(U)}{\xi(V)}\in\mathbb Q,
\]
and simultaneously
\[
\boxed{|\beta-k|<\varepsilon,\qquad
\beta-k\in M\mathbb Z_p\quad\text{for every }p\mid M.}
\tag{2}
\]

**Proof.** Choose the finite-index groups \(H_p\) as in Lemma 4.1 of Pasten–Vidaux, and let \(L\) be divisible by all their indices. The real Lie group \(E(\mathbb R)\) is compact. Since \(P\) is non-torsion, so is \(Q=[L]P\). In a compact group, the identity is an accumulation point of the positive multiples of any non-torsion element: there are integers \(m_j\to\infty\) with
\[
[m_j]Q\longrightarrow O
\quad\text{in }E(\mathbb R).
\]
Set \(V_j=[m_jL]P\). Each \(V_j\) belongs to every \(H_p\), is non-torsion, and tends to \(O\) in the real topology. For fixed \(k\ne0\), \(U_j=[k]V_j\) also tends to \(O\). For all sufficiently large \(j\), both points lie in \(\Gamma^\circ\), and \(\xi(V_j)\to0\) in \(\mathbb R\). Hence
\[
\frac{\xi(U_j)}{\xi(V_j)}-k
=\xi(V_j)F_k(\xi(V_j))\longrightarrow0
\]
in \(\mathbb R\). On the other hand \(\xi(V_j)\in M\mathbb Z_p\) and \(F_k\) has \(p\)-integral coefficients, so the same expression belongs to \(M\mathbb Z_p\) for every \(p\mid M\). Choose \(j\) large enough for the real inequality. \(\square\)

This is a strengthening of the paper's Lemma 4.1: it adds an arbitrarily small real error without sacrificing any of the prescribed local congruences. The new ingredient is recurrence of multiples on the compact real elliptic group.

## 3. The congruence is literally a sublattice condition

Write \(\beta=a/b\) in lowest terms with \(b>0\). Since \(b\) is an integer and (2) holds,
\[
v_p(a-kb)-v_p(b)\ge v_p(M)
\]
for every \(p\mid M\). Thus
\[
\boxed{a\equiv kb\pmod M.}
\tag{3}
\]
So the elliptic construction produces rational slopes \(a/b\) arbitrarily close in the real metric to \(k\), while their integer vectors \((b,a)\) lie in the index-\(M\) sublattice
\[
L_{k,M}=\{(q,p)\in\mathbb Z^2:p-kq\equiv0\pmod M\}.
\]
A basis is \((1,k),(0,M)\), so its covolume is \(M\).

More generally, start with any rational \(r=p_0/q_0\) in lowest terms. Apply the lemma with \(k=p_0\), then divide \(\beta\) by \(q_0\). The resulting rational \(\beta/q_0=a/(bq_0)\) can be made arbitrarily close to \(p_0/q_0\), and its vector \((Q,P)=(bq_0,a)\) lies in
\[
L_{p_0,q_0,M}
=\{(Q,P):Q\equiv0\pmod{q_0},\
P-p_0(Q/q_0)\equiv0\pmod M\}.
\]
This is a rank-two sublattice of covolume \(q_0M\), with basis \((q_0,p_0),(0,M)\).

Thus the elliptic formal group gives genuine points in congruence-restricted rational directions, with simultaneous real and \(p\)-adic control.

## 4. Exact obstruction: a fixed congruence sublattice is only an affine copy of the same approximation problem

For a fixed slope residue \(k\pmod M\), define the restricted approximation constant
\[
\kappa_{k,M}(\alpha)
=
\inf_{\substack{q\ge1,\ p\in\mathbb Z\\p-kq\in M\mathbb Z}}
q|p-\alpha q|.
\]
Writing \(p=kq+Mr\) gives
\[
\begin{aligned}
\kappa_{k,M}(\alpha)
&=\inf_{q\ge1,r\in\mathbb Z}
q|(k-\alpha)q+Mr|\\
&=M\inf_{q\ge1,r\in\mathbb Z}
q\left|r-q\frac{\alpha-k}{M}\right|.
\end{aligned}
\]
Hence
\[
\boxed{
\kappa_{k,M}(\alpha)
=M\,\kappa\!\left(\frac{\alpha-k}{M}\right),
}
\tag{4}
\]
where \(\kappa(\theta)=\inf_{q\ge1}q\|q\theta\|\).

Since rational affine transformations preserve bad approximability, the restricted sublattice has a positive lower approximation constant if and only if \(\alpha\) does. In particular, for \(\alpha=\pi\), studying one fixed congruence sublattice does not make the bad-approximability problem easier: it is exactly the same question for \((\pi-k)/M\).

A finite partition into residue classes cannot by itself help either. If all relevant primitive integer pairs are partitioned into finitely many classes, the global infimum is the minimum of the classwise infima (after treating non-unit denominator residues separately). A positive lower bound on every class would already imply the desired global lower bound.

## 5. The quantitative obstruction: the paper gives no height control

The strengthened lemma proves that \(\beta\) can be made arbitrarily close to \(k\) in \(\mathbb R\), while satisfying the local congruences. But it gives no upper bound for the reduced denominator \(b\) in terms of \(\varepsilon\), \(M\), or \(k\).

That missing height control is decisive for the lattice flow. If \(\beta=a/b\) is used to perturb a target rational \(p_0/q_0\), the new denominator is \(Q=bq_0\). Its approximation scale contains a factor \(b^2\):
\[
Q\,|P-\alpha Q|
=q_0b^2\,|\beta-\alpha q_0|.
\]
Thus merely making \(|\beta-p_0|\) small does not ensure a good value of \(Q|P-\alpha Q|\); the denominator growth can overwhelm the real improvement.

The paper's later height argument goes in a different direction: its Diophantine predicate propagates odd valuations to many primes, bounds a truncated counting function by the height of a conjugate difference, and then invokes the \(abc\)/Vojta inequality to force bounded height. That mechanism is designed to exclude non-rational elements of a number field, not to bound real rational approximation to the transcendental number \(\pi\).

## 6. A concrete hybrid object, and the exact point where it stalls

For \(\alpha=\pi\), define the congruence-restricted lattice flow
\[
\Lambda_{\pi;k,M}(t)
=
A_\pi(t)L_{k,M},
\qquad
A_\pi(t)=
\begin{pmatrix}e^{-t}&0\\-\pi e^t&e^t\end{pmatrix}.
\]
Its vectors are
\[
\left(qe^{-t},(p-\pi q)e^t\right),
\qquad p-kq\equiv0\pmod M.
\]
After normalizing covolume by \(M^{-1/2}\), this is a unimodular lattice flow. Its shortest-vector invariant is exactly (4), up to the factor \(2\) from the balance-time identity.

The elliptic lemma guarantees special rational vectors in \(L_{k,M}\) with slopes approaching the integer \(k\), and the strengthened lemma makes that real approximation simultaneous with the prescribed \(p\)-adic congruences. But unless one can bound their heights/denominators, these special vectors need not control the shortest-vector lower envelope.

**Conclusion:** the paper supplies a rigorous local-arithmetic generator; the continuous lattice flow supplies the real geometric optimizer. Their naive combination is equivalent to the original bad-approximation problem and does not yet improve it. The next viable research target is a quantitative formal-group construction with an explicit denominator-height bound, or a new global height inequality that controls the elliptic witnesses strongly enough to survive the \(b^2\) loss.

## 7. References and scope

- H. Pasten and X. Vidaux, *On Hilbert's tenth problem for number fields*, arXiv:2611.24748v1 (18 Nov 2026), especially Lemma 4.1, Sections 5–7, and the AI-use acknowledgement in Section 8.
- The real–\(p\)-adic strengthening in Section 2 above is a deduction from their construction plus compact recurrence on \(E(\mathbb R)\); it is not stated as their lemma.
- No result here proves \(\pi\) badly approximable, nor gives a positive lower bound for \(q\|q\pi\|\).
