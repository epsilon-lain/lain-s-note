# A continuous one-parameter lattice flow for the rational approximants of pi

Date: 2026-11-20  
Status: **exact geometric reformulation** of scalar rational approximation as a continuous geodesic flow on unimodular lattices. It is not a proof that pi is badly approximable.

## 1. Why the old clock should be replaced, not merely sped up

The triangular-number clock encodes every positive integer pair (q,p) once, then reads p/q. Its strengths are explicit surjectivity and a single integer input. Its weakness is structural: changing the speed only changes the order or timing of the same exhaustive list. It cannot, by itself, know which rational points are useful for approximating a fixed irrational target.

To remove redundant pairs, the new object should not enumerate pairs at all. It should carry the entire integer lattice continuously and expose only the vectors that are geometrically shortest at the current scale.

## 2. The explicit one-parameter curve

For alpha in R, define
\[
u_\alpha=\begin{pmatrix}1&0\\-\alpha&1\end{pmatrix},
\qquad
g_t=\begin{pmatrix}e^{-t}&0\\0&e^t\end{pmatrix},
\qquad t\in\mathbb R.
\]
Set
\[
\boxed{A_\alpha(t)=g_tu_\alpha
=\begin{pmatrix}e^{-t}&0\\-\alpha e^t&e^t\end{pmatrix}.}
\tag{1}
\]
This is a closed-form, floor-free, nonrecursive, continuous curve in SL_2(R), since det A_alpha(t)=1 for every t.

The associated lattice is
\[
\Lambda_\alpha(t)=A_\alpha(t)\mathbb Z^2.
\]
Changing an integer basis of the same lattice right-multiplies A_alpha(t) by an element of SL_2(Z). Thus the intrinsic trajectory is its projection to the three-dimensional moduli space of unimodular lattices.

For an integer pair (q,p), its lattice vector is exactly
\[
\boxed{A_\alpha(t)\binom qp
=\binom{qe^{-t}}{(p-\alpha q)e^t}.}
\tag{2}
\]
The two coordinates are competing costs: denominator/scale contracts, while approximation error expands. No rounding, floor, rational enumeration, or prescribed list appears in (1)–(2).

## 3. The exact optimization identity

For fixed integers q>=1 and p, write E=|p-alpha q|>0. The squared length of (2) is
\[
L_{p,q}(t)^2=q^2e^{-2t}+E^2e^{2t}.
\]
By AM-GM,
\[
L_{p,q}(t)^2\ge2qE.
\]
Equality occurs at the unique balance time
\[
\boxed{t_{p,q}=\frac12\log\frac qE.}
\tag{3}
\]
Therefore
\[
\boxed{\min_{t\in\mathbb R}L_{p,q}(t)^2
=2q|p-\alpha q|
=2q^2\left|\frac pq-\alpha\right|.}
\tag{4}
\]
This is the key bridge: the minimum length attained by the lattice trajectory by the vector (q,p) is exactly the square root of twice the classical bad-approximation quantity.

For alpha=pi, the relevant quantity is q^2|p/q-pi|. The curve has converted it into a shortest-vector event in a continuously deforming lattice.

## 4. The shortest-vector projection filters redundancy

Define the lattice systole
\[
\lambda_1(\Lambda_\alpha(t))
=\min_{(q,p)\in\mathbb Z^2\setminus\{(0,0)\}}
\left(q^2e^{-2t}+(p-\alpha q)^2e^{2t}\right)^{1/2}.
\tag{5}
\]
This is intrinsic to the current lattice, not a choice of rational enumeration.

At each time t, only the shortest nonzero lattice vectors are active. All other rational pairs are ignored at that time. As t changes, the active vector can switch. The resulting lower envelope is a nonredundant itinerary.

A dominance test removes many pairs. If two pairs satisfy q'<=q and |p'-alpha q'|<=|p-alpha q|, with at least one strict inequality, then the vector for (q',p') is shorter for every t. Hence (q,p) can never be a shortest vector. More generally, only points on the lower convex envelope of
\[
(q^2,\ |p-\alpha q|^2)
\]
can ever be active. Thus the path automatically suppresses dominated fractions, instead of visiting all p/q.

This active set is arithmetic information, not a prewritten schedule. It can be recovered by two-dimensional lattice reduction at each time. Do not claim that every continued-fraction convergent is necessarily active for this particular Euclidean norm without proving that additional statement; the exact safe statement is the shortest-vector/lower-envelope characterization above.

## 5. Bad approximability becomes boundedness of one continuous orbit

Let
\[
\kappa(\alpha)=\inf_{q\ge1}q\|q\alpha\|
=\inf_{\substack{q\ge1\\p\in\mathbb Z}}q|p-\alpha q|,
\]
where ||x|| is distance to the nearest integer. Then the exact identity is
\[
\boxed{
\inf_{t\ge0}\lambda_1(\Lambda_\alpha(t))^2
=2\kappa(\alpha).
}
\tag{6}
\]

**Proof.** For every q>=1 and p, AM-GM gives length squared at least 2q|p-alpha q|, hence at least 2kappa(alpha). Vectors with q=0 have length |p|e^t>=1 for t>=0. Since ||alpha||<=1/2, kappa(alpha)<=1/2, so these vertical vectors cannot lower the infimum below 2kappa.

Conversely, choose q and its nearest integer p, so E=|p-alpha q|=||q alpha||<=1/2<=q. Its balance time (3) is nonnegative. At that time its squared length is exactly 2q||q alpha||. Taking the infimum over q proves the reverse inequality. If alpha is rational, E=0 for some pair and the same identity holds with both sides zero by taking t to infinity.

Consequently
\[
\boxed{
\alpha\text{ is badly approximable}
\iff
\inf_{t\ge0}\lambda_1(\Lambda_\alpha(t))>0.
}
\tag{7}
\]
For alpha=pi, this says that pi is badly approximable exactly when the continuous lattice trajectory stays uniformly away from the cusp. This is a reformulation, not a resolution.

## 6. What this does and does not buy us

**It buys:**
1. A genuinely continuous one-parameter object with no floor or recurrence.
2. A two-dimensional lattice representation whose integer vectors encode every rational p/q.
3. Automatic filtering: only shortest/lower-envelope vectors matter; dominated fractions never appear in the active itinerary.
4. An exact identity between the best attainable vector length and q^2|p/q-alpha|.
5. A natural higher-dimensional generalization to simultaneous approximation and systems of linear forms via diagonal flows on spaces of unimodular lattices.

**It does not buy:**
1. A proof that pi is badly approximable.
2. A new lower bound for q||q pi||.
3. A proof that the active itinerary is finite or periodic.
4. A way to compute the exact shortest vector for pi without certified information about pi.
5. A proof that changing the speed schedule can change which rational approximants are active. Any strictly increasing reparameterization t->phi(t) with range [0,infinity) traces the same orbit and changes only timing.

## 7. The role of speed and the user's original idea

The original idea used an increasing scalar clock to visit integer pairs one by one. Here the speed is replaced by an exponential deformation:
\[
(q,p)\mapsto(qe^{-t},(p-\pi q)e^t).
\]
The first coordinate contracts while the second expands. Each rational pair has its own balance time (3). A useful approximant is one whose balanced vector becomes exceptionally short.

This is closer to a continuous geometric mechanism than a faster enumeration. The parameter does not say "next pair"; it changes the geometry of the entire lattice simultaneously.

If one wants a different speed law, replace t by any increasing u(t) with range [0,infinity):
\[
A_\alpha(t)=
\begin{pmatrix}e^{-u(t)}&0\\-\alpha e^{u(t)}&e^{u(t)}\end{pmatrix}.
\]
The set of lattice states is unchanged. A genuinely new itinerary requires changing the deformation path itself, not merely its speed.

## 8. Relation to the earlier one-parameter clock and infinite-dimensional work

This construction is deliberately different from the earlier one-parameter rational clock: it does not encode (p,q) by a triangular-number index and does not visit every pair. It is also different from the earlier infinite-dimensional kernel embedding: the two-dimensional lattice already captures the exact scalar approximation invariant, so extra coordinates are not needed for this reformulation.

The clock identity
\[
\liminf_{n\to\infty}n|\pi-f(n)|
=\frac{(\pi+1)^2}{2}\liminf_{q\to\infty}q\|q\pi\|
\]
and the lattice identity (4) are two distinct one-parameter views of the same classical quantity. The lattice flow is preferable when the goal is to filter redundant rational pairs; the triangular clock is preferable when the goal is explicit enumeration of every pair.

## 9. Next research questions

1. Give an exact characterization of the active shortest-vector itinerary in terms of continued fractions and possible intermediate convergents for the Euclidean norm in (5).
2. Derive the switching times between consecutive active vectors and determine how their gaps encode partial quotients.
3. Investigate whether the BBP 2-adic information already proved for the special family M=4N imposes any restriction on the geodesic orbit. Such a bridge must be proved; a denominator valuation on one rational family is not automatically a bound on all q.
4. Search for a pi-specific obstruction to cusp excursions. This is the genuine missing theorem.

## 10. References and attribution

The geodesic-flow criterion (7) is classical in Diophantine approximation and homogeneous dynamics; this note is a tailored derivation for the user's one-dimensional-to-lattice proposal, not a claim to have invented the criterion. Standard references include Cassels, *An Introduction to Diophantine Approximation*, and Kleinbock–Margulis, "Flows on homogeneous spaces and Diophantine approximation on manifolds," Annals of Mathematics 148 (1998), 339–360.

## 11. First finite lower-envelope experiment (2026-11-20)

For the finite candidate set \(1\le q\le200000\), the active vectors of the lower envelope were computed at 80 decimal digits. The active list begins
\[
(0,1),(1,3),(7,22),(113,355),(33215,104348),(66317,208341),(99532,312689).
\]
The finite experiment also shows that some continued-fraction convergents are skipped by this Euclidean lower envelope, including \(333/106\) and \(103993/33102\). This is expected: the active set is the lower convex hull of \((q^2,\|q\pi\|^2)\), not the full list of convergents. The last vector is active only relative to the cutoff \(q\le200000\); larger denominators may create later switches.

See [the experiment note](LOWER_ENVELOPE_EXPERIMENT.md), [the reproducible mpmath script](pi_lattice_flow_experiment.py), and [the JSON output](pi_lattice_flow_experiment.json). This is numerical exploration, not an infinite itinerary theorem.
